#!/bin/bash
# Exit on error
set -e

# ==== Input Parameters ====
sample_list="$1"
base_dir="$2"
genome="$3"
job_num="$4"
ref="$5"

# Validate input parameters
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <sample_list> <base_dir> <genome> <job_num> <fasta_file>"
    echo "  sample_list: path to sample list file"
    echo "  base_dir: base directory for processing"
    echo "  genome: hg38 or hg19"
    echo "  job_num: number of parallel jobs"
    echo "  ref: reference fasta file"
    exit 1
fi

# ==== Parallelism Control ====
# Counts running background sample jobs and waits until a slot is free
wait_for_free_samples() {
    while (( $(jobs -pr | wc -l) >= job_num )); do
        sleep 1
    done
}

# ==== Reference Configuration ====
if [[ "$genome" == "hg38" ]]; then
    echo "## Using hg38 reference files"
    ref_dir=".../reference/GRCh38/" # set reference path
    #ref="${ref_dir}/v0/Homo_sapiens_assembly38.fasta"
    known_sites1="${ref_dir}/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    known_sites2="${ref_dir}/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
    known_sites3="${ref_dir}/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

elif [[ "$genome" == "hg19" ]]; then
    echo "## Using hg19 reference files"
    ref_dir=".../reference/" # set reference path
    #ref="${ref_dir}/phase2_reference/hs37d5.fa"
    known_sites1="${ref_dir}/hg19/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    known_sites2="${ref_dir}/hg19/dbsnp_138.b37.vcf.gz"
    known_sites3="${ref_dir}/hg19/hapmap_3.3.b37.vcf.gz"

else
    echo "Error: genome must be 'hg38' or 'hg19'" >&2
    exit 1
fi

## Directories
cd "${base_dir}"
fastq_dir="${base_dir}/1_FQcombined"
clean_dir="${base_dir}/2_Clean"
alig_dir="${base_dir}/3_Alignment"
sort_dir="${base_dir}/4_Sort"
rmdup_dir="${base_dir}/5_Rmdup"
coclean_dir="${base_dir}/6_Coclean"
mutect2_dir="${base_dir}/7_MuTect2"
mutect2_log_dir="${base_dir}/logs_variant_calling"
mkdir -p "$fastq_dir" "$clean_dir" "$alig_dir" "$sort_dir" "$rmdup_dir" "$coclean_dir" \
         "$mutect2_dir" "$mutect2_log_dir"

## Helper functions
# Check BAM completeness
check_bam_completeness() {
    bam="$1"
    if ! samtools quickcheck -v "$bam"; then
        echo "Incomplete or corrupt BAM detected: $bam. Deleting..."
        rm -f "$bam"
        rm -f "${bam}.bai"
        exit 1
    fi
}

# Run MuTect2 (GATK 3.7, tumor-only) on the BQSR BAM
run_mutect2() {
    local id="$1"
    local BQSR_bam="$2"

    local out_vcf="${mutect2_dir}/${id}.mutect2.vcf"
    local tmp_vcf="${mutect2_dir}/${id}.mutect2.tmp.vcf"
    local log="${mutect2_log_dir}/${id}_mutect2.log"

    # Skip completed output
    if [[ -s "$out_vcf" ]] && grep -q '^#CHROM' "$out_vcf"; then
        echo "Skipping MuTect2 for ${id}: VCF already exists."
        return 0
    fi

    if [[ ! -s "$BQSR_bam" ]]; then
        echo "Error: BQSR BAM not found for ${id}: ${BQSR_bam}" >&2
        return 1
    fi

    # samtools and GATK may not be loaded if earlier steps were skipped
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate bwa_samtools-0.7.17_1.19.2
    module load gatk3/3.7

    if ! samtools quickcheck -v "$BQSR_bam"; then
        echo "Error: BQSR BAM is incomplete or corrupt for ${id}: ${BQSR_bam}" >&2
        return 1
    fi

    # Accept either BAM index naming convention
    if [[ ! -f "${BQSR_bam}.bai" && ! -f "${BQSR_bam%.bam}.bai" ]]; then
        echo "Creating BAM index for ${id}"
        samtools index -@ 20 "$BQSR_bam"
    fi

    echo "Running MuTect2 for ${id}"
    rm -f "$tmp_vcf"

    if java -Xmx16g -jar "$gatk3_37" \
        -T MuTect2 \
        -R "$ref" \
        -I:tumor "$BQSR_bam" \
        --dbsnp "$known_sites3" \
        -o "$tmp_vcf" \
        > "$log" 2>&1; then

        if [[ -s "$tmp_vcf" ]] && grep -q '^#CHROM' "$tmp_vcf"; then
            mv "$tmp_vcf" "$out_vcf"
            [[ -f "${tmp_vcf}.idx" ]] && mv "${tmp_vcf}.idx" "${out_vcf}.idx"
            echo "MuTect2 completed for ${id}: ${out_vcf}"
        else
            echo "Error: MuTect2 output missing or invalid for ${id}. Check ${log}" >&2
            rm -f "$tmp_vcf" "${tmp_vcf}.idx"
            return 1
        fi
    else
        echo "Error: MuTect2 failed for ${id}. Check ${log}" >&2
        rm -f "$tmp_vcf" "${tmp_vcf}.idx"
        return 1
    fi
}

## Process samples
tail -n +2 "${sample_list}" | tr -d '\r' | while IFS=$'\t' read -r Patient Tissue_Type Sample_ID FQ_1 FQ_2 Aliged_bam Sorted_bam Sorted_dedup_bam Realign_bam BQSR_bam || [ -n "$Patient" ]; do
    (
        # Clean variables of carriage returns
        Patient=$(echo "$Patient" | tr -d '\r')
        Tissue_Type=$(echo "$Tissue_Type" | tr -d '\r')
        Sample_ID=$(echo "$Sample_ID" | tr -d '\r')
        FQ_1=$(echo "$FQ_1" | tr -d '\r')
        FQ_2=$(echo "$FQ_2" | tr -d '\r')
        Sorted_bam=$(echo "$Sorted_bam" | tr -d '\r')
        Sorted_dedup_bam=$(echo "$Sorted_dedup_bam" | tr -d '\r')
        Realign_bam=$(echo "$Realign_bam" | tr -d '\r')
        BQSR_bam=$(echo "$BQSR_bam" | tr -d '\r')

        id="${Sample_ID}"
        echo "Running sample for ${id}"

        if [[ -s "${BQSR_bam}" ]]; then
            echo "BQSR BAM exists for ${id}: ${BQSR_bam}"
            echo "Skipping preprocessing and running MuTect2 directly."
        else
        echo "BQSR BAM does not exist for ${id}. Starting preprocessing."

        # 1. Trimmomatic and alignment
        if { [ -f "${Sorted_bam}.bai" ] || [ -f "${Sorted_dedup_bam}.bai" ] || [ -f "${Realign_bam}.bai" ] || [ -f "${BQSR_bam}.bai" ]; } && [ ! -f "${clean_dir}/${id}_R2_clean.fastq.gz" ]; then
            echo "Skipping trimming/alignment for ${id}, processed BAM already exists."
        else
            FQ_1=$(printf '%s' "$FQ_1" | tr -d '\r' | sed 's/^"//; s/"$//')
            FQ_2=$(printf '%s' "$FQ_2" | tr -d '\r' | sed 's/^"//; s/"$//')

            # Split on comma, clean each entry, load into arrays
            mapfile -t R1_FILES < <(printf '%s' "$FQ_1" | tr ',' '\n' | sed 's/^ *//; s/ *$//' | sed 's/^"//; s/"$//')
            mapfile -t R2_FILES < <(printf '%s' "$FQ_2" | tr ',' '\n' | sed 's/^ *//; s/ *$//' | sed 's/^"//; s/"$//')

            # Sanity checks
            if [ "${#R1_FILES[@]}" -ne "${#R2_FILES[@]}" ]; then
               echo "ERROR: R1/R2 count mismatch for ${id}"
               exit 1
            fi

            for f in "${R1_FILES[@]}" "${R2_FILES[@]}"; do
                if [ ! -f "$f" ]; then
                    echo "ERROR: FASTQ not found: [$f]"
                    exit 1
                fi
            done

            # Combine FASTQs
            cat "${R1_FILES[@]}" > "${fastq_dir}/${id}_combined_R1.fastq.gz"
            cat "${R2_FILES[@]}" > "${fastq_dir}/${id}_combined_R2.fastq.gz"

            eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
            conda activate trimmomatic-0.39

            trimmomatic PE -threads 20 -phred33 \
                "${fastq_dir}/${id}_combined_R1.fastq.gz" "${fastq_dir}/${id}_combined_R2.fastq.gz" \
                "${clean_dir}/${id}_R1_clean.fastq.gz" "${clean_dir}/${id}_R1_unpaired.fastq.gz" \
                "${clean_dir}/${id}_R2_clean.fastq.gz" "${clean_dir}/${id}_R2_unpaired.fastq.gz" \
                ILLUMINACLIP:/risapps/rhel7/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

            rm "${fastq_dir}/${id}_combined_R1.fastq.gz"
            rm "${fastq_dir}/${id}_combined_R2.fastq.gz"
            rm "${clean_dir}/${id}_R1_unpaired.fastq.gz"
            rm "${clean_dir}/${id}_R2_unpaired.fastq.gz"

            eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
            conda activate bwa_samtools-0.7.17_1.19.2

            bwa mem -t 20 -R "@RG\tPL:ILLUMINA\tID:${id}\tSM:${id}" "${ref}" \
                "${clean_dir}/${id}_R1_clean.fastq.gz" "${clean_dir}/${id}_R2_clean.fastq.gz" | \
            samtools view -S -b - | \
            samtools sort -O bam -@ 10 -o "${Sorted_bam}" -

            check_bam_completeness "${Sorted_bam}"
            samtools index -@ 20 "${Sorted_bam}"

            if [ -f "${Sorted_bam}.bai" ]; then
               echo "[INFO] BAM index created successfully. Removing original FASTQs for ${id}."
               rm "${clean_dir}/${id}_R1_clean.fastq.gz"
               rm "${clean_dir}/${id}_R2_clean.fastq.gz"

               for f in "${R1_FILES[@]}" "${R2_FILES[@]}"; do
                   if [ -f "$f" ]; then
                       #rm "$f"
                       echo "Keep $f"
                   else
                       echo "[WARN] FASTQ already missing: $f"
                   fi
               done
            else
               echo "[ERROR] BAM index missing for ${id}. FASTQs NOT removed."
               exit 1
            fi
        fi

        # 2. MARK DUPLICATES
        if { [ -f "${Sorted_dedup_bam}.bai" ] || [ -f "${Realign_bam}.bai" ] || [ -f "${BQSR_bam}.bai" ]; } && [ ! -f "${Sorted_bam}.bai" ]; then
            echo "Skipping duplicate marking for ${id}, deduplicated BAM already exists."
        else
            eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
            conda activate picard-2.27.4

            mkdir -p "${rmdup_dir}/${id}/tmp"

            picard -Xmx32g -Djava.io.tmpdir="${rmdup_dir}/${id}/tmp" MarkDuplicates \
                --INPUT "${Sorted_bam}" \
                --OUTPUT "${Sorted_dedup_bam}" \
                --METRICS_FILE "${rmdup_dir}/${id}.sorted_markdup_metrics.txt" \
                --REMOVE_DUPLICATES true \
                --VALIDATION_STRINGENCY LENIENT \
                > "${rmdup_dir}/${id}_rmdup.log" 2>&1

            eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
            conda activate bwa_samtools-0.7.17_1.19.2

            check_bam_completeness "${Sorted_dedup_bam}"
            samtools index -@ 20 "${Sorted_dedup_bam}"

            if [ -f "${Sorted_dedup_bam}.bai" ]; then
                rm "${Sorted_bam}" "${Sorted_bam}.bai"
            fi
        fi

        # 3. Co-cleaning Workflow
        if { [ -f "${Realign_bam}.bai" ] || [ -f "${BQSR_bam}.bai" ]; } && [ ! -f "${Sorted_dedup_bam}.bai" ]; then
            echo "Skipping realignment for ${id}, realigned BAM already exists."
        else
            module load gatk3/3.7

            java -jar "$gatk3_37" -T RealignerTargetCreator \
                                  -R "${ref}" \
                                  -I "${Sorted_dedup_bam}" \
                                  -known "${known_sites1}" \
                                  -known "${known_sites2}" \
                                  -known "${known_sites3}" \
                                  -o "${coclean_dir}/${id}.IndelRealigner.intervals" > "${coclean_dir}/${id}_REALIGNT.log" 2>&1

            java -jar "$gatk3_37" -T IndelRealigner \
                                  -R "${ref}" \
                                  -I "${Sorted_dedup_bam}" \
                                  -known "${known_sites1}" \
                                  -known "${known_sites2}" \
                                  -known "${known_sites3}" \
                                  --targetIntervals "${coclean_dir}/${id}.IndelRealigner.intervals" \
                                  -o "${Realign_bam}" > "${coclean_dir}/${id}_INDELREALIGNER.log" 2>&1

            eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
            conda activate bwa_samtools-0.7.17_1.19.2

            check_bam_completeness "${Realign_bam}"
            samtools index -@ 20 "${Realign_bam}"

            if [ -f "${Realign_bam}.bai" ]; then
                rm "${Sorted_dedup_bam}" "${Sorted_dedup_bam}.bai"
            fi
        fi

        # 4. Base Recalibrator
        if [ -f "${BQSR_bam}.bai" ] && [ ! -f "${Realign_bam}.bai" ]; then
            echo "Skipping BQSR for ${id}, BQSR BAM already exists."
        else
            module load gatk3/3.7

            java -jar "$gatk3_37" -T BaseRecalibrator -R "${ref}" -I "${Realign_bam}" \
                --knownSites "${known_sites1}" --knownSites "${known_sites2}" --knownSites "${known_sites3}" \
                -o "${coclean_dir}/${id}_recal_data.table" > "${coclean_dir}/${id}_BASERECALIBRATOR.log" 2>&1

            java -jar "$gatk3_37" -T PrintReads -R "${ref}" -I "${Realign_bam}" \
                --BQSR "${coclean_dir}/${id}_recal_data.table" -o "${BQSR_bam}" > "${coclean_dir}/${id}_PRINTREADS.log" 2>&1

            eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
            conda activate bwa_samtools-0.7.17_1.19.2

            check_bam_completeness "${BQSR_bam}"
            samtools index -@ 20 "${BQSR_bam}"

            if [ -f "${BQSR_bam}.bai" ]; then
                rm "${Realign_bam}" "${Realign_bam}.bai"
            fi
        fi

        fi   # end of preprocessing (skipped when BQSR BAM exists)

        # 5. MuTect2 variant calling on the BQSR BAM
        module load gatk3/3.7
        run_mutect2 "${id}" "${BQSR_bam}"

    ) &

    wait_for_free_samples

done

# Wait for all background jobs to complete
wait