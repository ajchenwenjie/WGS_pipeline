## Pipiline for mutation calling on cell line WES data
```sh
#!/bin/bash
# =============================================================================
# WGS_preprocess_cellline.sh
#
# WGS preprocessing + tumor-only MuTect2 calling (GATK 3.7) for cell lines.
#
# Pipeline per sample:
#   1. Combine FASTQs -> Trimmomatic -> BWA-MEM -> sort  (Sorted_bam)
#   2. Picard MarkDuplicates, duplicates removed          (Sorted_dedup_bam)
#   3. GATK IndelRealigner                                (Realign_bam)
#   4. GATK BaseRecalibrator + PrintReads                 (BQSR_bam)
#   5. GATK MuTect2, tumor-only                           (7_MuTect2/<id>.mutect2.vcf)
#
# If the BQSR BAM already exists, steps 1-4 are skipped and only MuTect2 runs.
# If the MuTect2 VCF already exists, MuTect2 is skipped too, so reruns are safe.
#
# Usage:
#   bash WGS_preprocess_cellline.sh <sample_list> <base_dir> <genome> <job_num> [fasta_file]
#
# Parameters:
#   sample_list  Tab-separated file with a header line. Columns, in this order:
#                  1  Patient           patient ID
#                  2  Tissue_Type       e.g. tumor / cell line
#                  3  Sample_ID         sample ID, used in output file names
#                  4  FQ_1              R1 FASTQ path(s); separate several with commas
#                  5  FQ_2              R2 FASTQ path(s), same order as FQ_1
#                  6  Aliged_bam        unused, but the column must be present
#                  7  Sorted_bam        output path: sorted BAM
#                  8  Sorted_dedup_bam  output path: deduplicated BAM
#                  9  Realign_bam       output path: indel-realigned BAM
#                 10  BQSR_bam          output path: recalibrated BAM (MuTect2 input)
#                Use absolute paths; the script changes into base_dir.
#   base_dir     Working directory. Subfolders 1_FQcombined ... 7_MuTect2 and
#                logs_variant_calling are created here.
#   genome       hg38 or hg19. Selects the reference FASTA and known-sites VCFs
#                (Mills indels, known indels / dbSNP, dbSNP / HapMap).
#   job_num      Number of samples processed in parallel. Each sample uses up to
#                20 threads and 32 GB RAM (MarkDuplicates), so size this to the node.
#   fasta_file   Optional. Reference FASTA overriding the genome default. It must
#                match the genome's known-sites VCFs and have .fai and .dict files.
#
# Example:
#   bash WGS_preprocess_cellline.sh samples.tsv /path/to/project hg38 2
#
# =============================================================================
```

# Reference files
```sh
# ==== Reference Configuration ====
# Genome:  hg38
ref_dir=".../reference/GRCh38/" 
ref="${ref_dir}/v0/Homo_sapiens_assembly38.fasta" # https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pli=1
known_sites1="${ref_dir}/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" # https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pli=1
known_sites2="${ref_dir}/v0/Homo_sapiens_assembly38.known_indels.vcf.gz" # https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pli=1
known_sites3="${ref_dir}/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz" # https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pli=1

# Genome: hg19
ref_dir=".../reference/"
ref="${ref_dir}/phase2_reference/hs37d5.fa" # https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/
known_sites1="${ref_dir}/hg19/Mills_and_1000G_gold_standard.indels.b37.vcf.gz" # ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
known_sites2="${ref_dir}/hg19/dbsnp_138.b37.vcf.gz" # ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/
known_sites3="${ref_dir}/hg19/hapmap_3.3.b37.vcf.gz" # ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/
```

# Exmaple run
```sh
path_code=".../code/" 
sample_list=".../sample_list_test.txt"
base_dir=".../"
genome="hg38"
job_num="20"
ref=".../Homo_sapiens_assembly38.fasta"
bash ${path_code}/WGS_preprocess_cellline.sh ${sample_list} \
                                             ${base_dir} \
                                             ${genome} \
                                             ${job_num} \
                                             ${ref}
```
