## Requried columns in sample_list.txt: Patient,Tissue_Type,Sample_ID,FQ_1,FQ_2,Aliged_bam,Sorted_bam,Sorted_dedup_bam,Realign_bam,BQSR_bam
## Novaseq
```sh
study="xxx"
path_OV="xxx"
path_base=${path_OV}/${study}/

  bsub <<EOF
#BSUB -W 240:00
#BSUB -q long
#BSUB -o ${path_base}/WGS_preprocess_output_nf_%J.log
#BSUB -e ${path_base}/WGS_preprocess_error_nf_%J.log
#BSUB -cwd ${path_base}
#BSUB -u wchen20@mdanderson.org
#BSUB -n 12
#BSUB -M 10
#BSUB -R "rusage[mem=10]"
#BSUB -P WGS_preprocess
#BSUB -J WGS_preprocess_nf

cd .../code/nextflow/WGS_preprocess
module load nextflow

nextflow clean -f

nextflow run main.nf \
  -w ${path_base}/work \
  -profile seadragon \
  --sample_list sample_list_newly_downloaded_novaseq20260331.txt \
  --cohort_dir ${path_base}/ \
  --data_dir ${path_base}/ \
  --genome hg38 \
  --ref ${path_base}/reference/GRCh38/v0/Homo_sapiens_assembly38.fasta

EOF
```
