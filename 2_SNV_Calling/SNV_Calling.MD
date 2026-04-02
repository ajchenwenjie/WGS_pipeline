## ----------------------------------------------------------
## 2. Mutations
## ----------------------------------------------------------
```sh
path_base="/xxx/"

bsub <<EOF
#BSUB -W 240:00
#BSUB -q long
#BSUB -J call_snv_nf
#BSUB -o ${path_base}/call_snv_nf_%J_out.log
#BSUB -e ${path_base}/call_snv_nf_%J_error.log
#BSUB -n 1
#BSUB -M 20
#BSUB -R "rusage[mem=20]"

cd .../code/nextflow/call_snv_hla
module load nextflow
#nextflow clean -f
nextflow run main.nf \
  -w ${path_base}/log \
  -profile seadragon \
  --cohort_dir ${path_base} \
  --parts paired_list_header_processed_all.txt \
  --genome hg38 \
  --job_num 20 \
  --ref .../reference/GRCh38/v0/Homo_sapiens_assembly38.fasta \

EOF
```
