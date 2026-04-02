## ----------------------------------------------------------
## 6. Consensus SV, ecDNA and MSI
## ----------------------------------------------------------
```sh
study="xxx"
path_OV="xxx"
path_base=${path_OV}/${study}/

bsub <<EOF
#BSUB -W 240:00
#BSUB -q long
#BSUB -J call_SV_nf
#BSUB -o ${path_base}/call_SV_nf_%J_out.log
#BSUB -e ${path_base}/call_SV_nf_%J_error.log
#BSUB -n 1
#BSUB -M 10
#BSUB -R "rusage[mem=10]"

cd .../code/nextflow/call_SV_consensus_new
module load nextflow
#nextflow clean -f
nextflow run main.nf \
  -profile seadragon \
  -w ${path_base}/log/ \
  --cohort_dir ${path_base} \
  --parts paired_list_header_processed_all.txt \
  --genome hg38 \
  --ref .../reference/GRCh38/v0/Homo_sapiens_assembly38.fasta

EOF
```
