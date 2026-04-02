## ----------------------------------------------------------
## 2. Mutations
## ----------------------------------------------------------
```sh
path_base="/rsrch5/home/hema_bio-Malignan/OvarianCancerWGS_data/2023_CancerCell_Lahtinen/"

bsub <<EOF
#BSUB -W 240:00
#BSUB -q long
#BSUB -J call_snv_nf
#BSUB -o ${path_base}/call_snv_nf_%J_out.log
#BSUB -e ${path_base}/call_snv_nf_%J_error.log
#BSUB -n 1
#BSUB -M 20
#BSUB -R "rusage[mem=20]"

cd /rsrch6/home/hema_bio-Malignan/wchen20/code/nextflow/call_snv_hla
module load nextflow
#nextflow clean -f
nextflow run main.nf \
  -w /rsrch5/home/hema_bio-Malignan/OvarianCancerWGS_data/2023_CancerCell_Lahtinen/log \
  -profile seadragon \
  --cohort_dir /rsrch5/home/hema_bio-Malignan/OvarianCancerWGS_data/2023_CancerCell_Lahtinen \
  --parts paired_list_header_processed_all.txt \
  --genome hg38 \
  --job_num 20 \
  --ref /rsrch6/home/hema_bio-Malignan/Gu_lab_projects/WenjieChen/reference/GRCh38/v0/Homo_sapiens_assembly38.fasta \

EOF
```
