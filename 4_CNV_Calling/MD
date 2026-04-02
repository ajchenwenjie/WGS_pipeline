## ----------------------------------------------------------
## 5. CNV
## ----------------------------------------------------------
### Battenberg - new version
```sh
## Prepare input
module load R/4.1.0

type="OV"

# Required header for paired_list: id tumour_id normal_id tumour_inbam normal_inbam tumour_bam normal_bam tumour_header normal_header type sex
paired_list="${path_base}/metadata/paired_list_header_processed_all.txt"
paired_list_bb=${path_base}/metadata/paired_list_header_BB.txt
snv_pre=${path_base}/Variantcalling/annotation/annovar/${type}/snv/
path_bam="non"
suffix="_snv.anno.hg38_multianno"

Rscript /rsrch6/home/hema_bio-Malignan/wchen20/code/source/BB_prepare.R -i ${paired_list} \
                                                                        -o ${paired_list_bb} \
                                                                        -b ${path_bam} \
                                                                        -p ${snv_pre} \
                                                                        -s ${suffix}

Rscript .../Battenberg/BB/Battenberg_pipeline.R ${paired_list_bb} ${path_base}/CNV/Battenberg/${type}

# Set base directory
#!/bin/bash
cd "${path_base}" || exit 1

for type in OV; do

  bsub <<EOF
#BSUB -W 240:00
#BSUB -q long
#BSUB -o ${path_base}/Battenberg_BB_%J_output.log
#BSUB -e ${path_base}/Battenberg_BB_%J_error.log
#BSUB -cwd ${path_base}
#BSUB -u xxx
#BSUB -N
#BSUB -n 1
#BSUB -M 5
#BSUB -R "rusage[mem=5]"
#BSUB -J OV_battenberg_${type}
#BSUB -P OV_battenberg_${type}

Dir="${path_base}/CNV/Battenberg/${type}"
cd "\$Dir" || exit 1

for i in */; do
  i=\${i%/}  # remove trailing slash
  cd "\$Dir/\$i" || continue

  # Skip if already completed
  if ls Refit_1/*_IDcard.pdf 1>/dev/null 2>&1; then
    echo "Skip \$i — Task Completed"
    continue
  fi

  joblist=\$(bjobs -u \$USER | tail -n +2 | wc -l)
  while [ \$joblist -gt 200 ]; do
    sleep 10m
    joblist=\$(bjobs -u \$USER | tail -n +2 | wc -l)
  done

  ./*_submission.sh
  echo "Submitted \$i"
done
EOF

done

### Summary
type="OV"
Dir="${path_base}/CNV/Battenberg/${type}"
find ${Dir} -type f -name "*IDcard.Rdata"|wc -l

module load R/4.1.0

type="OV"
path_batternbegout="${path_base}/CNV/Battenberg/${type}/"
path_griticin="${path_base}/GRITIC/input/${type}/"
metadata="${path_base}/metadata/paired_list_header_processed_all.txt"

Rscript .../code/source/BB_summary.R -p ${path_base} \
                                                                        -b ${path_batternbegout} \
                                                                        -g ${path_griticin} \
                                                                        -t ${type} \
                                                                        -f ${metadata}
```
