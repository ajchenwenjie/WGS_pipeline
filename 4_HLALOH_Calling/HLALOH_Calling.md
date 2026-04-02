## ----------------------------------------------------------
## 4. HLA LOH
## ----------------------------------------------------------
### Prepare input
```sh
module load R/4.1.0

path_code="xxx"
sample_list="${path_base}/metadata/purity_OV_clonal_BBall.txt"
path_battenberg="${path_base}/CNV/Battenberg/"
path_hlaloh="${path_base}/HLALOH/"

Rscript ${path_code}/HLALOH_input.R -l ${sample_list} -b ${path_battenberg} -o ${path_hlaloh}
```

### Run LOHHLA
```sh
#!/bin/bash
paired_base="${path_base}/metadata/paired_list_header_part"
HLALOH_dir="${path_base}/HLALOH/output/"
HLAmut_dir="${path_base}/HLAmutations/"
solutions_path="${path_base}/HLALOH/input/solutions/"
genome="hg38"
job_num="10"

path_base="xxx"
i=1

  bsub <<EOF
#BSUB -W 240:00
#BSUB -q long
#BSUB -o ${path_base}/HLALOH/HLAHLOH_output_${i}_%J.log
#BSUB -e ${path_base}/HLALOH/HLAHLOH_error_${i}_%J.log
#BSUB -cwd ${path_base}/HLALOH/
#BSUB -u wchen20@mdanderson.org
#BSUB -n 8
#BSUB -M 10
#BSUB -R "rusage[mem=10]"
#BSUB -P HLALOH
#BSUB -J HLALOH_${i}

export NXF_DISABLE_CHECK_LATEST=true
export NXF_OPTS='-Xms2g -Xmx16g'

cd xxx/code/nextflow/call_hlaloh_nf
module load nextflow

nextflow run main.nf \
  -profile seadragon \
  -w ${path_base}/log/ \
  --cohort_dir ${path_base} \
  --parts paired_list_header_processed_all.txt \
  --genome hg38 \
  --hlaloh_dir ${path_base}/HLALOH/output \
  --hlamut_dir ${path_base}/HLAmutations \
  --solutions_path ${path_base}/HLALOH/input/solutions
EOF

##
module load R/4.1.0

path_code="xxx/code/source"

type="OV"
path_output="${path_base}/HLALOH/output/Summary/${type}"
mkdir -p ${path_output}

echo "HLA LOH call done for ${type}"
find ${path_base}/HLALOH/output/${type} \
          -type f -name "*_LOHHLA.csv" \
          -exec grep -l "UnPairedPval_unique" {} + | \
          xargs -I {} cp {} "$path_output"

Rscript ${path_code}/HLALOHtoDF.R -p ${path_base} -t ${type} -o ${path_output}
```
