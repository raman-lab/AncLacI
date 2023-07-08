#!/bin/bash
pear=${PATH_TO_PEAR}

for i in $(cat sample_names)
{
  mkdir ${i}
  mv ${i}_*.fastq.gz ${i} 

  cd ${i}
  ${pear} -f *_R1_001.fastq.gz -r *_R2_001.fastq.gz -o combined
  cp combined.assembled.fastq combined.fastq
  sed -i '/@M07074/d' combined.fastq
  sed -i '2~3d' combined.fastq
  sed -i 'N;s/\n/ /' combined.fastq
  cd ..
}

