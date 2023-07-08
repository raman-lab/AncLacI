#!/bin/bash

for i in $(cat sample_names)
{
  cp quality_filter_and_translate.cpp ${i}/
  cd ${i}
  g++ -o qft.out quality_filter_and_translate.cpp
  ./qft.out
  
  screen=${i:0:3}

  if [[ ${screen} == "LGF" ]]
  then
    cp ../LGF_ref_lib.csv ./   
  elif [[ ${screen} == "DMS" ]]
  then
    cp ../DMS_ref_lib.csv ./
  fi

  for j in $(awk -F "," '{print $2}' *_ref_lib.csv)
  {
    sed -n "/${j}/p" translations.csv > temp
    c=$(wc -l temp | cut -d' ' -f1)
    echo ${j},${c} >> library_counts.csv 
  }
  cd ..
}

