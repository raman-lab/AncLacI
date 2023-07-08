These scripts are used in sequence to analyze the LGF and DMS NGS data.

Paired-end Illumina sequencing reads are first merged with PEAR (Paired-End Read Merger) using the default software parameters. This software is available online. Sequences are 'cleaned' by deleting the header, the third empty line, and putting the nucleotide sequence on the same line as the quality scores. This is done with a bash script:
	./clean_and_merge_fastq.sh

Poor quality reads are then removed with a C++ script (quality_filter_and_translate.cpp) which is exectuted with a bash script (process_data.sh)
	./process_data.sh

Single variants are then counted with enrichment_scores.R

