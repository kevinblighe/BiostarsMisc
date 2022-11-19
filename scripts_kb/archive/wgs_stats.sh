
  cd /shared/CLL/ ;
  mkdir -p qc_kb ;

  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/Raw/ --recursive |\
    grep -e "stats" | rev | cut -f1 -d" " | rev | \
      while read file ; do
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/Raw/Stats/ --dryrun ;
      done ;

  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/Raw/ --recursive |\
    grep -e "stats" | rev | cut -f1 -d" " | rev | \
      while read file ; do
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/Raw/Stats/ ;
      done ;

  # compile output
    # check alignment between R1 and R2 (should print nothing)
      paste \
        <(find WGS/Raw/Stats/ -name "*_1.fastq.stats" | grep -e "CD19-" | sort) \
        <(find WGS/Raw/Stats/ -name "*_2.fastq.stats" | grep -e "CD19-" | sort) |\
        sed 's/_[12].fastq.stats//g' | awk -F "\t" '$1!=$2' ;

    # unclipped
      paste\
        <(echo "Metric")\
        <(paste\
          <(find WGS/Raw/Stats/ -name "*_1.fastq.stats" | grep -e "CD19-" | sort | cut -f4 -d"/" | cut -f1 -d"_" | awk '{print $0" R1"}')\
          <(find WGS/Raw/Stats/ -name "*_2.fastq.stats" | grep -e "CD19-" | sort | cut -f4 -d"/" | cut -f1 -d"_" | awk '{print $0" R2"}') |\
          awk -F "\t" '{printf $1"\t"$2"\t"} END {printf "\n"}') \
          > qc_kb/wgs_stats.tsv ;
      paste\
        <(find WGS/Raw/Stats/ -name "*_1.fastq.stats" | grep -e "CD19-" | sort) \
        <(find WGS/Raw/Stats/ -name "*_2.fastq.stats" | grep -e "CD19-" | sort) |\
        awk -F "\t" 'NR==1 {printf "paste <(cut -f1 "$1")"}; {printf " <(cut -f2 "$1") <(cut -f2 "$2")"} END {printf "\n"}' \
        > qc_kb/.temp.sh ;
      bash qc_kb/.temp.sh >> qc_kb/wgs_stats.tsv ;
      rm qc_kb/.temp.sh ;

    # clipped
      paste\
        <(echo "Metric")\
        <(paste\
          <(find WGS/Raw/Stats/ -name "*_1.clipped.fastq.stats" | grep -e "CD19-" | sort | cut -f4 -d"/" | cut -f1 -d"_" | awk '{print $0" R1"}')\
          <(find WGS/Raw/Stats/ -name "*_2.clipped.fastq.stats" | grep -e "CD19-" | sort | cut -f4 -d"/" | cut -f1 -d"_" | awk '{print $0" R2"}') |\
          awk -F "\t" '{printf $1"\t"$2"\t"} END {printf "\n"}') \
          > qc_kb/wgs_stats.clipped.tsv ;
      paste\
        <(find WGS/Raw/Stats/ -name "*_1.clipped.fastq.stats" | grep -e "CD19-" | sort) \
        <(find WGS/Raw/Stats/ -name "*_2.clipped.fastq.stats" | grep -e "CD19-" | sort) |\
        awk -F "\t" 'NR==1 {printf "paste <(cut -f1 "$1")"}; {printf " <(cut -f2 "$1") <(cut -f2 "$2")"} END {printf "\n"}' \
        > qc_kb/.temp.sh ;
      bash qc_kb/.temp.sh >> qc_kb/wgs_stats.clipped.tsv ;
      rm qc_kb/.temp.sh ;
