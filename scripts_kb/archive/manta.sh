
  cd ;
  cd DA0000882/ ;

  export PATH=/home/blighek1/bcftools-1.15/bin:$PATH ;
  
  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/ --recursive |\
    grep -e "manta" | grep -e "candidateSV.vcf.gz" -e "somaticSV.vcf" | grep -e "tumor_only" -v |\
    rev | cut -f1 -d" " | rev |\
      while read file ; do
        outdir=$(echo -e "${file}" | cut -f6 -d"/") ;
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/manta/"${outdir}""/" --dryrun ;
      done ;

  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/ --recursive |\
    grep -e "manta" | grep -e "alignmentStatsSummary.txt" | grep -e "tumor_only" -v |\
    rev | cut -f1 -d" " | rev |\
      while read file ; do
        outdir=$(echo -e "${file}" | cut -f6 -d"/") ;
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/manta/"${outdir}""/" --dryrun ;
      done ;

  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/ --recursive |\
    grep -e "manta" | grep -e "candidateSV.vcf.gz" -e "somaticSV.vcf" | grep -e "tumor_only" -v |\
    rev | cut -f1 -d" " | rev |\
      while read file ; do
        outdir=$(echo -e "${file}" | cut -f6 -d"/") ;
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/manta/"${outdir}""/" ;
      done ;

  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/ --recursive |\
    grep -e "manta" | grep -e "alignmentStatsSummary.txt" | grep -e "tumor_only" -v |\
    rev | cut -f1 -d" " | rev |\
      while read file ; do
        outdir=$(echo -e "${file}" | cut -f6 -d"/") ;
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/manta/"${outdir}""/" ;
      done ;

  # compile output
    find WGS/manta/ -name "somaticSV.vcf.gz" | sort | while read vcf ;
    do
      paste -d "\t" \
        <(echo "${vcf}") \
        <(bcftools view -f 'PASS' "${vcf}" | bcftools stats -1 | grep -e "^SN" | grep -e "number of indels" | cut -f4 | awk '{printf $0"\t"} END {printf "\n"}' | sed 's/\t$//g')
    done |\
      sed 's/WGS\/manta\///g' | sed 's/\/somaticSV\.vcf\.gz//g' > manta_stats.tsv ;

    command=$(find WGS/manta/ -name "somaticSV.vcf.gz" | sort | while read vcf ;     do       paste -d "\t"         <(echo "<(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ""${vcf}"") \\") ; done ) ;
    cat "${command}" | awk '{print $1":"$2":"$3}' | sort | uniq -c | head ;
    cat "${command}" | awk '{print $1":"$2":"$3}' | sort | uniq -c | tail ;
    
