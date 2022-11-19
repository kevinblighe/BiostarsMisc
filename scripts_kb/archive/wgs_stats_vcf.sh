
  cd /shared/CLL/ ;
  mkdir -p qc_kb ;

  export PATH=/home/blighek1/bcftools-1.15/bin:$PATH ;

  echo -e "SNVs\tInDels\tpc > 40x" > qc_kb/pull_wgs_stats_vcf.tsv ;
  find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v | sort | while read vcf ;
  do
    bcftools view -h "${vcf}" > qc_kb/.test.vcf ;
    bcftools view -f 'PASS' "${vcf}" | grep -e "protein_coding" >> qc_kb/.test.vcf ;

    paste -d "\t" \
      <(echo "${vcf}") \
      <(bcftools view -f 'PASS' qc_kb/.test.vcf | bcftools stats -1 | grep -e "^SN" | grep -e "number of SNPs" -e "number of indels" | cut -f4 | awk '{printf $0"\t"} END {printf "\n"}' | sed 's/\t$//g') \
      <(bcftools view -f 'PASS' qc_kb/.test.vcf | bcftools stats -1 -d 0,39,39 | tail -1 | cut -f7) ;
done |\
    sed 's/WGS\/ProcessedData\/[0-9]*\/dbnsfp_annotation\///g' | sed 's/\.vcf\.gz//g' | sed 's/ \+/\t/g' >> qc_kb/pull_wgs_stats_vcf.tsv ;
  rm qc_kb/.test.vcf ;

  echo -e "SNVs\tInDels\tpc > 100x" > qc_kb/pull_wgs_stats_vcf_100X.tsv ;
  find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v | sort | while read vcf ;
  do
    bcftools view -h "${vcf}" > qc_kb/.test.vcf ;
    bcftools view -f 'PASS' "${vcf}" -i 'DP>100' | grep -e "protein_coding" >> qc_kb/.test.vcf ;

    paste -d "\t" \
      <(echo "${vcf}") \
      <(bcftools view -f 'PASS' qc_kb/.test.vcf | bcftools stats -1 | grep -e "^SN" | grep -e "number of SNPs" -e "number of indels" | cut -f4 | awk '{printf $0"\t"} END {printf "\n"}' | sed 's/\t$//g') \
      <(bcftools view -f 'PASS' qc_kb/.test.vcf | bcftools stats -1 -d 0,99,99 | tail -1 | cut -f7) ;
done |\
    sed 's/WGS\/ProcessedData\/[0-9]*\/dbnsfp_annotation\///g' | sed 's/\.vcf\.gz//g' | sed 's/ \+/\t/g' >> qc_kb/pull_wgs_stats_vcf_100X.tsv ;
  rm qc_kb/.test.vcf ;

  echo -e "SNVs\tInDels\tpc > 40x" > qc_kb/pull_wgs_stats_vcf_40X.tsv ;
  find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v | sort | while read vcf ;
  do
    bcftools view -h "${vcf}" > qc_kb/.test.vcf ;
    bcftools view -f 'PASS' "${vcf}" -i 'DP>40' | grep -e "protein_coding" >> qc_kb/.test.vcf ;

    paste -d "\t" \
      <(echo "${vcf}") \
      <(bcftools view -f 'PASS' qc_kb/.test.vcf | bcftools stats -1 | grep -e "^SN" | grep -e "number of SNPs" -e "number of indels" | cut -f4 | awk '{printf $0"\t"} END {printf "\n"}' | sed 's/\t$//g') \
      <(bcftools view -f 'PASS' qc_kb/.test.vcf | bcftools stats -1 -d 0,39,39 | tail -1 | cut -f7) ;
done |\
    sed 's/WGS\/ProcessedData\/[0-9]*\/dbnsfp_annotation\///g' | sed 's/ \+/\t/g' >> qc_kb/pull_wgs_stats_vcf_40X.tsv ;
  rm qc_kb/.test.vcf ;

