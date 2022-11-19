
  #!/bin/bash ;
  export PATH=/home/blighek1/bcftools-1.15/bin:$PATH ;
  export PATH=/home/blighek1/htslib-1.15.1/bin:$PATH ;
  export PATH=/home/blighek1/plink_linux_x86_64_20220402/:$PATH ;
  cd /shared/CLL/ ;

  # check for known oncogene drivers for del17p CLL
    #find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v | sort | while read vcf ;
    #do
    #  bcftools view -h "${vcf}" | grep -e "NOTCH1" -e "RPS15" -e "SF3B1" -e "DDX3X" ;
    #done ;

    bcftools view -h WGS/Merge.vcf.gz > WGS/Merge_DriverGenes.vcf ;
    bcftools view WGS/Merge.vcf.gz | grep -e "|NOTCH1|" -e "|RPS15|" -e "|SF3B1|" -e "|DDX3X|" >> WGS/Merge_DriverGenes.vcf ;
    grep -e "^#" -v WGS/Merge_DriverGenes.vcf | wc -l ;
    grep -e "^#" -v WGS/Merge_DriverGenes.vcf | cut -f7 | sort | uniq -c ;
    bcftools view --apply-filters 'PASS' WGS/Merge_DriverGenes.vcf > WGS/Merge_DriverGenes.Filt.vcf ;
    paste <(bcftools view WGS/Merge_DriverGenes.Filt.vcf |\
      awk -F"\t" 'BEGIN {print "Variant\tID"} \
        !/^#/ {print $1":"$2","$4">"$5"\t"$3}') \
      \
      <(cat <(echo "Gene") <(bcftools query -f '%ANN\n' WGS/Merge_DriverGenes.Filt.vcf | cut -f4 -d"|")) \
      \
      <(cat\
        <(echo "Impact")\
        <(bcftools query -f '%ANN\n' WGS/Merge_DriverGenes.Filt.vcf |\
          awk '{if ($0 ~/\|HIGH\|/ || $0 ~ /\|MODERATE\|/) {print "HIGH|MODERATE"} else {print "OTHER"}}')) \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge_DriverGenes.Filt.vcf |\
        awk 'BEGIN {print "nHet"} {print gsub(/0\/1|1\/0/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge_DriverGenes.Filt.vcf |\
        awk 'BEGIN {print "nHomAlt"} {print gsub(/1\/1/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge_DriverGenes.Filt.vcf |\
        awk 'BEGIN {print "nHomRef"} {print gsub(/0\/0/, "")}') \
      \
      <(bcftools view WGS/Merge_DriverGenes.Filt.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge_DriverGenes.Filt.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge_DriverGenes.Filt.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      | sed 's/,\t/\t/g' | sed 's/,$//g' | sed 's/\t/;/g' ;

    bcftools view -h WGS/Merge.vcf.gz > WGS/Merge_DriverGenes_TP53.vcf ;
    bcftools view WGS/Merge.vcf.gz | grep -e "|TP53|" >> WGS/Merge_DriverGenes_TP53.vcf ;
    grep -e "^#" -v WGS/Merge_DriverGenes_TP53.vcf | wc -l ;
    grep -e "^#" -v WGS/Merge_DriverGenes_TP53.vcf | cut -f7 | sort | uniq -c ;
    bcftools view --apply-filters 'PASS' WGS/Merge_DriverGenes_TP53.vcf > WGS/Merge_DriverGenes_TP53.Filt.vcf ;
    paste <(bcftools view WGS/Merge_DriverGenes_TP53.Filt.vcf |\
      awk -F"\t" 'BEGIN {print "Variant\tID"} \
        !/^#/ {print $1":"$2","$4">"$5"\t"$3}') \
      \
      <(cat <(echo "Gene") <(bcftools query -f '%ANN\n' WGS/Merge_DriverGenes_TP53.Filt.vcf | cut -f4 -d"|")) \
      \
      <(cat\
        <(echo "Impact")\
        <(bcftools query -f '%ANN\n' WGS/Merge_DriverGenes_TP53.Filt.vcf |\
          awk '{if ($0 ~/\|HIGH\|/ || $0 ~ /\|MODERATE\|/) {print "HIGH|MODERATE"} else {print "OTHER"}}')) \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge_DriverGenes_TP53.Filt.vcf |\
        awk 'BEGIN {print "nHet"} {print gsub(/0\/1|1\/0/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge_DriverGenes_TP53.Filt.vcf |\
        awk 'BEGIN {print "nHomAlt"} {print gsub(/1\/1/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge_DriverGenes_TP53.Filt.vcf |\
        awk 'BEGIN {print "nHomRef"} {print gsub(/0\/0/, "")}') \
      \
      <(bcftools view WGS/Merge_DriverGenes_TP53.Filt.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge_DriverGenes_TP53.Filt.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge_DriverGenes_TP53.Filt.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      | sed 's/,\t/\t/g' | sed 's/,$//g' | sed 's/\t/;/g' ;
