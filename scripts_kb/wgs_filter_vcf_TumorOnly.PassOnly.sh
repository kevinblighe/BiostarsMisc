
  #!/bin/bash ;
  export PATH=/home/blighek1/bcftools-1.15/bin:$PATH ;
  export PATH=/home/blighek1/htslib-1.15.1/bin:$PATH ;
  export PATH=/home/blighek1/plink_linux_x86_64_20220402/:$PATH ;
  cd /shared/CLL/ ;
  mkdir -p WGS/FilteredData.TumorOnly.PassOnly/ ;



  # determine which 'tumor-only' samples are matched to the T-N analysis
    # check that the tumour 'tumour-only' samples are matched to the normal 'tumour-only' samples
      # T
        find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
          sort | grep -e "-T" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort ;

      # N
        find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
          sort | grep -e "-N" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort ;

      # T but no N
        paste \
          <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
            sort | grep -e "-T" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort)\
          <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
            sort | grep -e "-N" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort) |\
        awk -F "\t" '{a[$1]=$1; b[$2]=$2} END {for (id in a) {if (!b[id]) print id}}' ;

      # N but no T
        paste \
          <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
            sort | grep -e "-T" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort)\
          <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
            sort | grep -e "-N" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort) |\
        awk -F "\t" '{a[$1]=$1; b[$2]=$2} END {for (id in b) {if (!a[id]) print id}}' ;

    # T-N paired samples
      paste \
        <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v |\
          sort | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort) \
        <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
          sort | grep -e "-T" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort) |\
        awk -F "\t" '{a[$1]=$1; b[$2]=$2} END {for (id in a) {if (!b[id]) print id}}' ;
      paste \
        <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v |\
          sort | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort) \
        <(find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" |\
          sort | grep -e "-N" | cut -f5 -d"/" | cut -f1 -d"-" | cut -f1 -d"_" | sort) |\
        awk -F "\t" '{a[$1]=$1; b[$2]=$2} END {for (id in a) {if (!b[id]) print id}}' ;



  # run the main filtration process: Tumor 'tumor-only'
    # exclude tumor samples 2019 3139 (from '# T but no N', above)
    find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" | sort | grep -e "-T" | grep -e "2019" -v | grep -e "3139" -v | while read vcf ;
    do
      echo -e "--Input VCF is:\t""${vcf}" ;

      out=$(echo "${vcf}" | sed 's/\.vcf\.gz/\.filt\.vcf\.gz/g' | sed 's/ProcessedData\/[0-9]*\/dbnsfp_annotation_tumor_only/FilteredData.TumorOnly.PassOnly/g') ;
      echo -e "--Output VCF is:\t""${out}" ;

      # Filter to main contigs only
        # bcftools filter -r "chr1...
      # Filter to variants with PASS in the FILTER column, and left-align indels / MNPs
        # --apply-filters 'PASS' # NOT USED
        # bcftools norm -f refgenome/hg38.fa
      # Remove any variants that are now wildtype, ./. or 0/0
        # --min-ac 1
      # Remove variants based on read depth
      # Remove variants based on AF in samples
        # TUMOR_depth >= 10
        # TUMOR_AF > 0.05
        # TUMOR_alt_depth >= 3
        # SUM(FORMAT/AD[0:*])>=10 && FORMAT/AF[0:*]>0.05 && FORMAT/AD[0:1]>=3
      # Remove common variants
        # as.numeric(dbNSFP_ExAC_AF) < 0.01 | is.na(as.numeric(dbNSFP_ExAC_AF))
        # as.numeric(dbNSFP_gnomAD_genomes_AF) < 0.01 | is.na(as.numeric(dbNSFP_gnomAD_genomes_AF))
        # as.numeric(dbNSFP_gnomAD_exomes_AF) < 0.01 | is.na(as.numeric(dbNSFP_gnomAD_exomes_AF))
        # INFO/dbNSFP_ExAC_AF<0.01 || INFO/dbNSFP_gnomAD_genomes_AF<0.01 || INFO/dbNSFP_gnomAD_exomes_AF<0.01
        bcftools filter -r "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chrX,chrY" "${vcf}" |\
          bcftools view --apply-filters 'PASS' --min-ac 1 | bcftools norm -f refgenome/hg38.fa |\
            bcftools view --include 'SUM(FORMAT/AD[0:*])>=10 && FORMAT/AF[0:*]>0.05 && FORMAT/AD[0:1]>=3' |\
              bcftools filter --include 'INFO/dbNSFP_ExAC_AF<0.01 || INFO/dbNSFP_ExAC_AF="."' |\
                bcftools filter --include 'INFO/dbNSFP_gnomAD_genomes_AF<0.01 || INFO/dbNSFP_gnomAD_genomes_AF="."' |\
                  bcftools filter --include 'INFO/dbNSFP_gnomAD_exomes_AF<0.01 || INFO/dbNSFP_gnomAD_exomes_AF="."' \
                    > "${vcf}".tmp ;

      # filter out unexpectedly high read depth variants
        # TUMOR_depth <= 3*mean(TUMOR_depth, na.rm=T)
        tumor_avg_3x=$(bcftools query -f '[%SAMPLE=%AD\t]\n' "${vcf}".tmp | cut -f1 | cut -f2 -d"=" | awk -F "," 'BEGIN {a=0} {a+=($1 + $2)} END {print 3*(a/NR)}') ;
        bcftools view --include "SUM(FORMAT/AD[0:*])<=$tumor_avg_3x" "${vcf}".tmp \
          > "${vcf}".tmp.tmp ;

      # include tumor sample only, and only include SNPs/SNVs
        tumor=$(bcftools view -h "${vcf}".tmp.tmp | grep -e "#tumor_sample" | cut -f2 -d"=") ;
        echo -e "--Tumor sample ID is:\t""${tumor}"". Filtering..." ;
        bcftools view --samples "${tumor}" --types snps -Oz "${vcf}".tmp.tmp > "${out}" ;
        tabix -p vcf "${out}" ;
        echo -e "Done." ;
        genotypes=$(bcftools view "${out}" | grep -e "^#" -v | cut -f10 | cut -f1 -d":" | sort | uniq) ;
        echo -e "--genotypes after filtering include:\t""${genotypes}" ;
        echo -e "--total variants=\t""$(bcftools view "${out}" | grep -e "^#" -v | wc -l)" ;
        echo -e "\n" ;

      rm "${vcf}".tmp "${vcf}".tmp.tmp ;
    done ;

  # run the main filtration process: Normal 'tumor-only'
    find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" | sort | grep -e "-N" | while read vcf ;
    do
      echo -e "--Input VCF is:\t""${vcf}" ;

      out=$(echo "${vcf}" | sed 's/\.vcf\.gz/\.filt\.vcf\.gz/g' | sed 's/ProcessedData\/[0-9]*\/dbnsfp_annotation_tumor_only/FilteredData.TumorOnly.PassOnly/g') ;
      echo -e "--Output VCF is:\t""${out}" ;

      # Filter to main contigs only
        # bcftools filter -r "chr1...
      # Filter to variants with PASS in the FILTER column, and left-align indels / MNPs
        # --apply-filters 'PASS' # NOT USED
        # bcftools norm -f refgenome/hg38.fa
      # Remove any variants that are now wildtype, ./. or 0/0
        # --min-ac 1
      # Remove variants based on read depth
      # Remove variants based on AF in samples
        # NORMAL_depth >= 10
        # NORMAL_AF < 0.05
        # SUM(FORMAT/AD[0:*])>=10 && FORMAT/AF[0:*]>0.25
      # Remove common variants
        # as.numeric(dbNSFP_ExAC_AF) < 0.01 | is.na(as.numeric(dbNSFP_ExAC_AF))
        # as.numeric(dbNSFP_gnomAD_genomes_AF) < 0.01 | is.na(as.numeric(dbNSFP_gnomAD_genomes_AF))
        # as.numeric(dbNSFP_gnomAD_exomes_AF) < 0.01 | is.na(as.numeric(dbNSFP_gnomAD_exomes_AF))
        # INFO/dbNSFP_ExAC_AF<0.01 || INFO/dbNSFP_gnomAD_genomes_AF<0.01 || INFO/dbNSFP_gnomAD_exomes_AF<0.01
        bcftools filter -r "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chrX,chrY" "${vcf}" |\
          bcftools view --apply-filters 'PASS' --min-ac 1 | bcftools norm -f refgenome/hg38.fa |\
            bcftools view --include 'SUM(FORMAT/AD[0:*])>=10 && FORMAT/AF[0:*]>0.25' |\
              bcftools filter --include 'INFO/dbNSFP_ExAC_AF<0.01 || INFO/dbNSFP_ExAC_AF="."' |\
                bcftools filter --include 'INFO/dbNSFP_gnomAD_genomes_AF<0.01 || INFO/dbNSFP_gnomAD_genomes_AF="."' |\
                  bcftools filter --include 'INFO/dbNSFP_gnomAD_exomes_AF<0.01 || INFO/dbNSFP_gnomAD_exomes_AF="."' \
                    > "${vcf}".tmp ;

      # filter out unexpectedly high read depth variants
        # NORMAL_depth <= 3*mean(NORMAL_depth,na.rm=T)
        normal_avg_3x=$(bcftools query -f '[%SAMPLE=%AD\t]\n' "${vcf}".tmp | cut -f1 | cut -f2 -d"=" | awk -F "," 'BEGIN {a=0} {a+=($1 + $2)} END {print 3*(a/NR)}') ;
        bcftools view --include "SUM(FORMAT/AD[0:*])<=$normal_avg_3x" "${vcf}".tmp \
          > "${vcf}".tmp.tmp ;

      # include normal sample only, and only include SNPs/SNVs
        normal=$(bcftools view -h "${vcf}".tmp.tmp | grep -e "#tumor_sample" | cut -f2 -d"=") ;
        echo -e "--Normal sample ID is:\t""${normal}"". Filtering..." ;
        bcftools view --samples "${normal}" --types snps -Oz "${vcf}".tmp.tmp > "${out}" ;
        tabix -p vcf "${out}" ;
        echo -e "Done." ;
        genotypes=$(bcftools view "${out}" | grep -e "^#" -v | cut -f10 | cut -f1 -d":" | sort | uniq) ;
        echo -e "--genotypes after filtering include:\t""${genotypes}" ;
        echo -e "--total variants=\t""$(bcftools view "${out}" | grep -e "^#" -v | wc -l)" ;
        echo -e "\n" ;

      rm "${vcf}".tmp "${vcf}".tmp.tmp ;
    done ;
