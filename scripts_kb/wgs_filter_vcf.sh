
  #!/bin/bash ;
  export PATH=/home/blighek1/bcftools-1.15/bin:$PATH ;
  export PATH=/home/blighek1/htslib-1.15.1/bin:$PATH ;
  export PATH=/home/blighek1/plink_linux_x86_64_20220402/:$PATH ;
  cd /shared/CLL/ ;
  mkdir -p WGS/FilteredData ;



  # check alignment stats
    find WGS/ProcessedData/wgsmetrics/ -name "*.qcstats" | sort | head -1 | while read qc ;
    do
      grep -e "## METRICS CLASS" -A 2 "${qc}" | sed '1d' > qc_kb/WGS_picard.qcstats.tmp ;
    done ;
    find WGS/ProcessedData/wgsmetrics/ -name "*.qcstats" | sort | sed '1d' | while read qc ;
    do
      grep -e "## METRICS CLASS" -A 2 "${qc}" | sed '1,2d' >> qc_kb/WGS_picard.qcstats.tmp ;
    done ;
    paste \
      <(cat \
        <(echo "SAMPLE") \
        <(find WGS/ProcessedData/wgsmetrics/ -name "*.qcstats" | sort | cut -f4 -d"/" | cut -f1 -d".")) \
      <(cat qc_kb/WGS_picard.qcstats.tmp) > qc_kb/WGS_picard.qcstats ;
    rm qc_kb/WGS_picard.qcstats.tmp ;



  # initial run to check position of tumor and normal in VCF
    find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v | sort | while read vcf ;
    do
      bcftools view -h "${vcf}" | grep -e "#tumor_sample" ;
      bcftools view -h "${vcf}" | grep -e "#normal_sample" ;
      bcftools view -h "${vcf}" | tail -1 ;
    done ;
  # In those samples with 'CD19' in the name, the tumor sample always comes first; in all others, its second.



  # run the main filtration process
    find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v | sort | while read vcf ;
    do
      echo -e "--Input VCF is:\t""${vcf}" ;

      out=$(echo "${vcf}" | sed 's/\.vcf\.gz/\.filt\.vcf\.gz/g' | sed 's/ProcessedData\/[0-9]*\/dbnsfp_annotation/FilteredData/g') ;
      echo -e "--Output VCF is:\t""${out}" ;

      # Filter to main contigs only
        # bcftools filter -r "chr1...
      # Filter to variants with PASS in the FILTER column, and left-align indels / MNPs
        # --apply-filters 'PASS'
        # bcftools norm -f refgenome/hg38.fa
      # Remove any variants that are now wildtype, ./. or 0/0
        # --min-ac 1
      # Remove variants based on read depth
      # Remove variants based on AF in samples
        # TUMOR_depth >= 10
        # NORMAL_depth >= 10
        # TUMOR_AF > 0.05
        # NORMAL_AF < 0.05
        # TUMOR_alt_depth >= 3
        # SUM(FORMAT/AD[0:*])>=10 && SUM(FORMAT/AD[1:*])>=10 && FORMAT/AF[0:*]>0.05 && FORMAT/AF[1:*]<0.05 && FORMAT/AD[0:1]>=3
      # Remove common variants
        # as.numeric(dbNSFP_ExAC_AF) < 0.01 | is.na(as.numeric(dbNSFP_ExAC_AF))
        # as.numeric(dbNSFP_gnomAD_genomes_AF) < 0.01 | is.na(as.numeric(dbNSFP_gnomAD_genomes_AF))
        # as.numeric(dbNSFP_gnomAD_exomes_AF) < 0.01 | is.na(as.numeric(dbNSFP_gnomAD_exomes_AF))
        # INFO/dbNSFP_ExAC_AF<0.01 || INFO/dbNSFP_gnomAD_genomes_AF<0.01 || INFO/dbNSFP_gnomAD_exomes_AF<0.01
        if [ $(echo "${vcf}" | grep -c "CD19") -eq 1 ]
        then
          bcftools filter -r "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chrX,chrY" "${vcf}" |\
            bcftools view --apply-filters 'PASS' --min-ac 1 | bcftools norm -f refgenome/hg38.fa |\
              bcftools view --include 'SUM(FORMAT/AD[0:*])>=10 && SUM(FORMAT/AD[1:*])>=10 && FORMAT/AF[0:*]>0.05 && FORMAT/AF[1:*]<0.05 && FORMAT/AD[0:1]>=3' |\
                bcftools filter --include 'INFO/dbNSFP_ExAC_AF<0.01 || INFO/dbNSFP_ExAC_AF="."' |\
                  bcftools filter --include 'INFO/dbNSFP_gnomAD_genomes_AF<0.01 || INFO/dbNSFP_gnomAD_genomes_AF="."' |\
                    bcftools filter --include 'INFO/dbNSFP_gnomAD_exomes_AF<0.01 || INFO/dbNSFP_gnomAD_exomes_AF="."' \
                      > "${vcf}".tmp ;
        else
          bcftools filter -r "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chrX,chrY" "${vcf}" |\
            bcftools view --apply-filters 'PASS' --min-ac 1 | bcftools norm -f refgenome/hg38.fa |\
              bcftools view --include 'SUM(FORMAT/AD[0:*])>=10 && SUM(FORMAT/AD[1:*])>=10 && FORMAT/AF[0:*]<0.05 && FORMAT/AF[1:*]>0.05 && FORMAT/AD[1:1]>=3' |\
                bcftools filter --include 'INFO/dbNSFP_ExAC_AF<0.01 || INFO/dbNSFP_ExAC_AF="."' |\
                  bcftools filter --include 'INFO/dbNSFP_gnomAD_genomes_AF<0.01 || INFO/dbNSFP_gnomAD_genomes_AF="."' |\
                    bcftools filter --include 'INFO/dbNSFP_gnomAD_exomes_AF<0.01 || INFO/dbNSFP_gnomAD_exomes_AF="."' \
                      > "${vcf}".tmp ;
        fi

      # filter out unexpectedly high read depth variants
      #   NORMAL_depth <= 3*mean(NORMAL_depth,na.rm=T)
      #   TUMOR_depth <= 3*mean(TUMOR_depth, na.rm=T),
        tumor_avg_3x=$(bcftools query -f '[%SAMPLE=%AD\t]\n' "${vcf}".tmp | cut -f1 | cut -f2 -d"=" | awk -F "," 'BEGIN {a=0} {a+=($1 + $2)} END {print 3*(a/NR)}') ;
        normal_avg_3x=$(bcftools query -f '[%SAMPLE=%AD\t]\n' "${vcf}".tmp | cut -f2 | cut -f2 -d"=" | awk -F "," 'BEGIN {a=0} {a+=($1 + $2)} END {print 3*(a/NR)}') ;
        bcftools view --include "SUM(FORMAT/AD[0:*])<=$tumor_avg_3x && SUM(FORMAT/AD[1:*])<=$normal_avg_3x" "${vcf}".tmp \
          > "${vcf}".tmp.tmp ;

      # include tumor sample only, and only include SNPs, InDels, MNPs
        tumor=$(bcftools view -h "${vcf}".tmp.tmp | grep -e "#tumor_sample" | cut -f2 -d"=") ;
        echo -e "--Tumor sample ID is:\t""${tumor}"". Filtering..." ;
        bcftools view --samples "${tumor}" --types snps,indels,mnps -Oz "${vcf}".tmp.tmp > "${out}" ;
        tabix -p vcf "${out}" ;
        echo -e "Done." ;
        genotypes=$(bcftools view "${out}" | grep -e "^#" -v | cut -f10 | cut -f1 -d":" | sort | uniq) ;
        echo -e "--genotypes after filtering include:\t""${genotypes}" ;
        echo -e "--total variants=\t""$(bcftools view "${out}" | grep -e "^#" -v | wc -l)" ;
        echo -e "\n" ;

      rm "${vcf}".tmp "${vcf}".tmp.tmp ;
    done ;



  # generate stats for >40X
    echo -e "Sample\tSNVs\tpc>40x\tProtein Coding Region\tMissense Variants\tHigh Impact\tMedium Impact\tHigh|Moderate Impact Genes\tTop 3 Most Mutated Genes (any)\tTop 3 Most Mutated Genes (missense)" > qc_kb/WGS_Filtered_Variant_Stats.tsv ;
    find WGS/FilteredData/ -name "*.vcf.gz" | sort | while read vcf ;
    do
      paste -d "\t" \
        <(echo "${vcf}") \
        <(bcftools stats -1 "${vcf}" | grep -e "^SN" |\
          grep -e "number of SNPs" | cut -f4 | awk '{printf $0"\t"} END {printf "\n"}' | sed 's/\t$//g') \
        <(bcftools stats -1 "${vcf}" -d 0,39,39 | tail -1 | cut -f7) \
        <(bcftools view "${vcf}" | grep -e "^#" -v | grep -e "|protein_coding|" | wc -l) \
        <(bcftools view "${vcf}" | grep -e "^#" -v | grep -e "|missense_variant|" | wc -l) \
        <(bcftools view "${vcf}" | grep -e "^#" -v | grep -e "|HIGH|" | wc -l) \
        <(bcftools view "${vcf}" | grep -e "^#" -v | grep -e "|MODERATE|" | wc -l) \
        <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\n' "${vcf}" | grep -e "|HIGH|" -e "|MODERATE|" |\
          cut -f4 -d"|" | sort | uniq | awk '{printf $0","} END {printf "\n"}' | sed 's/,$//g') \
        <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\n' "${vcf}" | cut -f4 -d"|" | sed '1d' |\
          sort | uniq -c | sed 's/^ \+//g' | sort -k1n | tail -3 | awk -F " " '{print $2"("$1")"}' |\
          awk '{printf $0","} END {printf "\n"}' | sed 's/,$//g') \
        <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\n' "${vcf}" | cut -f2,4 -d"|" | grep -e "missense" |\
          cut -f2 -d"|" | sed '1d' | sort | uniq -c | sed 's/^ \+//g' | sort -k1n | tail -3 | awk -F " " '{print $2"("$1")"}' |\
          awk '{printf $0","} END {printf "\n"}' | sed 's/,$//g') ;
    done | sed 's/WGS\/FilteredData\///g' | sed 's/\.vcf\.gz//g' | sed 's/ \+/\t/g' >> qc_kb/WGS_Filtered_Variant_Stats.tsv ;

    # Notes:
      # other categories for protein_coding are:
        # antisense
        # IG_D_gene
        # lincRNA
        # miRNA
        # misc_RNA
        # nonsense_mediated_decay
        # polymorphic_pseudogene
        # processed_pseudogene
        # processed_transcript
        # retained_intron
        # rRNA
        # scaRNA
        # sense_intronic
        # sense_overlapping
        # snoRNA
        # snRNA
        # TEC
        # transcribed_processed_pseudogene
        # transcribed_unprocessed_pseudogene
        # unitary_pseudogene
        # unprocessed_pseudogene
      # other categories for missense_variant consequences are:
        # 3_prime_UTR_variant
        # 5_prime_UTR_variant
        # downstream_gene_variant
        # intergenic_region
        # intron_variant
        # missense_variant&splice_region_variant
        # non_coding_transcript_exon_variant
        # sequence_feature
        # splice_acceptor_variant&intron_variant
        # splice_region_variant
        # splice_region_variant&intron_variant
        # splice_region_variant&non_coding_transcript_exon_variant
        # synonymous_variant
        # TF_binding_site_variant
        # upstream_gene_variant



  # generate stats for read depth
    find WGS/FilteredData/ -name "*.vcf.gz" | sort | while read vcf ;
    do
      out=$(echo "${vcf}" | sed 's/\.vcf\.gz/\.DP\.txt/g' | sed 's/WGS\/FilteredData/qc_kb/g') ;
      echo -e "Variant\tRead Depth\tMedian Base Quality" > "${out}" ;
      bcftools query -f '%CHROM:%POS,%REF>%ALT\t%DP[\t%MBQ]\n' "${vcf}" >> "${out}" ;
    done ;



  # merge all files and find most frequent variants, and frequently mutated genes (FILTERED VCF)
    mkdir -p .tmp ;
    find WGS/FilteredData/ -name "*.vcf.gz" | sort | while read vcf
    do
      out=$(echo "${vcf}" | sed 's/WGS\/FilteredData/\.tmp/g') ;
      bcftools annotate -x INFO/ECNT,INFO/IN_PON,INFO/NLOD,INFO/N_ART_LOD,INFO/POP_AF,INFO/P_GERMLINE,INFO/RPA,INFO/RU,INFO/STR,INFO/TLOD,INFO/PrimarySite,INFO/SiteSubtype_1,INFO/SiteSubtype_2,INFO/SiteSubtype_3,INFO/PrimaryHistology,INFO/HistologySubtype_1,INFO/HistologySubtype_2,INFO/HistologySubtype_3,INFO/Resistance,INFO/Status,INFO/SampleType,INFO/TumorOrigin,INFO/Occurences,INFO/LOF,INFO/NMD,INFO/CDA,INFO/OTH,INFO/S3D,INFO/WTD,INFO/dbSNPBuildID,INFO/SLO,INFO/NSF,INFO/R3,INFO/R5,INFO/NSN,INFO/NSM,INFO/G5A,INFO/COMMON,INFO/RS,INFO/RV,INFO/TPA,INFO/CFL,INFO/GNO,INFO/VLD,INFO/ASP,INFO/ASS,INFO/REF,INFO/U3,INFO/U5,INFO/WGT,INFO/MTP,INFO/LSD,INFO/NOC,INFO/DSS,INFO/SYN,INFO/KGPhase3,INFO/CAF,INFO/VC,INFO/MUT,INFO/KGPhase1,INFO/NOV,INFO/VP,INFO/SAO,INFO/GENEINFO,INFO/INT,INFO/G5,INFO/OM,INFO/PMC,INFO/SSR,INFO/RSPOS,INFO/HD,INFO/PM,INFO/CDS,INFO/AA,INFO/GENE,INFO/CNT,INFO/STRAND,INFO/CLNSIG,INFO/VARTYPE,INFO/SNP,INFO/MNP,INFO/INS,INFO/DEL,INFO/MIXED,INFO/HOM,INFO/HET,INFO/dbNSFP_ExAC_nonTCGA_SAS_AF,INFO/dbNSFP_ExAC_SAS_AF,INFO/dbNSFP_GERP___RS,INFO/dbNSFP_phastCons100way_vertebrate_rankscore,INFO/dbNSFP_ExAC_SAS_AC,INFO/dbNSFP_ExAC_nonpsych_AMR_AF,INFO/dbNSFP_ExAC_nonpsych_AMR_AC,INFO/dbNSFP_MutPred_protID,INFO/dbNSFP_gnomAD_genomes_AMR_AN,INFO/dbNSFP_MetaSVM_pred,INFO/dbNSFP_1000Gp3_EAS_AC,INFO/dbNSFP_FATHMM_pred,INFO/dbNSFP_ExAC_nonTCGA_SAS_AC,INFO/dbNSFP_MutPred_Top5features,INFO/dbNSFP_M_CAP_score,INFO/dbNSFP_MutationAssessor_score_rankscore,INFO/dbNSFP_ExAC_nonpsych_Adj_AF,INFO/dbNSFP_ExAC_nonpsych_NFE_AF,INFO/dbNSFP_gnomAD_genomes_AMR_AC,INFO/dbNSFP_gnomAD_exomes_SAS_AC,INFO/dbNSFP_gnomAD_genomes_NFE_AN,INFO/dbNSFP_ExAC_nonpsych_NFE_AC,INFO/dbNSFP_SIFT_score,INFO/dbNSFP_ExAC_nonpsych_Adj_AC,INFO/dbNSFP_phyloP100way_vertebrate_rankscore,INFO/dbNSFP_1000Gp3_EAS_AF,INFO/dbNSFP_gnomAD_exomes_SAS_AF,INFO/dbNSFP_ExAC_FIN_AC,INFO/dbNSFP_gnomAD_genomes_NFE_AF,INFO/dbNSFP_gnomAD_genomes_AMR_AF,INFO/dbNSFP_ExAC_FIN_AF,INFO/dbNSFP_gnomAD_exomes_SAS_AN,INFO/dbNSFP_ESP6500_EA_AC,INFO/dbNSFP_gnomAD_genomes_NFE_AC,INFO/dbNSFP_MetaSVM_rankscore,INFO/dbNSFP_ExAC_nonTCGA_EAS_AC,INFO/dbNSFP_FATHMM_score,INFO/dbNSFP_gnomAD_genomes_FIN_AC,INFO/dbNSFP_gnomAD_genomes_FIN_AF,INFO/dbNSFP_ExAC_nonTCGA_FIN_AC,INFO/dbNSFP_ExAC_nonpsych_FIN_AC,INFO/dbNSFP_LRT_score,INFO/dbNSFP_1000Gp3_SAS_AC,INFO/dbNSFP_ExAC_EAS_AC,INFO/dbNSFP_ExAC_nonTCGA_FIN_AF,INFO/dbNSFP_1000Gp3_SAS_AF,INFO/dbNSFP_ExAC_EAS_AF,INFO/dbNSFP_ESP6500_EA_AF,INFO/dbNSFP_gnomAD_genomes_FIN_AN,INFO/dbNSFP_MutPred_rankscore,INFO/dbNSFP_ExAC_nonTCGA_Adj_AF,INFO/dbNSFP_ExAC_nonpsych_FIN_AF,INFO/dbNSFP_GTEx_V6p_tissue,INFO/dbNSFP_ExAC_Adj_AC,INFO/dbNSFP_ExAC_Adj_AF,INFO/dbNSFP_ExAC_nonpsych_AFR_AF,INFO/dbNSFP_gnomAD_genomes_AC,INFO/dbNSFP_ExAC_nonpsych_AFR_AC,INFO/dbNSFP_gnomAD_genomes_AF,INFO/dbNSFP_SIFT_converted_rankscore,INFO/dbNSFP_gnomAD_genomes_OTH_AN,INFO/dbNSFP_ExAC_nonTCGA_Adj_AC,INFO/dbNSFP_gnomAD_genomes_OTH_AC,INFO/dbNSFP_gnomAD_genomes_OTH_AF,INFO/dbNSFP_Eigen_PC_raw_rankscore,INFO/dbNSFP_gnomAD_exomes_NFE_AC,INFO/dbNSFP_MutationTaster_converted_rankscore,INFO/dbNSFP_phyloP100way_vertebrate,INFO/dbNSFP_LRT_pred,INFO/dbNSFP_PROVEAN_pred,INFO/dbNSFP_phastCons100way_vertebrate,INFO/dbNSFP_MetaSVM_score,INFO/dbNSFP_1000Gp3_EUR_AC,INFO/dbNSFP_gnomAD_exomes_ASJ_AC,INFO/dbNSFP_Eigen_PC_phred,INFO/dbNSFP_1000Gp3_EUR_AF,INFO/dbNSFP_gnomAD_exomes_ASJ_AF,INFO/dbNSFP_gnomAD_genomes_AN,INFO/dbNSFP_ESP6500_AA_AF,INFO/dbNSFP_M_CAP_pred,INFO/dbNSFP_gnomAD_exomes_ASJ_AN,INFO/dbNSFP_ESP6500_AA_AC,INFO/dbNSFP_PROVEAN_score,INFO/dbNSFP_ExAC_nonTCGA_NFE_AF,INFO/dbNSFP_ALSPAC_AC,INFO/dbNSFP_SiPhy_29way_logOdds,INFO/dbNSFP_ExAC_NFE_AF,INFO/dbNSFP_gnomAD_genomes_EAS_AC,INFO/dbNSFP_GERP___NR,INFO/dbNSFP_gnomAD_genomes_EAS_AF,INFO/dbNSFP_FATHMM_converted_rankscore,INFO/dbNSFP_1000Gp3_AMR_AF,INFO/dbNSFP_ALSPAC_AF,INFO/dbNSFP_ExAC_nonTCGA_AMR_AC,INFO/dbNSFP_ExAC_nonTCGA_NFE_AC,INFO/dbNSFP_PROVEAN_converted_rankscore,INFO/dbNSFP_gnomAD_genomes_EAS_AN,INFO/dbNSFP_Eigen_PC_raw,INFO/dbNSFP_1000Gp3_AMR_AC,INFO/dbNSFP_gnomAD_exomes_NFE_AN,INFO/dbNSFP_ExAC_nonTCGA_AMR_AF,INFO/dbNSFP_ExAC_nonpsych_EAS_AC,INFO/dbNSFP_gnomAD_exomes_AMR_AF,INFO/dbNSFP_ExAC_nonpsych_EAS_AF,INFO/dbNSFP_gnomAD_exomes_AMR_AC,INFO/dbNSFP_SiPhy_29way_logOdds_rankscore,INFO/dbNSFP_gnomAD_exomes_NFE_AF,INFO/dbNSFP_Reliability_index,INFO/dbNSFP_ExAC_AFR_AF,INFO/dbNSFP_gnomAD_exomes_AMR_AN,INFO/dbNSFP_MutPred_AAchange,INFO/dbNSFP_gnomAD_genomes_ASJ_AF,INFO/dbNSFP_ExAC_AFR_AC,INFO/dbNSFP_gnomAD_exomes_AF,INFO/dbNSFP_gnomAD_genomes_ASJ_AC,INFO/dbNSFP_GERP___RS_rankscore,INFO/dbNSFP_gnomAD_exomes_AC,INFO/dbNSFP_MetaLR_score,INFO/dbNSFP_ExAC_nonpsych_AF,INFO/dbNSFP_gnomAD_genomes_ASJ_AN,INFO/dbNSFP_ExAC_nonpsych_AC,INFO/dbNSFP_gnomAD_exomes_AFR_AF,INFO/dbNSFP_1000Gp3_AFR_AC,INFO/dbNSFP_MetaLR_rankscore,INFO/dbNSFP_ExAC_AMR_AF,INFO/dbNSFP_1000Gp3_AFR_AF,INFO/dbNSFP_MutationTaster_pred,INFO/dbNSFP_gnomAD_exomes_AFR_AN,INFO/dbNSFP_MetaLR_pred,INFO/dbNSFP_ExAC_AMR_AC,INFO/dbNSFP_ExAC_NFE_AC,INFO/dbNSFP_gnomAD_exomes_AN,INFO/dbNSFP_gnomAD_exomes_FIN_AF,INFO/dbNSFP_gnomAD_genomes_AFR_AF,INFO/dbNSFP_gnomAD_exomes_OTH_AN,INFO/dbNSFP_MutationAssessor_score,INFO/dbNSFP_TWINSUK_AF,INFO/dbNSFP_TWINSUK_AC,INFO/dbNSFP_gnomAD_exomes_FIN_AN,INFO/dbNSFP_gnomAD_genomes_AFR_AN,INFO/dbNSFP_MutPred_score,INFO/dbNSFP_gnomAD_exomes_AFR_AC,INFO/dbNSFP_ExAC_nonTCGA_EAS_AF,INFO/dbNSFP_gnomAD_exomes_OTH_AC,INFO/dbNSFP_GTEx_V6p_gene,INFO/dbNSFP_GenoCanyon_score_rankscore,INFO/dbNSFP_gnomAD_genomes_AFR_AC,INFO/dbNSFP_gnomAD_exomes_OTH_AF,INFO/dbNSFP_ExAC_nonpsych_SAS_AF,INFO/dbNSFP_gnomAD_exomes_EAS_AN,INFO/dbNSFP_GenoCanyon_score,INFO/dbNSFP_ExAC_nonTCGA_AFR_AC,INFO/dbNSFP_ExAC_nonpsych_SAS_AC,INFO/dbNSFP_1000Gp3_AC,INFO/dbNSFP_ExAC_nonTCGA_AF,INFO/dbNSFP_1000Gp3_AF,INFO/dbNSFP_ExAC_AF,INFO/dbNSFP_ExAC_nonTCGA_AC,INFO/dbNSFP_ExAC_AC,INFO/dbNSFP_ExAC_nonTCGA_AFR_AF,INFO/dbNSFP_gnomAD_exomes_EAS_AF,INFO/dbNSFP_gnomAD_exomes_EAS_AC,INFO/dbNSFP_M_CAP_rankscore,INFO/dbNSFP_MutationTaster_score,INFO/dbNSFP_MutationAssessor_pred,INFO/dbNSFP_SIFT_pred,INFO/dbNSFP_gnomAD_exomes_FIN_AC "${vcf}" -Oz > "${out}" ;
      tabix -p vcf "${out}" ;
    done ;
    find .tmp/ -name "*.vcf.gz" | sort | awk 'BEGIN {printf "bcftools merge --info-rules DP:AVG "} {printf $0" "} END {printf "-Oz > WGS/Merge.Filtered.vcf.gz\n"}' ;
    rm -R .tmp ;

    paste <(bcftools view WGS/Merge.Filtered.vcf.gz |\
      awk -F"\t" 'BEGIN {print "Variant\tID"} \
        !/^#/ {print $1":"$2","$4">"$5"\t"$3}') \
      \
      <(cat <(echo "Gene") <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | cut -f4 -d"|")) \
      \
      <(cat\
        <(echo "Impact")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz |\
          awk '{if ($0 ~/\|HIGH\|/ || $0 ~ /\|MODERATE\|/) {print "HIGH|MODERATE"} else {print "OTHER"}}')) \
      \
      <(cat\
        <(echo "Function")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[2]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHet"} {print gsub(/0\/1|1\/0/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHomAlt"} {print gsub(/1\/1/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHomRef"} {print gsub(/0\/0/, "")}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      | sed 's/,\t/\t/g' | sed 's/,$//g' > qc_kb/WGS_Merged.Filtered_Variant_Stats.tsv ;

    paste <(bcftools view WGS/Merge.Filtered.vcf.gz |\
      awk -F"\t" 'BEGIN {print "Variant\tID"} \
        !/^#/ {print $1":"$2","$4">"$5"\t"$3}') \
      \
      <(cat <(echo "Gene") <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | cut -f4 -d"|")) \
      \
      <(bcftools query -f '[%AF,]\n' WGS/Merge.Filtered.vcf.gz | sed 's/\.,//g' | sed 's/,$//g' | awk 'BEGIN {print "AF"} {print $0}') \
      \
      <(cat\
        <(echo "Impact")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz |\
          awk '{if ($0 ~/\|HIGH\|/ || $0 ~ /\|MODERATE\|/) {print "HIGH|MODERATE"} else {print "OTHER"}}')) \
      \
      <(cat\
        <(echo "Function")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[2]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHet"} {print gsub(/0\/1|1\/0/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHomAlt"} {print gsub(/1\/1/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHomRef"} {print gsub(/0\/0/, "")}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      | sed 's/,\t/\t/g' | sed 's/,$//g' > qc_kb/WGS_Merged.Filtered_Variant_Stats.VAF.tsv ;



  # merge all files and find most frequent variants, and frequently mutated genes (UNFILTERED VCF)
    mkdir -p .tmp ;
    find WGS/ProcessedData/ -name "*.vcf.gz" | grep -e "dbnsfp_annotation" | grep -e "tumor_only" -v | sort | while read vcf ;
    do
      out=$(echo "${vcf}" | sed 's/WGS\/ProcessedData\/[0-9]*\/dbnsfp_annotation/.tmp/g') ;
      bcftools annotate -x INFO/ECNT,INFO/IN_PON,INFO/NLOD,INFO/N_ART_LOD,INFO/POP_AF,INFO/P_GERMLINE,INFO/RPA,INFO/RU,INFO/STR,INFO/TLOD,INFO/PrimarySite,INFO/SiteSubtype_1,INFO/SiteSubtype_2,INFO/SiteSubtype_3,INFO/PrimaryHistology,INFO/HistologySubtype_1,INFO/HistologySubtype_2,INFO/HistologySubtype_3,INFO/Resistance,INFO/Status,INFO/SampleType,INFO/TumorOrigin,INFO/Occurences,INFO/LOF,INFO/NMD,INFO/CDA,INFO/OTH,INFO/S3D,INFO/WTD,INFO/dbSNPBuildID,INFO/SLO,INFO/NSF,INFO/R3,INFO/R5,INFO/NSN,INFO/NSM,INFO/G5A,INFO/COMMON,INFO/RS,INFO/RV,INFO/TPA,INFO/CFL,INFO/GNO,INFO/VLD,INFO/ASP,INFO/ASS,INFO/REF,INFO/U3,INFO/U5,INFO/WGT,INFO/MTP,INFO/LSD,INFO/NOC,INFO/DSS,INFO/SYN,INFO/KGPhase3,INFO/CAF,INFO/VC,INFO/MUT,INFO/KGPhase1,INFO/NOV,INFO/VP,INFO/SAO,INFO/GENEINFO,INFO/INT,INFO/G5,INFO/OM,INFO/PMC,INFO/SSR,INFO/RSPOS,INFO/HD,INFO/PM,INFO/CDS,INFO/AA,INFO/GENE,INFO/CNT,INFO/STRAND,INFO/CLNSIG,INFO/VARTYPE,INFO/SNP,INFO/MNP,INFO/INS,INFO/DEL,INFO/MIXED,INFO/HOM,INFO/HET,INFO/dbNSFP_ExAC_nonTCGA_SAS_AF,INFO/dbNSFP_ExAC_SAS_AF,INFO/dbNSFP_GERP___RS,INFO/dbNSFP_phastCons100way_vertebrate_rankscore,INFO/dbNSFP_ExAC_SAS_AC,INFO/dbNSFP_ExAC_nonpsych_AMR_AF,INFO/dbNSFP_ExAC_nonpsych_AMR_AC,INFO/dbNSFP_MutPred_protID,INFO/dbNSFP_gnomAD_genomes_AMR_AN,INFO/dbNSFP_MetaSVM_pred,INFO/dbNSFP_1000Gp3_EAS_AC,INFO/dbNSFP_FATHMM_pred,INFO/dbNSFP_ExAC_nonTCGA_SAS_AC,INFO/dbNSFP_MutPred_Top5features,INFO/dbNSFP_M_CAP_score,INFO/dbNSFP_MutationAssessor_score_rankscore,INFO/dbNSFP_ExAC_nonpsych_Adj_AF,INFO/dbNSFP_ExAC_nonpsych_NFE_AF,INFO/dbNSFP_gnomAD_genomes_AMR_AC,INFO/dbNSFP_gnomAD_exomes_SAS_AC,INFO/dbNSFP_gnomAD_genomes_NFE_AN,INFO/dbNSFP_ExAC_nonpsych_NFE_AC,INFO/dbNSFP_SIFT_score,INFO/dbNSFP_ExAC_nonpsych_Adj_AC,INFO/dbNSFP_phyloP100way_vertebrate_rankscore,INFO/dbNSFP_1000Gp3_EAS_AF,INFO/dbNSFP_gnomAD_exomes_SAS_AF,INFO/dbNSFP_ExAC_FIN_AC,INFO/dbNSFP_gnomAD_genomes_NFE_AF,INFO/dbNSFP_gnomAD_genomes_AMR_AF,INFO/dbNSFP_ExAC_FIN_AF,INFO/dbNSFP_gnomAD_exomes_SAS_AN,INFO/dbNSFP_ESP6500_EA_AC,INFO/dbNSFP_gnomAD_genomes_NFE_AC,INFO/dbNSFP_MetaSVM_rankscore,INFO/dbNSFP_ExAC_nonTCGA_EAS_AC,INFO/dbNSFP_FATHMM_score,INFO/dbNSFP_gnomAD_genomes_FIN_AC,INFO/dbNSFP_gnomAD_genomes_FIN_AF,INFO/dbNSFP_ExAC_nonTCGA_FIN_AC,INFO/dbNSFP_ExAC_nonpsych_FIN_AC,INFO/dbNSFP_LRT_score,INFO/dbNSFP_1000Gp3_SAS_AC,INFO/dbNSFP_ExAC_EAS_AC,INFO/dbNSFP_ExAC_nonTCGA_FIN_AF,INFO/dbNSFP_1000Gp3_SAS_AF,INFO/dbNSFP_ExAC_EAS_AF,INFO/dbNSFP_ESP6500_EA_AF,INFO/dbNSFP_gnomAD_genomes_FIN_AN,INFO/dbNSFP_MutPred_rankscore,INFO/dbNSFP_ExAC_nonTCGA_Adj_AF,INFO/dbNSFP_ExAC_nonpsych_FIN_AF,INFO/dbNSFP_GTEx_V6p_tissue,INFO/dbNSFP_ExAC_Adj_AC,INFO/dbNSFP_ExAC_Adj_AF,INFO/dbNSFP_ExAC_nonpsych_AFR_AF,INFO/dbNSFP_gnomAD_genomes_AC,INFO/dbNSFP_ExAC_nonpsych_AFR_AC,INFO/dbNSFP_gnomAD_genomes_AF,INFO/dbNSFP_SIFT_converted_rankscore,INFO/dbNSFP_gnomAD_genomes_OTH_AN,INFO/dbNSFP_ExAC_nonTCGA_Adj_AC,INFO/dbNSFP_gnomAD_genomes_OTH_AC,INFO/dbNSFP_gnomAD_genomes_OTH_AF,INFO/dbNSFP_Eigen_PC_raw_rankscore,INFO/dbNSFP_gnomAD_exomes_NFE_AC,INFO/dbNSFP_MutationTaster_converted_rankscore,INFO/dbNSFP_phyloP100way_vertebrate,INFO/dbNSFP_LRT_pred,INFO/dbNSFP_PROVEAN_pred,INFO/dbNSFP_phastCons100way_vertebrate,INFO/dbNSFP_MetaSVM_score,INFO/dbNSFP_1000Gp3_EUR_AC,INFO/dbNSFP_gnomAD_exomes_ASJ_AC,INFO/dbNSFP_Eigen_PC_phred,INFO/dbNSFP_1000Gp3_EUR_AF,INFO/dbNSFP_gnomAD_exomes_ASJ_AF,INFO/dbNSFP_gnomAD_genomes_AN,INFO/dbNSFP_ESP6500_AA_AF,INFO/dbNSFP_M_CAP_pred,INFO/dbNSFP_gnomAD_exomes_ASJ_AN,INFO/dbNSFP_ESP6500_AA_AC,INFO/dbNSFP_PROVEAN_score,INFO/dbNSFP_ExAC_nonTCGA_NFE_AF,INFO/dbNSFP_ALSPAC_AC,INFO/dbNSFP_SiPhy_29way_logOdds,INFO/dbNSFP_ExAC_NFE_AF,INFO/dbNSFP_gnomAD_genomes_EAS_AC,INFO/dbNSFP_GERP___NR,INFO/dbNSFP_gnomAD_genomes_EAS_AF,INFO/dbNSFP_FATHMM_converted_rankscore,INFO/dbNSFP_1000Gp3_AMR_AF,INFO/dbNSFP_ALSPAC_AF,INFO/dbNSFP_ExAC_nonTCGA_AMR_AC,INFO/dbNSFP_ExAC_nonTCGA_NFE_AC,INFO/dbNSFP_PROVEAN_converted_rankscore,INFO/dbNSFP_gnomAD_genomes_EAS_AN,INFO/dbNSFP_Eigen_PC_raw,INFO/dbNSFP_1000Gp3_AMR_AC,INFO/dbNSFP_gnomAD_exomes_NFE_AN,INFO/dbNSFP_ExAC_nonTCGA_AMR_AF,INFO/dbNSFP_ExAC_nonpsych_EAS_AC,INFO/dbNSFP_gnomAD_exomes_AMR_AF,INFO/dbNSFP_ExAC_nonpsych_EAS_AF,INFO/dbNSFP_gnomAD_exomes_AMR_AC,INFO/dbNSFP_SiPhy_29way_logOdds_rankscore,INFO/dbNSFP_gnomAD_exomes_NFE_AF,INFO/dbNSFP_Reliability_index,INFO/dbNSFP_ExAC_AFR_AF,INFO/dbNSFP_gnomAD_exomes_AMR_AN,INFO/dbNSFP_MutPred_AAchange,INFO/dbNSFP_gnomAD_genomes_ASJ_AF,INFO/dbNSFP_ExAC_AFR_AC,INFO/dbNSFP_gnomAD_exomes_AF,INFO/dbNSFP_gnomAD_genomes_ASJ_AC,INFO/dbNSFP_GERP___RS_rankscore,INFO/dbNSFP_gnomAD_exomes_AC,INFO/dbNSFP_MetaLR_score,INFO/dbNSFP_ExAC_nonpsych_AF,INFO/dbNSFP_gnomAD_genomes_ASJ_AN,INFO/dbNSFP_ExAC_nonpsych_AC,INFO/dbNSFP_gnomAD_exomes_AFR_AF,INFO/dbNSFP_1000Gp3_AFR_AC,INFO/dbNSFP_MetaLR_rankscore,INFO/dbNSFP_ExAC_AMR_AF,INFO/dbNSFP_1000Gp3_AFR_AF,INFO/dbNSFP_MutationTaster_pred,INFO/dbNSFP_gnomAD_exomes_AFR_AN,INFO/dbNSFP_MetaLR_pred,INFO/dbNSFP_ExAC_AMR_AC,INFO/dbNSFP_ExAC_NFE_AC,INFO/dbNSFP_gnomAD_exomes_AN,INFO/dbNSFP_gnomAD_exomes_FIN_AF,INFO/dbNSFP_gnomAD_genomes_AFR_AF,INFO/dbNSFP_gnomAD_exomes_OTH_AN,INFO/dbNSFP_MutationAssessor_score,INFO/dbNSFP_TWINSUK_AF,INFO/dbNSFP_TWINSUK_AC,INFO/dbNSFP_gnomAD_exomes_FIN_AN,INFO/dbNSFP_gnomAD_genomes_AFR_AN,INFO/dbNSFP_MutPred_score,INFO/dbNSFP_gnomAD_exomes_AFR_AC,INFO/dbNSFP_ExAC_nonTCGA_EAS_AF,INFO/dbNSFP_gnomAD_exomes_OTH_AC,INFO/dbNSFP_GTEx_V6p_gene,INFO/dbNSFP_GenoCanyon_score_rankscore,INFO/dbNSFP_gnomAD_genomes_AFR_AC,INFO/dbNSFP_gnomAD_exomes_OTH_AF,INFO/dbNSFP_ExAC_nonpsych_SAS_AF,INFO/dbNSFP_gnomAD_exomes_EAS_AN,INFO/dbNSFP_GenoCanyon_score,INFO/dbNSFP_ExAC_nonTCGA_AFR_AC,INFO/dbNSFP_ExAC_nonpsych_SAS_AC,INFO/dbNSFP_1000Gp3_AC,INFO/dbNSFP_ExAC_nonTCGA_AF,INFO/dbNSFP_1000Gp3_AF,INFO/dbNSFP_ExAC_AF,INFO/dbNSFP_ExAC_nonTCGA_AC,INFO/dbNSFP_ExAC_AC,INFO/dbNSFP_ExAC_nonTCGA_AFR_AF,INFO/dbNSFP_gnomAD_exomes_EAS_AF,INFO/dbNSFP_gnomAD_exomes_EAS_AC,INFO/dbNSFP_M_CAP_rankscore,INFO/dbNSFP_MutationTaster_score,INFO/dbNSFP_MutationAssessor_pred,INFO/dbNSFP_SIFT_pred,INFO/dbNSFP_gnomAD_exomes_FIN_AC "${vcf}" -Oz > "${out}" ;
      tabix -p vcf "${out}" ;
    done ;
    find .tmp/ -name "*.vcf.gz" | sort | awk 'BEGIN {printf "bcftools merge --info-rules DP:AVG "} {printf $0" "} END {printf "-Oz > WGS/Merge.UnFiltered.vcf.gz\n"}' ;
    rm -R .tmp ;

    paste <(bcftools view WGS/Merge.UnFiltered.vcf.gz |\
      awk -F"\t" 'BEGIN {print "Variant\tID"} \
        !/^#/ {print $1":"$2","$4">"$5"\t"$3}') \
      \
      <(cat <(echo "Gene") <(bcftools query -f '%ANN\n' WGS/Merge.UnFiltered.vcf.gz | cut -f4 -d"|")) \
      \
      <(bcftools query -f '[%AF,]\n' WGS/Merge.UnFiltered.vcf.gz | sed 's/\.,//g' | sed 's/,$//g' | awk 'BEGIN {print "AF"} {print $0}') \
      \
      <(cat\
        <(echo "Impact")\
        <(bcftools query -f '%ANN\n' WGS/Merge.UnFiltered.vcf.gz |\
          awk '{if ($0 ~/\|HIGH\|/ || $0 ~ /\|MODERATE\|/) {print "HIGH|MODERATE"} else {print "OTHER"}}')) \
      \
      <(cat\
        <(echo "Function")\
        <(bcftools query -f '%ANN\n' WGS/Merge.UnFiltered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[2]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.UnFiltered.vcf.gz |\
        awk 'BEGIN {print "nHet"} {print gsub(/0\/1|1\/0/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.UnFiltered.vcf.gz |\
        awk 'BEGIN {print "nHomAlt"} {print gsub(/1\/1/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.UnFiltered.vcf.gz |\
        awk 'BEGIN {print "nHomRef"} {print gsub(/0\/0/, "")}') \
      \
      <(bcftools view WGS/Merge.UnFiltered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.UnFiltered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.UnFiltered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      | sed 's/,\t/\t/g' | sed 's/,$//g' > qc_kb/WGS_Merged.UnFiltered_Variant_Stats.tsv ;



  # extract TRIM44 variants
    cat \
      <(head -1 qc_kb/WGS_Merged.UnFiltered_Variant_Stats.tsv) \
      <(awk -F "\t" '$3 ~ /TRIM44/' qc_kb/WGS_Merged.UnFiltered_Variant_Stats.tsv) \
      > qc_kb/WGS_Merged.UnFiltered_Variant_Stats.TRIM44.tsv ;



  # extract more information at transcript isoform level for TP53 and MED12
    paste <(bcftools view WGS/Merge.Filtered.vcf.gz |\
      awk -F"\t" 'BEGIN {print "Variant\tID"} \
        !/^#/ {print $1":"$2","$4">"$5"\t"$3}') \
      \
      <(cat <(echo "Gene") <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | cut -f4 -d"|")) \
      \
      <(bcftools query -f '[%AF,]\n' WGS/Merge.Filtered.vcf.gz | sed 's/\.,//g' | sed 's/,$//g' | awk 'BEGIN {print "AF"} {print $0}') \
      \
      <(cat\
        <(echo "Impact1")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz |\
          awk '{if ($0 ~/\|HIGH\|/ || $0 ~ /\|MODERATE\|/) {print "HIGH|MODERATE"} else {print "OTHER"}}')) \
      \
      <(cat\
        <(echo "Function")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[2]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(cat\
        <(echo "Impact2")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[3]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(cat\
        <(echo "Symbol")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[4]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(cat\
        <(echo "Ensembl Gene")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[5]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(cat\
        <(echo "Ensembl Transcript")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[7]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(cat\
        <(echo "Biotype")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[8]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(cat\
        <(echo "HGVS")\
        <(bcftools query -f '%ANN\n' WGS/Merge.Filtered.vcf.gz | sed 's/,/\t/g' |\
          awk -F "\t" '{for (i=1; i<=NF; i++) {split($(i), arr, /\|/); printf arr[10]";"; if (i==NF) {printf "\n"}}}' | sed 's/;$//g')) \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHet"} {print gsub(/0\/1|1\/0/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHomAlt"} {print gsub(/1\/1/, "")}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' WGS/Merge.Filtered.vcf.gz |\
        awk 'BEGIN {print "nHomRef"} {print gsub(/0\/0/, "")}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      <(bcftools view WGS/Merge.Filtered.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
      \
      | sed 's/,\t/\t/g' | sed 's/,$//g' |\
        awk -F "\t" '{if (NR==1) print}; {if ($3 == "TP53" || $3 == "MED12") print}' > qc_kb/WGS_Merged.Filtered_Variant_Stats.TP53.MED12.tsv ;
