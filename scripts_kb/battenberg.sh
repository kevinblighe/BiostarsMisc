
  cd /shared/aws-backup-restore_2022-02-11T14-34-30-948Z/CLL ;

  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/ --recursive |\
    grep -e "battenberg" | rev | cut -f1 -d" " | rev |\
      while read file ; do
        outdir=$(echo -e "${file}" | cut -f6 -d"/" | cut -f1 -d"_") ;
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/Battenberg/"${outdir}""/" --dryrun ;
      done ;

  aws s3 ls s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/ --recursive |\
    grep -e "battenberg" | rev | cut -f1 -d" " | rev |\
      while read file ; do
        outdir=$(echo -e "${file}" | cut -f6 -d"/" | cut -f1 -d"_") ;
        aws s3 cp s3://celgene-rnd-riku-lgp/"${file}" WGS/Battenberg/"${outdir}""/" ;
      done ;

