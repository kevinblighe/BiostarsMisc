
  cd /shared/CLL/ ;

  aws s3 sync s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/ WGS/ProcessedData/ ;

  aws s3 sync s3://celgene-rnd-riku-lgp/DA0000882/WGS/ProcessedData/wgsmetrics/ WGS/ProcessedData/wgsmetrics/ ;

  aws s3 sync s3://celgene-rnd-riku-lgp/DA0000882/WGS/AnalyzedData/ WGS/AnalyzedData/ ;
