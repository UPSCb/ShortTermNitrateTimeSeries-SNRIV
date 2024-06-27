files<-list.files("data/raw", pattern="R1_001.fastq.gz", full.names=TRUE)
readr::write_csv(tibble::tibble(sample=sub("_S.*","",basename(files)),
                                fastq_1=files,
                                fastq_2=sub("R1","R2",files),
                                strandedness="auto"),file="doc/sample_sheet.csv")
