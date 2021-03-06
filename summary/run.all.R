
source("source.all.r")
library(stringr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("run.R polished_assembly_NCI_57.fasta NCI_57.list.txt NCI_57.report.txt", call.=FALSE)
} 


#dir <- "~/Desktop/HpGP/IMS_help/script/";
dir <- getwd();
#list.files(dir)

fastaf  <- paste(dir,"/", args[1], sep="");
listf   <- paste(dir,"/", args[2], sep="");
reportf <- paste(dir,"/", args[3], sep="");

if (length(reportf) == 1) {
  reportf <- scan(reportf, what="character", sep="\n", quiet=TRUE)
}

#fastaf  <- paste(dir, "polished_assembly_NCI_57.fasta", sep="");
#listf   <- paste(dir, "NCI_57.list.txt", sep="");
#reportf <- paste(dir, "NCI_57.report.txt", sep="");
#genelist <- read.csv('~/Desktop/HpGP/IMS_help/script/genelist.txt',header=FALSE);

genelist <- read.csv('genelist.all.txt',header=FALSE);
genelist1 <- as.vector(genelist[,"V1"])


##out.file <- "~/Desktop/HpGP/IMS_help/script/output.txt.xls"
for (i in 1 : length(genelist1))
{
   gene     <- genelist1[i];
##      gene  <- "HP0024";
    ret      <- seq_gene_sub(gene, listf, fastaf, reportf,NULL)
    ##print(ret);
   summary <- c(sum(ret$Count_forward),sum(ret$Count_reverse))
   
    if(sum(ret$Count_forward)>=0 )
    { 
     cat(gene,summary," ")
    } else {
     cat(gene,"-1 -1 ")
   }
   
    forward_reverse_count(gene,reportf)
}

