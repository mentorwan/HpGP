


seq_gene_loc <- function(gene, report.file, id.col=1, start.col=7, end.col=8) {

  ret <- NULL
  if (length(report.file) == 1) {
    report.file <- scan(report.file, what="character", sep="\n", quiet=TRUE)
  }
  tmp <- grep(gene, report.file, fixed=TRUE)
  if (!length(tmp)) return(NULL)

  row <- tmp[1]
  y   <- report.file[row+1]
  if (substr(y, 1, 2) != "##") return(NULL)
  vec <- strsplit(y, "\t", fixed=TRUE)[[1]]
  if (length(vec) < end.col) return(NULL)

  id    <- vec[id.col]
  id    <- gsub("#", "", id, fixed=TRUE)
  start <- as.numeric(vec[start.col])-1
  end   <- as.numeric(vec[end.col])-1

  list(id=id, start=start, end=end)

} # END: seq_gene_loc


forward_reverse_count <- function(gene, report.file) {
  ret <- NULL
  if (length(report.file) == 1) {
    report.file <- scan(report.file, what="character", sep="\n", quiet=TRUE)
  }
  tmp <- grep(gene, report.file, fixed=TRUE)
   if (!length(tmp)) return(NULL);
   
   row <- tmp[1]
  y   <- report.file[row+1]
  vec <- strsplit(y, "\t", fixed=TRUE)[[1]]
  
  start <- as.numeric(vec[7])-1
  end   <- as.numeric(vec[8])-1
  
  

  forwardcount =0
  reversecount =0
  row <- tmp[length(tmp)]
  y <- report.file[row+1]
  
  while (!is.na(y) & substr(y,1,2) != "HP") 
  {
   vec <- strsplit(y, "\t", fixed=TRUE)[[1]]
  if(grepl("\\+",y)){
      pos <- as.numeric(vec[2])
      
      if((pos+length(vec[5]))<=end) forwardcount=forwardcount+1
  } else if (grepl("\\-",y)) { 
      
      pos <- as.numeric(vec[2])
      if( pos <= end ) reversecount=reversecount+1
  } 
      row = row+1
      y <- report.file[row+1]
      
  }
  cat(forwardcount,reversecount,"\n")
}   #END: forward_reverse_count



# Function to get subset of fasta file
seq_fasta_subset <- function(fasta, id, start, end) {

  if (length(fasta) == 1) {
    fasta <- scan(fasta, what="character", sep="\n", quiet=TRUE)
  }

  tmp <- grep(id, fasta, fixed=TRUE)
  m   <- length(tmp)
  if (m > 1) {
    str <- paste("ERROR in seq_fasta_subset with id = ", id, sep="")
    print(str)
    stop(str)
  }
  if (!m) return(NULL)
  n <- length(fasta)
  x <- fasta[tmp:n]
  n <- length(x)
  if (n < 2) stop("ERROR2: in seq_fasta_subset")
  x <- x[-1] # Remove id row

  tmp <- grep(">", x, fixed=TRUE)
  if (length(tmp)) x <- x[1:(tmp[1]-1)]
  ret <- paste(x, collapse="", sep="")
  len <- nchar(ret)
  if (len < end) stop("ERROR3: in seq_fasta_subset")
  ret <- substr(ret, start, end)
  ret <- seq_norm(ret)

  ret  

} # END: seq_fasta_subset

# Function to normalize sequences
seq_norm <- function(seq) {

  seq <- toupper(seq)
  seq <- gsub("U", "T", seq, fixed=TRUE)

  seq

} # END: seq_norm

# Function to get reverse strand
seq_reverse_strand <- function(seq) {

  for (i in 1:length(seq)) {
    vec     <- strsplit(seq[i], "", fixed=TRUE)[[1]]
    tA      <- vec %in% "A"
    tC      <- vec %in% "C"
    tG      <- vec %in% "G"
    tT      <- vec %in% "T"
    vec[tA] <- "T"
    vec[tC] <- "G"
    vec[tG] <- "C"
    vec[tT] <- "A"
    vec     <- rev(vec)

    seq[i]  <- paste(vec, collapse="", sep="")
  }

  seq

} # END: seq_reverse_strand

# Function to get all possible sequences given the map and seuence
seq_allPossible <- function(seq, map) {

  vec <- strsplit(seq, "", fixed=TRUE)[[1]]
  n   <- length(vec)
  tmp <- !(vec %in% c("A", "C", "G", "T"))
  m   <- sum(tmp)
  if (!m) return(seq)
  miss    <- (1:n)[tmp]
  letters <- vec[tmp]
  ncomb   <- 1     
  for (j in 1:m) {
    ncomb <- ncomb*length(map[[letters[j], exact=TRUE]])
  }
  if (!ncomb) stop("ERROR in seq_allPossible")
    
  # From a matrix of all combinations
  mat <- matrix(vec, byrow=TRUE, nrow=ncomb, ncol=n)
  for (j in 1:m) {
     vec  <- map[[letters[j]]]
     nvec <- length(vec)
     if (j == 1) {
       m <- ncomb/nvec
     } else {
       m <- m/nvec
     }
     mat[, miss[j]] <- rep(vec, each=m)
  }

  # Convert to character strings
  ret <- rep(seq, ncomb)
  for (j in 1:ncomb) ret[j] <- paste(mat[j, ], collapse="", sep="")

  ret

} # END: seq_allPossible

# Function to count the number of occurrences
seq_count <- function(seq, pat, map) {

  # Get all possible sequences
  allseq <- seq_allPossible(pat, map)
  sum <- 0
  for (i in 1:length(allseq)) {
    sum <- sum + str_count(seq, paste0("(?=",pattern=allseq[i],")"))
  } 
   
  sum

} # END: seq_count

seq_gene_sub <- function(gene, list.file, fasta.file, report.file, out.file) {

  # gene           Gene name to search for in the report file
  # list.file      File containing 1 column of (sequences) patterns to match
  # fasta.file     fasta file
  # report.file
  # out.file       NULL or path to output file to save results

  # Map for characters
  map <- list(R=c("A", "G"), Y=c("C", "T"), S=c("G", "C"), W=c("A", "T"), K=c("G", "T"), M=c("A", "C"),
              B=c("C", "G", "T"), D=c("A", "G", "T"), H=c("A", "C", "T"), V=c("A", "C", "G"),
              N=c("A", "C", "G", "T"))

  # Read in the patterns to match
  patterns  <- scan(list.file, what="character", quiet=TRUE)
  npatterns <- length(patterns)
  if (!npatterns) stop("ERROR: no patterns to match")
  patterns  <- seq_norm(patterns)

  # Get the id and genomic range
  tmp   <- seq_gene_loc(gene, report.file)
  if (is.null(tmp)) {
    ret <- data.frame(patterns, -1, -1, stringsAsFactors=FALSE)
    colnames(ret) <- c("Pattern", "Count_forward", "Count_reverse")
    if (!is.null(out.file)) write.table(ret, file=out.file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    return(ret)
  }
  id    <- tmp$id
  start <- tmp$start
  end   <- tmp$end

  # Read in fasta file and subset
  x <- seq_fasta_subset(fasta.file, id, start, end)

  cntfow <- rep(0, npatterns)
  cntrev <- rep(0, npatterns)

  # Get counts on forward strand
  for (i in 1:npatterns) {
    cntfow[i] <- seq_count(x, patterns[i], map)
  }  

  # Get reverse strand
  x <- seq_reverse_strand(x)

  # Get counts on revese strand
  for (i in 1:npatterns) {
    cntrev[i] <- seq_count(x, patterns[i], map)
  }  

  ret <- data.frame(patterns, cntfow, cntrev, stringsAsFactors=FALSE)
  colnames(ret) <- c("Pattern", "Count_forward", "Count_reverse")

  if (!is.null(out.file)) write.table(ret, file=out.file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

  ret

} # END: seq_gene_sub


