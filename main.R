library(tercen)
library(dplyr)
library(tidyr)

is.POSIXct <- function(x) inherits(x, "POSIXct")

straf2freqtab <- function(df) {
  
  filename <- tempfile()
  writeBin(ctx$client$fileService$download(df$documentId[1]), filename)
  on.exit(unlink(filename))
  
  dat <- readLines(filename)
  dat <- strsplit(dat, "[\t]")
  
  mat <- matrix(
    unlist(dat),
    nrow = length(dat),
    ncol = length(dat[[1]]),
    byrow = TRUE
  )
  mat[mat == "0"] <- NA ###
  colnames(mat) <- mat[1, ]
  rownames(mat) <- mat[, 1]
  mat <- mat[-1, ]
  loci <- unique(colnames(mat[, -1:-2]))
  freqTAB <- NULL
  mat2 <- sub("[.]", "-", mat)
  
  for(i in 1:length(loci)) {
    
    ids <- which(colnames(mat) == loci[i])
    alleles <- unique(c(mat[, ids]))
    alleles <- sub("[.]", "-", alleles)
    alleles <- alleles[!is.na(alleles)] ###
    nameCol <- paste(loci[i], ".", alleles, sep = "")
    
    newmat <- matrix(
      NA,
      ncol = length(nameCol),
      nrow = dim(mat)[1]
    )
    
    for(ii in 1:length(alleles)) {
      newmat[,ii] <- apply(mat2[,ids]==alleles[ii],1,sum)
      colnames(newmat) <- nameCol
    }
    
    freqTAB <- cbind(freqTAB, newmat)
  }
  
  rownames(freqTAB) <- mat[, 1]
  colnames(freqTAB) <- sub(" ", "", colnames(freqTAB))
  freqTAB <- as.data.frame(freqTAB)
  freqTAB$population <- mat[, "pop"]
  freqTAB$sample <- mat[, 1]
  
  freqTAB <- freqTAB %>% 
    as.data.frame() %>%
    gather(key = "allele", value = "frequency", -sample, -population) %>%
    mutate(loc_allele = allele) %>%
    separate(allele, c("locus", "allele"), "\\.", extra = "merge", remove = FALSE) %>%
    mutate_if(is.POSIXct, as.character) %>%
    mutate_if(is.logical, as.character) %>%
    mutate_if(is.integer, as.double) %>%
    mutate(.ci= rep_len(df$.ci[1], nrow(.)))
  
  return(freqTAB)
}

ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required") 

df <- ctx$cselect() %>% 
  mutate(.ci= 1:nrow(.)-1) %>%
  split(.$.ci) %>%
  lapply(straf2freqtab) %>%
  bind_rows() %>%
  ctx$addNamespace() %>%
  ctx$save()
