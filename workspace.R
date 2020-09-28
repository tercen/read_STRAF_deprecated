library(tercen)
library(dplyr)

options("tercen.workflowId" = "b5eee6b1ed83d50347d04d6ba20a0f29")
options("tercen.stepId"     = "53ebc609-5f25-4cdd-9930-4ef5efb29c1b")
 
is.POSIXct <- function(x) inherits(x, "POSIXct")

straf_to_tercen = function(df){
  filename = tempfile()
  writeBin(ctx$client$fileService$download(df$documentId[1]), filename)
  on.exit(unlink(filename))
  
  straf_in <- readLines(filename)
  straf_in <- strsplit(straf_in, "\t")
  straf <- do.call("rbind", straf_in)

  headr <- straf[1, -1:-2]
  clsid <- straf[-1, 1:2]
  geno <- straf[-1, -1:-2]

  df_o <- list()
  for(i in 1:nrow(clsid)) {
    df_o[[i]] <- data.frame(
      sample = clsid[i, 1],
      population = clsid[i, 2],
      locus = headr,
      genotype = geno[i, ]
    )
  }
  straf_long <- do.call(rbind, df_o) %>%
    mutate_if(is.POSIXct, as.character) %>%
    mutate_if(is.logical, as.character) %>%
    mutate_if(is.integer, as.double) %>%
    mutate(.ci= rep_len(df$.ci[1], nrow(.)))
  
  return(straf_long)
}
 
ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required") 
 
ctx$cselect() %>% 
  mutate(.ci= 1:nrow(.)-1) %>%
  split(.$.ci) %>%
  lapply(straf_to_tercen) %>%
  bind_rows() %>%
  ctx$addNamespace() %>%
  ctx$save()
