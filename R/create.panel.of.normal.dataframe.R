create.panel.of.normal.dataframe <- function(files=NULL) {
  if (is.null(files)) stop("Mandatory argument files is missing")

  panel.data <- data.frame()
  for(i in 1:length(files)) {
    data <- load(files[i])
    panel.data <- rbind(panel.data, get(data))
  }
  panel.data <- panel.data %>% mutate(snvid=paste(
    sep='_',
    annovar_chr,
    annovar_start,
    annovar_end,
    annovar_ref,
    annovar_alt
    ))
  return(panel.data)
  }
