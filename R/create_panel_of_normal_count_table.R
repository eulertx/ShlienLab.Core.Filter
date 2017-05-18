create_panel_of_normal_count_table <- function(path=NULL, recurrence=2, pattern='annotated_filtered.rda$') {
  if (is.null(path)) stop("Mandatory argument path is missing")

  files <- list.files(
    path=pon,
    pattern=pattern,
    recursive=FALSE,
    full.names=TRUE
  )

  # create the panel of normal dataframe using available normals
  # this also adds a column representing the SNV ID which will be
  # used to identify recurrent mutations among the samples
  data <- create.panel.of.normal.dataframe(files=files)
  count <- get.recurrent.snvids(data=pon_data)
  count <- count %>% filter(count >= recurrence)

  return(count)
  }
