getOrdrAndGenera <- function(nid) {
  order <- genus <- NA
  lng <- getNdLng(tree, nid)
  pssbls <- orders[orders %in% lng]
  if(length(pssbls) > 0) {
    order <- pssbls[length(pssbls)]
  }
  pssbls <- genera[genera %in% lng]
  if(length(pssbls) > 0) {
    genus <- pssbls[length(pssbls)]
  }
  data.frame(nid, order, genus)
}

makeMdlData <- function(ed_files) {
  extrct <- function(j) {
    t0 <- ed_slice[j, ]
    t1 <- ed_slice[j+1, ]
    tmsplt <- paste0(rownames(ed_slice)[j:(j+1)], collapse='-')
    tmp <- data.frame(t0=log(as.numeric(t0)),
                      t1=log(as.numeric(t1)),
                      id=names(t0), cnt=1,
                      tmsplt=tmsplt, stringsAsFactors=FALSE)
    tmp <- na.omit(tmp)
    pull <- !tmp$id %in% mdl_data$id[mdl_data$tmsplt == tmsplt]
    if(sum(pull) > 0) {
      mdl_data <- rbind(mdl_data, tmp[pull, ])
    }
    if(sum(!pull) > 0) {
      tmp <- tmp[!pull, ]
      tmsplt_i <- which(mdl_data$tmsplt == tmsplt)
      mtchng <- match(tmp$id, mdl_data$id[tmsplt_i])
      mtchng <- tmsplt_i[mtchng]
      mdl_data[mtchng, 't0'] <-
        (mdl_data[mtchng, 't0'] + tmp[['t0']])/2
      mdl_data[mtchng, 't1'] <-
        (mdl_data[mtchng, 't1'] + tmp[['t1']])/2
      mdl_data[mtchng, 'cnt'] <- mdl_data[mtchng, 'cnt'] + 1
    }
    mdl_data <<- mdl_data
    NULL
  }
  mdl_data <- data.frame(t0=NA, t1=NA, tmsplt=NA, id=NA, cnt=NA)
  for(ed_file in ed_files) {
    i <- which(ed_files == ed_file)
    cat('....[', i, '/', length(ed_files),
        ']\n', sep='')
    if(!file.exists(ed_file)) {
      next
    }
    load(ed_file)
    tmp <- sapply(1:(nrow(ed_slice) - 1), extrct)
    rm(ed_slice)
  }
  mdl_data <- mdl_data[-1, ]
  mdl_data
}