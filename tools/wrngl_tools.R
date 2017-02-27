getOrdrAndGenera <- function(nid) {
  order <- genus <- NA
  lng <- getNdLng(tree, nid)
  pssbls <- lng[lng %in% orders]
  if(length(pssbls) > 0) {
    order <- pssbls[length(pssbls)]
  }
  pssbls <- lng[lng %in% genera]
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
    age <- as.numeric(rownames(ed_slice)[j])
    tmp <- data.frame(t0=log(as.numeric(t0)),
                      t1=log(as.numeric(t1)),
                      id=names(t0), cnt=1,
                      tmsplt=tmsplt, age=age,
                      stringsAsFactors=FALSE)
    tmp <- na.omit(tmp)
    tmp[['n']] <- length(unique(tmp[['id']]))
    tmp
  }
  t0t1s <- data.frame(t0=NA, t1=NA, tmsplt=NA, id=NA, cnt=NA,
                         n=NA, age=NA)
  for(ed_file in ed_files) {
    i <- which(ed_files == ed_file)
    cat('....[', i, '/', length(ed_files),
        ']\n', sep='')
    if(!file.exists(ed_file)) {
      next
    }
    load(ed_file)
    tmp <- plyr::mdply(1:(nrow(ed_slice) - 1), extrct)[ ,-1]
    t0t1s <- rbind(t0t1s, tmp)
    rm(ed_slice)
  }
  t0t1s <- t0t1s[-1, ]
  t0t1s$ed <- (t0t1s$t0 + t0t1s$t1)/2
  mdl_data <- plyr::ddply(t0t1s, c('id', 'tmsplt', 'age'), plyr::summarise,
                          t0=mean(t0), t1=mean(t1), mean_ed=mean(ed),
                          sd_ed=sd(ed), cnt=sum(cnt), n=mean(n))
  mdl_data
}