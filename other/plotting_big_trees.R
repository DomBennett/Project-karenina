
countPtid <- function(n) {
  if(!is.null(n[['ptid']])) {
    n <- length(n[['ptid']])
  } else {
    n <- 0
  }
  n
}
nptids <- unlist(lapply(tree@nodelist, countPtid))
hist(nptids)
which(nptids == 1)

tree@nodelist[["n1806"]]


data("mammals")
nptids_m <- unlist(lapply(mammals@nodelist, countPtid))
hist(nptids_m)
which(nptids_m > 100)





library(ggplot2)

plotKidsAge <- function(tree1, tree2) {
  kids1 <- getNodesKids(tree1, tree1@nodes)
  nkids1 <- log(unlist(lapply(kids1, function(k) length(k))))
  ages1 <- getNodesAge(tree1, tree1@nodes)
  kids2 <- getNodesKids(tree2, tree2@nodes)
  nkids2 <- log(unlist(lapply(kids2, function(k) length(k))))
  ages2 <- getNodesAge(tree2, tree2@nodes)
  type <- c(rep(1, length(nkids1)), rep(2, length(nkids2)))
  p_data <- data.frame(nkids=c(nkids1, nkids2), ages=c(ages1, ages2),
                       type)
  p <- ggplot(p_data, aes(x=ages, y=nkids, colour=factor(type))) +
    stat_smooth(aes(fill=factor(type))) + geom_point(alpha=0.05) + scale_x_reverse()
  p + theme_bw()
}


plotKidsAge(mammals, tree)

data(mammals)
tree1 <- mammals
tree2 <- tree

plotKidsAge(mammals, tree)

hist(log(nkids))
sum(log(nkids) == log(2))

tree@nodelist[[10922]]
log(nkids)[10922]


load('1_pin/Mammalia_real/1.RData')

load('big_tree.RData')
writeTree(tree, file='big_tree.tre')


writeTree <- function(tree, file, nodeLabels=function(n){NULL}) {
  tipBytip <- function(i) {
    ids <- c(ndlst[[prid]][['kids']], prid,
             ndlst[[prid]][['prid']])
    id <<- ids[!ids %in% deja_vues][1]
    deja_vues[i] <<- id
    spn <- ndlst[[id]][['span']]
    if(id %in% tids) {
      dpth <- which(ndlst[[id]][['prid']] == prid) - 1
      prid <<- ndlst[[id]][['prid']][[1]]
      tpstr <- paste0(id, ':', spn)
      if(dpth > 0) {
        brckts <- paste0(rep('(', dpth), collapse='')
        trstr <<- paste0(trstr, ',', brckts, tpstr)
      } else {
        trstr <<- paste0(trstr, ',', tpstr)
      }
    } else {
      prid <<- ndlst[[id]][['prid']][[1]]
      ndlbl <- nodeLabels(ndlst[[id]])
      trstr <<- paste0(trstr, ')', ndlbl,':', spn)
    }
    NULL
  }
  # start with first tip
  # loop through tree structure adding tip by tip to string
  # unpack
  ndlst <- tree@nodelist
  tids <- tree@tips
  nids <- tree@nodes
  rid <- tree@root
  # add first tip
  id <- tids[1]
  trstr <-  ''
  deja_vues <- rep(NA, length(ndlst))
  deja_vues[1] <- id
  spn <- ndlst[[id]][['span']]
  dpth <- length(ndlst[[id]][['prid']])
  prid <- ndlst[[id]][['prid']][[1]]
  tpstr <- paste0(id, ':', spn)
  trstr <- paste0(rep('(', dpth), collapse='')
  trstr <- paste0(trstr, tpstr)
  # loop through nodes
  m_ply(2:(length(ndlst) - 1), .fun=tipBytip)
  trstr <- paste0(trstr, ');')
  write.table(x=trstr, file=file, quote=FALSE, row.names=FALSE,
              col.names=FALSE)
}



ndlst[["n36"]]
