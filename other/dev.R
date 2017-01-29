tree <- rmTips(tree, 'Hylobates_sp')
.getTxnyms <- function(txnym, ...) {
  txnym
}
txnyms <- plyr::mlply(tree@ndlst, .fun=.getTxnyms)
txnyms <- txnyms[1:length(txnyms)]
names(txnyms) <- names(tree@ndlst)
tree_age <- getAge(tree)

tid <- "Hylobates_sp"
lng <- c("Chordata", "Mammalia", "Primates", "Hylobatidae", "Hylobates", "sp")
end <- 0.017
for(j in length(lng):1) {
  spns <- names(txnyms)[which(txnyms %in% lng[j])]
  if(length(spns) == 0) {
    next
  }
  spns <- unique(c(spns, unlist(getNdsPtids(tree, spns))))
  spns <- spns[spns != tree@root]
  rngs <- getSpnsAge(tree, ids=spns, tree_age=tree_age)
  bool <- rngs[ ,'start'] > end
  if(any(bool)) {
    rngs <- rngs[bool, ]
    print(rngs)
  }
}