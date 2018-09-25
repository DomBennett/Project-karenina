
calcFrPrp <- function(tree, tids, progress = "none") {
  .calc <- function(i) {
    id <- tree@all[i]
    spn <- getNdSlt(tree, "spn", id)
    kids <- getNdKids(tree, id)
    if (length(kids) == 0) {
      spn_shres[i, id] <<- spn
    } else {
      spn_shre <- spn/length(kids)
      spn_shres[i, kids] <<- spn_shre
    }
  }
  spn_shres <- bigmemory::big.matrix(init = 0, ncol = tree@ntips,
                                     nrow = tree@nall,
                                     shared = FALSE)
  options(bigmemory.allow.dimnames = TRUE)
  colnames(spn_shres) <- tree@tips
  plyr::m_ply(.data = data.frame(i = 1:tree@nall), .fun = .calc,
              .progress = progress)
  colSums(spn_shres[, tids])
}

timeslice <- function(tree, tm_ct, nd_spns) {
  # Internals
  rmPtnds <- function(nd) {
    ptids <- nd[['ptid']]
    nd[['ptid']] <- ptids[!ptids %in% to_rmv]
    nd
  }
  updateSpns <- function(i) {
    id <- nd_spns[['spn']][i]
    end <- nd_spns[['end']][i]
    ndlst[[id]][['spn']] <- ndlst[[id]][['spn']] - (tm_ct - end)
    ndlst[[id]]
  }
  # find all spans after the time split or that pass through it
  to_rmv <- nd_spns[['spn']][nd_spns[['start']] <= tm_ct]
  nd_spns <- nd_spns[!nd_spns[['spn']] %in% to_rmv, ]
  bool <- !tree@all %in% to_rmv
  # remove them
  ndlst <- tree@ndlst[bool]
  # remove all references to remvoed nodes
  ndlst <- lapply(ndlst, rmPtnds)
  # recalculate lengths of spans that pass through
  spns <- which(nd_spns[['end']] < tm_ct)
  ndlst[spns] <- lapply(spns, updateSpns)
  # update tree
  newtree <- tree
  newtree@ndlst <- ndlst
  newtree <- pstMnp(newtree)
  newtree <- updateSlts(newtree)
  if(!is.null(tree@ndmtrx)) {
    tree@ndmtrx <- bigmemory::as.big.matrix(tree@ndmtrx[bool, bool])
  }
  newtree
}

calcEDBySlice <- function(tree, time_cuts) {
  # Return ED values for clades at different time slices
  # Internals
  slcd <- fp_vals <- NULL
  getNdFP <- function(id) {
    kids <- getNdKids(slcd, id)
    mean(fp_vals[kids])
  }
  # for time callibrated tree
  tree_age <- getAge(tree)
  time_cuts <- time_cuts[time_cuts < tree_age]
  time_cuts <- sort(time_cuts, decreasing = TRUE)
  # gen res mtrx
  res <- matrix(NA, nrow = length(time_cuts),
                ncol = tree['nall'])
  rownames(res) <- time_cuts
  colnames(res) <- tree['all']
  nd_spns <- getSpnsAge(tree, tree@all, tree_age = tree_age)
  nd_spns[['spn']] <- as.character(nd_spns[['spn']])
  for (i in 1:length(time_cuts)) {
    # drop tips extinct by time cut
    tp_ages <- nd_spns[nd_spns[['spn']] %in% tree@tips, c('spn', 'end')]
    bool <- tp_ages[['end']] > time_cuts[i] + 0.00001 # zero tolerance control
    if (any(bool)) {
      to_drp <- tp_ages[['spn']][bool]
      tree <- rmTips(tree, tids = to_drp)
      nd_spns <- nd_spns[nd_spns[['spn']] %in% tree@all, ]
    }
    # slice tree at interval
    slcd <- timeslice(tree = tree, tm_ct = time_cuts[i],
                      nd_spns = nd_spns)
    # get ed vals (tips and mean tips for nodes)
    fp_vals <- calcFrPrp(slcd, slcd['tips'])
    kids <- getNdsKids(slcd, slcd['nds'])
    fp_vals <- c(fp_vals, sapply(slcd['nds'], getNdFP))
    # add to res
    res[i, names(fp_vals)] <- fp_vals
  }
  return(res)
}