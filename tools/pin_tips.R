# Correctable pinTips functions
pinTips2 <- function(tree, tids, lngs, end_ages, tree_age) {
  .getPtntls <- function(lng, end) {
    sccs <- FALSE
    for (i in length(lng):1) {
      pull <- txnyms %in% lng[i]
      if (sum(pull) == 0) {
        next
      }
      ptntls <- names(txnyms)[pull]
      if (length(ptntls) == 1) {
        prnt <- ptntls
      }
      else {
        prnt <- ptntls[which.max(spn_data[ptntls, "end"])]
      }
      ptntls <- c(prnt, getNdPtids(tree, prnt))
      ptntls <- ptntls[ptntls != rid]
      pull <- spn_data[ptntls, "start"] > end
      if (any(pull)) {
        ptntls <- ptntls[pull]
        sccs <- TRUE
        break
      }
    }
    if (sccs) {
      return(ptntls)
    }
    NULL
  }
  .getPTxnym <- function(tip_txnym, sid) {
    gp_txnym <- txnyms[[getNdSlt(tree, "prid", sid)]]
    s_txnym <- txnyms[[sid]]
    if (s_txnym == tip_txnym) {
      pid_txnym <- tip_txnym
    }
    else {
      pid_txnym <- gp_txnym
    }
    pid_txnym
  }
  .pin <- function(i) {
    tid <- tids[i]
    end <- end_ages[i]
    lng <- lngs[[i]]
    ptntls <- .getPtntls(lng, end)
    if (is.null(ptntls)) {
      message(paste0("[", tid, "] could not be added"))
      return(NULL)
    }
    rngs <- spn_data[ptntls, , drop = FALSE]
    rngs[rngs[, "end"] <= end, "end"] <- end
    prbs <- rngs[, "start"] - rngs[, "end"]
    if (sum(prbs) < 1E-8) {
      message(paste0("[", tid, "] could not be added"))
      return(NULL)
    }
    sid <- as.vector(sample(ptntls, prob = prbs, size = 1))
    start <- runif(min = rngs[sid, "end"], max = rngs[sid, 
                                                      "start"], n = 1)
    tip_txnym <- lng[length(lng)]
    pid_txnym <- .getPTxnym(tip_txnym, sid)
    pid <- paste0("p_", tid, sep = "")
    tree <- addTip(tree, tid = tid, sid = sid, strt_age = start, 
                   end_age = end, pid = pid, tree_age = tree_age)
    tree@ndlst[[tid]][["txnym"]] <- tip_txnym
    tree@ndlst[[pid]][["txnym"]] <- pid_txnym
    tid_spn <- getSpnAge(tree, tid, tree_age)
    spn_data[tid, "start"] <<- tid_spn[, "start"]
    spn_data[tid, "end"] <<- tid_spn[, "end"]
    pid_spn <- getSpnAge(tree, pid, tree_age)
    spn_data[pid, "start"] <<- pid_spn[, "start"]
    spn_data[pid, "end"] <<- pid_spn[, "end"]
    sid_spn <- getSpnAge(tree, sid, tree_age)
    spn_data[sid, "start"] <<- sid_spn[, "start"]
    spn_data[sid, "end"] <<- sid_spn[, "end"]
    txnyms[[tid]] <<- tip_txnym
    txnyms[[pid]] <<- pid_txnym
    tree <<- tree
  }
  .testLngs <- function(lng) {
    for (l in lng) {
      if (grepl("[^a-zA-Z_0-9]", l)) {
        stop(paste0("Unsuitable characters in [", l, 
                    "]"))
      }
    }
    NULL
  }
  if (!tree@wtxnyms) {
    stop("tree has no txnyms")
  }
  if (any(end_ages < 0)) {
    warning("One or more end ages are less than zero, this will change the age of tree.")
  }
  mapply(.testLngs, lngs)
  txnyms <- getNdsSlt(tree, "txnym", tree@all)
  txnyms <- c(txnyms, rep(NA, length(tids) * 2))
  names(txnyms) <- c(names(tree@ndlst), tids, paste0("p_", 
                                                     tids))
  spn_data <- matrix(NA, nrow = (length(tids) * 2) + tree@nall, 
                     ncol = 2)
  colnames(spn_data) <- c("start", "end")
  tmp_spn_data <- getSpnsAge(tree, tree@all, tree_age)
  rownames(spn_data) <- c(tree@all, tids, paste0("p_", tids))
  spn_data[tree@all, "start"] <- tmp_spn_data[["start"]]
  spn_data[tree@all, "end"] <- tmp_spn_data[["end"]]
  rm(tmp_spn_data)
  rid <- tree@root
  ordrd <- order(end_ages, decreasing = TRUE)
  plyr::m_ply(ordrd, .pin)
  tree
}