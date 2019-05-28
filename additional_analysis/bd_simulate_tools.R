# Libs ----
library(treeman)
makeMdlData <- calcEDBySlice <- NULL
source(file.path('tools', 'slice_tools.R'))
source(file.path('tools', 'wrngl_tools.R'))

# Functions ----
bd_simulate <- function(ntips, b = 3, b_true = 1, d = 1, burnin = 10,
                        type = c('bd', 'pan', 'de', 'pf', 'eph')) {
  type <- match.arg(type)
  tree <- as(ape::rtree(round(ntips/10)), 'TreeMan')
  tree_age <- getAge(tree)
  extnt <- tree["tips"]
  extnct <- NULL
  i <- 0L
  ntip_total <- length(extnt) + length(extnct)
  while (ntip_total <= ntips) {
    i <- i + 1L
    if (length(extnt) < 3) {
      message('trying again....')
      return(bd_simulate(ntips, b, b_true, d, burnin, type))
    }
    if (i > burnin) {
      b <- b_true
    }
    if (type != 'bd') {
      # calculate fair proportion
      fps <- calcPrtFrPrp(tree, tids = extnt, ignr = extnct)
      sp_prob <- switch(type, pan = 1/fps, de = 1/fps, eph = fps, pf = fps)
      ex_prob <- switch(type, pan = 1/fps, de = fps, eph = fps, pf = 1/fps)
      
    } else {
      ex_prob <- sp_prob <- NULL
    }
    # add/remove based on b and d
    to_add <- sample(c(TRUE, FALSE), size = 1, prob = c(b,d))
    if (to_add) {
      sid <- sample(extnt, prob = sp_prob, size = 1)  # sister ID of the new tip
      tid <- paste0('t', ntip_total + 1)  # new tip ID
      tree <- addTip(tree, tid = tid, sid = sid, strt_age = 0, end_age = 0,
                     tree_age = tree_age)
      extnt <- c(extnt, tid)
    } else {
      tid <- sample(extnt, prob = ex_prob, size = 1)
      extnct <- c(extnct, tid)
      extnt <- extnt[extnt != tid]
    }
    # grow tree
    spns <- getNdsSlt(tree, slt_nm = "spn", ids = extnt)
    tree <- setNdsSpn(tree, ids = extnt, vals = spns + 1)
    tree_age <- tree_age + 1
    # update while logic
    ntip_total <- length(extnt) + length(extnct)
  }
  tree
}

slice_and_save <- function(tree, i, flpth) {
  # age = 1
  spns <- getNdsSlt(tree = tree, slt_nm = 'spn', ids = tree['all'])
  spns <- spns/getAge(tree)
  tree <- setNdsSpn(tree, ids = tree['all'], vals = spns)
  # print(getAge(tree))
  # unique IDs
  ids <- tree@all
  new_ids <- paste0(i, '_', ids)
  tree <- setNdsID(tree, ids = ids, vals = new_ids)
  # slice, only consider second half of tree to ignore burnin
  intrvls <- 0.5/ncuts
  time_cuts <- seq(intrvls, 0.5, intrvls)
  ed_slice <- calcEDBySlice(tree, time_cuts)
  # save
  save(ed_slice, file = file.path(flpth, paste0(i, '.RData')))
}

tree_simulate <- function(type, flpth, overwrite = FALSE) {
  flpth <- folder_gen(file.path(flpth, 'trees'))
  itrns <- seq_len(niterations)
  if (overwrite) {
    itrns <- itrns
  } else {
    itrns <- itrns[!file.exists(file.path(flpth, paste0(itrns, '.RData')))]
    if (length(itrns) == 0) {
      return(invisible(NULL))
    }
  }
  foreach(i = itrns) %dopar% {
    print(i)
    tree <- bd_simulate(ntips = ntips, type = type)
    saveRDS(object = tree, file = file.path(flpth, paste0(i, '.RData')))
    invisible(NULL)
  }
}

slice <- function(flpth) {
  tree_files <- list.files(path = file.path(flpth, 'trees'), pattern = '.RData')
  for (i in seq_along(tree_files)) {
    tree <- readRDS(file = file.path(file.path(flpth, 'trees'),
                                     tree_files[[i]]))
    slice_and_save(tree = tree, i = i, flpth = file.path(flpth, 'slices'))
  }
}

makeMdlData <- function(flpth) {
  ed_files <- list.files(path = file.path(flpth, 'slices'), pattern = '.RData')
  ed_slice <- NULL
  extrct <- function(j) {
    t0 <- ed_slice[j, ]
    t1 <- ed_slice[j + 1, ]
    tmsplt <- paste0(rownames(ed_slice)[j:(j + 1)], collapse = '-')
    age <- as.numeric(rownames(ed_slice)[j])
    tmp <- data.frame(t0 = log(as.numeric(t0)), t1 = log(as.numeric(t1)),
                      id = names(t0), cnt = 1, tmsplt = tmsplt, age = age,
                      stringsAsFactors = FALSE)
    tmp <- na.omit(tmp)
    tmp[['n']] <- length(unique(tmp[['id']]))
    tmp
  }
  t0t1s <- data.frame(t0 = NA, t1 = NA, tmsplt = NA, id = NA, cnt = NA, n = NA,
                      age = NA)
  for (ed_file in ed_files) {
    i <- which(ed_files == ed_file)
    cat('....[', i, '/', length(ed_files),
        ']\n', sep = '')
    if (!file.exists(ed_file)) {
      next
    }
    load(ed_file)
    tmp <- try(plyr::mdply(1:(nrow(ed_slice) - 1), extrct)[ ,-1], silent = TRUE)
    if (inherits(tmp, 'try-error')) {
      next
    }
    t0t1s <- rbind(t0t1s, tmp)
    rm(ed_slice)
  }
  t0t1s <- t0t1s[-1, ]
  t0t1s$ed <- (t0t1s$t0 + t0t1s$t1)/2
  t0 <- t1 <- ed <- cnt <- n <- NULL
  mdl_data <- plyr::ddply(t0t1s, c('id', 'tmsplt', 'age'), plyr::summarise,
                          t0 = mean(t0), t1 = mean(t1), mean_ed = mean(ed),
                          sd_ed = sd(ed), cnt = sum(cnt), n = mean(n))
  mdl_data
}

folder_gen <- function(flpth) {
  if (!dir.exists(flpth)) {
    dir.create(flpth)
  }
  flpth
}