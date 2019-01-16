# libs ----
library(treeman)
makeMdlData <- calcEDBySlice <- NULL
source(file.path('tools', 'slice_tools.R'))
source(file.path('tools', 'wrngl_tools.R'))

# functions ----
slice_and_save <- function(tree, i, flpth) {
  # age = 1
  spns <- getNdsSlt(tree = tree, slt_nm = 'spn', ids = tree['all'])
  spns <- spns/getAge(tree)
  tree <- setNdsSpn(tree, ids = tree['all'], vals = spns)
  # unique IDs
  ids <- tree@all
  new_ids <- paste0(i, '_', ids)
  tree <- setNdsID(tree, ids = ids, vals = new_ids)
  # slice
  ed_slice <- calcEDBySlice(tree, c(0.05, 0.15, 0.25, 0.35))
  # save
  save(ed_slice, file = file.path(flpth, paste0(i, '.RData')))
}
bd_simulate <- function(ntips, b = 3, b_true = 1, d = 1, burnin = 10) {
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
      return(bd_simulate(ntips, b, b_true, d, burnin))
    }
    if (i > burnin) {
      b <- b_true
    }
    # calculate fair proportion
    #fps <- calcPrtFrPrp(tree, tids=extnt, ignr=extnct)
    # add/remove based on b and d
    to_add <- sample(c(TRUE, FALSE), size = 1, prob = c(b,d))
    if (to_add) {
      #sid <- sample(extnt, prob=1/fps, size = 1)  # sister ID of the new tip
      sid <- sample(extnt, size = 1)
      tid <- paste0('t', ntip_total + 1)  # new tip ID
      tree <- addTip(tree, tid = tid, sid = sid, strt_age = 0, end_age = 0,
                     tree_age = tree_age)
      extnt <- c(extnt, tid)
    } else {
      #tid <- sample(extnt, prob=1/fps, size=1)
      tid <- sample(extnt, size = 1)
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
edbmm_simulate <- function(ntips, b = 3, b_true = 1, d = 1, burnin = 10) {
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
      return(bd_simulate(ntips, b, b_true, d, burnin))
    }
    if (i > burnin) {
      b <- b_true
    }
    # calculate fair proportion
    fps <- calcPrtFrPrp(tree, tids = extnt, ignr = extnct)
    # add/remove based on b and d
    to_add <- sample(c(TRUE, FALSE), size = 1, prob = c(b,d))
    if (to_add) {
      sid <- sample(extnt, prob = 1/fps, size = 1)  # sister ID of the new tip
      tid <- paste0('t', ntip_total + 1)  # new tip ID
      tree <- addTip(tree, tid = tid, sid = sid, strt_age = 0, end_age = 0,
                     tree_age = tree_age)
      extnt <- c(extnt, tid)
    } else {
      tid <- sample(extnt, prob = 1/fps, size = 1)
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
makeMdlData <- function(ed_files) {
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

# vars ----
ncuts <- 10
niterations <- 20
ntips <- 1000

# dirs ----
outdir <- 'expected'
if (!dir.exists(outdir)) {
  dir.create(outdir)
}
bs_dir <- file.path(outdir, 'broken_stick')
if (!dir.exists(bs_dir)) {
  dir.create(bs_dir)
}
bd_dir <- file.path(outdir, 'birth_death')
if (!dir.exists(bd_dir)) {
  dir.create(bd_dir)
}
edbmm_dir <- file.path(outdir, 'edbmm')
if (!dir.exists(edbmm_dir)) {
  dir.create(edbmm_dir)
}

# broken stick ----
for (i in 21:100) {
  # gen random tree with broken stick
  print(i)
  success <- FALSE
  while (!success) {
    success <- tryCatch(expr = {
      tree <- as(ape::rtree(ntips), 'TreeMan')
      slice_and_save(tree = tree, i = i, flpth = bs_dir)
      TRUE
    }, error = function(e) {
      FALSE
    })
  }
}
ed_files <- list.files(path = bs_dir, pattern = '.RData')
bs_data <- makeMdlData(ed_files = file.path(bs_dir, ed_files))

# birth-death ----
for (i in 11:20) {
  print(i)
  tree <- bd_simulate(ntips = ntips)
  slice_and_save(tree = tree, i = i, flpth = bd_dir)
}
ed_files <- list.files(path = bd_dir, pattern = '.RData')
bd_data <- makeMdlData(ed_files = file.path(bd_dir, ed_files))

# EDBMM ----
for (i in 11:20) {
  print(i)
  tree <- edbmm_simulate(ntips = ntips)
  slice_and_save(tree = tree, i = i, flpth = edbmm_dir)
}
ed_files <- list.files(path = edbmm_dir, pattern = '.RData')
edbmm_data <- makeMdlData(ed_files = file.path(edbmm_dir, ed_files))


all_data <- rbind(bs_data, bd_data, edbmm_data)
all_data$type <- c(rep('BS', nrow(bs_data)),
                   rep('BD', nrow(bd_data)),
                   rep('EDBMM', nrow(edbmm_data)))
# drop multiples
all_data <- all_data[duplicated(all_data$id), ]

library(ggplot2)
ggplot(all_data, aes(x = t0, y = t1)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1) +
  facet_grid(~type) + theme(legend.position = 'none')
