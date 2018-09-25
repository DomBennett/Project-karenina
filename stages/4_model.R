# D.J. Bennett
# 18/07/2015
# Model EDt1 ~ EDt0

# PARAMETERS
source('parameters.R')

# LIBS
library(ggplot2)
library(gridExtra)
library(lme4)

# DIRS
input_dir <- '3_wrngl'
output_dir <- '4_model'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

# INPUT
load(file=file.path(input_dir, paste0(parent, '.RData')))

# QUICK LOOK
p1 <- ggplot(mdl_data, aes(x=t0, y=t1, colour=epoch)) +
  geom_point(alpha=0.05) + theme_bw() +
  xlab(expression('ED'['t0'])) +
  ylab(expression('ED'['t1'])) +
  theme(legend.position='none')
p2 <- ggplot(mdl_data, aes(x=t0, y=t1, colour=epoch)) +
  geom_smooth(se=TRUE, formula=y~x) + theme_bw() +
  xlab(expression('ED'['t0'])) +
  ylab('') +
  theme(legend.title=element_blank())
text_tt <- theme(text=element_text(size=6))
tiff(file.path('4_model', 'gam_and_points.tiff'), width=14, height=9, units="cm",
     res=1200)
grid.arrange (p1+text_tt, p2+text_tt, ncol=2)
dev.off()
p1 <- ggplot(rnd_data, aes(x=t0, y=t1, colour=epoch)) +
  geom_point(alpha=0.05) + theme_bw() +
  xlab(expression('ED'['t0'])) +
  ylab(expression('ED'['t1'])) +
  theme(legend.position='none')
p2 <- ggplot(rnd_data, aes(x=t0, y=t1, colour=epoch)) +
  geom_smooth(se=TRUE, formula=y~x) + theme_bw() +
  xlab(expression('ED'['t0'])) +
  ylab('') +
  theme(legend.title=element_blank())
text_tt <- theme(text=element_text(size=6))
tiff(file.path('4_model', 'gam_and_points_rand.tiff'), width=14, height=9, units="cm",
     res=1200)
grid.arrange (p1+text_tt,
              p2+text_tt, ncol=2)
dev.off()

# BASIC STATS
nrow(mdl_data)
quantile(mdl_data$cnt)
nrow(rnd_data)
quantile(rnd_data$cnt)
p_data <- data.frame(cnt=c(mdl_data$cnt, rnd_data$cnt),
                     real=c(rep('Real', nrow(mdl_data)),
                            rep('Random', nrow(rnd_data))))
ggplot(p_data, aes(cnt, colour=real, fill=real)) + geom_density(alpha=0.5)
# mean number of species in t0 by epoch
tapply(mdl_data$n, mdl_data$epoch, mean)

# SORT DATA
# add tax info to rnd
ints <- match(rnd_data$id, mdl_data$id)
rnd_data$genus <- mdl_data$genus[ints]
rnd_data$order <- mdl_data$order[ints]
# drop last two epochal transitions, seem to be different
mdl_data <- mdl_data[!mdl_data$epoch %in% c("JU-CL", "CL-CU"), ]
rnd_data <- rnd_data[!rnd_data$epoch %in% c("JU-CL", "CL-CU"), ]
# log time difference
ages <- tapply(mdl_data$age, mdl_data$epoch, mean)
tm <- ages - c(0, ages[-(length(ages))])
mdl_data$tm <- tm[mdl_data$epoch]
order_data <- mdl_data[!is.na(mdl_data[['order']]), ]
genus_data <- order_data[!is.na(order_data[['genus']]), ]
genus_data$age <- log(genus_data$age)
genus_data$n <- log(genus_data$n)
genus_data$tm <- log(genus_data$tm)
# t0 dummy, rounded t0 between 0 and 1
# I found this method better than other methods because it does not over-represent
# the extreme values
# bins <- seq(min(genus_data$t0)-0.00001, max(genus_data$t0), length.out=5)
# genus_data$t0_dummy <- as.numeric(cut(genus_data$t0, bins))
# genus_data$t0_dummy <- as.numeric(genus_data$t0 > mean(genus_data$t0))
genus_data$t0_dummy <- round(genus_data$t0)
genus_data$t0_dummy <- genus_data$t0_dummy/max(genus_data$t0_dummy)
hist(genus_data$t0_dummy)
plot(genus_data$t0_dummy, genus_data$t0)
# clean up
rm(order_data)

# ASSESS WHETHER RND AND REAL DIFFER
tree_ids <- mdl_data$id[!as.logical(mdl_data$fssl_nd)]
all_data <- data.frame(ids=c(mdl_data$id, rnd_data$id),
                       t0=c(mdl_data$t0, rnd_data$t0),
                       t1=c(mdl_data$t1, rnd_data$t1),
                       sd=c(mdl_data$sd_ed, rnd_data$sd_ed),
                       epoch=c(mdl_data$epoch, rnd_data$epoch),
                       real=c(rep('Real', nrow(mdl_data)),
                              rep('Random', nrow(rnd_data))))
all_data$fssl_nd <- !all_data$ids %in% tree_ids
all_data$tdff <- all_data$t0 - all_data$t1
xlbl <- expression(paste('ED'['t0'], ' - ED'['t1']))
# tdff
p1 <- ggplot(all_data, aes(tdff, colour=real, fill=real)) +
  geom_density(alpha=0.5) + theme_bw() + xlab(xlbl) + ylab('') +
  theme(legend.position='none')
# tdff without tree ids
ggplot(all_data[all_data$fssl_nd, ], aes(tdff, colour=real, fill=real)) +
  geom_density(alpha=0.5) + theme_bw() + xlab(xlbl) + ylab('') +
  theme(legend.title = element_blank())
# F-test
rl_tdffs <- all_data[all_data$real == 'Real', 'tdff']
rnd_tdffs <- all_data[all_data$real == 'Random', 'tdff']
var.test(rnd_tdffs, rl_tdffs)
var(rnd_tdffs, na.rm=TRUE)
var(rl_tdffs, na.rm=TRUE)
# check SD
p2 <- ggplot(all_data[all_data$fssl_nd, ], aes(sd, colour=real, fill=real)) +
  geom_density(alpha=0.5) + theme_bw() + xlab('Std. Dv. of ED') + ylab('') +
  theme(legend.title = element_blank())
# F-test
rl_sds <- all_data[all_data$real == 'Real' & all_data$fssl_nd, 'sd']
rnd_sds <- all_data[all_data$real == 'Random' & all_data$fssl_nd, 'sd']
t.test(rnd_sds, rl_sds)
mean(rnd_sds, na.rm=TRUE)
mean(rl_sds, na.rm=TRUE)
# save plots
text_tt <- theme(text=element_text(size=6))
tiff(file.path('4_model', 'diff_rand_real.tiff'), width=14, height=9, units="cm",
     res=1200)
grid.arrange(p1+text_tt, p2+text_tt, ncol=2)
dev.off()

# CHECK PINNING/FOSSILS
# is there a relationship between ED and number of times a clade appears?
genus_data$mean_ed <- (genus_data$t0+genus_data$t1)/2
m <- lm(mean_ed~cnt, data=genus_data)
summary(m) # VERY small decrease, little evidence that random clades are steering the results
p1 <- ggplot(genus_data, aes(y=mean_ed, x=cnt)) + geom_point(alpha=0.1) +
  geom_smooth(method='lm', se = TRUE) + xlab('No. clade counts') +
  ylab('ED') + theme_bw()
# is there a difference between fossil and non-fossil nodes?
genus_data$Fossil <- genus_data$fssl_nd
p2 <- ggplot(genus_data, aes(mean_ed, colour=fssl_nd, fill=Fossil)) +
  geom_density(alpha=0.5) + ylab('') + xlab('ED') + theme_bw()
t.test(genus_data$mean_ed[as.logical(genus_data$fssl_nd)],
       genus_data$mean_ed[!as.logical(genus_data$fssl_nd)])
t.test(genus_data$tm[as.logical(genus_data$fssl_nd)],
       genus_data$tm[!as.logical(genus_data$fssl_nd)])
# yes, presumably because fossil nodes went extinct, while ones in the tree survived


# MODEL SELECTION
# exp linear model
m0 <- lm(t1~t0, data=genus_data)
m1a <- lmer(t1~t0 + (1|epoch), data=genus_data, REML=FALSE)
m1b <- lmer(t1~t0 + (t0|epoch), data=genus_data, REML=FALSE)
m2a <- lmer(t1~t0 + (t0|epoch) + (1|order), data=genus_data, REML=FALSE)
m2b <- lmer(t1~t0 + (t0|epoch) + (1|genus), data=genus_data, REML=FALSE)
m2c <- lmer(t1~t0 + (t0|epoch) + (1|order/genus), data=genus_data, REML=FALSE)
m2d <- lmer(t1~t0 + (t0|epoch) + (1|id), data=genus_data, REML=FALSE)
m2e <- lmer(t1~t0 + (t0|epoch) + (1|order/id), data=genus_data, REML=FALSE)
m2f <- lmer(t1~t0 + (t0|epoch) + (t0|order), data=genus_data, REML=FALSE)
m2g <- lmer(t1~t0 + (t0|epoch) + (t0|genus), data=genus_data, REML=FALSE)
m2h <- lmer(t1~t0 + (t0|epoch) + (t0|id), data=genus_data, REML=FALSE)
AIC(m0)
anova(m1a, m1b)
anova(m1b, m2a)
anova(m2a, m2b)
anova(m2b, m2c)
anova(m2c, m2d)
anova(m2c, m2e)
anova(m2c, m2f)
anova(m2c, m2g)
anova(m2g, m2h)
# m2g is best

# fitting polynomials
m3a <- lmer(t1~poly(t0, 2) + (t0|epoch) + (t0|genus), data=genus_data, REML=FALSE)
anova(m2g, m3a)
m3b <- lmer(t1~poly(t0, 3) + (t0|epoch) + (t0|genus), data=genus_data, REML=FALSE)
anova(m3a, m3b)
m3c <- lmer(t1~poly(t0, 4) + (t0|epoch) + (t0|genus), data=genus_data, REML=FALSE)
anova(m3b, m3c)
m3d <- lmer(t1~poly(t0, 5) + (t0|epoch) + (t0|genus), data=genus_data, REML=FALSE)
anova(m3c, m3d)
m3e <- lmer(t1~poly(t0, 6) + (t0|epoch) + (t0|genus), data=genus_data, REML=FALSE)
anova(m3d, m3e)
# opt for 3c
obs_mdl <- m3c
save(m3c, file=file.path('4_model', 'obs_mdl.RData'))

# fitting exp ln
n0a <- lm(t1~tm, data=genus_data)
n0b <- lm(t1~n, data=genus_data)
n0c <- lm(t1~tm+n, data=genus_data)
n0d <- lm(t1~t0_dummy, data=genus_data)
n0e <- lm(t1~t0_dummy+tm+n, data=genus_data)
n0f <- lm(t1~t0_dummy+tm+n+age, data=genus_data)
anova(n0a, n0b)
anova(n0a, n0c)
anova(n0c, n0d)
anova(n0d, n0e)
AIC(n0a, n0b, n0c, n0d, n0e)
n1a <- lmer(t1~t0_dummy+tm+n+(1|order), data=genus_data, REML=FALSE)
n1b <- lmer(t1~t0_dummy+tm+n+(1|genus), data=genus_data, REML=FALSE)
n1c <- lmer(t1~t0_dummy+tm+n+(1|id), data=genus_data, REML=FALSE)
anova(n1a, n1b)
anova(n1b, n1c)
n1d <- lmer(t1~t0_dummy+tm+n+(t0_dummy|genus), data=genus_data, REML=FALSE)
n1e <- lmer(t1~t0_dummy+tm+n+(tm|genus), data=genus_data, REML=FALSE)
n1f <- lmer(t1~t0_dummy+tm+n+(n|genus), data=genus_data, REML=FALSE)
anova(n1b, n1d)
anova(n1e, n1d)
anova(n1f, n1e)
n1g <- lmer(t1~t0_dummy+tm+n+(t0_dummy|genus)+(tm|genus), data=genus_data, REML=FALSE)
anova(n1e, n1g)
n1h <- lmer(t1~t0_dummy+tm+n+(1|order/genus), data=genus_data, REML=FALSE)
anova(n1h, n1g)
# opt for n1g
exp_mdl <- n1g
save(exp_mdl, file=file.path('4_model', 'exp_mdl.RData'))

# OPTED MODELS
load(file.path('4_model', 'exp_mdl.RData'))
load(file.path('4_model', 'obs_mdl.RData'))

# PLOTTING
# create representative dataset of equal numbers epoch and sample genera
t0 <- seq(min(genus_data$t0), max(genus_data$t0), length.out=100)
rpsnttv <- expand.grid(t0=t0,
                       genus=sample(unique(genus_data$genus), 100),
                       epoch=unique(genus_data$epoch))
rpsnttv$order <- genus_data$order[match(rpsnttv$genus, genus_data$genus)]
ns <- tapply(genus_data$n, genus_data$epoch, mean)
rpsnttv$n <- ns[rpsnttv$epoch]
tm <- tapply(genus_data$tm, genus_data$epoch, mean)
rpsnttv$tm <- tm[rpsnttv$epoch]
rpsnttv$t0_dummy <- rpsnttv$t0/max(rpsnttv$t0)
rpsnttv$pfit <- predict(obs_mdl, rpsnttv)
rpsnttv$nfit <- predict(exp_mdl, rpsnttv)
# overall
p_data <- plyr::ddply(rpsnttv, c('t0'), plyr::summarise,
                      pfit=median(pfit), nfit=median(nfit))
p <- ggplot(p_data, aes(x=t0, y=pfit)) + 
  geom_line(aes(y=pfit), lwd=1) +
  geom_line(aes(y=nfit), lwd=1, lty=2) +
  xlab(expression('ED'['t0'])) +
  ylab(expression('ED'['t1'])) +
  theme_bw()
tiff(file.path('4_model', 'overall.tiff'), width=9, height=9, units="cm",
     res=1200)
print(p)
dev.off()
# by epoch
p_data <- plyr::ddply(rpsnttv, c('t0', 'epoch'), plyr::summarise,
                      pfit=median(pfit), nfit=median(nfit))
p <- ggplot(p_data, aes(x=t0, y=pfit)) + 
  geom_line(aes(y=pfit), lwd=1) +
  geom_line(aes(y=nfit), lwd=1, lty=2) +
  xlab(expression('ED'['t0'])) +
  ylab(expression('ED'['t1'])) +
  theme_bw()
p <- p + facet_grid(epoch ~ .)
tiff(file.path('4_model', 'by_epoch.tiff'), width=9, height=14, units="cm",
     res=1200)
print(p + theme(text=element_text(size=6)))
dev.off()


# DIAGNOSTICS
library(sjPlot)
library(sjmisc)
sjp.lmer(exp_mdl, type = "fe")
sjp.lmer(obs_mdl, type = "fe")
sjp.lmer(obs_mdl, type = "re.qq")


# COMPARE REAL AND RANDOM POLY 3
# use (1|genus) because (t0|genus) is not computable in real time
rndmdl <- lmer(t1~poly(t0, 3)+(t0|epoch)+(1|genus), data=rnd_data, REML=FALSE)
rlmdl <- lmer(t1~poly(t0, 3)+(t0|epoch)+(1|genus), data=genus_data, REML=FALSE)
# create representative p_data
t0 <- seq(min(rnd_data$t0), max(rnd_data$t0), length.out=100)
rnd_rep <- rl_rep <- expand.grid(t0=t0, genus=sample(genus_data$genus, 100),
                                 epoch=as.character(unique(rnd_data$epoch)))
rnd_rep$fit <- predict(rndmdl, rnd_rep)
rl_rep$fit <- predict(rlmdl, rl_rep)
rpsnttv <- rbind(rnd_rep, rl_rep)
rpsnttv$real <- c(rep('Random', nrow(rnd_rep)),
                  rep('Real', nrow(rl_rep)))
p_data <- plyr::ddply(rpsnttv, c('t0', 'real'), plyr::summarise,
                      fit=mean(fit))
# plot
p <- ggplot(p_data, aes(x=t0, y=fit, colour=real)) + 
  geom_line(aes(y=fit), lwd=1) +
  xlab(expression('ED'['t0'])) +
  ylab(expression('ED'['t1'])) +
  theme_bw() + theme(legend.title=element_blank())
tiff(file.path('4_model', 'real_rndm_p3.tiff'), width=9, height=9, units="cm",
     res=1200)
print(p + theme(text=element_text(size=6)))
dev.off()
