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

# INPUT ----
load(file=file.path(input_dir, paste0(parent, '.RData')))

# QUICK LOOK ----
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
grid.arrange (p1+text_tt, p2+text_tt, ncol=2)
dev.off()

# BASIC STATS
nrow(mdl_data)
quantile(mdl_data$cnt)
nrow(rnd_data)
quantile(rnd_data$cnt)
# mean number of species in t0 by epoch
tapply(mdl_data$n, mdl_data$epoch, mean)
# n data points per epoch
data_per_epoch <- table(mdl_data$epoch)
mean(data_per_epoch[!names(data_per_epoch) %in% c('CL-CU', 'JU-CL')])

# SORT DATA ----
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
# sort data
# t0 dummy, rounded t0 between 0 and 1
# I found this method better than other methods because it does not over-represent
# the extreme values
# bins <- seq(min(genus_data$t0)-0.00001, max(genus_data$t0), length.out=5)
# genus_data$t0_dummy <- as.numeric(cut(genus_data$t0, bins))
# genus_data$t0_dummy <- as.numeric(genus_data$t0 > mean(genus_data$t0))
genus_data$t0_dummy <- round(x = genus_data$t0, digits = 0)
genus_data$t0_dummy <- genus_data$t0_dummy/max(genus_data$t0_dummy)
hist(genus_data$t0_dummy)
plot(genus_data$t0_dummy, genus_data$t0)
# clean up
rm(order_data)

# ASSESS WHETHER RND AND REAL DIFFER ----
tree_ids <- mdl_data$id[!as.logical(mdl_data$fssl_nd)]
all_data <- data.frame(ids=c(mdl_data$id, rnd_data$id),
                       t0=c(mdl_data$t0, rnd_data$t0),
                       t1=c(mdl_data$t1, rnd_data$t1),
                       sd=c(mdl_data$sd_ed, rnd_data$sd_ed),
                       epoch=c(mdl_data$epoch, rnd_data$epoch),
                       real=c(rep('Real', nrow(mdl_data)),
                              rep('Random', nrow(rnd_data))),
                       cnt=c(mdl_data$cnt, rnd_data$cnt))
all_data$fssl_nd <- !all_data$ids %in% tree_ids
all_data$tdff <- all_data$t0 - all_data$t1
# mean diffs
tapply(X = all_data$tdff, INDEX = all_data$real, FUN = mean)
# tdff
xlbl <- expression(paste('ED'['t0'], ' - ED'['t1']))
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
t.test(rnd_tdffs, rl_tdffs)
var.test(rnd_tdffs, rl_tdffs)
var(rnd_tdffs, na.rm=TRUE)
var(rl_tdffs, na.rm=TRUE)
# cnt
p2 <- ggplot(all_data, aes(cnt, colour=real, fill=real)) +
  geom_density(alpha=0.5) + theme_bw() + xlab('No. occurrences across iterations') + ylab('') +
  theme(legend.position='none')
# check SD
p3 <- ggplot(all_data[all_data$fssl_nd, ], aes(sd, colour=real, fill=real)) +
  geom_density(alpha=0.5) + theme_bw() + xlab('Std. Dv. of ED') + ylab('') +
  theme(legend.position='none')
ggplot(all_data, aes(sd, colour=real, fill=real)) +
  geom_density(alpha=0.5) + theme_bw() + xlab('Std. Dv. of ED of shared nodes') + ylab('') +
  theme(legend.title = element_blank())
# SD of Shared Nodes
rl_sds <- all_data[all_data$real == 'Real' & all_data$fssl_nd, 'sd']
rnd_sds <- all_data[all_data$real == 'Random' & all_data$fssl_nd, 'sd']
t.test(rnd_sds, rl_sds)
mean(rnd_sds, na.rm=TRUE)
mean(rl_sds, na.rm=TRUE)
sum(!is.na(rnd_sds))
sum(!is.na(rl_sds))
# save plots
text_tt <- theme(text=element_text(size=6))
# get legend
p_null <- ggplot(all_data, aes(cnt, colour=real, fill=real)) +
  geom_density(alpha=0.5) + theme_bw() + theme(legend.title = element_blank())
tmp <- ggplot_gtable(ggplot_build(p_null)) 
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
legend <- tmp$grobs[[leg]]
tiff(file.path('4_model', 'diff_rand_real.tiff'), width=14, height=14, units="cm",
     res=1200)
grid.arrange(p1+text_tt, p2+text_tt, p3+text_tt, legend, ncol=2)
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

# MODEL SELECTION ----
# obs linear model
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
m2i <- lmer(t1~t0 + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
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
anova(m2g, m2i)
AIC(m2g, m2i)
# m2g is best
obs_mdl_1 <- m2i
save(obs_mdl_1, file = file.path('4_model', 'obs_mdl_1.RData'))

# fitting polynomials
m3a <- lmer(t1~poly(t0, 2) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
anova(obs_mdl_1, m3a)
m3b <- lmer(t1~poly(t0, 3) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
anova(m3a, m3b)
m3c <- lmer(t1~poly(t0, 4) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
anova(m3b, m3c)
m3d <- lmer(t1~poly(t0, 5) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
anova(m3c, m3d)
m3e <- lmer(t1~poly(t0, 6) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
anova(m3d, m3e)
m3f <- lmer(t1~poly(t0, 7) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
anova(m3e, m3f)
# too little improvement, too complex
# m3g <- lmer(t1~poly(t0, 8) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
# anova(m3e, m3g)
# m3h <- lmer(t1~poly(t0, 9) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
# anova(m3g, m3h)
# m3i <- lmer(t1~poly(t0, 10) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
# anova(m3h, m3i)
# m3j <- lmer(t1~poly(t0, 11) + (t0|epoch) + (t0|order/genus), data=genus_data, REML=FALSE)
# anova(m3i, m3j)
# opt for 3b
obs_mdl_2 <- m3b
save(obs_mdl_2, file = file.path('4_model', 'obs_mdl_2.RData'))

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
AIC(n1a, n1b, n1c)
n1d <- lmer(t1~t0_dummy+tm+n+(t0_dummy|genus), data=genus_data, REML=FALSE)
n1e <- lmer(t1~t0_dummy+tm+n+(tm|genus), data=genus_data, REML=FALSE)
n1f <- lmer(t1~t0_dummy+tm+n+(n|genus), data=genus_data, REML=FALSE)
anova(n1b, n1d)
anova(n1e, n1d)
anova(n1f, n1e)
AIC(n1d, n1e, n1f)
# optimisation error
n1g <- lmer(t1~t0_dummy+tm+n+(t0_dummy|genus)+(1|genus), data=genus_data,
            REML=FALSE)
# polynomial
n2a <- lmer(t1~poly(t0_dummy, 2)+tm+n+(t0_dummy|genus),
            data=genus_data, REML=FALSE)
n2b <- lmer(t1~poly(t0_dummy, 3)+tm+n+(t0_dummy|genus),
            data=genus_data, REML=FALSE)
n2c <- lmer(t1~poly(t0_dummy, 4)+tm+n+(t0_dummy|genus),
            data=genus_data, REML=FALSE)
anova(n2a, n2b, n2c)
AIC(n2a, n2b, n2c)
# significant, but minor gains
exp_mdl <- n2b
save(exp_mdl, file=file.path('4_model', 'exp_mdl.RData'))

# OPTED MODELS ----
load(file.path('4_model', 'exp_mdl.RData'))
load(file.path('4_model', 'obs_mdl_1.RData'))
load(file.path('4_model', 'obs_mdl_2.RData'))

# PLOTTING ----
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
# https://www.rdocumentation.org/packages/lme4/versions/1.1-19/topics/bootMer
rpsnttv$pfit <- predict(obs_mdl_2, rpsnttv)
rpsnttv$nfit <- predict(exp_mdl, rpsnttv)
rpsnttv$lfit <- predict(obs_mdl_1, rpsnttv)
# overall
p_data <- plyr::ddply(rpsnttv, c('t0'), plyr::summarise,
                      pmed = median(pfit),
                      nmed = median(nfit),
                      lmed = median(lfit),
                      pupper = quantile(pfit, probs = .95),
                      nupper = quantile(nfit, probs = .95),
                      lupper = quantile(lfit, probs = .95),
                      plower = quantile(pfit, probs = .05),
                      nlower = quantile(nfit, probs = .05),
                      llower = quantile(lfit, probs = .05))
med_ed <- c(p_data$pmed, p_data$lmed, p_data$nmed)
upper_ed <- c(p_data$pupper, p_data$lupper, p_data$nupper)
lower_ed <- c(p_data$plower, p_data$llower, p_data$nlower)
fitted_type <- c(rep('Polynomial', nrow(p_data)),
                 rep('Linear', nrow(p_data)),
                 rep('Expected', nrow(p_data)))
poly_pdata <- data.frame(t0 = p_data$t0, type = 'Polynomial',
                         med_ed = p_data$pmed, upper_ed = p_data$pupper,
                         lower_ed = p_data$plower, expected_ed = p_data$nmed)
line_pdata <- data.frame(t0 = p_data$t0, type = 'Linear',
                         med_ed = p_data$lmed, upper_ed = p_data$lupper,
                         lower_ed = p_data$llower, expected_ed = p_data$nmed)
theme_settings <- theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size = 18),
        title = element_text(size = 11))
p1 <- ggplot(poly_pdata, aes(x = t0, y = med_ed, ymin = lower_ed,
                             ymax = upper_ed)) + 
  geom_line(mapping = aes(x = t0, y = expected_ed)) +
  geom_abline(slope = 1, lty = 3, alpha = 0.75, lwd = 0.75) +
  geom_ribbon(alpha = .25) +
  geom_line(lwd = 2, colour = '#f4663f') +
  xlab(expression('log(ED'['t0']~')')) +
  ylab(expression('log(ED'['t1']~')')) +
  theme_settings + ylim(0.9, 5.4)
p2 <- ggplot(line_pdata, aes(x = t0, y = med_ed, ymin = lower_ed,
                             ymax = upper_ed)) + 
  geom_line(mapping = aes(x = t0, y = expected_ed)) +
  geom_abline(slope = 1, lty = 3, alpha = 0.75, lwd = 0.75) +
  geom_ribbon(alpha = .25) +
  geom_line(lwd = 2, colour = '#f4663f') +
  xlab(expression('log(ED'['t0']~')')) +
  ylab(expression('log(ED'['t1']~')')) +
  theme_settings + ylim(0.9, 5.4)
tiff(file.path('4_model', 'overall.tiff'), width=18, height=9, units="cm",
     res=1200)
grid.arrange(p2 + ggtitle(label = paste0('Best observed linear (m2i)')),
             p1 + ggtitle(label = 'Best observed nonlinear (m3b)') + ylab(''),
             nrow = 1, ncol = 2)
dev.off()
# by epoch
p_data <- plyr::ddply(rpsnttv, c('t0', 'epoch'), plyr::summarise,
                      pfit=median(pfit), nfit=median(nfit))
p <- ggplot(p_data, aes(x=t0, y=pfit)) + 
  geom_line(mapping = aes(x = t0, y = nfit), lwd = 0.75) +
  geom_abline(slope = 1, lty = 3, alpha = 0.75, lwd = 0.25) +
  geom_line(colour = '#f4663f') +
  xlab(expression('log(ED'['t0']~')')) +
  ylab(expression('log(ED'['t1']~')')) +
  theme_bw()
p <- p + facet_wrap(epoch ~ ., ncol = 2)
tiff(file.path('4_model', 'by_epoch.tiff'), width=9, height=14, units="cm",
     res=1200)
print(p + theme(text=element_text(size=6)))
dev.off()


# DIAGNOSTICS
# library(sjPlot)
# library(sjmisc)
# sjPlot::plot_model(exp_mdl, )
# sjPlot::sjplot(exp_mdl)
# sj(exp_mdl, type = "fe")
# sjp.lmer(obs_mdl, type = "fe")
# sjp.lmer(obs_mdl, type = "re.qq")


# COMPARE REAL AND RANDOM POLY 3
# use (1|genus) because (t0|genus) is not computable in real time
rndmdl <- lmer(t1~poly(t0, 3)+(t0|epoch), data=rnd_data, REML=FALSE)
rlmdl <- lmer(t1~poly(t0, 3)+(t0|epoch), data=genus_data, REML=FALSE)
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

# PLOTTING CONFIDENT POINTS
mdl_data$c50 <- '< 50'
mdl_data$c50[mdl_data$cnt > 50] <- '> 50'
mdl_data$c50[mdl_data$cnt < 50] <- '< 50'
p <- ggplot(mdl_data, aes(x=t0, y=t1, colour=c50)) +
  geom_smooth(se=TRUE, formula=y~x) + theme_bw() +
  xlab(expression('ED'['t0'])) +
  ylab(expression('ED'['t1'])) +
  theme(legend.title=element_blank())
tiff(file.path('4_model', 'confident_vs_unconfident.tiff'), width=9, height=9, units="cm",
     res=1200)
print(p + theme(text=element_text(size=6)))
dev.off()
