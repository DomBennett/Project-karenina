# D.J. Bennett
# 18/07/2015
# Model EDt=1 ~ EDt=0

# LIBS
library(ggplot2)

# DIRS
input_dir <- '3_wrngl'
output_dir <- '4_model'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

# INPUT
load(file.path(input_dir, 'data.RData'))






ggplot(mdl_data, aes(x=t0, y=t1, colour=tmsplt)) +
  geom_smooth(method='lm', se=TRUE, formula=y~x)

mdl_data <- mdl_data[!mdl_data$tmsplt %in% c("122.75-83.25", "154.25-122.75"), ]


library(lme4)

# fit t0 and tmpslt
m1 <- lmer(t1~1 + (1|tmsplt), data=mdl_data, REML=FALSE)
m2 <- lmer(t1~t0 + (1|tmsplt), data=mdl_data, REML=FALSE)
m3 <- lmer(t1~t0 + (t0|tmsplt), data=mdl_data, REML=FALSE)
anova(m1, m2, m3)
anova(m2, m3)  # definitely need random slopes
AIC(m2, m3)
# fit ^2
m4 <- lmer(t1~I(t0^2)+(t0|tmsplt), data=mdl_data, REML=FALSE)
m5 <- lmer(t1~t0+I(t0^2)+(t0|tmsplt), data=mdl_data, REML=FALSE)
anova(m3, m4)
anova(m4, m5)
# fit ^3
m6 <- lmer(t1~I(t0^3)+(t0|tmsplt), data=mdl_data, REML=FALSE)
m7 <- lmer(t1~t0+I(t0^3)+(t0|tmsplt), data=mdl_data, REML=FALSE)
anova(m3, m6)
anova(m6, m7)
anova(m5, m7) # third polynomial is worse fit than second
# fit polys
m8 <- lmer(t1~poly(t0, 2)+(t0|tmsplt), data=mdl_data, REML=FALSE)
m9 <- lmer(t1~poly(t0, 3)+(t0|tmsplt), data=mdl_data, REML=FALSE)
m10 <- lmer(t1~poly(t0, 4)+(t0|tmsplt), data=mdl_data, REML=FALSE)
m11 <- lmer(t1~poly(t0, 5)+(t0|tmsplt), data=mdl_data, REML=FALSE)
anova(m8, m9, m10, m11)
# compare to non-orthos
anova(m5, m8)
anova(m7, m9)
# m5 is best with t0, m9 is best poly
anova(m5, m9)
anova(m9, m10)

# Experimental
m12 <- lmer(t1~I(t0^.5)+(t0|tmsplt), data=mdl_data, REML=FALSE)
anova(m12, m10)
anova(m12, m5)
anova(m9, m10, m12)
AIC(m9, m10)

# PLOT
mdl_data$fit <- predict(m9)
m3
ggplot(mdl_data, aes(x=t0, y=t1, group=tmsplt, col=tmsplt)) + 
  geom_abline(intercept=1.4225, slope=0.6537, lty=3) +
  geom_point(alpha=0.01) +
  geom_line(aes(y=fit), lwd=2) +
  theme_bw()

ggplot(mdl_data, aes(x=t0, y=t1)) +
  geom_point(aes(colour=tmsplt), alpha=.01) +
  stat_smooth(method='lm', formula=y~x, se=FALSE,
              colour='black', lty=3) +
  stat_smooth(method='lm', formula=y~poly(x, 3), se=FALSE,
              colour='black', lwd=2) +
  theme_bw()
  


# fit t0 and tmpslt
m1 <- lm(t1~1, data=mdl_data)
m2 <- lm(t1~1+factor(tmsplt)-1, data=mdl_data)
m3 <- lm(t1~t0+factor(tmsplt)-1, data=mdl_data)
anova(m1, m2, m3)  # adding tmsplt and t0 massively improves model
# fit ^2
m4 <- lm(t1~I(t0^2)+factor(tmsplt)-1, data=mdl_data)
m5 <- lm(t1~t0+I(t0^2)+factor(tmsplt)-1, data=mdl_data)
m6 <- lm(t1~poly(t0, 2)+factor(tmsplt)-1, data=mdl_data)
m7 <- lm(t1~t0+poly(t0, 2)+factor(tmsplt)-1, data=mdl_data)
anova(m4, m5, m6, m7) # raw ^2 with t0 is best fit, little improvement with complexity
# fit ^3
m8 <- lm(t1~I(t0^3)+factor(tmsplt)-1, data=mdl_data)
m9 <- lm(t1~t0+I(t0^3)+factor(tmsplt)-1, data=mdl_data)
m10 <- lm(t1~poly(t0, 3)+factor(tmsplt)-1, data=mdl_data)
m11 <- lm(t1~t0+poly(t0, 3)+factor(tmsplt)-1, data=mdl_data)
anova(m8, m9, m10, m11) # slightly better fit for orthogonal poly 3
anova(m5, m9, m10)  # marginally better fit for orthogonal poly 3
anova(m5, m9)  # no difference between poly 2 and 3
# fit ^4
m12 <- lm(t1~I(t0^4)+factor(tmsplt)-1, data=mdl_data)
m13 <- lm(t1~t0+I(t0^4)+factor(tmsplt)-1, data=mdl_data)
m14 <- lm(t1~poly(t0, 4)+factor(tmsplt)-1, data=mdl_data)
m15 <- lm(t1~t0+poly(t0, 4)+factor(tmsplt)-1, data=mdl_data)
anova(m12, m13, m14, m15) # best fit for orthogonal poly 4
anova(m5, m9, m10, m13)  # orthogonal out-perform
anova(m5, m9, m13) # no difference between each of the polys
# which is best?
AIC(m1, m2, m3, m5, m9, m10, m13)  # M5 is best performing non-orthogonal
anova(m10, m14)

# PLOT
# M5
ggplot(mdl_data, aes(x=t0, y=t1, colour=factor(tmsplt))) +
  geom_point(alpha=.5) +
  stat_smooth(method='lm', se=TRUE, formula=y~x,
              colour='blue', fill='blue') +
  stat_smooth(method='lm', se=TRUE, formula=y~I(x^2),
              colour='red', fill='red')
# M10
ggplot(mdl_data, aes(x=t0, y=t1, colour=factor(tmsplt))) +
  geom_point(alpha=.5) +
  stat_smooth(method='lm', se=TRUE, formula=y~poly(x, 3),
              colour='red', fill='red')
# M14
ggplot(mdl_data, aes(x=t0, y=t1, colour=factor(tmsplt))) +
  geom_point(alpha=.5) +
  stat_smooth(method='lm', se=TRUE, formula=y~poly(x, 4),
              colour='red', fill='red')

# OPTING FOR M5:
# - we want to find the model on top of the linear regression with t0, only m5 has the linear regression and a polynomial
# - m5 is simplest
# - there is a step change in model improvement up to m5, additional parameters lead to only marginal AIC drop