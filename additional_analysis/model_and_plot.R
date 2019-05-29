# Libs ----
library(ggplot2)
library(gridExtra)
library(ghibli)

# Functions ----
predict_t1 <- function(p_data, pull, mdl) {
  predicted <- predict(mdl, p_data[pull,], se.fit = TRUE, level = 0.95)
  p_data$t1[pull] <- predicted$fit
  p_data$t1_upper[pull] <- predicted$fit + predicted$se.fit
  p_data$t1_lower[pull] <- predicted$fit - predicted$se.fit
  p_data
}

# Read ----
load(file = file.path('additional_analysis', 'sim_results.RData'))

# All plot ----
all_data <- rbind(bd_data, pan_data, de_data, pf_data)
all_data$type <- c(rep('Null', nrow(bd_data)),
                   rep('Pan.', nrow(pan_data)),
                   rep('Rel.', nrow(de_data)),
                   rep('P.F.', nrow(pf_data)))
# control level order
all_data$type <- factor(all_data$type,
                        levels = c("Null", "Pan.", "Rel.", 'P.F.'))
all_data$id <- paste0(all_data$type, '_', all_data$id)
# drop multiples
#all_data <- all_data[duplicated(all_data$id), ]
# individual
p1 <- ggplot(all_data, aes(x = t0, y = t1)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"),
              aes(colour = type), se = TRUE, level = 0.95) +
  geom_abline(slope = 1) +
  facet_wrap(~type, ncol = 2) +
  theme_bw() +
  scale_colour_ghibli_d("MarnieMedium1") +
  theme(legend.position = 'none') +
  xlab(expression('log(ED'['t0']~')')) +
  ylab(expression('log(ED'['t1']~')'))
# linear
p2 <- ggplot(all_data, aes(x = t0, y = t1, colour = type)) +
  geom_smooth() +
  geom_abline(slope = 1)
tiff(file.path('additional_analysis', 'simulations_points.tiff'), width = 18,
     height = 18, units = "cm", res = 1200)
print(p1 + theme(legend.title = element_blank(), text = element_text(size = 18),
                 title = element_text(size = 11)))
dev.off()

# Linear model ---
# bd
m1 <- lm(t1 ~ t0, data = bd_data)
m2 <- lm(t1 ~ poly(t0, 2), data = bd_data)
m3 <- lm(t1 ~ poly(t0, 3), data = bd_data)
anova(m1, m2, m3)
bd_mdl <- m2
# pan
m1 <- lm(t1 ~ t0, data = pan_data)
m2 <- lm(t1 ~ poly(t0, 2), data = pan_data)
m3 <- lm(t1 ~ poly(t0, 3), data = pan_data)
m4 <- lm(t1 ~ poly(t0, 4), data = pan_data)
anova(m1, m2, m3, m4)
pan_mdl <- m3
# de
m1 <- lm(t1 ~ t0, data = de_data)
m2 <- lm(t1 ~ poly(t0, 2), data = de_data)
m3 <- lm(t1 ~ poly(t0, 3), data = de_data)
m4 <- lm(t1 ~ poly(t0, 4), data = de_data)
anova(m1, m2, m3, m4)
de_mdl <- m3
# pf
m1 <- lm(t1 ~ t0, data = pf_data)
m2 <- lm(t1 ~ poly(t0, 2), data = pf_data)
m3 <- lm(t1 ~ poly(t0, 3), data = pf_data)
m4 <- lm(t1 ~ poly(t0, 4), data = pf_data)
anova(m1, m2, m3, m4)
pf_mdl <- m3

# Build plot data ----
# all based on bd_data
t0 <- bd_data$t0
p_data <- data.frame(t0 = rep(t0, 4), type = c(rep('Null', length(t0)),
                                               rep('Pan.', length(t0)),
                                               rep('Rel.', length(t0)),
                                               rep('P.F.', length(t0))),
                     t1 = NA, t1_upper = NA, t1_lower = NA)
# control level order
p_data$type <- factor(p_data$type, levels = c("Null", "Pan.", "Rel.", 'P.F.'))
# predict
p_data <- predict_t1(p_data = p_data, pull = p_data$type == 'Null',
                     mdl = bd_mdl)
p_data <- predict_t1(p_data = p_data, pull = p_data$type == 'Pan.',
                     mdl = pan_mdl)
p_data <- predict_t1(p_data = p_data, pull = p_data$type == 'Rel.',
                     mdl = de_mdl)
p_data <- predict_t1(p_data = p_data, pull = p_data$type == 'P.F.',
                     mdl = pf_mdl)

# Plot ----
p <- ggplot(p_data, aes(x = t0, y = t1, ymin = t1_lower, ymax = t1_upper)) +
  geom_abline(slope = 1, linetype = 2, lwd = 1) +
  geom_line(lwd = 1.5, aes(colour = type)) + facet_wrap(~ type, ncol = 2) +
  geom_ribbon(alpha = .1, aes(fill = type)) + theme_bw() +
  scale_colour_ghibli_d("MarnieMedium1") + theme(legend.position = 'none') +
  xlab(expression('log(ED'['t0']~')')) +
  ylab(expression('log(ED'['t1']~')'))
tiff(file.path('additional_analysis', 'simulations_models.tiff'), width = 18,
     height = 18, units = "cm", res = 1200)
print(p + theme(legend.title = element_blank(), text = element_text(size = 18),
                title = element_text(size = 11)))
dev.off()

# Dev ----
# scenario_plot <- function(p_data, scn_mdl, bd_mdl, colour) {
#   p_data$bd_fit <- predict(bd_mdl, p_data)
#   p_data$scn_fit <- predict(scn_mdl, p_data)
#   ggplot(p_data, aes(x = t0, y = t1)) +
#     geom_abline(slope = 1, linetype = 2) +
#     geom_line(mapping = aes(x = t0, y = scn_fit), colour = colour, lwd = 1.1) +
#     geom_line(mapping = aes(x = t0, y = bd_fit), colour = '#4286f4') + 
#     xlab('') + ylab('') + theme_bw()
# }
# 
# p1 <- ggplot(bd_data, aes(x = t0, y = t1)) + 
#   #geom_point(alpha = 0.1) +
#   geom_abline(slope = 1, linetype = 2) +
#   stat_smooth(formula = y ~ poly(x, 2), method = 'lm', se = FALSE,
#               colour = '#4286f4') +
#   xlab('') + ylab('') + theme_bw()
# p2 <- scenario_plot(p_data = bd_data, scn_mdl = pan_mdl, bd_mdl = bd_mdl,
#                     colour = '#990000')
# p3 <- scenario_plot(p_data = bd_data, scn_mdl = de_mdl, bd_mdl = bd_mdl,
#                     colour = '#38761d')
# p4 <- scenario_plot(p_data = bd_data, scn_mdl = pf_mdl, bd_mdl = bd_mdl,
#                     colour = '#741b47')
# gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
# 


