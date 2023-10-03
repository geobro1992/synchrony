# required libraries
library(lme4)
library(lubridate)
library(vegan)
library(rstatix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggeffects)
library(gridExtra)
library(DescTools)
library(report)
library(circular)
library(bpnreg)

# read in data
bc = read.csv("df_dat.csv")

bc = bc %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

# filter most abundant amphibians
spp = c("AMBBIS", "EURQUA", "LITUTR", "ACRGRY", "PSEORN", "GASCAR", "BUFTER")

bc = bc %>%
  filter((Species %in% spp)) %>%
  mutate(Species = factor(Species, levels = c("EURQUA", "AMBBIS", 
                                              "LITUTR", "PSEORN", "GASCAR", "BUFTER", "ACRGRY"))) %>%
  droplevels()

levels(bc$Species) = c("E. quadridigitata", "A. bishopi", "L. sphenocephalus", "P. ornata", "G. carolinensis", "A. terrestris", "A. gryllus")


pdf(file = "fig1.pdf", width = 8, height = 8)

ggplot(bc, aes(Species, Date)) +
  geom_sina(aes(color = Species), size = 0.7) +
  scale_y_continuous(limits = c(15200, 17300), breaks = c(15340, 15706, 16071, 16436, 16801, 17167), labels = c("2011-2012", "2012-2013", "2013-2014", "2014-2015", "2015-2016", "2016-2017")) +
  coord_flip() +
  labs(x = "") +
  scale_color_brewer(palette="Dark2") +
  theme_pubclean()+ 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_hline(yintercept = 2) + theme(axis.text.y = element_text(face = "italic")) 

dev.off()


#bc$mon = month(bc$Date, label = T, abbr = T)

# calculate abundance for each species in each month in each year 
#mon.counts = aggregate(Date ~ Season + mon + Species, bc, length)
#mons = month.abb

#mon.counts = mon.counts %>%
#  complete(mon = month.abb, nesting(Season, Species), fill = list(Date = 0)) 

#mon.counts$mon = factor(mon.counts$mon, levels = month.abb)

#ggplot(mon.counts, aes(x = mon, y = Date, fill = Species)) +
#  geom_col(position = "dodge") +
#  coord_polar() + 
#  ylim(-100, 450) +
#  facet_wrap(~Season) +
# theme_pubclean() +
#  scale_fill_brewer(palette="Dark2")

bc$degrees <- as.numeric(strftime(bc$Date, format = "%j"))/366*360 # convert year to 360 degrees

bc$circ <- circular(bc$degrees, units="degrees", template="geographics")

# Plot the results
plot(bc$circ, stack=TRUE, axes=FALSE)
axis <- circular(c(0, 90, 180, 270), units="degrees", template="geographics")
axis.circular(at=axis, labels=c("Jan", "April", "July", "October"),
              template="geographics", zero=0, tcl.text=.15)
arrows.circular(quantile(bc$circ, prob=c(.10, .90)), lwd=2)

# circular regression
#preds = cbind(factor(bc$Season-2010), factor(bc$Species))
#mcl <- lm.circular(y = bc$circ, x = preds, init=c(0,0), type="c-l", tol = 1e-4, verbose=T)
#mcl

# perform regression comparison
m1 = bpnr(circ ~ factor(Season) * factor(Species), data = bc, its = 1000, burn = 500)
m2 = bpnr(circ ~ factor(Season), data = bc, its = 1000, burn = 500)
m3 = bpnr(circ ~ factor(Species), data = bc, its = 1000, burn = 500)
m4 = bpnr(circ ~ 1, data = bc, its = 1000, burn = 500)
m5 = bpnr(circ ~ factor(Season) + factor(Species), data = bc, its = 1000, burn = 500)


# compare model fit
fit(m1) # full model is the best
fit(m2)
fit(m3)
fit(m4)
fit(m5)

# extract model coefficients
m1
preds = m1$circ.coef.cat


# calculate abundance for each species in each month in each year 
day.counts = aggregate(Date ~ Species + Season + degrees, bc, length)

g1 = ggplot(day.counts, aes(x = degrees, y = log(Date), fill = Species)) +
  geom_col(position = "dodge") +
  xlab("") + ylab("") +
  coord_polar() + 
  ylim(-5, 5) +
  facet_wrap(~Season) +
  theme_pubclean() +
  scale_fill_brewer(palette="Dark2") +
  labs(caption = "Figure S1. Radial plot of species arrival times over the six years of drift-fence monitoring. \nDay of year has been converted to degrees such that January 1st equals 1 and December 31st equals 360. \nHeights of bars are scaled to the number of individuals captured per day.") +
  theme(plot.caption = element_text(hjust=0, size=rel(1.4)))


ggsave("sal_circles.pdf", g1, device = cairo_pdf, width = 16, height = 12)


###################################################
# calculate abundance for each species in each year 
abs = aggregate(Date ~ Season + Species, bc, length)

abs = abs %>%
  spread(Species, Date)

abundances = as.matrix(abs[,-1])

# scale by sampling days
samp.days = c(150, 135, 255, 195, 165, 100)
BCI = abundances / samp.days
rowSums(abundances)/samp.days
# convert abundances to biomass using mass estimates from literature
mass = rev(c(ACRGRY = 0.45, BUFTER = 19, GASCAR = 1, PSEORN = 5, LITUTR = 22, AMBBIS = 6, EURQUA = 0.75))
BCI = t(t(BCI) * mass)
biomass = rowSums(BCI)

# dominance proportions
BCI/biomass


##################
# Evenness metrics
##################

abundances = abundances[,c(2,5,7,10,13,14,29,33,36,37,38,39,42,44,51,52,53,54,55,56,59,62)]

## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(abundances, 2) - 1
## Fisher alpha
alpha <- fisher.alpha(abundances)

# Simpson's diversity index 
simp <- diversity(BCI, "simpson")
invsimp <- diversity(BCI, "inv")

## Species richness (S), Sannon's diversity index (H), and Pielou's evenness (J):
BCI[is.na(BCI)] = 0
S <- specnumber(BCI)
H <- diversity(BCI)
J <- H/log(S)



######################################
# plot biomass and evenness across time

add.label <- function(label, ...) legend("topleft", legend=" ", title=label,
                                         bty='n', ...)


pdf(file = "fig3.pdf", width = 9, height = 8)

par(mfrow = c(2,1), cex.main = 1.5, mar = c(4, 6, 1, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1, bty = "n", las = 1)

digitsize <- 1.2
x <- c(1:6)
seas = c("2011-2012", "2012-2013", "2013-2014", "2014-2015", "2015-2016", "2016-2017")

plot(x, rep(-10000, 6), type = "p", ylab = " ", 
     xlab = " ", cex = 1.5, ylim = c(0, 140), xlim = c(1, 7), lwd = 2, pch = 5, 
     axes = F, main = " ")

axis(1, at = seq(1.5, 6.5, by = 1), labels = as.character(seas))
mtext(expression(paste("Season")), side = 1, line = 3, cex = 1.5, font = 2)
axis(2, at = c(0, 20, 40, 60, 80, 100, 120, 140))
mtext(expression(paste("Biomass (g ", day^-1,")")), side = 2, line = 4, cex = 1.5, font = 2, las = 0)

x = seq(1.5, 6.5, by = 1)
points(x, biomass, cex = 1.5, lwd = 2, pch = 19)
lines(seq(1.5, 6.5, by = 1), biomass, lwd = 2, type = "c")
add.label(label = "A", cex = 2)

# evenness
x <- c(1:6)
plot(x, rep(-10000, 6), type = "p", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(0, 1.1), xlim = c(1, 7), lwd = 2, axes = FALSE, main = " ")
axis(1, at = seq(1.5, 6.5, by = 1), labels = as.character(seas))
mtext(expression(paste("Season")), side = 1, line = 3, cex = 1.5, font = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), las = 1)
mtext(expression(paste("Pielou's Evenness")), side = 2, line = 4, cex = 1.5, font = 2, las = 0)

x = seq(1.5, 6.5, by = 1)
points(x, J, cex = 1.5, lwd = 2, pch = 8)
lines(seq(1.5, 6.5, by = 1), J, lwd = 2, lty = "dashed")
add.label(label = "B", cex = 2)

dev.off()


#############################################################
# read in dates for the 75%  threshold of population movement
dat = read.csv("75_arrivals.csv")
#dat = read.csv("50_arrivals.csv") # same result with 50% threshold
#dat = read.csv("90_arrivals.csv") # same result with 90% threshold

dat = dat %>%
  mutate(Species = factor(Species, levels = c("EURQUA", "AMBBIS", 
                                              "LITUTR", "PSEORN", "GASCAR", "BUFTER", "ACRGRY"))) %>%
  droplevels()

levels(dat$Species) = c("E. quadridigitata", "A. bishopi", "R. sphenocephalus", "P. ornata", "G. carolinensis", "A. terrestris", "A. gryllus")
# extract Julian day
dat$day.75 = as.numeric(format(as.Date(dat$Date.75)-180, "%j"))

# friedman test (non-parametric equivalent of repeated measure anova)
res.fried <- dat %>% friedman_test(day.75 ~ Season | Species)
res.fried

dat %>% friedman_effsize(day.75 ~ Season | Species)

# pairwise comparisons
ConoverTest(day.75 ~ Season, data = dat, method = "fdr")

# plot arrival dates
bxp = ggboxplot(dat, x = "Season", y = "day.75") +
  geom_point(aes(color = Species)) +
  ylab("") +
  labs(subtitle = get_test_label(res.fried,  detailed = TRUE)) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_x_discrete(breaks = c("2", "3", "4", "5", "6", "7"), labels = c("2011-2012", "2012-2013", "2013-2014", "2014-2015", "2015-2016", "2016-2017")) +
  scale_y_continuous(breaks = c(125, 155, 186, 217, 245, 276, 306), labels = c("Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1", "Apr 1", "May 1")) +
  color_palette("Dark2") + guides(colour = guide_legend(nrow = 1, reverse = T))

ggsave(filename = "fig2.pdf", bxp, device = cairo_pdf, width = 8, height = 6, units = "in", dpi = 600) 


# Libraries
library(GGally)

# Data set is provided by R natively
data = spread(dat[,-c(1,4)], key = Season, value = day.75)

# Plot
p1 = ggparcoord(data, columns = 2:7, groupColumn = 1, 
                scale="globalminmax", showPoints = TRUE, alphaLines = 0.7)

p1 = p1 + xlab("Season") + ylab("") + 
  labs(subtitle = get_test_label(res.fried,  detailed = TRUE)) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_x_discrete(breaks = c("2", "3", "4", "5", "6", "7"), labels = c("2011-2012", "2012-2013", "2013-2014", "2014-2015", "2015-2016", "2016-2017")) +
  scale_y_continuous(breaks = c(125, 155, 186, 217, 245, 276, 306), labels = c("Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1", "Apr 1", "May 1")) +
  color_palette("Dark2") + 
  guides(colour = guide_legend(reverse = T)) + 
  theme_Publication() +
  theme(legend.position=c(.15,.85))

ggsave(filename = "fig2.pdf", p1, device = cairo_pdf, width = 8, height = 6, units = "in", dpi = 600) 

####################
# hydroperiod analysis
####################

# ggplot publication theme
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(0.1,0.1,0.1,0.1, "cm"),
            legend.title = element_blank(),
            legend.text=element_text(size=12, face = "italic"),
            plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#######################
# hydroperiods

hy = c(7, 33, 230, 57, 252, 127)

#######################
# biomass ~ hydroperiod 

# merge hydroperiods with biomass data
ab.dat = data.frame(hy, BCI)
ab.dat = ab.dat[order(ab.dat$hy),]
ab.dat = data.frame(hy = rep(ab.dat$hy, 7), gather(ab.dat[,2:8]))

ab.dat = ab.dat %>%
  mutate(key = factor(key, levels = c("E..quadridigitata", "A..bishopi", 
                                      "L..sphenocephalus", "P..ornata", "G..carolinensis", "A..terrestris", "A..gryllus")))

levels(ab.dat$key) = c("E. quadridigitata", "A. bishopi", "R. sphenocephalus", "P. ornata", "G. carolinensis", "A. terrestris", "A. gryllus")


mod1 = lmer(log(value) ~ hy + (1|key), data = ab.dat)

p1 = ggpredict(mod1, terms = c("hy", "key"), type = "re", ci.lvl = 0) 

g1 = plot(p1) +
  theme_Publication() + 
  xlab("") + 
  ggtitle("") +
  ylab(expression(paste("Biomass (g ", day^-1,")"))) +
  scale_color_brewer(palette="Dark2")+ guides(colour = guide_legend(reverse = T)) +
  theme(legend.position=c(.2,.85))


#########################
# intraspecific synchrony
dat = readRDS("Movement_Thresholds.rds")

dat$day.90 = as.numeric(format(dat$Date.90-180, "%j"))
dat$day.75 = as.numeric(format(dat$Date.75-180, "%j"))
dat$day.50 = as.numeric(format(dat$Date.50-180, "%j"))
dat$delta.day = (dat$day.90 - dat$day.50) + 1

hy.dat = data.frame(hy, Season = 1:6)
dat.delta = merge(dat, hy.dat, by = "Season")

dat.delta = dat.delta %>%
  mutate(Species = factor(Species, levels = c("EURQUA", "AMBBIS", 
                                              "RANUTR", "PSEORN", "GASCAR", "BUFTER", "ACRGRY")))

levels(dat.delta$Species) = c("E. quadridigitata", "A. bishopi", "R. sphenocephalus", "P. ornata", "G. carolinensis", "A. terrestris", "A. gryllus")


mod2 = lmer(log(delta.day) ~ log(hy) + (1|Species), data = dat.delta)

p2 = ggpredict(mod2, terms = c("hy", "Species"), type = "re", ci.lvl = 0) 

g2 = plot(p2) +
  theme_Publication() + 
  xlab("Hydroperiod (days)") + 
  ggtitle("") +
  ylab(expression(paste(Delta, " days (intraspecific)"))) +
  scale_color_brewer(palette="Dark2")+ guides(colour = guide_legend(reverse = T)) +
  theme(legend.position="none")


gg = ggarrange(g1, g2, 
              labels = c("A", "B"),
              ncol = 1, nrow = 2,
              common.legend = F)


ggsave(filename = "ggfig4.pdf", plot = gg, device = cairo_pdf, width = 6, height = 8, units = "in")

#########################
# interspecific synchrony
int = aggregate(day.90 ~ Season, data = dat, FUN = function(x) max(x)-min(x))
int.dat = merge(int, hy.dat, by = "Season")

mod3 = glm(log(day.90) ~ hy, data = int.dat)
summary(mod3) # differences in arrival times between species is not related to hydroperiod

