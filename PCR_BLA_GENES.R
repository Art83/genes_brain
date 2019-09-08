#Second order conditioning BLA/PRC genes 

library(dplyr)
library(ggplot2)
library(outliers)
library(xlsx)

dat <- read.xlsx('camk_four.xlsx', sheetIndex = 1, header = F)
dat <- read.xlsx("camk.xlsx", sheetIndex = 1, header = F)


fold_change <- as.numeric(readClipboard())
groups <- c(rep("2nd order", 3), rep("control", 3))
time <- c(rep(c("30", "90", "270"), 2))
new_dat <- data.frame(groups, time, fold_change)
new_dat$fold_change <- c(mean(dat$X1[1:6]), mean(dat$X1[7:13]), mean(dat$X1[14:21]), mean(dat$X1[22:28]), mean(dat$X1[29:34]), mean(dat$X1[35:41]))
new_dat$time <- factor(new_dat$time, levels = c(30, 90, 270))
new_dat$groups <- factor(new_dat$groups, levels = c("2nd order", "control"))

breaks = c(seq(0.4, 2, by = 0.2), 1)
tiff(file = "camkIV.tiff", width = 4000, height = 3200, units = "px", res = 800)
camkIV_plot <- ggplot(new_dat, aes(x = time, y = fold_change, shape = groups, group = groups, colour = groups)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.4, 2), breaks = breaks, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("CamkIV")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))

camkIV_plot + theme(axis.text=element_text(size=18),
                    legend.title = element_text(size = 18),
                    legend.text = element_text(size = 16),
                    axis.title=element_text(size=18),
                    plot.title=element_text(size=22,hjust=0.5))

dev.off()
# Dnmt3a

f_ch <- as.numeric(readClipboard())
f_ch <- unlist(read.table("dnmt", header = F))
dnmt3a_tab <- new_dat
dnmt3a_tab$fold_change <- f_ch

tiff(file = "dnmt3a.tiff", width = 4000, height = 3200, units = "px", res = 800)
dnmt3a_plot <- ggplot(dnmt3a_tab, aes(x = time, y = fold_change, shape = groups, group = groups, colour = groups)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.4, 2), breaks = breaks, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("Dnmt3a")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))

dnmt3a_plot + theme(axis.text=element_text(size=18),
                    legend.title = element_text(size = 18),
                    legend.text = element_text(size = 16),
                    axis.title=element_text(size=18),
                    plot.title=element_text(size=22,hjust=0.5))
dev.off()


# mapk1
f_ch <- as.numeric(readClipboard())
f_ch <- unlist(read.table("mapk1", header = F))
mapk1_tab <- new_dat
mapk1_tab$fold_change <- f_ch

breaks = c(seq(0.4, 1.4, by = 0.2), 1)
tiff(file = "mapk1.tiff", width = 4000, height = 3200, units = "px", res = 800)
mapk1_plot <- ggplot(mapk1_tab, aes(x = time, y = fold_change, shape = groups, group = groups, colour = groups)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.4, 2), breaks = breaks, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("Mapk1")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))

mapk1_plot+theme(axis.text=element_text(size=18),
                 legend.title = element_text(size = 18),
                 legend.text = element_text(size = 16),
                 axis.title=element_text(size=18),
                 plot.title=element_text(size=22,hjust=0.5))

dev.off()
# Egr-1
f_ch <- as.numeric(readClipboard())
f_ch <- unlist(read.table("egr1", header = F))
egr1_tab <- new_dat
egr1_tab$fold_change <- f_ch

breaks = c(seq(0.4, 2, by = 0.2), 1)
tiff(file = "egr1.tiff", width = 4000, height = 3200, units = 'px', res = 800)
egr1_plot <- ggplot(egr1_tab, aes(x = time, y = fold_change, shape = groups, group = groups, colour = groups)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.4, 2), breaks = breaks, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("Egr1")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))

egr1_plot + theme(axis.text=element_text(size=18),
                  legend.title = element_text(size = 18),
                  legend.text = element_text(size = 16),
                  axis.title=element_text(size=18),
                  plot.title=element_text(size=22,hjust=0.5))
dev.off()

#PRC 
# CamkIV
groups_prc <- c(rep("second-order", 2), rep("spc", 2))
time_prc <- c(rep(c("30", "90"), 2))
f_ch <- c(0.646, 1.125, 2.581, 0.934)
camkIV_prc_tab <- data.frame(groups_prc, time_prc)
camkIV_prc_tab$fold_change <- f_ch

breaks_new = c(seq(0.4, 2.6, by = 0.2), 1)
ggplot(camkIV_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.6, 2.6), breaks = breaks_new, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("CamkIV")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))


# Dnmt3a
f_ch <- c(0.9946457126, 1.1008632222, 1.279597269, 1.0018791032)
dnmt3a_prc_tab <- camkIV_prc_tab
dnmt3a_prc_tab$fold_change <- f_ch

ggplot(dnmt3a_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.8, 1.4), breaks = breaks, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("Dnmt3a")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))


#PRC second-order group vs first-order 
# CamkIV
groups_prc <- c(rep("second-order", 2), rep("first-order", 2))
time_prc <- c(rep(c("30", "90"), 2))
f_ch <- c(0.871096841, 1.124798484,0.713892078,1.073762016)

camkIV_prc_tab <- data.frame(groups_prc, time_prc)
camkIV_prc_tab$fold_change <- f_ch

breaks_new = c(seq(0.4, 2.0, by = 0.2), 1)
ggplot(camkIV_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.6, 2.0), breaks = breaks_new, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("CamkIV PRC")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))

# EGR1
groups_prc <- c(rep("second-order", 2), rep("first-order", 2))
time_prc <- c(rep(c("30", "90"), 2))
f_ch <- c(2.119831912,1.147580525,1.498898043, 1.150104519)

EGR1_prc_tab <- data.frame(groups_prc, time_prc)
EGR1_prc_tab$fold_change <- f_ch

breaks_new = c(seq(0.4, 2.2, by = 0.2), 1)
ggplot(EGR1_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.6, 2.2), breaks = breaks_new, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("EGR1 PRC")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))


# Dnmt3a
groups_prc <- c(rep("second-order", 2), rep("first-order", 2))
time_prc <- c(rep(c("30", "90"), 2))
f_ch <- c(0.994645713, 1.100863222, 0.998911186, 1.304869378)

Dnmt3a_prc_tab <- data.frame(groups_prc, time_prc)
Dnmt3a_prc_tab$fold_change <- f_ch

breaks_new = c(seq(0.4, 2.0, by = 0.2), 1)
ggplot(Dnmt3a_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.6, 2.0), breaks = breaks_new, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("Dnmt3a PRC")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))

# PKCb1
groups_prc <- c(rep("second-order", 2), rep("first-order", 2))
time_prc <- c(rep(c("30", "90"), 2))
f_ch <- c(1.110822658, 1.091477295, 0.860903931, 0.954321414)

PKCb1_prc_tab <- data.frame(groups_prc, time_prc)
PKCb1_prc_tab$fold_change <- f_ch

breaks_new = c(seq(0.4, 2.0, by = 0.2), 1)
ggplot(PKCb1_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.6, 2.0), breaks = breaks_new, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("PKCb1 PRC")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))

# PRC second-order vs sensory preconditioning

# CamkIV
groups_prc <- c(rep("spc", 2), rep("second-order", 2))
time_prc <- c(rep(c("30", "90"), 2))
f_ch <- c(1.63106653, 0.968369468,1.200864517,1.162955522)

CamkIV_prc_tab <- data.frame(groups_prc, time_prc)
CamkIV_prc_tab$fold_change <- f_ch

breaks_new = c(seq(0.4, 2.2, by = 0.2), 1)
ggplot(CamkIV_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.6, 2.2), breaks = breaks_new, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("CamkIIb BLA")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))


#BLA
groups_prc <- c(rep("spc", 2), rep("second-order", 2))
time_prc <- c(rep(c("30", "90"), 2))
f_ch <- c(1.279597269, 1.0055901,0.994645713,1.100863222)

CamkIV_prc_tab <- data.frame(groups_prc, time_prc)
CamkIV_prc_tab$fold_change <- f_ch

breaks_new = c(seq(0.3, 2.0, by = 0.2), 1)
ggplot(CamkIV_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
  geom_line(size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0.3, 2.0), breaks = breaks_new, name = "Fold change")+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
  ggtitle("Dnmt3a PRC")+
  xlab("Timepoint")+
  scale_color_manual(values=c("steelblue3", "red1"))


# for preliminary proccessing

pcr_analysis_30 <- function(so, fo, sp){
  dat <- as.numeric(readClipboard())
  groups <- c(rep("second_order", so), rep("first_order", fo), rep("sensory_preconditioning",sp))
  dat_table <- data.frame(groups, as.numeric(dat))
  colnames(dat_table) <- c("groups", "ddCt")
  p.val_SW <- numeric()
  for (i in 1:length(levels(dat_table$groups))){
    p.val_SW <- c(p.val_SW, shapiro.test(dat_table[dat_table$groups == levels(dat_table$groups)[i], 2])$p.value)
  }
  if (any(p.val_SW <= 0.05)){
    print("Nonnormal distribution")
    so_fo = wilcox.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "first_order", 2])
    so_sp = wilcox.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    fo_sp = wilcox.test(dat_table[dat_table$groups == "first_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    cat("sec_order_vs_fir-order",so_fo$p.value,' ')
    cat("sec_order_vs_spc",so_sp$p.value, ' ')
    cat("fir_order_vs_spc",fo_sp$p.value)
  } else {
    print("Normal distribution ")
    so_fo = t.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "first_order", 2])
    so_sp = t.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    fo_sp = t.test(dat_table[dat_table$groups == "first_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    cat("sec_order_vs_fir-order",so_fo$p.value,' ')
    cat("sec_order_vs_spc",so_sp$p.value, ' ')
    cat("fir_order_vs_spc",fo_sp$p.value)
  }
}
pcr_analysis_30(7,8,6)

pcr_analysis_90 <- function(sp, fo, so){
  dat <- as.numeric(readClipboard())
  groups <- c(rep("sensory_preconditioning", sp), rep("first_order", fo), rep("second_order",so))
  dat_table <- data.frame(groups, as.numeric(dat))
  colnames(dat_table) <- c("groups", "ddCt")
  p.val_SW <- numeric()
  for (i in 1:length(levels(dat_table$groups))){
    p.val_SW <- c(p.val_SW, shapiro.test(dat_table[dat_table$groups == levels(dat_table$groups)[i], 2])$p.value)
  }
  if (any(p.val_SW <= 0.05)){
    print("Nonnormal distribution")
    so_fo = wilcox.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "first_order", 2])
    so_sp = wilcox.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    fo_sp = wilcox.test(dat_table[dat_table$groups == "first_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    cat("sec_order_vs_fir-order",so_fo$p.value,' ')
    cat("sec_order_vs_spc",so_sp$p.value, ' ')
    cat("fir_order_vs_spc",fo_sp$p.value)
  } else {
    print("Normal distribution ")
    so_fo = t.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "first_order", 2])
    so_sp = t.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    fo_sp = t.test(dat_table[dat_table$groups == "first_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
    cat("sec_order_vs_fir-order",so_fo$p.value,' ')
    cat("sec_order_vs_spc",so_sp$p.value, ' ')
    cat("fir_order_vs_spc",fo_sp$p.value)
  }
}

pcr_analysis_90(6,7,7)





dat <- as.numeric(readClipboard())

# Arc
f_ch <- c(1.611018693, 1.178432691, 1.96312974, 2.074641378, 1.244378343, 1.384323151)
groups_prc <- c(rep("second-order", 3), rep("spc", 3))
time_prc <- c(rep(c("30", "90", "270"), 2))
arc_prc_tab <- data.frame(groups_prc, time_prc)
arc_prc_tab$fold_change <- f_ch
arc_prc_tab$time_prc <- factor(arc_prc_tab$time_prc, levels = c(30,90,270))
breaks_new = c(seq(0.4, 2.2, by = 0.2), 1)
ggplot(arc_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
geom_line(size = 1)+
geom_point(size = 3)+
scale_y_continuous(limits = c(0.4, 2.2), breaks = breaks_new, name = "Fold change")+
geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
ggtitle("Arc PRC")+
xlab("Timepoint")+
scale_color_manual(values=c("steelblue3", "red1"))


#CamkIIb
f_ch <- c(1.436186037, , 1.009887408, 1.397267825, 0.751312395)
camkIIb_prc_tab <- data.frame(groups_prc, time_prc)
camkIIb_prc_tab$fold_change <- f_ch
camkIIb_prc_tab$time_prc <- factor(camkIIb_prc_tab$time_prc, levels = c(30,270))
breaks_new = c(seq(0.4, 1.6, by = 0.2), 1)
ggplot(camkIIb_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
geom_line(size = 1)+
geom_point(size = 3)+
scale_y_continuous(limits = c(0.4, 1.6), breaks = breaks_new, name = "Fold change")+
geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
ggtitle("CamkIIb PRC")+
xlab("Timepoint")+
scale_color_manual(values=c("steelblue3", "red1"))


#Prkaca
f_ch <- c(1.375424151, 1.210942874, 1.014587108, 1.572183543, 1.230686256, 0.762353125)
prkaca_prc_tab <- data.frame(groups_prc, time_prc)
prkaca_prc_tab$fold_change <- f_ch
prkaca_prc_tab$time_prc <- factor(prkaca_prc_tab$time_prc, levels = c(30,90,270))

breaks_new = c(seq(0.4, 1.7, by = 0.2), 1)
breaks_new = c(seq(0.4, 1.7, by = 0.2), 1)
ggplot(prkaca_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
geom_line(size = 1)+
geom_point(size = 3)+
scale_y_continuous(limits = c(0.4, 1.7), breaks = breaks_new, name = "Fold change")+
geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
ggtitle("Prkaca PRC")+
xlab("Timepoint")+
scale_color_manual(values=c("steelblue3", "red1"))


#CamkIIb
f_ch <- c(1.436186037, 1.003587673, 1.009887408, 1.397267825, 1.05627492, 0.751312395)
camkIIb_prc_tab <- data.frame(groups_prc, time_prc)
camkIIb_prc_tab$fold_change <- f_ch
camkIIb_prc_tab$time_prc <- factor(camkIIb_prc_tab$time_prc, levels = c(30,90, 270))
breaks_new = c(seq(0.4, 1.6, by = 0.2), 1)
ggplot(camkIIb_prc_tab, aes(x = time_prc, y = fold_change, shape = groups_prc, group = groups_prc, colour = groups_prc)) +
geom_line(size = 1)+
geom_point(size = 3)+
scale_y_continuous(limits = c(0.4, 1.6), breaks = breaks_new, name = "Fold change")+
geom_hline(yintercept = 1, linetype = 'dashed', color = 'black')+
ggtitle("CamkIIb PRC")+
xlab("Timepoint")+
scale_color_manual(values=c("steelblue3", "red1"))


outlier <- as.numeric(readClipboard())
dixon.test(dat)
di.test(dat)

dat <- dat[complete.cases(dat)]

hist(as.numeric(dat), breaks = 20)
plot(rep(1, 8), dat, col="blue")
boxplot(dat)
grubbs.test(dat, opposite = F, type = 20)



groups <- c(rep("second_order", 6), rep("first_order", 5), rep("sensory_preconditioning",7))
groups <- c(rep("sensory_preconditioning",7), rep("first_order", 6), rep("second_order", 8))
groups <- c(rep("first_order", 8), rep("second_order", 7))

groups <- c(rep("second_order", 8), rep("first_order", 8), rep("sensory_preconditioning",8))
groups <- c(rep("second_order", 8),rep("sensory_preconditioning",8), rep("context", 8))
groups <- c(rep("second_order", 7),rep("sensory_preconditioning",7), rep("context", 7))
groups <- c(rep("second_order", 6),rep("first_order",8), rep("context", 7))
groups <- c(rep("second_order", 7), rep("first_order", 7), rep("sensory_preconditioning",7))
groups <- c(rep("second_order", 6), rep("sensory_preconditioning",8),rep("context", 7) )


groups <- c(rep("second_order", 6), rep("first_order", 6), rep("sensory_preconditioning",8),rep("context", 7))
groups <- c(rep("second_order", 6), rep("first_order", 7),rep("context", 8))

dat_table <- data.frame(groups, as.numeric(dat))
colnames(dat_table) <- c("groups", "ddCt")



shapiro.test(dat_table[dat_table$groups == "second_order", 2])
shapiro.test(dat_table[dat_table$groups == "first_order", 2])
shapiro.test(dat_table[dat_table$groups == "sensory_preconditioning", 2])
shapiro.test(dat_table[dat_table$groups == "context", 2])

dat_table <- dat_table[-c(7,11),]
dat_table <- dat_table[-10,]
dat_table <- dat_table[-21,]
dat_table <- dat_table[-c(13,18),]


t.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "first_order", 2],var.equal = T)
t.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
t.test(dat_table[dat_table$groups == "first_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
t.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "context", 2])
t.test(dat_table[dat_table$groups == "first_order", 2], dat_table[dat_table$groups == "context", 2])

wilcox.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "first_order", 2])
wilcox.test(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])
wilcox.test(dat_table[dat_table$groups == "first_order", 2], dat_table[dat_table$groups == "sensory_preconditioning", 2])


pdf(file = "calc.pdf")
print(everybodyFree)
dev.off()




boxplot(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "first_order", 2],dat_table[dat_table$groups == "sensory_preconditioning", 2])
boxplot(dat_table[dat_table$groups == "second_order", 2], dat_table[dat_table$groups == "context", 2],dat_table[dat_table$groups == "sensory_preconditioning", 2])



pairwise.t.test(dat_table$ddCt, dat_table$groups)




library(ggplot2)
ggplot(dat_table, aes(x = groups, y = ddCt, col = groups)) +
  geom_boxplot()


res.aov <- aov(ddCt ~ groups, data = dat_table)
summary(res.aov)
kruskal.test(ddCt ~ groups, data = dat_table)


TukeyHSD(res.aov)


plot(res.aov, 1)
plot(res.aov,2)
library(car)

leveneTest(ddCt ~ groups, data = dat_table)


# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )



library(ggpubr)
compare_means(as.numeric.dat.~ groups,method = "t.test", data = dat_table)



p <- ggboxplot(dat_table, x = "groups", y = "as.numeric.dat.",
               color = "groups", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

# Change method
p + stat_compare_means(method = "anova")

p



pairwise.wilcox.test(dat_table$ddCt, groups)

groups <- factor(dat_table[,1])

dat_table[,1]


t.test(dat_table[dat_table$groups == "second_order",2], dat_table[dat_table$groups == "context",2])
t.test(dat_table[dat_table$groups == "first_order",2], dat_table[dat_table$groups == "context",2])


a <- dat_table[dat_table$groups == "first_order",2]
a[-2]
t.test(a[-2], dat_table[dat_table$groups == "context",2])


#Bayesian stats

fit <- bayes.t.test(dat_table[dat_table$groups == "second_order",2], dat_table[dat_table$groups == "sensory_preconditioning",2])
plot(fit)


bayes.t.test(dat_table[dat_table$groups == "second_order",2], dat_table[dat_table$groups == "sensory_preconditioning",2])
bayes.t.test(dat_table[dat_table$groups == "second_order",2], dat_table[dat_table$groups == "context",2])
bayes.t.test(dat_table[dat_table$groups == "second_order",2], dat_table[dat_table$groups == "first_order",2])


summary(dat_table[dat_table$groups == "second_order", 2])
sd(dat_table[dat_table$groups == "second_order", 2])
summary(dat_table[dat_table$groups == "first_order", 2])
sd(dat_table[dat_table$groups == "first_order", 2])
summary(dat_table[dat_table$groups == "sensory_preconditioning", 2])
sd(dat_table[dat_table$groups == "sensory_preconditioning", 2])
summary(dat_table[dat_table$groups == "context", 2])
sd(dat_table[dat_table$groups == "context", 2])



tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length=mean(len))
ToothGrowth


#Behaviour
# S2
freez_S2_30 <- as.numeric(readClipboard())
freez_S2_SEM_30 <- as.numeric(readClipboard())
# 30 minutes S2
#36
#72
#70
#66

# SEM
#14
#12
#12
#14





freez_S2_90 <- as.numeric(readClipboard())
freez_S2_SEM_90 <- as.numeric(readClipboard())

trials <- c(rep(seq(1:4), 3))
time <- c(rep("30 min", 4), rep("90 min", 4), rep("270 min", 4))
time <- factor(time, levels = c("30 min", "90 min", "270 min"))

freez_S2_270 <- as.numeric(readClipboard())
freez_S2_SEM_270 <- as.numeric(readClipboard())



tab_behaviour_S2 <- data.frame(time, trials, c(freez_S2_30,freez_S2_90, freez_S2_270), c(freez_S2_SEM_30, freez_S2_SEM_90, freez_S2_SEM_270))
colnames(tab_behaviour_S2) <- c("Timepoint", "Trial", "Freezing", "SEM")


tiff(file = "s2.tiff", width = 4000, height = 3200, units = "px", res = 800)
behav_plot_S2 <- ggplot(tab_behaviour_S2, aes(x = Trial, y = Freezing, group = Timepoint,colour = Timepoint))+
         geom_errorbar(aes(ymin=Freezing-SEM, ymax=Freezing+SEM), width =.1, size=1, position = position_dodge(0.1))+
         geom_line(size=1, position = position_dodge(0.1))+
         geom_point(size=2,position = position_dodge(0.1))+
         geom_point(aes(x = 0.5, y = 14), color = "#00BA38", size = 2)+
         geom_errorbar(aes(x=0.5, ymin=14-5, ymax=14+5), width = .1, size = 1, color = "#00BA38")+
         geom_point(aes(x = 0.5, y = 1), color = "#F8766D", size = 2)+
         geom_errorbar(aes(x=0.5, ymin=1-1, ymax=1+1), width = .1, size = 1, color = "#F8766D")+
         geom_point(aes(x = 0.6, y = 2), color = "#619CFF", size = 2)+
         geom_errorbar(aes(x=0.6, ymin=2-1, ymax=2+1), width = .1, size = 1, color = "#619CFF")+
         scale_y_continuous(limits=c(0,100))+
         ylab("Freezing (%)")+
         ggtitle("Stage 2")
       
behav_plot_S2 + theme(axis.text=element_text(size=18),
                   legend.title = element_text(size = 18),
                   legend.text = element_text(size = 16),
                   axis.title=element_text(size=18),
                   plot.title=element_text(size=22,hjust=0.5))
dev.off()




#1
freez_S1_30 <- as.numeric(readClipboard())
freez_S1_SEM_30 <- as.numeric(readClipboard())

freez_S1_90 <- as.numeric(readClipboard())
freez_S1_SEM_90 <- as.numeric(readClipboard())

freez_S1_270 <- as.numeric(readClipboard())
freez_S1_SEM_270 <- as.numeric(readClipboard())


tab_behaviour_S1 <- data.frame(time, trials, c(freez_S1_30,freez_S1_90, freez_S1_270), c(freez_S1_SEM_30, freez_S1_SEM_90, freez_S1_SEM_270))
colnames(tab_behaviour_S1) <- c("Timepoint", "Trial", "Freezing", "SEM")



tiff(file = "s1.tiff", width = 5000, height = 3200, units = "px", res = 800)
behav_plot_S1 <- ggplot(tab_behaviour_S1, aes(x = Trial, y = Freezing, group = Timepoint,colour = Timepoint))+
  geom_errorbar(aes(ymin=Freezing-SEM, ymax=Freezing+SEM), width =.1, size=1, position = position_dodge(0.1))+
  geom_line(size=1, position = position_dodge(0.1))+
  geom_point(size=2,position = position_dodge(0.1))+
  scale_y_continuous(limits=c(0,102))+
  ylab("Freezing (%)")+
  ggtitle("Stage 2 (Light)")

behav_plot_S1 + theme(axis.text=element_text(size=18),
                      legend.title = element_text(size = 18),
                      legend.text = element_text(size = 16),
                      axis.title=element_text(size=18),
                      plot.title=element_text(size=22,hjust=0.5))

dev.off()



#619CFF
#00BA38
#F8766D


tiff("test.tiff", units="in", width=5, height=5, res=300)
behav_plot_S1 <- ggplot(tab_behaviour_S1, aes(x = Trial, y = Freezing, group = Timepoint,colour = Timepoint))+
  geom_errorbar(aes(ymin=Freezing-SEM, ymax=Freezing+SEM), width =.1, size=1, position = position_dodge(0.1))+
  geom_line(size=1, position = position_dodge(0.1))+
  geom_point(size=2,position = position_dodge(0.1))+
  scale_y_continuous(limits=c(0,100))+
  ylab("Freezing,%")+
  ggtitle("S1 test (light)")

behav_plot_S1 + theme(axis.text=element_text(size=14),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 14),
                      axis.title=element_text(size=14),
                      plot.title=element_text(size=18,hjust=0.5))
dev.off()



for (i in 1:length(levels(dat_table$groups))){
  cat(unlist(levels(dat_table$groups)[i]))
}




