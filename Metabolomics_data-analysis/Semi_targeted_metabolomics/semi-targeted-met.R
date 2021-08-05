library(tibble)
library(dplyr)

# Load data 
# For calculating t0 variation use file = GC_MS_table_resT0.txt

data1 <- read.table("GC_MS_table_res.txt", check.names = FALSE, h = T)

# Normalize to internal standard

Internal_std <- data1$L_Norvaline

data <- sweep(data1[,3:112], 1, Internal_std, FUN = '/')

# Finish df
data <- cbind(Group = data1$Group, data)

data <- cbind(Sample = data1$Sample, data)

# Calculate log2FC, produce volcano plots ----
# This is done per time point
# t0 

T0 <-c("183_T0","184_T0","185_T0","186_T0")

T1 <-c("183_T1","184_T1","185_T1","186_T1")

# The following function will perform a t-test and will originate 1 file per strain
# with log2FC and SE results as well as p-value adjusted. 

t_test <- function(data, group1, group2) {
  group1_data <- subset(data, data$Group==group1)
  group2_data <- subset(data, data$Group==group2)
  sample_data <- rbind(group1_data,group2_data)
  sample_data <- droplevels(sample_data)
  nb_ions <- length(colnames(sample_data))
  t.test.list <- lapply(sample_data[,c(3:nb_ions)],function(x) t.test(x~sample_data$Group))
  group1_log2 <- log2(group1_data[,c(3:nb_ions)])
  group2_log2 <- log2(group2_data[,c(3:nb_ions)])
  group1_sd <- apply(group1_log2,2,function(x){sd(x)})
  group2_sd <-  apply(group2_log2,2,function(x){sd(x)})
  final_se <- sqrt((group1_sd^2) + (group2_sd^2))
  t.test.list <- list(t.test.list,final_se)
  return(t.test.list)
}

nb_samples <- length(T0)
index <- c(1:nb_samples)
pdf("volcano_plot0.pdf")
par(mfrow=c(2,2))
for (i in index) {
  T0_sample <- T0[i]
  T0_sample <- T0[i]
  strain <- substr(T0_sample,1,6)
  table_name <- paste0(strain,"_FC.txt")
  t.test <- t_test(data, T0_sample, T0_sample)
  t.test.list <- t.test[[1]]
  SE <- t.test[[2]]
  mean.values.list <- lapply(t.test.list, function(x) x$estimate)
  mean.values.T0 <- unlist(lapply(mean.values.list, function(x) x[1]))
  mean.values.T0 <- unlist(lapply(mean.values.list, function(x) x[2]))
  p.values.list <- lapply(t.test.list,function(x) x$p.value)
  p.values <- unlist(p.values.list)
  p.values.adjust <- p.adjust(p.values, method="BH")
  p.values.log <- -log10(p.values.adjust)
  fold_change <- mean.values.T0/mean.values.T0
  FC_log2 <- log2(fold_change)
  ions <- c(1:length(p.values))
  out.table <- data.frame(ions, mean.values.T0, mean.values.T0, fold_change, FC_log2, SE, p.values.adjust)
  names(out.table) <- c("ions",T0_sample,T0_sample,"FC","log2_FC","SE","p.value.adjusted")
  write.table(out.table,file=table_name,row.names=FALSE,quote=FALSE,sep="\t")
  small_index <- which(FC_log2 > -1 & FC_log2 < 1)
  large_index <- which(FC_log2 < -1 | FC_log2 > 1)
  small_x <- FC_log2[small_index]
  small_y <- p.values.log[small_index]
  large_x <- FC_log2[large_index]
  large_y <- p.values.log[large_index]
  plot(large_x,large_y, pch=19, col="blue", main=strain, ylab="-log10 p-value", xlab="log2(Fold-change)", ylim=c(0,10), xlim=c(-6,6), cex=0.7)
  points(small_x, small_y, pch=19, col="grey", cex=0.7)
  abline(h=2,col="red",lty=4)
  grid()
} 
dev.off()

# t8 

T0 <-c("183_T0","184_T0","185_T0","186_T0")

T8 <-c("183_T8","184_T8","185_T8","186_T8")

nb_samples <- length(T0)
index <- c(1:nb_samples)
pdf("volcano_plot0.pdf")
par(mfrow=c(2,2))
for (i in index) {
  T0_sample <- T0[i]
  T8_sample <- T8[i]
  strain <- substr(T0_sample,1,6)
  table_name <- paste0(strain,"_FC.txt")
  t.test <- t_test(data, T8_sample, T0_sample)
  t.test.list <- t.test[[1]]
  SE <- t.test[[2]]
  mean.values.list <- lapply(t.test.list, function(x) x$estimate)
  mean.values.T0 <- unlist(lapply(mean.values.list, function(x) x[1]))
  mean.values.T8 <- unlist(lapply(mean.values.list, function(x) x[2]))
  p.values.list <- lapply(t.test.list,function(x) x$p.value)
  p.values <- unlist(p.values.list)
  p.values.adjust <- p.adjust(p.values, method="BH")
  p.values.log <- -log10(p.values.adjust)
  fold_change <- mean.values.T8/mean.values.T0
  FC_log2 <- log2(fold_change)
  ions <- c(1:length(p.values))
  out.table <- data.frame(ions, mean.values.T0, mean.values.T8, fold_change, FC_log2, SE, p.values.adjust)
  names(out.table) <- c("ions",T0_sample,T8_sample,"FC","log2_FC","SE","p.value.adjusted")
  write.table(out.table,file=table_name,row.names=FALSE,quote=FALSE,sep="\t")
  small_index <- which(FC_log2 > -1 & FC_log2 < 1)
  large_index <- which(FC_log2 < -1 | FC_log2 > 1)
  small_x <- FC_log2[small_index]
  small_y <- p.values.log[small_index]
  large_x <- FC_log2[large_index]
  large_y <- p.values.log[large_index]
  plot(large_x,large_y, pch=19, col="blue", main=strain, ylab="-log10 p-value", xlab="log2(Fold-change)", ylim=c(0,10), xlim=c(-6,6), cex=0.7)
  points(small_x, small_y, pch=19, col="grey", cex=0.7)
  abline(h=2,col="red",lty=4)
  grid()
} 
dev.off()

# t16

T0 <-c("183_T0","184_T0","185_T0","186_T0")

T16 <-c("183_T16","184_T16","185_T16","186_T16")

nb_samples <- length(T0)
index <- c(1:nb_samples)
pdf("volcano_plot0.pdf")
par(mfrow=c(2,2))
for (i in index) {
  T0_sample <- T0[i]
  T16_sample <- T16[i]
  strain <- substr(T0_sample,1,6)
  table_name <- paste0(strain,"_FC.txt")
  t.test <- t_test(data, T16_sample, T0_sample)
  t.test.list <- t.test[[1]]
  SE <- t.test[[2]]
  mean.values.list <- lapply(t.test.list, function(x) x$estimate)
  mean.values.T0 <- unlist(lapply(mean.values.list, function(x) x[1]))
  mean.values.T16 <- unlist(lapply(mean.values.list, function(x) x[2]))
  p.values.list <- lapply(t.test.list,function(x) x$p.value)
  p.values <- unlist(p.values.list)
  p.values.adjust <- p.adjust(p.values, method="BH")
  p.values.log <- -log10(p.values.adjust)
  fold_change <- mean.values.T16/mean.values.T0
  FC_log2 <- log2(fold_change)
  ions <- c(1:length(p.values))
  out.table <- data.frame(ions, mean.values.T0, mean.values.T16, fold_change, FC_log2, SE, p.values.adjust)
  names(out.table) <- c("ions",T0_sample,T16_sample,"FC","log2_FC","SE","p.value.adjusted")
  write.table(out.table,file=table_name,row.names=FALSE,quote=FALSE,sep="\t")
  small_index <- which(FC_log2 > -1 & FC_log2 < 1)
  large_index <- which(FC_log2 < -1 | FC_log2 > 1)
  small_x <- FC_log2[small_index]
  small_y <- p.values.log[small_index]
  large_x <- FC_log2[large_index]
  large_y <- p.values.log[large_index]
  plot(large_x,large_y, pch=19, col="blue", main=strain, ylab="-log10 p-value", xlab="log2(Fold-change)", ylim=c(0,10), xlim=c(-6,6), cex=0.7)
  points(small_x, small_y, pch=19, col="grey", cex=0.7)
  abline(h=2,col="red",lty=4)
  grid()
} 
dev.off()

# t24

T0 <-c("183_T0","184_T0","185_T0","186_T0")

T24 <-c("183_T24","184_T24","185_T24","186_T24")

nb_samples <- length(T0)
index <- c(1:nb_samples)
pdf("volcano_plot0.pdf")
par(mfrow=c(2,2))
for (i in index) {
  T0_sample <- T0[i]
  T24_sample <- T24[i]
  strain <- substr(T0_sample,1,6)
  table_name <- paste0(strain,"_FC.txt")
  t.test <- t_test(data, T24_sample, T0_sample)
  t.test.list <- t.test[[1]]
  SE <- t.test[[2]]
  mean.values.list <- lapply(t.test.list, function(x) x$estimate)
  mean.values.T0 <- unlist(lapply(mean.values.list, function(x) x[1]))
  mean.values.T24 <- unlist(lapply(mean.values.list, function(x) x[2]))
  p.values.list <- lapply(t.test.list,function(x) x$p.value)
  p.values <- unlist(p.values.list)
  p.values.adjust <- p.adjust(p.values, method="BH")
  p.values.log <- -log10(p.values.adjust)
  fold_change <- mean.values.T24/mean.values.T0
  FC_log2 <- log2(fold_change)
  ions <- c(1:length(p.values))
  out.table <- data.frame(ions, mean.values.T0, mean.values.T24, fold_change, FC_log2, SE, p.values.adjust)
  names(out.table) <- c("ions",T0_sample,T24_sample,"FC","log2_FC","SE","p.value.adjusted")
  write.table(out.table,file=table_name,row.names=FALSE,quote=FALSE,sep="\t")
  small_index <- which(FC_log2 > -1 & FC_log2 < 1)
  large_index <- which(FC_log2 < -1 | FC_log2 > 1)
  small_x <- FC_log2[small_index]
  small_y <- p.values.log[small_index]
  large_x <- FC_log2[large_index]
  large_y <- p.values.log[large_index]
  plot(large_x,large_y, pch=19, col="blue", main=strain, ylab="-log10 p-value", xlab="log2(Fold-change)", ylim=c(0,10), xlim=c(-6,6), cex=0.7)
  points(small_x, small_y, pch=19, col="grey", cex=0.7)
  abline(h=2,col="red",lty=4)
  grid()
} 
dev.off()

# Log2FC values reported in Supplementary Dataset 8, calculate how many compounds across
# strains have significant log2FC. The significant ones (log2FC <=-1 and log2FC >=1 and p-value <=0.01)
# will be used for continuing the analysis and produce lineplots. 

# Load packages
library(ggplot2)
library(dplyr)


# Load and parse data
line_data <- read.table("line_plots_all.txt", header = T, check.names = FALSE) %>%
  mutate(
    time = factor(time, c("t0", "t8", "t16", "t24")),
    metabolite = factor(metabolite, sort(unique(metabolite)))
  )

# Open new pdf
pdf("line_plots_all.pdf")

# Iterate over metabolites
for (metaboliteName in levels(line_data$metabolite)) {
  # Get metabolite specific data
  metaboliteData <- line_data %>% filter(metabolite == metaboliteName)
  
  # Create metabolite plot
  p <- ggplot(metaboliteData, aes(x = time, y = log2FC, group = species, color = species)) +
    # Add error bar
    geom_errorbar(aes(ymin = log2FC - SE, ymax = log2FC + SE), width = 0.1) +
    # Add line and points
    geom_line(size=1) + geom_point() +
    # Set custom color scale and theme
    scale_color_manual(values = c("Lapi" = "#df2e5c", "Lhel" = "#f7c655", "Lmel" = "#43a980", "Lkul" = "#2476a1")) + theme_minimal() +
    # Set plot title
    ggtitle(metaboliteName) +
    scale_x_discrete(labels=c("0", "8", "16", "24")) +
    # To edit plot title and axis legends and much more
    # https://ggplot2.tidyverse.org/reference/theme.html
    # https://ggplot2.tidyverse.org/reference/element.html
    theme(
      plot.title = element_text(size = 20, face = 'bold'),
      axis.title = element_text(size = 15),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 15)),
      axis.text.x = element_text(size=16),
      axis.text.y = element_text(size=16),
      legend.position = 'none'
    )
  
  # Print plot on pdf
  print(p)
}

# Close pdf
dev.off()

# Make logistic growth curves 


data <- read.table("copy_numb.txt", check.names = FALSE, h = T)
library(growthcurver)

pdf("lapi_lc.pdf")

data_lapi <- dplyr::filter(data, !grepl('Lmel|hel|Lkul', DV))
gc_fit_lapi <- SummarizeGrowth(data_lapi$group, data_lapi$Copy_number)
plot(gc_fit_lapi, xlab = "time", ylab = "Copy number", main = "Lapi")

dev.off()

pdf("lhel_lc.pdf")

data_lhel <- dplyr::filter(data, !grepl('Lmel|Lapi|Lkul', DV))
gc_fit_lhel <- SummarizeGrowth(data_lhel$group, data_lhel$Copy_number)
plot(gc_fit_lhel, xlab = "time", ylab = "Copy number", main = "Lhel")

dev.off()

pdf("lmel_lc.pdf")

data_lmel <- dplyr::filter(data, !grepl('Lhel|Lapi|Lkul', DV))
gc_fit_lmel <- SummarizeGrowth(data_lmel$group, data_lmel$Copy_number)
plot(gc_fit_lmel, xlab = "time", ylab = "Copy number", main = "Lmel")

dev.off()

pdf("lkul_lc.pdf")

data_lkul <- dplyr::filter(data, !grepl('Lhel|Lapi|Lmel', DV))
gc_fit_lkul <- SummarizeGrowth(data_lkul$group, data_lkul$Copy_number)
plot(gc_fit_lkul, xlab = "time", ylab = "Copy number", main = "Lkul")

dev.off()
