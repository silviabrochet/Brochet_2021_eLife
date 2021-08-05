# In vivo

mydata_invivo <- read.csv("CFUs_stats_invivo.csv")

# Test for normality data-distribution 

shapiro.test(mydata_invivo$value)
mydata_invivo$value <- log10(mydata_invivo$value)

qqnorm(mydata$value)

# Perform ANOVA

anova_CFUs <- aov(value~condition*group, data=mydata_invivo)
summary(anova_CFUs)

# Perform Tuckey post-hoc test to identify each comparison's p-value

TukeyHSD <- TukeyHSD(anova_CFUs)

# In vitro

mydata_invitro <- read.csv("CFUs_stats_invitro.csv")

# Test for normality data-distribution 

shapiro.test(mydata_invitro$value)
mydata_invitro$value <- log10(mydata_invitro$value)

qqnorm(mydata$value)

# Perform ANOVA

anova_CFUs <- aov(value~condition*group, data=mydata_invitro)
summary(anova_CFUs)

# Perform Tuckey post-hoc test to identify each comparison's p-value

TukeyHSD <- TukeyHSD(anova_CFUs)

