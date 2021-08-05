# load data
cts <- read.csv("ct_TC.csv")
cts$X <- NULL

######## CALCULATIONS ########

## Calculate CopyNum for each target
library(readxl)
StdCurves <- read_excel("StdCurves.xlsx", sheet=1) # export intercept, slope and threshold of detection values for each target
StdCurves <- as.data.frame(StdCurves) # table from excel is loaded as 'tibble' -> needs to be fixed
StdCurves$LOD_copies <- NULL
StdCurves$num16Sor18S <- NULL

n <- merge(cts, StdCurves) 

# calculate CopyNums of each target
n$CopyNum <- (n$Efficiency^(n$Intercept - n$Ct))*50 #*elution volume

# push all CopyNums to '1' when the Ct value is higher than its threshold or NA -> "1" (undetected), better for statistical analysis
n$CopyNum[n$Ct >= n$Threshold] <- 1 # later do the same for CellNums
n$CopyNum[is.na(n$Ct)] <- 1

# calculate number of each bacterium in gut - divide by the number of 16S/18S in genome of the organism
n$Target <- factor(n$Target) # to get rid of old factors

n$Intercept <- NULL # get rid of extra columns
n$Slope <- NULL
n$Efficiency <- NULL

write.table(n, "copy_numb_TC.txt")


