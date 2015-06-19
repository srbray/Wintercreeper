################################################################################
#                                                                              #
#  Wintercreeper: Bacterial Community Alpha Diversity                          #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Sarah Bray                                                       #
#                                                                              #
#	Initiated: 16 June 2015                                                      #
#                                                                              #
################################################################################
#
# Notes: Adapted from source code from Mario Muscarella                        #
#
################################################################################

# Set working environment
setwd("C:/Users/sbray/GitHub/Wintercreeper")
source("./DiversityFunctions.r")
require("vegan")

#Definition of inputs and parameters
design<- "./Data/WC.design"
shared<- "./Data/WC.final.shared"
level <- "0.03"
size <- 60000 #SGN3 had only 67390 counts
design<-read.delim(design, header=T, row.names=1)
se<-function(x, ...){sd(x,...)/sqrt(length(x))}

# Alpha Diversity with Resampling
rich <- round(richness.iter(input=shared, cutoff=level, size=size, iters=1000), 3)
even<- round(evenness.iter(input=shared, cutoff=level, size=size, iters=100),3) #produced something crazy: 20x20005 matrix
shan<- round(diversity.iter(input=shared, cutoff=level, size=size, iters=100),3)

# Calculation of mean richness and SE
rich_data<- merge(design, rich, by="row.names")
rich_data$mean<-round(apply(rich, 1, mean, na.rm = TRUE),3) #1 signifies calculating mean across rows)
rich_data$SE <- round(apply(rich, 1, se, na.rm = TRUE), 3)
rich_data$Design<-c("ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbWC",
                    "ArbWC","ArbWC","ArbWC","ArbWC","SGNative","SGNative","SGNative","SGNative",
                    "SGNative","SGWC","SGWC","SGWC","SGWC","SGWC")
rich_data.means<-rich_data[,1004:1007]


# ttest of differences in richness within site by cover
rich_Scott<-rich_data[11:20, 1004:1007]
t.test(mean~Design, data=rich_Scott)
rich_Arb<-rich_data[1:10, 1004:1007]
t.test(mean~Design, data=rich_Arb)

# Richness Plot
attach(rich_data.means)
Design<-factor(Design)
Richness<-plot(Design,mean,xlab="Location and Cover", ylab="OTU Richness")
pdf("./Plots/rich.pdf")
plot(Design,mean,xlab="Location and Cover", ylab="OTU Richness")
dev.off()

# Calculation of mean shannon and SE
shan_data<- merge(design, shan, by="row.names")
shan$mean<-round(apply(shan, 1, mean, na.rm = TRUE),3) #1 signifies calculating mean across rows)
shan_data$SE <- round(apply(shan, 1, se, na.rm = TRUE), 3)
shan_data$Design<-c("ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbWC",
                    "ArbWC","ArbWC","ArbWC","ArbWC","SGNative","SGNative","SGNative","SGNative",
                    "SGNative","SGWC","SGWC","SGWC","SGWC","SGWC")
shan_data.means<-shan_data[,1005:1007]


# ttest of differences in shannon within site by cover
shan_Scott<-shan_data[11:20, 1005:1007]
t.test(mean~Design, data=shan_Scott)
shan_Arb<-shan_data[1:10, 1005:1007]
t.test(mean~Design, data=shan_Arb)

# shannon Plot
attach(shan_data.means)
Design<-factor(Design)
shannon<-plot(Design,mean,xlab="Location and Cover", ylab="Shannon Index")
pdf("./Plots/even.pdf")
plot(Design,mean,xlab="Location and Cover", ylab="OTU evenness")
dev.off()
