################################################################################
#                                                                              #
#  Wintercreeper: Bacterial Community Alpha Diversity                          #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Sarah Bray                                                       #
#                                                                              #
#	Initiated: 16 June 2015   
# Revised: 22 June 2015
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
even<- round(evenness.iter(input=shared, cutoff=level, size=size, iters=1000),3) #produced something crazy: 20x20005 matrix
shan<- round(diversity.iter(input=shared, cutoff=level, size=size, iters=1000),3)

# Calculation of mean richness and SE
rich_data<- merge(design, rich, by="row.names")
rich_data$mean<-round(apply(rich, 1, mean, na.rm = TRUE),3) #1 signifies calculating mean across rows)
rich_data$SE <- round(apply(rich, 1, se, na.rm = TRUE), 3)
rich_data$Design<-c("ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbWC",
                    "ArbWC","ArbWC","ArbWC","ArbWC","SGNative","SGNative","SGNative","SGNative",
                    "SGNative","SGWC","SGWC","SGWC","SGWC","SGWC")
rich_data.means<-rich_data[,1004:1007]
rich_data.trt.mean<-tapply(rich_data$mean,list(rich_data$Design), mean)
rich_data.trt.se<-tapply(rich_data$mean,list(rich_data$Design), se)

# ttest of differences in richness within site by cover
rich_Scott<-rich_data[11:20, 1004:1007]
t.test(mean~Design, data=rich_Scott)
rich_Arb<-rich_data[1:10, 1004:1007]
t.test(mean~Design, data=rich_Arb)

# Richness Plot
pdf("./Plots/rich.pdf")
Richness<-plot(as.factor(rich_data.means$Design),rich_data.means$mean,xlab="Location and Cover", ylab="OTU Richness")
dev.off()

# Calculation of mean evenness and SE
even_data<- merge(design, even, by="row.names")
even_data$mean<-round(apply(even, 1, mean, na.rm = TRUE),3) #1 signifies calculating mean across rows)
even_data$SE <- round(apply(even, 1, se, na.rm = TRUE), 3)
even_data$Design<-c("ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbWC",
                    "ArbWC","ArbWC","ArbWC","ArbWC","SGNative","SGNative","SGNative","SGNative",
                    "SGNative","SGWC","SGWC","SGWC","SGWC","SGWC")
even_data.anova<-even_data[,c(2:4,1005:1007)]
even_data.trt.mean<-tapply(even_data$mean,list(even_data$Design), mean)
even_data.trt.se<-tapply(even_data$mean,list(even_data$Design), se)


# ttest of differences in evenness within site by cover
even_Scott<-even_data[11:20, 105:107]
t.test(mean~Design, data=even_Scott)
even_Arb<-even_data[1:10, 105:107]
t.test(mean~Design, data=even_Arb)

#t-test of site effect evenness
t.test(mean~site, data=even_data.anova)


# evenness Plot
pdf("./Plots/even.pdf")
par(mfrow=c(1,1), mar=c(5,5,1,1)) 
evenness<-plot(as.factor(even_data.anova$Design),even_data.anova$mean,xlab="Location and Cover", ylab="Evenness")
dev.off()

# Calculation of mean shannon and SE
shan_data<- merge(design, shan, by="row.names")
shan_data$mean<-round(apply(shan, 1, mean, na.rm = TRUE),3) #1 signifies calculating mean across rows)
shan_data$SE <- round(apply(shan, 1, se, na.rm = TRUE), 3)
shan_data$Design<-c("ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbNative", "ArbWC",
                    "ArbWC","ArbWC","ArbWC","ArbWC","SGNative","SGNative","SGNative","SGNative",
                    "SGNative","SGWC","SGWC","SGWC","SGWC","SGWC")
shan_data.anova<-shan_data[,c(2:4,1005:1007)]
shan_data.trt.mean<-tapply(shan_data$mean,list(shan_data$Design), mean)
shan_data.trt.se<-tapply(shan_data$mean,list(shan_data$Design), se)

#t-test of site effect shannon
t.test(mean~site, data=shan_data.anova)

# ttest of differences in shannon within site by cover
shan_Scott<-shan_data[11:20, 1005:1007]
t.test(mean~Design, data=shan_Scott)
shan_Arb<-shan_data[1:10, 1005:1007]
t.test(mean~Design, data=shan_Arb)

# Shannon Plots
pdf("./Plots/shan.pdf")
par(mfrow=c(1,1), mar=c(5,5,1,1)) 
shannon<-plot(as.factor(shan_data.anova$Design),shan_data.anova$mean,xlab="Location and Cover", ylab="Shannon Index")
dev.off()
