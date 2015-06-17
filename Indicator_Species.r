################################################################################
#                                                                              #
#  Wintercreeper: Indicator Species Analysis                                   #
#                                                                              #
################################################################################
#                                                                              #                                                                          #
#	Written by: Sarah Bray
#
# First draft: 17 June 2015
#
#	Updated: 
#                                                                              #
################################################################################

# Set your directory
WCdir<-getwd()
setwd("C:/Users/sbray/Dropbox/Sabbatical/2014-Bray_WC_Analysis")

source("./R_files/Diversity_Functions/DiversityFunctions.r")  
require(vegan)  
require(BiodiversityR)
require(labdsv)
require(reshape)

WC <- read.otu(shared ="./WC.final.shared", "0.03")
design <- read.delim(file="./WC.design", header=T, row.names=1)
Arbdesign<-design[1:10,]
Scottdesign<-design[11:20,]

# Make separate matrices for Arboretum and Scotts Grove
Arb<- WC[1:10,]
Scotts<- WC[11:20,]

# Remove Zero Sum OTUs
Arb <- Arb[,!(colSums(abs(Arb)) ==0)]
Scotts<-Scotts[,!(colSums(abs(Scotts))==0)]

# Calculate Presense Absence
dataPA <- (WC_data > 0)*1 

# Calculating Relative Abundance
ArbREL <- Arb
for(i in 1:nrow(Arb)){
  ArbREL[i,] = Arb[i,]/sum(Arb[i,])
}  

ScottsREL <- Scotts
for(i in 1:nrow(Scotts)){
  ScottsREL[,i] = Scotts[,i]/sum(Scotts[,i])
} 

# Species Indicator Values, see Legendre and Legendre 2012, Numerical Ecology
# Code modified from p-meso.indval.r by Mario Muscarella
# Groups: Invaded = 1, Native = 2
Scotts.iva <- indval(ScottsREL, design$invasion, numitr=1000)
gr <- Scotts.iva$maxcls[Scotts.iva$pval <= 0.05] #mxcls = the class each species has maximum value for
iv <- Scotts.iva$indcls[Scotts.iva$pval <= 0.05] #indcls = indicator value for each species to its maximum class 
pv <- Scotts.iva$pval[Scotts.iva$pval   <= 0.05] #pval = p-value 
fr <- apply(ScottsREL>0, 2, sum)[Scotts.iva$pval <= 0.05] #fr is # of species that are indicators? 
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr) #making new data frame of significant outputs of indval 
fidg <- fidg[order(fidg$group, -fidg$indval), ] #ordering data frame by category

# Combine Indicator Values with Taxonomy
tax.raw3 <- (read.delim("./DataFiles/WC.taxonomy.03.tab")) #taxonomy tab delimited
include<-row.names(fidg) #creates vector of OTUs from significant species indicator values
include<-include[order(include)]#put OTUs in numerical order
tax3<-tax.raw3[include,] #attemping to subset only significant INDVal OTUs from taxonomy file; yeilds matrix of correct dimensions, but all NAs
Scotts.iva.data <- merge(fidg, tax2, by = "row.names", all.x=T ) #resulting dataframe has NAs in columns coming from taxonomy file
