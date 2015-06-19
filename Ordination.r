################################################################################
#                                                                              #
#	Wintercreeper: Microbial Community PcoA and PERMANOVA             #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Modified by: Sarah Bray
#
#	Updated: 16 June 2014
#                                                                              #
################################################################################

  # Set your directory
  WCdir<-getwd()
  setwd("C:/Users/sbray/GitHub/Wintercreeper")
  
  source("./DiversityFunctions.r")  
  require(vegan)  
  require(BiodiversityR)

  WC <- read.otu(shared ="./Data/WC.final.shared", "0.03")
  design <- read.delim(file="./Data/WC.design", header=T, row.names=1)
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

  # Create Distance Matrix
  ArbREL.dist <- vegdist(decostand((ArbREL),method="log", logbase=10),method="bray")
  ScottsREL.dist <- vegdist(decostand((ScottsREL),method="log", logbase=10),method="bray")


  # Principal Coordinates Analysis
  Arbpcoa <- cmdscale(ArbREL.dist,k=3,eig=TRUE,add=FALSE) 
  ArbpcoaS<-add.spec.scores(Arbpcoa,(ArbREL),method="cor.scores",Rscale=TRUE,scaling=1,multi=1) 
  Scottspcoa <- cmdscale(ScottsREL.dist,k=3,eig=TRUE,add=FALSE)
  ScottspcoaS<-add.spec.scores(Scottspcoa,(ScottsREL),method="cor.scores",Rscale=TRUE,scaling=1,multi=1) 

    # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
    # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

  # OTU correlations with PCoA axes (Arboretum)
  Arbcor<-as.data.frame(ArbpcoaS$cproj)
  colnames(Arbcor)[1]<-"Corr1"
  colnames(Arbcor)[2]<-"Corr2"
  Arbcor<-Arbcor[order(Arbcor[,1],decreasing=TRUE),] # sorts based on r of PCoA axis 1
  Arbcor<-Arbcor[!is.na(Arbcor[,1]) | !(is.na(Arbcor[,2])),] # remove rows where correlations = NA
  Arbresponders_a1<-Arbcor[Arbcor[,1]>=0.7|Arbcor[,2]<=-0.7,] # subset PCoA axis 1 based on r > |0.7|; see Fowler J, Cohen L, Jarvis P
  Arbresponders_a2<-Arbcor[Arbcor[,2]>=0.7|Arbcor[,2]<=-0.7,] # subset PCoA axis  2based on r > |0.7|;

  # OTU correlations with PCoA axes (Scotts Grove)
  Scottscor<-as.data.frame(ScottspcoaS$cproj)
  colnames(Scottscor)[1]<-"Corr1"
  colnames(Scottscor)[2]<-"Corr2"
  Scottscor<-Scottscor[order(Scottscor[,1],decreasing=TRUE),] # sorts based on r of PCoA axis 1
  Scottscor<-Scottscor[!is.na(Scottscor[,1]) | !(is.na(Scottscor[,2])),] # remove rows where correlations = NA
  Scottsresponders_a1<-Scottscor[Scottscor[,1]>=0.7|Scottscor[,2]<=-0.7,] # subset PCoA axis 1 based on r > |0.7|; see Fowler J, Cohen L, Jarvis P
  Scottsresponders_a2<-Scottscor[Scottscor[,2]>=0.7|Scottscor[,2]<=-0.7,] # subset PCoA axis  2based on r > |0.7|;

  # Creating Vectors for soil variables
  soilarb<-read.csv("./DataFiles/soilarb.csv", header=T, row.names=1)
  Arbsoil<-envfit(Arbpcoa, soilarb)

  soil<-read.csv("./DataFiles/soil.csv", header=T, row.names=1)
  Scottssoil<- envfit(Scottspcoa, soil, na.rm=T)



  # Percent Variance Explained Using PCoA (Axis 1,2,3)- Arboretum
  explainvar1 <- round(Arbpcoa$eig[1]/sum(Arbpcoa$eig)*100,2) 
  explainvar2 <- round(Arbpcoa$eig[2]/sum(Arbpcoa$eig)*100,2)
  explainvar3 <- round(Arbpcoa$eig[3]/sum(Arbpcoa$eig)*100,2)

  # Percent Variance Explained Using PCoA (Axis 1,2,3)- Scotts Grove
  explainvar1 <- round(Scottspcoa$eig[1]/sum(Scottspcoa$eig)*100,2) 
  explainvar2 <- round(Scottspcoa$eig[2]/sum(Scottspcoa$eig)*100,2)
  explainvar3 <- round(Scottspcoa$eig[3]/sum(Scottspcoa$eig)*100,2)
  
  pcoap <- merge(as.data.frame(Arbpcoa$points),Arbdesign,by=0,all.x=T)[,-1]
  rownames(pcoap) <- rownames(Arbpcoa$points)

  pcoap <- merge(as.data.frame(Scottspcoa$points),Scottdesign,by=0,all.x=T)[,-1]
  rownames(pcoap) <- rownames(Scottspcoa$points)
  
  # Plot Parameters Arb
par(mfrow=c(1,1), mar=c(5,5,1,1)) 
layout(rbind(1, 2), height=c(7, 1)) 
y.dim <- c(min(pcoap$V2)+min(pcoap$V2)*0.2,max(pcoap$V2)+max(pcoap$V2)*0.2)
x.dim <- c(min(pcoap$V1)+min(pcoap$V1)*0.2,max(pcoap$V1)+max(pcoap$V1)*2.2)

# Initiate Plot Arb
plot(pcoap$V1, pcoap$V2, xlab=paste("PCoA Axis 1 (",explainvar1, "%)", sep="")
     , ylab=paste("PCoA Axis 2 (",explainvar2, "%)", sep=""), 
     xlim=x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
     yaxt="n", cex.lab=1.5, cex.axis=1.2)  
axis(side=1, las=1)   
axis(side=2, las=1)    
abline(h=0, lty="dotted")  
abline(v=0, lty="dotted")

invasion.color <- rep(NA, dim(pcoap)[1])
for (i in 1:length(invasion.color)){
  if (pcoap$invasion[i] == "Invaded") {invasion.color[i] = "purple"}
  else {invasion.color[i] = "white"}
} 
points(pcoap$V1, pcoap$V2,pch=22, cex=2.0, col="black", bg=invasion.color, lwd=2)   
text(pcoap$V1, pcoap$V2, labels=pcoap$sample.1, pos=4, offset=0.5, cex=0.7)
plot(Scottssoil, p.max=0.1, cex = 0.75)
legend(-0.5,0.35,c("Invaded","Native"), pch=c(22,22),pt.lwd=2, col="black", pt.bg=c("purple", "white"), bty='y')



# Plot Parameters: Scotts Grove
par(mfrow=c(1,1), mar=c(5,5,1,1)) 
layout(rbind(1, 2), height=c(7, 1)) 
y.dim <- c(min(pcoap$V2)+min(pcoap$V2)*0.2,max(pcoap$V2)+max(pcoap$V2)*0.2)
x.dim <- c(min(pcoap$V1)+min(pcoap$V1)*0.2,max(pcoap$V1)+max(pcoap$V1)*2.2)

# Initiate Plot
plot(pcoap$V1, pcoap$V2, xlab=paste("PCoA Axis 1 (",explainvar1, "%)", sep="")
     , ylab=paste("PCoA Axis 2 (",explainvar2, "%)", sep=""), 
     xlim=x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
     yaxt="n", cex.lab=1.5, cex.axis=1.2)  
axis(side=1, las=1)   
axis(side=2, las=1)    
abline(h=0, lty="dotted")  
abline(v=0, lty="dotted")

invasion.color <- rep(NA, dim(pcoap)[1])
for (i in 1:length(invasion.color)){
  if (pcoap$invasion[i] == "Invaded") {invasion.color[i] = "purple"}
  else {invasion.color[i] = "white"}
} 
points(pcoap$V1, pcoap$V2,pch=22, cex=2.0, col="black", bg=invasion.color, lwd=2)   
text(pcoap$V1, pcoap$V2, labels=pcoap$sample.1, pos=4, offset=0.5, cex=0.7)
plot(Scottssoil, p.max=0.1, cex = 0.75)
#box(lwd=2)
#par(mar=c(0, 3, 0, 0)) #marios
#par(mar=c(4,4,4,4))
#plot.new()
legend(0.45,0.35,c("Invaded","Native"), pch=c(22,22),pt.lwd=2, col="black", pt.bg=c("purple", "white"), bty='y')
# Adonis (PERMANOVA)
# Adonis runs a PERMANOVA (Created by Marti J. Anderson) 
  # this is very similar to ANOVA but for multivariate data
  # You can make very complex experimental designs with it
  # The default distance measure is bray-curtis, but other measures 
  # (Chao, Jaccard, Euclidean) can be used when specified  
  Adonis <- adonis(ArbREL.dist ~ Arbdesign$invasion, method="bray", 
    permutations=1000)
  AdonisScott<-adonis(ScottsREL.dist ~ Scottdesign$invasion, method="bray", permutations=1000)

  


  
  
