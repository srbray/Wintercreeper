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
#	Date: 9 Sept 2014
#                                                                              #
################################################################################

  # Set your directory
  WCdir<-getwd()
  setwd("C:/Users/sbray/Dropbox/Sabbatical/2014-Bray_WC_Analysis")
  
#wc.bysite <- function(shared = " ", design = " ", plot.title = "PCoA"){

  source("./R_files/Diversity_Functions/DiversityFunctions.r")  
  require(vegan)                                   

  WC <- read.otu(shared ="./WC.final.shared", "0.03")
  design <- read.delim(file="./WC.design", header=T)
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
  ScottsREL.dist <- vegdist(decostand(t(ScottsREL),method="log", logbase=10),method="bray")


  # Principal Coordinates Analysis
  Arbpcoa <- cmdscale(ArbREL.dist,k=3,eig=TRUE,add=FALSE) 
  Scottspcoa <- cmdscale(ScottsREL.dist,k=3,eig=TRUE,add=FALSE)
    # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
    # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(Arbpcoa$eig[1]/sum(Arbpcoa$eig)*100,2) 
  explainvar2 <- round(Arbpcoa$eig[2]/sum(Arbpcoa$eig)*100,2)
  explainvar3 <- round(Arbpcoa$eig[3]/sum(Arbpcoa$eig)*100,2)
  
  pcoap <- merge(as.data.frame(Arbpcoa$points),design,by=0,all.x=T)[,-1]
  rownames(pcoap) <- rownames(Arbpcoa$points)
  
  # Plot Parameters
  par(mfrow=c(1,1), mar=c(5,5,1,1)) 
  layout(rbind(1, 2), height=c(7, 1)) 
  y.dim <- c(min(pcoap$V2)+min(pcoap$V2)*0.2,max(pcoap$V2)+max(pcoap$V2)*0.2)
  x.dim <- c(min(pcoap$V1)+min(pcoap$V1)*0.2,max(pcoap$V1)+max(pcoap$V1)*0.2)
  
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
  text(pcoap$V1, pcoap$V2, labels=pcoap$sample, pos=1, offset=0.5, cex=0.7)
  box(lwd=2)
  par(mar=c(0, 3, 0, 0))
  plot.new()
  #legend("center", c(paste("Scotts Grove; Invaded"),
    #paste("Scotts Grove; Native"), 
    #paste("Arboretum; Invaded"),
    #paste("Arboretum; Native")), 
    #pt.lwd=2, col="black", pt.bg=c("purple", "white", "purple", 
    #"white"), pch=c(21,21,22,22), bty='n', ncol=2, cex=1.5, pt.cex=2)
  dev.copy2pdf(file=paste("./plots/",plot.title,".pdf",sep=""))
  dev.copy(png, file=paste("./plots/",plot.title,".png",sep=""), width=72*(7*4), 
    height=72*(8*4), res=72*4)
  dev.off()
  
  # Adonis (PERMANOVA)
  # Adonis runs a PERMANOVA (Created by Marti J. Anderson) 
  # this is very similar to ANOVA but for multivariate data
  # You can make very complex experimental designs with it
  # The default distance measure is bray-curtis, but other measures 
  # (Chao, Jaccard, Euclidean) can be used when specified  
  Adonis <- adonis(t(ArbREL.dist) ~ Arbdesign$invasion, method="bray", 
    permutations=1000)
  AdonisScott<-adonis(t(ScottsREL.dist) ~ Scottdesign$invasion, method="bray", 
                      permutations=1000)