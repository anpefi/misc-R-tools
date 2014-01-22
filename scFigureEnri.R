# sciFigureEnri.R v.1.0
# Script for plotting results from Go Enrichment (fisher_analysis) obtained by Blast2go software
# by Andrés Pérez-Figueroa Jan'14

# HOW TO USE
#----------------------------
#
# 1. Obtain a txt file with the results for the Enrichment Analysis in Blast2Go by using the Save as TXT option
#   Note1: Blast2Go would ask wich value of FDR or p-value youwant to use as threshold to perform the analysis and save the txt.
#          I suggest to use a large value as you can trim it later with this script
#   Note2: If Blast2Go was executed in a system with a Spanish locale (or other using a different character from . to separate decimals)
#          Then the commas should be replaced by dots. (in unix-based systems: cat fisher_reults_raw.txt |  sed 's/,/./g' > enrichment.txt )
#
# 2. From the directory where the file is saved, open an R session
#
# 3. source this script: source("[PATH/TO/]sciFigureEnri.R")
#
# 4. For the basic plotting, using the default parameters, just type:
#           plotEnriGO("name_of_the_file")
#     This will produce a pdf file called outPlotEnri.pdf with the barplot for those terms with a FDR<0.05
#
# 5. For further refinements you can set any of the following parameters by adding them into the brackets (comma-separated) of the calling to plotEnriGO
# 
#   fdr=[value] - Threshold value for FDR to be included in the plot (default: 0.05)
#
#   by.cats=[TRUE/FALSE] - 
#
#   colores=[Array of 2 strings] - array with the colors to be used as Observed and Expected bars (default: c("red","blue"))
#
#   show.pvals=[TRUE/FALSE] - Add the p-values (FDR corrected) to the right of the bars (default: FALSE)
#
#   type=["pdf"/"eps"/"screen"] - Set the way to generate the output. 
#                           "screen" will use the R default device to show the plot. 
#                           "pdf" (default) will save the plot in a pdf file.
#                           "eps" will save the plot in a eps file suitable for publication.
#
#   outfile=[string] - Name of the outputfile (no extension required) if type is not "screen".
#
#   fig.dim=[c(width, height)] - Dimension, in inches, for the pdf/eps figure.
#
#   margins=[c(bottom,left,top,right)] - inner margins to be passed to the graphical parameter par(mar)
#
#   lab.legend=[c("O","E")/"none"] - array with the legend labels for the Observed and Expected bars. If "none" is provided, then no legend is shown
#
#############################################################################
plotEnriGO <- function(file, fdr=0.05, by.cats=T, colores=c("red","blue"), lab.legend=c("O","E"),show.pvals=F, type="pdf", outfile="outPlotEnri", 
                       fig.dim=c(6,9), margins = c(5,19,5,3), cex=0.7, xlab="% of sequences"){
  
  # Read data
  clases <- c("character","character","factor","numeric","numeric","numeric","numeric","numeric","numeric","character","character","character")
  data <- read.table(file, , header=F, skip=1, sep="\t", colClasses=clases, quote="\"")[,1:10]
  colnames(data)<-c("GO-ID","Term","Category","FDR","P-Value","n.Test","n.Ref","n.notAnnotTest","n.notAnnotRef","Over/Under")
  
  #Filter data by category and threshold FDR
  CC <- data[which((data$Category=="C")&(data$FDR<fdr)),] #Cellular Component
  MF <- data[which((data$Category=="F")&(data$FDR<fdr)),] #Molecular Function
  BP <- data[which((data$Category=="P")&(data$FDR<fdr)),] #Biological Process
 
  #Initialize some variables
  pval<-NULL
  allData <- NULL
  
  #Get the values to be plotted, as percentage, and taking care of missing Go categories.
  if(nrow(CC)>0){
    CCTest <- c(rev(100*(CC$n.Test/(CC$n.Test+CC$n.notAnnotTest))),NA)
    CCRef <- c(rev(100*(CC$n.Ref/(CC$n.Ref+CC$n.notAnnotRef))),NA)
    pval<-c(pval,rev(CC$FDR), NA) #Reverse order for p-values
  } else {
    CCTest <- NULL
    CCRef <- NULL
  }
  if(nrow(MF)>0){
    MFTest <- c(rev(100*(MF$n.Test/(MF$n.Test+MF$n.notAnnotTest))),NA)
    MFRef <- c(rev(100*(MF$n.Ref/(MF$n.Ref+MF$n.notAnnotRef))),NA)
    pval<-c(pval,rev(MF$FDR), NA)
  } else {
    MFTest <- NULL
    MFRef <- NULL
  }
  if(nrow(BP)>0){
    BPTest <- c(rev(100*(BP$n.Test/(BP$n.Test+BP$n.notAnnotTest))))
    BPRef <- c(rev(100*(BP$n.Ref/(BP$n.Ref+BP$n.notAnnotRef))))
    pval<-c(pval,rev(BP$FDR))
  }  else {
    BPTest <- NULL
    BPRef <- NULL
  }

  
  if(nrow(CC)>0) allData <- cbind(allData,matrix(c(CCTest,CCRef), nrow=2, byrow=T))[2:1,] #Reverse order to get grom most to less significant
  if(nrow(MF)>0) allData <- cbind(allData,matrix(c(MFTest,MFRef), nrow=2, byrow=T))[2:1,]
  if(nrow(BP)>0) allData <- cbind(allData,matrix(c(BPTest,BPRef), nrow=2, byrow=T))[2:1,]
 
  Terms <- c(rev(CC$Term),"",rev(MF$Term),"",rev(BP$Term))
  #Terms <- sapply(Terms, splitstring, simplify="array")
  
  # define graphic device
  if (type=="pdf") pdf(file=paste(outfile,"pdf",sep="."), fig.dim[1], fig.dim[2])  
  if (type=="eps"){
    #setEPS()
    postscript(file=paste(outfile,"eps",sep="."), onefile=F, horizontal=F, paper="special",height=fig.dim[2],width=fig.dim[1] )
  }
  
  par( mar = margins, las=1, cex=cex)
  
  barplot(allData, xlim=c(0,ceiling(max(allData, na.rm=T))), xlab=xlab, names.arg= Terms, beside=T, space=c(0,0.5), col=colores, horiz=T, xpd=T, width=0.4)
  
  if (show.pvals) text( maxVal+0.5,1:dim(allData)[2]-0.4,labels=bquote(italic("p")~"="~.(sprintf("%.1E", pval))), cex=0.8)
  if(by.cats){
    abline(h=c(nrow(CC)-1, nrow(CC)+nrow(MF))+1.5)
    if(nrow(CC)>0) mtext("CC", side=4, line=0.5, at=(nrow(CC)/2)+0.2)
    if(nrow(MF)>0) mtext("MF", side=4, line=0.5, at=nrow(CC)+(nrow(MF)/2)+0.2)
    if(nrow(BP)>0) mtext("BP", side=4, line=0.5, at=nrow(CC)+nrow(MF)+(nrow(BP)/2)+0.2)
  }
  box()
  if(lab.legend[1]!="none") legend("top",legend=lab.legend, inset=c(0,-0.07),  horiz=T, fill=colores[2:1], xpd=T, bty="n" )
  if((type=="pdf")||(type=="eps")) dev.off()  
}

