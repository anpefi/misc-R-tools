# sciFigureGOdist.R v.1.0
# Script for plotting GO Terms distributions obtained by Blast2go software
# by Andrés Pérez-Figueroa Jan'14

# HOW TO USE
#----------------------------
#
# 1. Obtain a txt file with the results from the combined graphs in Blast2Go by using the Save as TXT option. The file should contain the three ontology categories and no headers
#  
#   Note2: If Blast2Go was executed in a system with a Spanish locale (or other using a different character from . to separate decimals)
#          Then the commas should be replaced by dots. (in unix-based systems: cat GOdist.txt |  sed 's/,/./g' > enrichment.txt )
#
# 2. From the directory where the file is saved, open an R session
#
# 3. source this script: source("[PATH/TO/]sciFigureGOdist.R")
#
# 4. For the basic plotting, using the default parameters, just type:
#           plotdistGO("name_of_the_file")
#     This will produce a pdf file called outPlotGOdist.pdf with the barplot for those terms with a FDR<0.05
#
# 5. For further refinements you can set any of the following parameters by adding them into the brackets (comma-separated) of the calling to plotdistGO
# 
#   fil.seq=[value] - Threshold value number of sequences in a Term to be represented (default: 10)
#
#   fil.level=[value] - Ontology level to be represented (default: 3)
#
#   by.cats=[TRUE/FALSE] - Box by categories
#
#   colores=[Array of 3 strings] - array with the colors to be used as BP, MF and CC
#
#  type=["pdf"/"eps"/"screen"] - Set the way to generate the output. 
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
#
#############################################################################
plotdistGO <- function(file, fil.seq=10, fil.level=3, by.cats=T, type="pdf", outfile="outPlotGOdist", fig.dim=c(6,9), margins = c(5,19,5,3), cex=1, xlab="# Loci", colores=c("gray","gray","gray")){
  
  # Read data
  clases <- c("factor","character","character","factor","numeric","numeric","character")
  data <- read.table(file, header=F, sep="\t", colClasses=clases, quote="\"")[,1:6]
  colnames(data)<-c("Level","GO-ID","Term","Category","Seqs","Score")
  levels(data$Category) <- c("BP","CC","MF")
  
  #Filter data by category, level and #Seq FDR
  CC <- data[which((data$Category=="CC")&(data$Seqs>=fil.seq)&(data$Level==fil.level)),c(3,5)] #Cellular Component
  MF <- data[which((data$Category=="MF")&(data$Seqs>=fil.seq)&(data$Level==fil.level)),c(3,5)] #Molecular Function
  BP <- data[which((data$Category=="BP")&(data$Seqs>=fil.seq)&(data$Level==fil.level)),c(3,5)] #Biological Process
  
  #Sorting by Seq frequency
  CC <- CC[order(CC$Seqs),]
  MF <- MF[order(MF$Seqs),]
  BP <- BP[order(BP$Seqs),]
 
  
  Terms <- c(CC$Term,"",MF$Term,"",BP$Term)
  Seqs <- c(CC$Seqs,NA,MF$Seqs,NA,BP$Seqs)
  colors <- c( rep(colores[1], nrow(CC)+1), rep(colores[2], nrow(MF)+1), rep(colores[3], nrow(BP)+1) )
  
  # define graphic device
  if (type=="pdf") pdf(file=paste(outfile,"pdf",sep="."), fig.dim[1], fig.dim[2])  
  if (type=="eps") postscript(file=paste(outfile,"eps",sep="."), onefile=F, horizontal=F, paper="special",height=fig.dim[2],width=fig.dim[1] )
    
  par( mar = margins, las=1, cex=cex)
  
  barplot(Seqs,  names.arg= Terms, xlim=c(0,1.1*max(Seqs, na.rm=T)),  col=colors, xlab=xlab, main=paste("Level",fil.level),space=0.35, horiz=T, xpd=F, offset=0, width=0.7)
  abline(v=0)
  legend("top",legend=c("BP", "MF", "CC"),  horiz=T, fill=colores[3:1], xpd=T, bty="n" )
  if(by.cats){
    abline(h=c(nrow(CC), nrow(CC)+nrow(MF)))
    box()
  }
  
  if((type=="pdf")||(type=="eps")) dev.off()  
}

