###################################################################
## make a barplot of KLD vs. number of clusters, with bars divided in individual group contributions
## use as: R-2.9 --vanilla --args [inputfile] [output_graph] < R-KLDvsClust_barplot.R

cmd_args = commandArgs(TRUE);

table = cmd_args[1];
ograph = cmd_args[2]; 

options(warn=-1)

## move data into a matrix
frame=read.delim(table,header=T)

heights=t(as.matrix(frame))

## row names
rowlab=as.numeric(row.names(frame))

## number of clusters
ming=min(rowlab)
maxg=max(rowlab)

nr=maxg
nc=maxg-ming+1

##labels
Glab=mat.or.vec(nr, nc)

if (nc==1 && nr==1) { 
  Glab="1"
} else if (nc==1) {
      for (x in c(1:nr)) {
         if (heights[x,1]!=0) {
	     Glab[x]=x
	  } else {
	     Glab[x]=""
	  }
      }
} else {
     for (x in c(1:nr)) {
        for (y in c(1:nc)) {
	  if (heights[x,y]!=0) {
	      Glab[x,y]=x
	   } else {
	      Glab[x,y]=""
	   }
        }
     } 
}
Glab=t(Glab)

png(ograph,height=600,width=600)
#pdf(ograph)

##make barplot
b=barplot(heights,horiz=FALSE, names.arg=rowlab, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1, col="black", border="white", xlab="Number of clusters", ylab="Kullbach Leibler distance")

ypos <- apply(heights, 2, cumsum)
ypos <- ypos - heights/2
ypos=t(ypos)
text(b,ypos,Glab,col="white",cex=1.5)

dev.off()


