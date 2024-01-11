#code for creating plots from STRUCTURE output files for both SNP and microsatellite loci 

options(rlib_downstream_check = FALSE)
library(pophelper)
library(ggplot2)
library(gridExtra)

#STRUCTURE RUNS
#code from http://www.royfrancis.com/pophelper/articles/index.html

#input "results" files interactively as list 
slist=readQ(files=choose.files(multi=TRUE),filetype="structure")
# check class of ouput
class(slist)
# view head of first converted file
head(slist[[1]])

#check attributes
attributes(slist)

#check tabulated list
head(tabulateQ(slist))

#summarized table of structure runs
tr1 <- tabulateQ(qlist=slist)
sr1 <- summariseQ(tr1)

#best k, evanno
evannoMethodStructure(data=sr1)
em <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=F)
evannoMethodStructure(data=sr1,exportplot=T,writetable=T,na.rm=T,exportpath=getwd())
sr1 <- summariseQ(tabulateQ(slist))
p <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p)

#add in labels 
#import population file 
snp_pops=read.table("pop_relate_substructure.txt", head=T)
onelabset <- snp_pops[,2,drop=FALSE]
head(onelabset)

# length of labels equal to number of individuals?
#nrow(snp_pops)
nrow(snp_pops)

# check if labels are a character data type
#sapply(snp_pops, is.character)
sapply(snp_pops, is.character)

#FULL FIGURE
slist=readQ(files=choose.files(multi=TRUE),filetype="structure", indlabfromfile=T)
slist_1 <- alignK(slist)
slist_2 <- mergeQ(slist_1)

slist2=readQ(files=choose.files(multi=TRUE),filetype="structure", indlabfromfile=T)
slist2_1 <- alignK(slist2)
slist2_2 <- mergeQ(slist2_1)

slist3=readQ(files=choose.files(multi=TRUE),filetype="structure", indlabfromfile=T)
slist3_1 <- alignK(slist3)
slist3_2 <- mergeQ(slist3_1)

slist4=readQ(files=choose.files(multi=TRUE),filetype="structure", indlabfromfile=T)
slist4_1 <- alignK(slist4)
slist4_2 <- mergeQ(slist4_1)

p2 <- plotQ(slist_2, returnplot=T,exportplot=F, basesize=11,
            grplab=onelabset,subsetgrp=c("IB", "IT","ALP", "BAV", "BW","RH", "SAX", "RU"),ordergrp=TRUE, grplabsize=0,
            linesize=0.8,pointsize=4, splabsize=12, 
            clustercol=c("mediumseagreen",  "lightsteelblue",  "wheat", "darkolivegreen", "tomato", "rosybrown3", "gray53", "lavender", "tan"), showtitle=F, titlesize=20,
            showsubtitle=F,  subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1, barbordercolour="white",barbordersize=0)


p3 <- plotQ(slist2_2, returnplot=T,exportplot=F, basesize=11,
            grplab=onelabset,subsetgrp=c("IB", "IT","ALP", "BAV", "BW","RH", "SAX", "RU"),ordergrp=TRUE, grplabsize=0,
            linesize=0.8,pointsize=4, splabsize=12, 
            clustercol=c("mediumseagreen",  "lightsteelblue",  "wheat", "darkolivegreen", "tomato", "rosybrown3", "gray53", "lavender", "sienna3"), showtitle=F, titlesize=20,
            showsubtitle=F,  subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1, barbordercolour="white",barbordersize=0)

p4 <- plotQ(slist3_2, returnplot=T,exportplot=F, basesize=11,
            grplab=onelabset,subsetgrp=c("IB", "IT","ALP", "BAV", "BW","RH", "SAX", "RU"),ordergrp=TRUE, grplabsize=5,
            linesize=0.8,pointsize=4, splabsize=12, 
            clustercol=c("mediumseagreen",  "lightsteelblue",  "wheat", "darkolivegreen", "tomato", "rosybrown3", "gray53", "lavender", "sienna3"), showtitle=F, titlesize=20,
            showsubtitle=F,  subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1,  barbordercolour="white",barbordersize=0)


grid.arrange(p2$plot[[1]],p3$plot[[1]], p4$plot[[1]],nrow=3)


p5 <- plotQ(slist4_2, returnplot=T,exportplot=F, basesize=11,
            grplab=onelabset,subsetgrp=c("IB", "IT","ALP", "BAV", "BW","RH", "SAX", "RU"),ordergrp=TRUE, grplabsize=5,
            linesize=0.8,pointsize=4, splabsize=12, 
            clustercol=c("mediumseagreen",  "lightsteelblue",  "wheat", "darkolivegreen", "tomato", "rosybrown3", "gray53", "lavender", "sienna3"), showtitle=F, titlesize=20,
            showsubtitle=F,  subtitlesize=16, height=1.6,indlabsize=2.3,indlabheight=0.08,indlabspacer=-1,  barbordercolour="white",barbordersize=0)

grid.arrange(p2$plot[[1]],p3$plot[[1]], p4$plot[[1]],p5$plot[[1]],nrow=4)
grid.arrange(p5$plot[[1]],nrow=1)