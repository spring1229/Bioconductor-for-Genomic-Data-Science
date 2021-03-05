plotRanges<-function(x, xlim=x, main=deparse(substitue(x)),
+ col='black',sep=0.5, ...) {
+ height<-1
+ if (is(xlim, 'Ranges'))
+ xlim <-c(min(start(xlim)), max(end(xlim)))
+ bins<-disjointBins(IRanges(start(x), end(x)+1))
+ plot.new()
+ plot.window(xlim, c(0, max(bins)*(height+sep)))
+ ybottom<-bins *(sep+height)-height
+ rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom+height, col=col, ...)
+ title(main)
+ axis(1)
+ }
plotRanges(ir)
par(mfrow=c(2,1))
countOverlaps(ir1,ir2)
library(GenomicRanges)
library(GenomeInfoDb)
gr=GRanges(seqnames=c('chr1'), strand=c('+','-','+'), ranges=IRanges(start=c(1,3,5),
                                                                     width=3))
gr
seqinfo(gr)
seqlength(gr)=c('chr1'=10)
seqlengths(gr)=c('chr1'=10)
seqlevels(gr)
gaps(gr)
sort(gr)
genome(gr)='hg19'
gr
gf2=gr
genome(gr2)='hg18'
findOverlaps(gr,gf2)
ir=IRanges(start=1:3, width=2)
df=DataFrame(ir=ir, score=rnorm(3))
df2=DataFrame(ir=ir)
df2
gr<-GRanges(seqnames='chr1', strand=c('+','-','+'))
gr
subsetByOverlaps(gr2,gr)
gr2<-DataFrame(seqnames='chr1', start=1:3, end=4:6, score=rnorm(3))
gr2
library(GenomicRanges)
gr=GRanges(seqnames=c('chr1','chr2'), ranges=IRanges(start=1:2, end=4:5))
gr
dropSeqlevels(gr,'chr2')
keepStandardChromosomes(gr)
unique(ah$dataprovider)
library(AnnotationHub)
ah=AnnotationHub()
ah
ah[1]
unique(ah$dataprovider)
unique(ah$species)
ah=subset(ah,species='Homo sapiens')
ah
query(ah, 'H3K4me3')
ah2=display(ah)
ah2
ahub=AnnotationHub()
ahub=subset(ahub, species=='Homo sapiens')
qhs=query(ahub, c('H3K4me3','Gm12878'))
qhs
gr1=qhs[[2]]
gr2=qhs[[4]]
summary(width(gr1))
summary(width(gr2))
table(width(gr2))
qhs[4]
qhs=query(ahub,'RefSeq')
qhs$genome
genes=qhs[[1]]
table(table(genes$name))
prom=promoters(genes)
table(width(prom))
prom
args(promoters)
findOverlaps(prom, peaks)
length(unique(queryHits(prom)))
length(unique(queryHits(ov)))
length(subsetByOverlaps(peaks, prom, ignore.strand=True))/length(prom)
sum(width(reduce(peaks, ignore.strand=True))
    +)
inOut=matrix(0,ncol=2,nrow=2)
colnames(inOut)=c('in','out')
rownames(inOut)=c('in','out')
inOut
inOut[1,1]=sum(width(intersect(peaks, prom, ignore.strand=True)))
inOut[1,2]=sum(width(setdiff(peaks, prom, ignore.strand=True)))
inOut[2,1]=sum(width(setdiff(prom, peaks, ignore.strand=True)))
inOut[2,2]=3*10^9-sum(inOut)
inOut
OddsRatio=inOut[1,1]*inOut[2,2]/(inOut[2,1]*inOut[1,2])
OddsRatio
ah=AnnotationHub()
ah_human_CpG <- query(ah, c('CpG Islands', 'hg19'))
ah_human_CpG
ah_human_CpG_data<- ah_human_CpG[['AH5086']]
ah_human_CpG_reduce<-reduce(ah_human_CpG_data)
autosome<-c(paste('chr', 1:22, sep=''))
> split_data_by_chr<-split(ah_human_CpG_reduce, seqnames(ah_human_CpG_reduce))
> autosome_CpG_data<-split_data_by_chr[autosome]
> unlist(autosome_CpG_data)
#question 3
autosome_CpG_data[4]
#question 4
ah_H3K4me<- query(ah, c('H3K4me3','E003'))
ah_H3K4me_data <- ah_H3K4me[['AH29884']]
ah_H3K4me_autosome_data <-subset(ah_H3K4me_data, seqnames %in% autosome)
unlist(ah_H3K4me_autosome_data)

library(GenomicRanges)
ir1<-IRanges(start=c(1,3,5), end=c(3,5,7), width=3)
ir1
plotRanges(ir1)
resize(ir1,width=1,fix='start')
gr=GRanges(seqnames=c('chr1'), strand=c('+','-','+'), 
           ranges=IRanges(start=c(1,3,5), width=3))
gr
df=DataFrame(ir=ir1, score=rnorm(3))
df
df[1,1]
df$ir
df[1]
values(gr)=DataFrame(score = rnorm(3))
values(gr)
gr
gr$score
gr$score2=gr$score/3
gr2 <- GRanges(seqnames=c('chr1','chr2','chr1'), strand='*', ranges=IRanges(start=c(1,3,5), width=3))
findOverlaps(gr,gr2)
subsetByOverlaps(gr2,gr)
ir3<-makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
