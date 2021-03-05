library(AnnotationHub)
source("http://www.bioconductor.org/biocLite.R")
biocLite("AnnotationHub")
biocLite(c("rtracklayer","AnnotationHub"))

library(GenomicRanges)
library(rtracklayer)
ahub=AnnotationHub()
ah<-subset(ahub, species=='Homo sapiens')
ah_1<-query(ah,'H3K4me3')
ah_2<-query(ah_1, 'Gm12878')
ah_2
#ah<-query(ahub, c('H3K4me3','H3K27me3','E003'))
ah_2$title
ah_2$dataprovider
gr1<-subset(ah_2, title == 'wgEncodeUwHistoneGm12878H3k4me3StdPkRep1.narrowPeak.gz')
gr1
ah_2<-query(ah,'refseq')
ah_2$genome
genome(gr1)
refseq<-ah_2[ah_2$genome=='hg19'&title=='Refseq Genes']
genome()

data_CpG<-query(ah, 'CpG Islands')
data_CpG[['AH5086']] #hg19
data_CpG[['AH5204']]
data_CpG[['AH5344']]
data_CpG[['AH5463']]
data_CpG<-data_CpG[['AH5086']]
summary(width(data_CpG))
seqinfo(data_CpG)
gaps(data_CpG)
sum(as.numeric(seqlengths(reduce(data_CpG))))

reduce_data<-reduce(data_CpG)
(coverage(reduce_data))[1:5]    
4925+3377+2327+2063+2459+2507+3157+1793+2075+2453+2287+2735+2443+1211+1577+1585+2983+3269+1015+1603+363+5083+1439+731
(coverage(reduce_data))[6:10]
(coverage(reduce_data))[11:15]
(coverage(reduce_data))[16:20]
(coverage(reduce_data))[21:25]
reduce_data

#split the data as the chromosome title
split_data<-split(reduce_data,seqnames(reduce_data))

#subset the data to autosome chr1-22
chromosome<-c(paste("chr",1:22,sep=""))
autosomoe_data<-split_data[chromosome]
#Q1:How many islands exists on the autosomes?
unlist(autosomoe_data)

#Q2:How many CpG Islands exists on chromosome 4.
autosomoe_data[4]

#Q3:How many bases does these regions cover?
ah_H3K4me_autosome_data<-subset(ah_H3K4me_data, seqnames %in% autosome)
sum(width(unlist(ah_H3K4me_autosome_data)))

#Q4: H3K27me3 histone modification for the H1 cell line from Epigenomics Roadmap
ah_H3K27me3<-query(ah, c('H3K4mo3','narrowPeak','E003'))
ah_H3K27me3_data <- ah_H3K27me3[['AH29892']]
ah_H3K27me3_autosome_data <- subset(ah_H3K27me3_data, seqnames %in% autosome)
ah_H3K27me3_autosome_data_mean <- mean(ah_H3K27me3_autosome_data$signalValue)
ah_H3K27me3_autosome_data_mean

#Q5:Bivalent regions are bound by both H3K4me3 and H3K27me3.
bivalent_data <- intersect(unlist(ah_H3K4me_autosome_data), unlist(ah_H3K27me3_autosome_data))
sum(width(reduce(bivalent_data)))

#Q6:how big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap 
# one or more CpG Islands?
CpG_bivalent_data <- findOverlaps(bivalent_data, unlist(autosome_CpG_data))
fraction_bi <- length(unique(queryHits(CpG_bivalent_data)))/length(bivalent_data)
fraction_bi

#Q7:How big a fraction (expressed as a number between 0 and 1) of the bases which
#are part of CpG Islands, are also bivalent marked
ov_CpG_bivalent <- intersect(bivalent_data, unlist(autosome_CpG_data))
fraction_CpG <- sum(width(reduce(ov_CpG_bivalent))/sum(width(unlist(autosome_CpG_data))))

#Q8: How many bases are bivalently marked within 10kb of CpG Islands?
autosome_CpG_data
CpG_10k <- resize(unlist(autosome_CpG_data), 
                  width=20000+ width(unlist(autosome_CpG_data)),
                  fix='center')
CpG_10k_bivalent <- intersect(CpG_10k, bivalent_data)
sum(width(CpG_10k_bivalent))

#Q9: How big a fraction (expressed as a number between 0 and 1) of the human genome 
#is contained in a CpG Island?
genome <- ah[['AH5081']]
genome <- keepSeqlevels(genome, c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                  'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                                  'chr20','chr21','chr22'))
genome_size <- sum(as.numeric(seqlengths(genome)))
sum(as.numeric(width(unlist(autosome_CpG_data))))/genome_size

#Q10:Compute an odds-ratio for the overlap of bivalent marks with CpG islands.
inOut = matrix(0, ncol=2, nrow=2)
colnames(inOut)=c('in','out')
rownames(inOut)=c('in','out')
inOut[1,1] = sum(width(intersect(bivalent_data,
                                 unlist(autosome_CpG_data),
                                 ignore.strand=True)))
inOut[1,2] = sum(width(setdiff(bivalent_data,
                               unlist(autosome_CpG_data),
                               ignore.strand=True)))
inOut[2,1] = sum(width(setdiff(unlist(autosome_CpG_data),
                               bivalent_data, ignore.strand=True)))
inOut[2,2] = genome_size - sum(inOut)

odd_ratio <- inOut[1,1]*inOut[2,2]/(inOut[1,2]*inOut[2,1])
odd_ratio