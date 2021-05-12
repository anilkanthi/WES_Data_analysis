library(ExomeDepth)
data(exons.hg19)
print(head(exons.hg19))
data(exons.hg19)

my.counts <- getBamCounts(bed.frame = exons.hg19,bam.files = "sudhiksha.0_bwamem.sort.rmdup.readfiltered.bam",include.chr = FALSE,referenceFasta = "human_g1k_v37.fasta")

#OR (To include multiple samples)
my.counts <- getBamCounts(bed.frame = exons.hg19,bam.files = c("",""),include.chr = FALSE,referenceFasta = "/mnt/exome/Softwares/anil/human_g1k_v37.fasta")

ExomeCount.dafr <- as(my.counts[,  colnames(my.counts)],'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space), pattern = 'chr', replacement ='')
print(head(ExomeCount.dafr))
my.test <- ExomeCount.dafr$X62295713_S21_realigned.bam
my.ref.samples <- c('X62295699_S13_realigned.bam','X62295711_S20_realigned.bam','X62295724_S16_realigned.bam','X62295727_S17_realigned.bam','X62295794_S35_realigned.bam')
my.reference.set <- as.matrix(ExomeCount.dafr[, my.ref.samples])
my.choice <- select.reference.set (test.counts = my.test,reference.counts = my.reference.set,bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,n.bins.reduced = 10000)
print(my.choice[[1]])

my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])
my.reference.selected <- apply(X = my.matrix,MAR = 1,FUN = sum)
#OR
my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = TRUE])
my.reference.selected <- apply(X = my.matrix,MAR = 1,FUN = sum)

all.exons <- new('ExomeDepth',test = my.test,reference = my.reference.selected,formula = 'cbind(test, reference) ~ 1')

all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4, chromosome = ExomeCount.dafr$space, start = ExomeCount.dafr$start, end = ExomeCount.dafr$end, name = ExomeCount.dafr$names)
print(head(all.exons@CNV.calls))
output.file <- 'exome_calls.csv'
write.csv(file = output.file,x = all.exons@CNV.calls,row.names = FALSE)

#To Annotate
data(Conrad.hg19)
head(Conrad.hg19.common.CNVs)
levels(GenomicRanges::seqnames(Conrad.hg19.common.CNVs))
all.exons <- AnnotateExtra(x = all.exons, reference.annotation = Conrad.hg19.common.CNVs,min.overlap = 0.5,column.name = 'Conrad.hg19')
print(head(all.exons@CNV.calls))
exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),names = exons.hg19$name)
all.exons <- AnnotateExtra(x = all.exons,reference.annotation = exons.hg19.GRanges,min.overlap = 0.0001,column.name = 'exons.hg19')
all.exons@CNV.calls[3:6,]

plot (all.exons, sequence = '15', xlim = c(22318837 - 100000, 22483140 + 100000),count.threshold = 20,main = 'C15orf2 gene',cex.lab = 0.8,with.gene = TRUE)
#OR
plot (all.exons, sequence = '14', xlim = c(74951119 - 10000, 74953139 + 10000),count.threshold = 20,main = 'NPC2 gene',cex.lab = 0.8,with.gene = TRUE, col = 'black')
dev.off()
