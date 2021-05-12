#To change chr in BAM files-->
samtools view -H lohith_seqmule.bam | sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | sed -e 's/SN:MT/SN:chrM/' | samtools reheader - lohith_seqmule.bam > test_lohith_seqmule.bam
samtools index .bam

R
BAMFiles <- list.files(pattern=".bam$")
segments <- read.table("nexterarapidcapture_exome_targetedregions_v1.2.bed",sep="\t",as.is=TRUE)

library(GenomeInfoDb)
seqlevels(gr) <- gsub("chr","",seqlevels(gr))
gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
X <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr,mode="paired")
X
resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
(resCNMOPS)
cnvs(resCNMOPS)[1:5]
plot(resCNMOPS,which=1)
segm <- as.data.frame(segmentation(resCNMOPS))
CNVs <- as.data.frame(cnvs(resCNMOPS))
CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(segm,file="segmentation.csv")
write.csv(CNVs,file="cnvs.csv")
write.csv(CNVRegions,file="cnvr.csv")
