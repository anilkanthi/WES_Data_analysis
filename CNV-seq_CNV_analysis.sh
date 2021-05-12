#For affected patient-->
samtools view -F 4 ../test_rajeshwari_seqmule.bam | perl -lane'print "$F[2]\t$F[3]"'> rajeshwari.hits

#For control -->
samtools view -F 4 /mnt/exome/WES_data/Backup_21_12_2016/Diagnostic_exome/2096734/Lava/62271272_S50.bsrt.bam | perl -lane'print "$F[2]\t$F[3]"'> ref.hits

./cnv-seq.pl --test rajeshwari.hits --ref ref.hits --genome human --annotate --minimum-windows-required 4 --bigger-window 1.5 --p-value 0.001 --log2-threshold 0.6 --Rexe /usr/bin/R

#To plot CNVs -->

R
library(cnv)
data <- read.delim("rajeshwari.hits-vs-ref.hits.log2-0.6.pvalue-0.001.minw-4.cnv")
cnv.summary(data)
plot.cnv(data, CNV=4, upstream=4e+6, downstream=4e+6)
ggsave("sample.pdf")
