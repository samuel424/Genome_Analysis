library("DESeq2")
directory <- "/home/samuel50/Genom_Analysis/htseq"

sampleFiles <- c('counts-BH-1.txt','counts-BH-2.txt','counts-BH-3.txt','counts-Serum-1.txt','counts-Serum-2.txt','counts-Serum-3.txt')
sampleCondition <- c('BH','BH','BH','Serum','Serum','Serum')
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)

sampleTable$condition <- factor(sampleTable$condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

dds <- DESeq(ddsHTSeq)

res <- results(dds)

summary(res)

plotMA(res)

library(ReportingTools)
report <- HTMLReport(shortName = 'Differential expression analysis BH vs Serum', title = 'Differential expression analysis BH vs Serum', reportDirectory = '.')
publish(dds, report, pvalueCutoff=0.05, make.plots = TRUE, factor = sampleTable$condition, reportDir = ".")
finish(report) 
