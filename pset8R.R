library('DESeq2')
library('tximport')
files=c("C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458493\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458500\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458507\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458514\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458521\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458528\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458878\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458885\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458892\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\Pset8\\ERR458899\\abundance.tsv")

names(files)=c("ERR458493","ERR458500","ERR458507","ERR458514","ERR458521","ERR458528","ERR458878","ERR458885","ERR458892","ERR458899")

txdat = tximport(files, type ="kallisto", txOut=TRUE)

coldata=data.frame(condition=c("ERR458493","ERR458500","ERR458507","ERR458514","ERR458521","ERR458528","ERR458878","ERR458885","ERR458892","ERR458899"))

rownames(coldata) =names(files)

coldata

dds = DESeqDataSetFromTximport(txdat, colData=coldata, design=~ condition)

des=DESeq(dds)

res = results(des)
plotMA(res)

#sizeFactorEst = estimateSizeFactors(des)
#dispEst = estimateDispersions(des)
plotDispEsts(des, 1)

res$padj
res$pvalue
sum(res$padj<0.05, na.rm=TRUE)
sum(res$pvalue<0.05, na.rm=TRUE)


write.table(res,file="pset8_result.txt", sep="\t", quote=FALSE)
