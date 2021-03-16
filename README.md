
**Load functions**
'''
source('BHOT.R')
'''

**Load files (rccFiles is a vector of file paths)**
```
ns.data <- parseRCC(rccFiles)

ns.counts <- ns.data$counts
ns.att <- ns.data$attributes
```

**Comparisons for differential expressoin analysis**
```
exp_design <- data.frame(ID=ref.df$ID,
                         rejection=ifelse(ref.df$Dx %in% c("Chronic Active ABMR", "Active ABMR", "Chronic active TCMR", "Acute TCMR"), 1, 0),
                         abmr=ifelse(ref.df$Dx %in% c("Chronic Active ABMR", "Active ABMR"), 1, 0),
                         tcmr=ifelse(ref.df$Dx %in% c("Chronic active TCMR", "Acute TCMR"), 1, 0),
                         g_score=ifelse(ref.df$g > 0, 1, 0),
                         cg_score=ifelse(ref.df$cg > 0, 1, 0),
                         ptc_score=ifelse(ref.df$ptc > 0, 1, 0),
                         i_score=ifelse(ref.df$i > 1, 1, 0),
                         t_score=ifelse(ref.df$t > 1, 1, 0))
factor_vars <- c("rejection", "abmr", "tcmr", "g_score", "cg_score", "ptc_score", "i_score", "t_score")
exp_design[,factor_vars] <- lapply(exp_design[,factor_vars], as.factor)
```

**QC**

QC attributes: imaging and binding density
```
att.qc <- qcAttributes(ns.att)
```
QC counts
```
count.qc <- qcCounts(ns.counts)
```
per sample endogenous gene stats
```
sampleStats(ns.counts)
```

**Select housekeeping genes for normalization**
Include group labels (optional) to perform NB regression between hk_genes ~ biological conditions/outcomes of interest
```
group_labels <- ifelse(colnames(ns.counts[,!colnames(ns.counts) %in% c("CodeClass", "Name", "Accession")]) %in% exp_design[exp_design$abmr==1,"ID"], 1, 0)
group_labels <- ifelse(colnames(ns.counts[,!colnames(ns.counts) %in% c("CodeClass", "Name", "Accession")]) %in% exp_design[exp_design$tcmr==1,"ID"], 1, 0)

hk_qc <- hkQC(ns.counts, group_labels)
```

Remove unstable housekeeping genes
```
ns.counts <- ns.counts[!ns.counts$Name %in% hk_qc$hk_genes_remove,]
rownames(ns.counts) <- ns.counts$Name
```

**Normalization**

create expression set (count column names must match phenoData rows)
```
mat <- as.matrix(ns.counts[,-c(1:3)])

pdat <- data.frame(ID=colnames(mat))
pdat <- merge(pdat, exp_design, by="ID", all.x=TRUE, sort=FALSE)
rownames(pdat) <- pdat$ID
pdat <- AnnotatedDataFrame(data=pdat)

ns.anno <- ns.counts[,c(1:3)]
rownames(ns.anno) <- ns.anno$Name
fdat <- AnnotatedDataFrame(data=ns.anno)

eset <- newSeqExpressionSet(mat, phenoData=pdat, featureData=fdat)
```

nSolver normalization
```
ns_norm <- NS_norm(eset, background_correction=FALSE, take_log=FALSE)
```

RUV normalization [Bhattacharya et al. 2020](https://doi.org/10.1093/bib/bbaa163)
```
ruv_norm <- RUV_norm(eset, k=1)
```
raw counts
```
counts(ruv_norm$eset)
```
normalized counts
```
assay(ruv_norm$vsd)
```

RCR normalization
```
RCR_norm(ns.counts)
```

**RLE plots**
raw counts
```
raw_counts <- as.matrix(ns.counts[ns.counts$CodeClass=="Endogenous",-c(1:3)])
plotRLE(raw_counts, is_logged=FALSE, main="Raw counts")
```
RUV normalized counts
```
endo_genes <- ns.counts[ns.counts$CodeClass=="Endogenous","Name"]
ruv.norm.endo <- assay(ruv_norm$vsd)[rownames(assay(ruv_norm$vsd)) %in% endo_genes,]
plotRLE(ruv.norm.endo, is_logged=TRUE, main="RUV norm")
```

**Differential Expression Analysis (negative binomial model)**

Input raw counts to DESeq2
```
ruv_endo <- counts(ruv_norm$eset)[rownames(counts(ruv_norm$eset)) %in% endo_genes,]
```

ABMR v. all
```
ruv.abmr <- DESeqDataSetFromMatrix(countData = ruv_endo, colData = pData(ruv_norm$eset), design = ~ W_1 + abmr)
ruv.abmr <- DESeq(ruv.abmr, test="Wald")
ruv.abmr.sig <- data.frame(results(ruv.abmr))
hist(ruv.abmr.sig$pvalue)
```

TCMR v. all
```
ruv.tcmr <- DESeqDataSetFromMatrix(countData = ruv_endo, colData = pData(ruv_norm$eset), design = ~ W_1 + tcmr)
ruv.tcmr <- DESeq(ruv.tcmr, test="Wald")
ruv.tcmr.sig <- data.frame(results(ruv.tcmr))
hist(ruv.tcmr.sig$pvalue)
```

Volcano Plot
```
plotVolcano(ruv.abmr.sig, p_cutoff=0.05, plot_title="ABMR v. all")
```
