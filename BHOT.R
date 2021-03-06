
## Load dependencies
library_list <- c("DESeq2", "dplyr", "ggplot2", "ggrepel", "knitr", "MASS", "RUVSeq", "RCRnorm", "NormqPCR")
missing_libraries <- library_list[!(library_list %in% installed.packages()[,"Package"])]

#BiocManager::install("RUVSeq")
#BiocManager::install("NormqPCR")

if(length(missing_libraries) > 0) {
	cat("The following libraries are missing:\n", missing_libraries)
}

## Load all packages
lapply(library_list, library, character.only=TRUE)

##------------------------------
## geometric mean
geoMean <- function(x) {
    x[x<1] <-1;
    exp(mean(log(na.omit(x))));
}

##------------------------------
## parse a single RCC file
importRCC <- function(rccFile) {
    
    rccFile <- readLines(rccFile)
    rccFile <- data.frame(line = rccFile, tag = rep("", length(rccFile)), stringsAsFactors = FALSE)

    ## lane attributes contain Fov, BD
    tags <- c("Header", "Sample_Attributes", "Lane_Attributes", "Code_Summary", "Messages")

    ## assign tag to each line
    for (i in 1:length(tags)) {
        tag_start <- which(rccFile$line == paste0("<", tags[i], ">"))
        tag_end <- which(rccFile$line == paste0("</", tags[i], ">"))
        rccFile$tag[tag_start:tag_end] <- tags[i]
    }

    rccFile <- rccFile[rccFile$tag != "",] ## exclude lines with empty tags
    rccFile <- rccFile[!grepl("<", rccFile$line),] ## remove lines with tag names start
    rccFile <- rccFile[!grepl("</", rccFile$line),] ## remove linear with tag names end

    ## get sample ID
    samp_attr <- rccFile[rccFile$tag=="Sample_Attributes",]
    id_line <- samp_attr[grep("^ID", samp_attr$line),]
    sampID <- do.call('rbind', strsplit(as.character(id_line$line), ','))[2]

    all_attr <- rccFile[rccFile$tag!="Code_Summary",]

    parse_attr <- data.frame(do.call('rbind', strsplit(as.character(all_attr$line), ',')), stringsAsFactors=FALSE)
    colnames(parse_attr) <- c("variable", sampID)

    ID_lnum <- grep("^ID", parse_attr[,1])
    parse_attr[,1][ID_lnum] <- c("sampID", "laneID")

    ## get raw counts: split columns and remove header row
    parse_counts <- rccFile[rccFile$tag=="Code_Summary",]
    parse_counts <- data.frame(do.call('rbind', strsplit(as.character(parse_counts$line), ',', fixed=TRUE)), stringsAsFactors=FALSE)
    parse_counts <- parse_counts[grep("^Code", parse_counts[,1], invert=TRUE),]
    colnames(parse_counts) <- c("CodeClass", "Name", "Accession", sampID)
    parse_counts[,sampID] <- sapply(parse_counts[,sampID], as.numeric)
    
    ## return a nested list of raw counts and sample attributes
    rcc_summary <- list(counts=parse_counts, 
                        attributes=parse_attr)
    return(rcc_summary)
}

## call importRCC function to import and parse data, combine all into count table and attribute table
parseRCC <- function(rccFiles) {
    
    count_attr_list <- lapply(rccFiles, importRCC)
    
    allCounts <- do.call(rbind, count_attr_list)[,1]
    allAttributes <- do.call(rbind, count_attr_list)[,2]
    
    ## collapse raw counts to single table
    countTable <- Reduce(function(x,y) merge(x, y, all=TRUE, by=c("CodeClass", "Name", "Accession")), allCounts)
    attTable <- Reduce(function(x,y) merge(x, y, all=TRUE, by=c("variable")), allAttributes)
    
    ## Reporter Library File is used during image processing to assign target identities to barcodes, based on CodeSet name
    cat('----------------------------------', '\n');
    cat("CodeSet Name:", unique(as.character(attTable[attTable$variable=="GeneRLF",][-1])));
    cat('\n');
    cat("Successfully imported data for\n");
    cat(dim(countTable[,-c(1:3)])[2], "samples", '\n');
    cat(nrow(countTable[countTable$CodeClass=="Endogenous",]), "endogenous genes", '\n');
    cat('----------------------------------\n');
    
    return(list(counts=countTable,
                attributes=attTable))
}

##------------------------------
## QC RCC attributes: output binding density plot and per sample qc table
qcAttributes <- function(attributes_table) {
    
    ## Imaging/Field of View measure of % of requested fields of view successfully scanned in each cartridge lane
    qcTab <- attributes_table[attributes_table$variable %in% c("BindingDensity", "FovCount", "FovCounted"),]
    
    pct_fov <- as.numeric(qcTab[qcTab$variable=="FovCounted",][-1]) / as.numeric(qcTab[qcTab$variable=="FovCount",][-1])
    qcTab <- rbind(qcTab, c("pct_fov", pct_fov))
    
    if (any(qcTab[qcTab$variable=="pct_fov",][-1] < 0.75)) {
        cat("Imaging flag (FoV < 0.75):\n")
        cat(names(qcTab[qcTab$variable=="pct_fov",][-1][which(qcTab[qcTab$variable=="pct_fov",][-1] < 0.75)]))
    } else {
        cat(">>All samples pass Imaging QC\n")
    }
    
    ## Binding density
    if (any(qcTab[qcTab$variable=="BindingDensity",][-1] < 0.05 | qcTab[qcTab$variable=="BindingDensity",][-1] > 2.25)) {
        cat("Binding Density flag (0.05 < BD < 2.25):\n")
        cat(names(qcTab[qcTab$variable=="BindingDensity",][-1][which(qcTab[qcTab$variable=="BindingDensity",][-1] < 0.05)]))
        cat(names(qcTab[qcTab$variable=="BindingDensity",][-1][which(qcTab[qcTab$variable=="BindingDensity",][-1] > 2.25)]))
    } else {
        cat(">>All samples pass Binding Density QC\n")
    }
    
    ## per sample binding density plot
    bd <- data.frame(t(qcTab[qcTab$variable=="BindingDensity",][-1]), stringsAsFactors=FALSE)
    colnames(bd) <- "bd"
    bd$bd <- as.numeric(bd$bd)
    bd$id <- rownames(bd)
    bd <- bd[order(bd$bd, decreasing=TRUE),]
    bd$id <- factor(bd$id, levels=as.character(bd[order(bd$bd, decreasing=TRUE),"id"]))
    
    bd_plot <- ggplot(bd, aes(x=id, y=bd)) +
        geom_bar(stat='identity', show.legend = FALSE) +
        xlab("") + ylab("binding density") + 
        geom_hline(yintercept=2.25, col="dodgerblue", linetype="dashed") +
        geom_hline(yintercept=0.05, col="tomato", linetype="dashed") +
        theme(panel.grid.major=element_blank(),
              panel.background=element_blank(), 
              axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="black", fill=NA, size=1),
              plot.title = element_text(size=12, face = "bold"),
              axis.text=element_text(size=6, angle=270, family="sans", colour="black"), 
              axis.title.x=element_text(size=8, family="sans", colour="black"), 
              axis.title.y=element_text(size=8, family="sans", colour="black"))
    
    return(list(density_plot=bd_plot, 
                qcTable=qcTab))
}

##------------------------------
## Per sample control probes quality control checks
## TODO: output QC summary table+plots
qcCounts <- function(count_table) {
    
    cat('----------------------------------\n');
    cat('Quality Control metrics:\n');
    cat('----------------------------------\n');
    
    samp_ids <- colnames(count_table[,!colnames(count_table) %in% c("CodeClass", "Name", "Accession")])
    endo_tab <- count_table[count_table$CodeClass=="Endogenous",c("Name", samp_ids)]
    rownames(count_table) <- count_table$Name

    cat('\nOverall Assay Efficiency:\n')
    cat("(low efficiency: geometric mean of positive controls > 3x the geometric mean of endogenous genes)\n")
    
    geo_means <- data.frame(ID=samp_ids, stringsAsFactors=FALSE)
    geo_means$hk_mean <- rep("NA", nrow(geo_means))
    geo_means$pos_mean <- rep("NA", nrow(geo_means))
    geo_means$endo_mean <- rep("NA", nrow(geo_means))
    geo_means$endo_sd3 <- rep("NA", nrow(geo_means))
    
    for (i in 1:length(samp_ids) ) {
        sampID <- samp_ids[i]
        hk_counts <- count_table[count_table$CodeClass=="Housekeeping",][sampID][,1]
        pos_counts <- count_table[count_table$CodeClass=="Positive",][sampID][,1]
        endo_counts <- count_table[count_table$CodeClass=="Endogenous",][sampID][,1]
        geo_means[geo_means$ID==sampID,"hk_mean"] <- round(geoMean(hk_counts), 3) # geometric mean of HK genes
        geo_means[geo_means$ID==sampID,"pos_mean"] <- round(geoMean(pos_counts), 3) # geometric mean of positive controls
        geo_means[geo_means$ID==sampID,"endo_mean"] <- round(geoMean(endo_counts), 3) # geometric mean of endogenous genes
        geo_means[geo_means$ID==sampID,"endo_sd3"] <- round(3*geoMean(endo_counts), 3)
    }
    
    if (any(geo_means$pos_mean < geo_means$sd3)) {
        cat("Assay Efficiency is low for the following samples:\n");
        cat(geo_means[which(geo_means$pos_mean < geo_means$endo_sd3),"ID"], sep="\n")
    } else {
        cat(">>OK for all samples\n");
    }
    
    cat('----------------------------------\n',
        'Positive Control Linearity:\n',
        '(correlation between the observed counts and concentrations of Positive ERCC probes)\n')
    
    pos_tab <- count_table[count_table$CodeClass=="Positive",c("Name", samp_ids)]
    
    ## POS_F is considered below the limit of detection: remove when calculating linearity
    pos_tab <- pos_tab[grep("POS_F", pos_tab$Name, invert=TRUE),]
    
    if (!all(grepl("[[:digit:]]", pos_tab$Name))) {
        stop("Positive controls must have concentrations, eg. POS_A(128)")
    }
    
    ## extract pos control concentration
    pos_tab$concentration <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", pos_tab$Name))
    
    ## calculate positive control linearity for each sample
    positiveLinearityQC <- apply(pos_tab[,samp_ids], 2, function(x) {
        round(summary(lm(x ~ pos_tab$concentration))$r.squared, 3)
    })
    
    if (any(positiveLinearityQC < 0.95)) {
        cat(">>Linearity performance is low for:\n",
            names(which(positiveLinearityQC < 0.95)))
    } else {
        cat(">>Positive control linearity is > 0.95 for all samples\n")
    }
    
    tab <- pos_tab
    tab$Name <- NULL
    ## supress annoying data table warning message
    suppressWarnings({
    tab <- data.table::melt(tab, measure.vars=samp_ids, verbose=FALSE)
    colnames(tab) <- c("concentration", "ID", "count")
    })
    
    linearity_plot <- ggplot(tab, aes(x=concentration, y=count, fill=ID)) +
        geom_point(shape = 21, size=3) + 
        xlab("raw counts") + ylab("concentration (fM)") +
        geom_smooth(method = "lm", fullrange=TRUE, se=TRUE, size=1, 
                    color="slategray", formula = y ~ x, linetype="dashed")
    
    cat('----------------------------------\n',
        'Positive Control Limit of Detection:\n',
        '(POS_E control probe counts > mean negative control + 2*SD)\n')
    
    pos_e = count_table[grep("POS_E", count_table$Name),]
    neg_raw <- count_table[count_table$CodeClass=="Negative",]
    
    ncgMean = apply(neg_raw[,samp_ids], 2, mean)
    ncgSD = apply(neg_raw[,samp_ids], 2, sd)
    lod = ncgMean + 2*ncgSD
    pos_e_counts = pos_e[,samp_ids]
    
    ## POS_E counts should be > lod
    if (any(pos_e_counts < lod)) {
        cat(">>Postive control limit of detection is low for the following samples:\n",
            names(pos_e_counts[which(pos_e_counts < lod)]), '\n')
    } else {
        cat(">>Positive control E counts are > LoD for all samples\n")
    }
    
    ## % endogenous genes above LoD per sample
    pctLOD <- data.frame(ID=samp_ids, stringsAsFactors=FALSE)
    pctLOD$pctLOD <- rep(NA, nrow(pctLOD))
    
    for (i in 1:length(samp_ids)) {
        i_counts <- endo_tab[,samp_ids][i]
        i_lod <- lod[i]
        pctLOD[i,"pctLOD"] <- round(length(i_counts[which(i_counts[,1] > i_lod),]) / nrow(i_counts) * 100, 3)
    }
    
    cat('----------------------------------\n');
    print(kable(pctLOD[pctLOD$pctLOD < 20,], caption="<20% endogenous genes above the Limit of Detection"), 
          warnings=FALSE, row.names=FALSE)
    
    return(linearity_plot=linearity_plot)
    
}

##------------------------------
## QC table and plot functions
sampleStats <- function(tab) {
    
    samp_ids <- colnames(tab[,!colnames(tab) %in% c("CodeClass", "Name", "Accession")])
    tab <- tab[tab$CodeClass=="Endogenous",]
    tab <- tab[,!colnames(tab) %in% c("CodeClass", "Name", "Accession")]
    
    endoMissing <- apply(tab, MARGIN = 2, FUN = function(x) { sum(x <= 0, na.rm=TRUE) / length(x) *100; });
    endoMedian <- apply(tab, 2, FUN=median, na.rm=TRUE);
    endoMean <- apply(tab, 2, FUN=mean, na.rm=TRUE);
    endoSD <- apply(tab, 2, FUN=sd, na.rm=TRUE);
    
    stats_df <- data.frame(row.names=names(tab[,samp_ids]),
                              Median=endoMedian, 
                              Mean=endoMean, 
                              SD=endoSD,
                              Missing=endoMissing);
    return(stats_df)
}

##------------------
## batch effects/outliers
#ns.norm$batch.effects #output form nanostringnorm

##------------------
## Housekeeping gene QC and selection
## candidate housekeeping genes have a high mean and very low variance

hkQC <- function(raw_counts, group_labels) {
    
    samp_ids <- colnames(ns.counts[,!colnames(ns.counts) %in% c("CodeClass", "Name", "Accession")])
    
    hk_counts <- ns.counts[ns.counts$CodeClass=="Housekeeping",c("Name", samp_ids)]
    rownames(hk_counts) <- hk_counts$Name
    hk_counts$Name <- NULL
    
    ## TODO: exlcude HK genes below limit of detection
    
    ## caclulate CV for each HK gene
    #hk_cv = apply(hk_counts, 1, sd, na.rm=TRUE) / apply(hk_counts, 1, mean, na.rm=TRUE) * 100
    
    ## if available, NB regression between hk_genes ~ biological conditions/outcomes of interest
    if (length(group_labels)>0) { 
        
        hk_pvals = data.frame(gene=rownames(hk_counts), pval=rep(NA, nrow(hk_counts)))
        
        for (i in 1:nrow(hk_counts)){
            hk_nb <- MASS::glm.nb(as.numeric(hk_counts[i,]) ~ as.factor(group_labels))
            hk_pvals$'pval'[i] = coef(summary(hk_nb))[2,4]
        }
    }
    
    ## 01: geNorm method, Vandesompele et al. 2002
    ## rank of genes from most to least stable (2 most stable cannot be ranked)
    ## calculate average pairwise variation (sd of the log transformed expression ratios) with all other genes to determine
    ## gene-stability value M: the average pairwise variation of a particular gene with all other control genes. 
    ## Genes with the lowest M values have the most stable expression
    
    hk_gn <- NormqPCR::selectHKs(t(hk_counts), method = "geNorm", Symbols=colnames(t(hk_counts)), na.rm=TRUE, minNrHK=2, log=FALSE)
    
    hk_rank <- data.frame(rank=names(hk_gn$rank), gene=hk_gn$rank, stringsAsFactors=FALSE)
    
    hk_m <- data.frame(n_genes=names(hk_gn$meanM), meanM=hk_gn$meanM, stringsAsFactors=FALSE) #stability values
    hk_m <- hk_m[order(hk_m$meanM, decreasing=TRUE),]
    hk_m$n_genes <- factor(hk_m$n_genes, levels=hk_m$n_genes)
    
    ## rank least to most stable
    gp_hk <- ggplot(data=hk_m, aes(x=n_genes, y=meanM)) + geom_point(col="dodgerblue", size=2, pch=19) + 
        ggtitle("Mean HK gene stability") +  xlab("# HK genes") +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_line(color="gray"),
              panel.background=element_blank(), axis.line=element_line(color="white"),
              legend.title=element_blank(),
              axis.text=element_text(size=10), 
              axis.text.x=element_text(face="bold", size=10, angle=0),
              axis.title.x=element_text(face="bold", size=10, angle=0), 
              axis.title.y=element_text(size=10),  
              panel.border=element_rect(colour="black", fill=NA, size=1))
    
    ## select HK genes to be used in normalization
    ## manually or automatically (M<0.7): experimental value shown to eliminate most variable genes in microarray data
    num_hk_genes=as.numeric(as.character(hk_m[hk_m$meanM<0.7,][1,"n_genes"]))
    #num_hk_genes=5
    
    hk_genes_keep <- hk_rank[1:num_hk_genes,"gene"]
    hk_genes_remove <- rownames(hk_counts)[!rownames(hk_counts) %in% hk_genes_keep]
    
    ## 02: NormFinder method, Andersen et al. 2004; use if groups available
    #hk_nf <- NormqPCR::selectHKs(t(hk_counts), method = "NormFinder", group=group_labels, Symbols=colnames(mat), minNrHK = 2, log = FALSE)
    ## ideally: start with 3 control genes > stepwise inclusion of more control genes until the (n + 1)th gene has 
    ## no significant contribution to the newly calculated normalization factor
    
    ## output staibilty plot, regression pvals, hk genes to exclude (if any)
    hk_list <- list(hk_plot=gp_hk, nb_pvals=hk_pvals, exclude=hk_genes_remove)
    return(hk_list[lengths(hk_list) != 0])
    
}

##------------------
## Normalize data
## normalization: geNorm algorithm (Vandesompele 2002) to choose best genes for normalization 
## (exc. ref genes to minimize pairwise variation)

## OPTION 1: nSolver
NS_norm <- function(eset, background_correction, take_log) {
    
    ## 01 CodeCount normalization (adjust each sample based on relative value to all samples)
    #https://github.com/chrisrgc/NanoStringNorm/blob/master/R/code.count.normalization.R
    rownames(ns.counts) <- ns.counts$Name
    
    posTab <- ns.counts[ns.counts$CodeClass=="Positive",!colnames(ns.counts) %in% c("CodeClass", "Name", "Accession")]
    ## exclude POS_F? (considered below lod)
    
    ## arithmetic mean of the geometric means of positive controls
    pos.sample <- apply(posTab, MARGIN=2, FUN=geoMean);
    pos.norm.factors <- mean(pos.sample) / pos.sample;
    
    if ( any(pos.norm.factors < 0.3) | any(pos.norm.factors > 3) ) {
        cat("Normalization factors out of range:\n");
        cat(names(pos.norm.factors[pos.norm.factors<0.3 | pos.norm.factors>3]))
    }
    
    # multiply normalization factor by raw counts
    tab <- ns.counts[,!colnames(ns.counts) %in% c("CodeClass", "Name", "Accession")]
    rownames(tab) <- ns.counts$Name
    x <- t(apply(tab, MARGIN = 1, FUN = '*', pos.norm.factors));
    
    ## 02 Background correction (subtract background from each sample)
    ## recommended for experiments in which low expressing targets common
    ## mean is least conservative; max and mean.2sd are most robust to FP
    #https://github.com/chrisrgc/NanoStringNorm/blob/master/R/background.normalization.R
    
    if (background_correction == TRUE) {
        x_neg <- x[rownames(x) %in% ns.anno[ns.anno$CodeClass=="Negative","Name"],]
        background='mean.2sd'
        
        ## calculate max or mean of negative control probes
        background.level <- if( background != 'none') {
            if (background == 'mean.2sd') {
                mean.plus.2sd <- function(x) mean(x, na.rm=TRUE) + 2 * sd(x, na.rm=TRUE);
                background.level <- apply(x_neg, MARGIN=2, FUN=mean.plus.2sd)
            } else if (background == 'max') {
                background.level <- apply(x_neg, MARGIN=2, FUN=max, na.rm=TRUE);
            }
        }
        
        background.zscores <- data.frame(background.zscore = (background.level - mean(background.level)) / sd(background.level));
        
        if ( any(background.zscores > 3)) {
            cat("\nBackground score out of range\n")
        }
        
        ## subtract background from all genes
        x <- t(apply(x, MARGIN=1, FUN='-', background.level))
        x[x<0] <- 0 ## convert negative values to zero
        
    } else {
        cat("\nBackground correction skipped\n")
        x <- x
    }
    
    ## 03 SampleContent (normalize to HK genes to account for sample or RNA content ie. pipetting fluctuations)
    ## Normalize by substracting geometric mean of housekeeping genes from each endogenous gene
    #https://github.com/chrisrgc/NanoStringNorm/blob/master/R/sample.content.normalization.R
    ## another option is to take top endo or HK genes with lowest CV (low.cv.geo.mean)
    
    hk_genes=ns.anno[ns.anno$CodeClass=="Housekeeping","Name"]
    
    # calculate the normalization factor: arithmetic mean of the geometric means of positive controls
    rna.content <- apply(x[rownames(x) %in% hk_genes,], MARGIN=2, FUN=geoMean);
    hk.norm.factor <- mean(rna.content) / rna.content
    
    ## flag outlier samples
    rna.content.zscore <- data.frame(rna.zscore = (rna.content - mean(rna.content)) / sd(rna.content));
    
    if (any(abs(rna.content.zscore) > 3)) {
        cat('SampleContent: The following samples have sample/rna content greater than \n\t3 standard deviations from the mean:\n');
        print(signif(subset(rna.content.zscore, abs(rna.content.zscore) > 3),3));
        cat('\n');
    }
    
    ## adjust data based on normalization factors
    x.norm <- t(apply(x, MARGIN = 1, FUN = '*', hk.norm.factor));
    
    if (take_log == TRUE) {
        cat('Counts are log transformed')
        x.norm <- log2(x.norm);
    }
    
    ## TODO:additional norm
    #https://github.com/cran/NanoStringNorm/blob/master/R/other.normalization.quantile.R
    #https://github.com/cran/NanoStringNorm/blob/master/R/other.normalization.vsn.R
    
    return(x.norm)
}

## OPTION 2: RUV
RUV_norm <- function(eset, k) {
    
    cat("k=", k, '\n');
    ## upper quantile normalization (Bullard 2010)
    eset <- betweenLaneNormalization(eset, which="upper")
    
    ## RUVg using housekeeping genes
    hk_idx <- rownames(eset)[fData(eset)$CodeClass == "Housekeeping"]
    ## estimated factors of unwanted variation (W)
    eset <- RUVg(eset, hk_idx, k=k)
    
    dds <- DESeqDataSetFromMatrix(counts(eset), colData=pData(eset), design=~1)
    #rowData(dds) <- fData(eset)
    
    ## estimate size factors
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersionsGeneEst(dds)
    
    cts <- counts(dds, normalized=TRUE)
    disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
    mcols(dds)$dispGeneEst <- disp
    dds <- estimateDispersionsFit(dds, fitType="mean")
    
    ## variance stabilizing log transformation
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    mat <- assay(vsd)
    
    ## remove unwanted variation as estimated by RUVg
    covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))), drop=FALSE])
    mat <- removeBatchEffect(mat, covariates=covars)
    assay(vsd) <- mat
    
    ## return expression set with unwanted variation added and log, normalized counts
    return(list(eset = eset, 
                vsd = vsd))
}

## OPTION 3: RCR
RCR_norm <- function(counts) {
    
    ns.anno <- ns.counts[,c(1:3)];
    rownames(ns.anno) <- ns.anno$Name
    
    ns.data <- ns.counts[,-c(1:3)];
    rownames(ns.data) <- ns.anno$Name
    
    ## input is a list of data frames for each code class
    ns.data.rcr <- list(pos_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Positive","Name"],],
                        neg_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Negative","Name"],],
                        hk_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Housekeeping","Name"],],
                        reg_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Endogenous","Name"],])
    
    #dput(as.numeric(gsub(".*\\((.*)\\).*", "\\1", ns.anno[ns.anno$CodeClass=="Positive","Name"])))
    
    rcr.norm <- RCRnorm(ns.data.rcr, pos_conc=log10(c(128, 32, 8, 2, 0.5, 0.125)),
                        fast_method=FALSE, iter=10000, warmup=5000)
    
    return(rcr.norm)
    
}

##------------------
## RLE (log2(expression/median across all assays))
plotRLE <- function(counts, is_logged, main) {
    
    x_endo <- counts[rownames(counts) %in% ns.counts[ns.counts$CodeClass=="Endogenous","Name"],]
    
    endo_median <- apply(x_endo, 1, median) #per gene median
    
    if (is_logged == TRUE) {
        endo_rle = x_endo - endo_median
    } else {
        endo_rle <- log2(x_endo/endo_median)
    }
    
    suppressWarnings({
        endo_melt <- data.table::melt(endo_rle, measure.vars=c(1:dim(endo_rle)[2]))
    })
    endo_melt <- na.omit(endo_melt)
    
    ggplot(endo_melt, aes(x=reorder(Var2, value, median, order=TRUE), y=value, fill=Var2)) + 
        ggtitle(main) + ylab("log2(counts/median)") +
        geom_hline(yintercept=0, linetype="dashed", color="firebrick") +
        geom_boxplot(outlier.size=1, outlier.shape=21, outlier.fill="whitesmoke") +
        theme(legend.position="none",
              axis.text.x=element_text(size=6, angle=90),
              axis.title.x=element_blank())
}

##--------------------------------------------------------------------------------
## Differential expression analysis:
## TODO: mean-diff plot (with and w/o control genes), pval histogram

## 1 Mixture model: stats4::MLE()
## 2 simplified model: MASS::glm.nb()
## 3 loglinear model: lm()

## feature selection: 10-50 genes

## Volcano plot
plotVolcano <- function(df.fit, p_cutoff, plot_title) {
    
    df.fit$gene <- rownames(df.fit)
    
    df.fit$logP <- -log10(df.fit$pvalue)
    df.fit$logPadj <- -log10(df.fit$padj)
    
    #df.fit <- plyr::mutate(df.fit, sig=ifelse(df.fit$padj < p_cutoff & abs(df.fit$log2FoldChange)>0.5, "sig", "not"))
    df.fit <- plyr::mutate(df.fit, sig=ifelse(df.fit$padj < p_cutoff, "sig", "not"))
    df.fit <- df.fit[order(df.fit$pvalue, decreasing=FALSE),]
    
    ## label top genes by pvalue or logFC
    gp_volcano <- ggplot() + ylab("-log10(pval)") + xlab("log2FC") + ggtitle(plot_title) +
        #geom_point(data=df.fit, aes( x=log2FoldChange, y=logP), colour="slategray", size=3, alpha=0.7) +
        geom_point(data=df.fit, aes(x=log2FoldChange, y=logP, color=sig), shape=19, size=3, alpha=0.7) +
        geom_text_repel(data=head(df.fit, 20), aes(x=log2FoldChange, y=logP, label=gene), colour="gray", size=3) +
        scale_color_manual(values=c("dodgerblue", "salmon")) + 
        geom_vline(xintercept=0, linetype="dashed", size=0.4) +
        geom_vline(xintercept=1, linetype="dashed", size=0.2) + geom_vline(xintercept=-1, size=0.2, linetype="dashed") +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.background=element_blank(), axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="black", fill=NA, size=1),
              plot.title = element_text(size=12, face = "bold"), legend.position="none",
              axis.text=element_text(size=12, family="sans", colour="black"), 
              axis.title.x=element_text(size=12, family="sans", colour="black"), 
              axis.title.y=element_text(size=12, family="sans", colour="black"))
    
    df.fit[df.fit$gene %in% c("CCL4", "CXCL11", "CXCL10", "PLA1A", "GNLY", "ROBO4", "FGFBP2",
                              "SH2D1B", "CD160", "DARC", "ROBO4", "CDH5", "CDH13", "SOST"),] #ABMR
    
    df.fit[df.fit$gene %in% c("ADAMDEC1", "ANKRD22", "CTLA4", "IFNG", "ICOS", "BTLA", 
                              "CD96", "LAG3", "SIRPG", "LCP2", "DUSP2", "CD8A"),] #TCMR
    
    return(gp_volcano)
    
}

##------------------
## TODO: QC plots
## create plotting function to call each plot individually OR to generate r markdown report


