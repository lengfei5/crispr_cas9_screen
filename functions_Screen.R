Import_Request = function()
{
  if(Import_1st_Request)
  {
    aa = read.table(paste('../output_5_RKi0MM_In_Ab.txt', sep=''), header = TRUE, sep='\t')
    dim(aa)
    
    jj = grep('_In', colnames(aa))
    kk = grep('_Ab', colnames(aa))
    
    aa = data.frame(aa[, 1], aa[,kk], aa[,jj], stringsAsFactors = FALSE)
    colnames(aa)[1] = 'guidename'
    
    Changecolnames = function(x)
    {
      #x = colnames(aa)[2] 
      test = unlist(strsplit(as.character(x), '_'))
      return(paste(c(paste0('F', test[2]), test[3]), collapse = '_'))
    }
    names = sapply(colnames(aa)[-1], Changecolnames)
    
    colnames(aa)[-1] = names
    
    samples = paste0('F', c(1, 2, 4, 5, 7, 8, 28, 29, 31, 10, 11, 13, 19, 20, 22, 14, 16, 17, 32, 34, 35, 23, 25, 26, 37, 38, 40))
    kk = match(paste(samples, 'Ab', sep='_'), colnames(aa))
    jj = match(paste(samples, 'In', sep='_'), colnames(aa))
    
    bb = data.frame(aa[, 1], aa[, c(kk, jj)], stringsAsFactors = FALSE)
    colnames(bb)[1] = 'sgRNA'
    
    find.gene.target = function(x)
    {
      test = unlist(strsplit(as.character(x), '_'))
      return(as.character(test[1]))
    }
    gene = sapply(bb$sgRNA, find.gene.target)
    
    bb = data.frame(bb$sgRNA, gene, bb[, -1], stringsAsFactors = FALSE)
    colnames(bb)[1] = 'sgRNA'
    
    jj = c(grep('F1_', colnames(bb)), grep('F2_', colnames(bb)), grep('F4_', colnames(bb)))
    bg0 = bb[, c(1, jj)]
    colnames(bg0)[c(2:ncol(bg0))] = c('Ab_bg0_F1', "In_bg0_F1", "Ab_bg0_F2",  "In_bg0_F2",  "Ab_bg0_F4",  "In_bg0_F4")
    
    jj = c(grep('F28_', colnames(bb)), grep('F29_', colnames(bb)), grep('F31_', colnames(bb)))
    bg_s3 = bb[, c(1, jj)]
    colnames(bg_s3)[c(2:ncol(bg_s3))] = c('Ab_bg_s3_F28', "In_bg_s3_F28", "Ab_bg_s3_F29",  "In_bg_s3_F29",  "Ab_bg_s3_F31",  "In_bg_s3_F31")
    
    ##### remove samples of low quality
    #jj = c(grep('F11', colnames(bb)), grep('F13', colnames(bb)))
    #if(length(jj)>0) bb = bb[, -jj]
  }
}


find.gene.target = function(x)
{
  test = unlist(strsplit(as.character(x), '_'))
  return(as.character(test[1]))
}


Check.RNAseq.Quality = function(read.count, design.matrix)
{
  #read.count = raw[, sels]; design.matrix = design.matrix[sels,]
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr"); 
  library("ggplot2")
  #library("dplyr")
  #library("ggplot2")
  #load(file=paste0('Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  # kk = grep('Ab', colnames(bb))
  if(ncol(design.matrix)>2){cc = apply(design.matrix[, -1], 1, paste0, collapse="_")
  }else{cc = design.matrix[, -1]}
  #o1 = order(cc)
  #read.count = read.count[o1,]
  #cc = cc[o1]
  read.count = as.matrix(read.count)
  #xx = raw
  #dim(raw)
  read.count[which(is.na(read.count))] = 0
  #xx = raw;
  
  par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,16,2,0.8)+0.1, tcl = -0.3)
  par(mfrow=c(1,1))
  
  total = apply(read.count, 2, sum)
  cc.uniq = unique(cc);
  cols = match(cc, cc.uniq)
  #cols = (c(1:length(unique(cc)))-1)%/%3+1
  barplot(total/10^6, horiz = TRUE, names.arg = colnames(read.count), las=1, col = cols, 
          main='Total nb of reads quantified for features', xlab='number of reads (Million)')
  abline(v=c(1, 2, seq(5, 20, by=5)), col='red', lty=1, lwd=2.0);#abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
  
  #xx = raw
  #xx[which(xx==0)] = 0.1
  
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
  boxplot(log10(read.count+0.1), las=3, col=cols, ylab='log10(nb of reads)', main='Read distribution for features')
  abline(v=-1, col='gray', lwd=1.5)
  
  library(vioplot)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
  read.norm = read.count+0.1;
  for(n in 1:ncol(read.norm)) read.norm[, n] = read.norm[, n]/sum(read.norm[,n])*100;
  plot(c(0, (ncol(read.count)+1)), c(0.01, (ncol(read.count))), ylim=range(0, 10), log='', xlab=NA, ylab="fractions of total counts", type='n')
  for(n in 1:ncol(read.count)){
   vioplot(read.norm[,n], names=rownames(read.count)[n], 
                    col=cols[n], at=n, add = TRUE) 
    mtext(colnames(read.count)[n], side = 1, at=n, las=2)
  }
  
  #x1 <- mtcars$mpg[mtcars$cyl==4]
  #x2 <- mtcars$mpg[mtcars$cyl==6]
  #x3 <- mtcars$mpg[mtcars$cyl==8]
  
  #title("Violin Plots of Miles Per Gallon")
  #p <- ggplot(mtcars, aes(factor(cyl), mpg))
  #p + geom_violin()
  
  ### make DESeq object using read counts and design matrix
  countData = ceiling(read.count)
  conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
  eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
  #dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
  dds <- dds[ rowSums(counts(dds)) > 0, ]
  
  dds <- estimateSizeFactors(dds, type="poscounts")
  #dds1 <- estimateSizeFactors(dds, type="poscounts")
  #dds2 <- estimateSizeFactors(dds, type="iterate")
  
  plot(sizeFactors(dds), total/10^6, xlab='Size factor', ylab="Library size (millions)", log='xy')
  #plot(sizeFactors(dds), sizeFactors(dds1), xlab='Size factor', ylab="Library size (millions)", log='xy')
  
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  ### boxplot (distributions of fpm) for all samples
  fpm = as.matrix(fpm)
  par(mfrow=c(1,1))
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
  for(n in 1:ncol(xx)){
    #kk = which(xx[,n]>0);
    if(n==1) boxplot(log2(fpm[, n]+2^-6), horizontal = FALSE, las=2, at=(n), ylim=c(-6, 23), xlim=c(1, (ncol(xx)+1)), names=as.character(colnames(xx)[n]),
                     las=1, width=0.6, ylab='log2(fpm)', col=cols[n], main="Distribution of normalized signals (cpm)")
    else boxplot(log2(fpm[, n]+2^-6), horizontal = FALSE, las=1, add=TRUE, at=(n), names=colnames(xx)[n], width=0.6, col=cols[n])
    mtext(colnames(fpm)[n], side = 1, at=n, las=2)
  }
  
  par(mfrow=c(1,1))
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
  plot(c(1:(ncol(fpm))), c(1:(ncol(fpm))), ylim=range(log2(fpm+2^-6)), log='', xlab=NA, ylab="fractions of total counts", type='n')
  for(n in 1:ncol(fpm)){
    vioplot(log2(fpm[,n]+2^-6), names=rownames(fpm)[n], 
            col=cols[n], at=n, add = TRUE) 
    mtext(colnames(read.count)[n], side = 1, at=n, las=2)
  }
  abline(v=-6, col='gray', lwd=1.5)
  
  ### pairwise correlations for fpm
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  library(corrplot)
  col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                             "cyan", "#007FFF", "blue","#00007F"))
  xx = as.matrix(fpm)
  xx[which(xx==0)] = NA
  M <- cor(xx, use = "na.or.complete")
  #corrplot(M, method="circle", type = 'upper', order="hclust")
  corrplot(M, method="ellipse", order="hclust", tl.cex=1.2, cl.cex=0.7, tl.col="black", addrect=ceiling(ncol(xx)/2), col=col1(100), rect.col=c('green'), rect.lwd=2.0)
  
  ### count transformation using vsd
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    #as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  ## illustration of count transformation using two samples
  colnames(df)[1:2] <- c("x", "y")  
  vsd.transform=ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  #print(vsd.transform)
  
  ### heatmap clustering samples using sample distance  
  sampleDists <- dist(t(assay(vsd)))
  #sampleDists
  #rld <- rlog(dds, blind=FALSE)
  sampleDistMatrix <- as.matrix( sampleDists )
  #rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  
  ## heatmap clustering samples using poisson distance 
  library("PoiClaClu")
  poisd <- PoissonDistance(t(counts(dds)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- colnames(dds)
  colnames(samplePoisDistMatrix) <- NULL
  pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors)
  
  ### clustering samples using PCA analysis
  pca=plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = FALSE)
  print(pca)
  
  Show.sample.names.PCA.Clusters = FALSE 
  if(Show.sample.names.PCA.Clusters){
    pca2save = as.data.frame(plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = TRUE))
    #p = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape=time)) + geom_point(size=3)
    #p + geom_text(hjust = 0.5, nudge_y = 0.1, size=2.5) 
    ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape=time)) + geom_point(size=3) +
      geom_text(hjust = 0.7, nudge_y = 1, size=2.5)  
    plot(ggp);
  }
  
  ###  pairwise correlation and fitting (to chek if there is batch effect)
  if(ncol(fpm)<20)
  {
    yy = as.matrix(fpm)
    yy[which(yy==0)] = NA;
    yy = log2(yy)
    pairs(yy, lower.panel=NULL, upper.panel=panel.fitting)
  }
  
  #ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  #  geom_point(size=3) +
  #  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #  coord_fixed()
  #write.table(processed, file=paste0('Table_quanseq_WT_tbx_negative_positive', version.analysis, '.txt'), sep='\t', quote=FALSE, col.names = TRUE, row.names = FALSE)
  #kk = grep('fpm', colnames(process
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.fitting = function (x, y, bg = NA, pch = par("pch"), cex = 0.01, col='black') 
{
  #x = yy[,1];y=yy[,2];
  #kk = which(x>0 & y>0); x=x[kk];y=y[kk]
  lims = range(c(x, y), na.rm = TRUE)
  points(x, y, pch = 1, col = col, cex = cex, xlim=lims, ylim=lims)
  abline(0, 1, lwd=1.5, col='red')
  R = cor(x, y, use="na.or.complete", method='pearson')
  text(lims[2]*0.2, lims[2]*0.9, paste0('R = ', signif(R, d=2)), cex=1., col='red')
  jj = which(!is.na(x) & !is.na(y))
  fit = lm(y[jj] ~ x[jj])
  #slope=summary(fit)$coefficients[1]
  slope = fit$coefficients[2]
  intercept = fit$coefficients[1]
  pval=summary(fit)$coefficients[4]
  abline(intercept, slope, lwd=1.2, col='darkblue', lty=3)
  #text(lims[2]*0.1, lims[2]*0.7, paste0('slop = ', signif(slope, d=2)), cex=1., col='blue')
  #text(lims[2]*0.1, lims[2]*0.6, paste0('pval = ', signif(pval, d=2)), cex=1., col='blue')
  #ok <- is.finite(x) & is.finite(y)
  #if (any(ok)) 
  #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
  #        col = col.smooth, ...)
}

### 
### This function is to check the normalized signals for control genes across all samples
Check.positve.controls = function(read.count, design.matrix, control.list=NULL, norm = "median", control.guides=NULL)
{
  if(is.null(control.list)){
    cat("need list of control genes---\n")
  }else{
    # read.count=raw; design.matrix=design.matrix; control.list=control.list;
    require(lattice);
    require(ggplot2)
    require('DESeq2');
    library("vsn");
    library("pheatmap");
    library("RColorBrewer");
    library("dplyr"); 
    library("ggplot2")
    #load(file=paste0('Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
    # kk = grep('Ab', colnames(bb))
    if(ncol(design.matrix)>2){cc = apply(design.matrix[, -1], 1, paste0, collapse="_")
    }else{cc = design.matrix[, -1]}
    #o1 = order(cc)
    #read.count = read.count[o1,]
    #cc = cc[o1]
    raw = as.matrix(read.count)
    #xx = raw
    dim(raw)
    raw[which(is.na(raw))] = 0
    xx = raw;
    
    ## total number of reads for each samples
    #par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,16,2,0.8)+0.1, tcl = -0.3)
    #par(mfrow=c(1,1))
    #total = apply(raw, 2, sum)
    #cc.uniq = unique(cc);
    #cols = match(cc, cc.uniq)
    #cols = (c(1:length(unique(cc)))-1)%/%3+1
    #barplot(total/10^6, horiz = TRUE, names.arg = colnames(raw), las=1, col = cols, main='Total nb of reads quantified for features', xlab='number of reads (Million)')
    #abline(v=c(1, 2, seq(5, 20, by=5)), col='red', lty=1, lwd=2.0);#abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
    #xx = raw
    #xx[which(xx==0)] = NA
    
    ## read distribution for each sample
    #par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
    #boxplot(log10(xx), las=3, col=cols, ylab='log10(nb of reads)', main='Read distribution for features')
    
    ### make DESeq object using read counts and design matrix
    countData = ceiling(raw)
    conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
    eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
    #dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
    dds <- dds[ rowSums(counts(dds)) > 0, ]
    dds <- estimateSizeFactors(dds, type="poscounts")
    #dds = estimateSizeFactors(dds, type="ratio")
    if(norm=="median") cpm = fpm(dds, robust = TRUE)
    if(norm=="total") cpm = fpm(dds, robust = FALSE)
    ### boxplot of normalization read counts (cpm)
    xx = as.matrix(cpm)
    par(mfrow=c(1,1))
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
    #for(n in 1:ncol(xx)){
    #  kk = which(xx[,n]>0);
    #  if(n==1) boxplot(log2(xx[kk, n]), horizontal = FALSE, las=2, at=(n), ylim=c(-6, 23), xlim=c(1, (ncol(xx)+1)), names=as.character(colnames(xx)[n]),
    #                   las=1, width=0.6, ylab='log2(fpm)', col=cols[n], main="Distribution of normalized signals (cpm)")
    #  else boxplot(log2(xx[kk, n]), horizontal = FALSE, las=1, add=TRUE, at=(n), names=colnames(xx)[n], width=0.6, col=cols[n])
    #  mtext(colnames(xx)[n], side = 1, at=n, las=2)
    #}
    
    ggs = sapply(rownames(cpm), function(x) unlist(strsplit(as.character(x), '_'))[1])
    
    sample.names = colnames(cpm)
    ranks = apply(cpm, 2, rank)
    #jj = order(sample.names)
    #sample.names = sample.names[jj]
    for(n in 1:length(control.list))
    {
      # n = 3
      cat(control.list[n], "\n")
      kk = which(ggs==control.list[n])
      if(length(kk)>0)
      {
        par(mfrow=c(1,2))
        ## counts
        matplot(t(log2(cpm[kk, ]+2^(-6))), type='b', lwd=2.0, cex=1.5, main=control.list[n], ylab="log2(normalized signals)", cex.main=1.5)
        mtext(colnames(cpm), side = 1, at=c(1:ncol(cpm)), las=2, cex=1.2)
        abline(h=-6, col='gray', lwd=1.5)
        ## ranks
        matplot(t((ranks[kk, ])), type='b', lwd=2.0, cex=1.5, main=control.list[n], ylab="ranks", cex.main=1.5)
        mtext(colnames(ranks), side = 1, at=c(1:ncol(ranks)), las=2, cex=1.2)
        abline(h=-6, col='gray', lwd=1.5)
      }
    }
    
    if(!is.null(control.guides))
    {
      kk = grep(control.guides, ggs)
      if(length(kk)>0)
      {
        par(mfrow=c(1,2))
        ## counts
        matplot(t(log2(cpm[kk, ]+2^(-6))), type='l', lwd=0.7, cex=1.5, main="control guides", ylab="log2(normalized signals)", cex.main=1.5, col='gray', lty=1)
        mtext(colnames(cpm), side = 1, at=c(1:ncol(cpm)), las=2, cex=1.2)
        points(c(1:ncol(cpm)), apply((log2(cpm[kk, ]+2^(-6))), 2, median), col='blue', lwd=2.0, type='b')
        abline(h=-6, col='gray', lwd=1.5)
        ## ranks
        matplot(t((ranks[kk, ])), type='l', lwd=0.7, cex=1.0, main="control guides", ylab="ranks", cex.main=1.5, col='gray', lty = 1)
        mtext(colnames(ranks), side = 1, at=c(1:ncol(ranks)), las=2, cex=1.2)
        abline(h=-6, col='gray', lwd=1.5)
      }
    }
  }
}

Check.positve.controls.sgRNA.UMI = function(read.count, design.matrix, control.list=NULL)
{
  if(is.null(control.list)){
    cat("need list of control genes---\n")
  }else{
    # read.count=raw; design.matrix=design.matrix; control.list=control.list;
    require(lattice);
    require(ggplot2)
    require('DESeq2');
    library("vsn");
    library("pheatmap");
    library("RColorBrewer");
    library("dplyr"); 
    library("ggplot2")
    #load(file=paste0('Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
    # kk = grep('Ab', colnames(bb))
    if(ncol(design.matrix)>2){cc = apply(design.matrix[, -1], 1, paste0, collapse="_")
    }else{cc = design.matrix[, -1]}
    #o1 = order(cc)
    #read.count = read.count[o1,]
    #cc = cc[o1]
    raw = as.matrix(read.count)
    #xx = raw
    dim(raw)
    raw[which(is.na(raw))] = 0
    #xx = raw;
    
    ## total number of reads for each samples
    par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,16,2,0.8)+0.1, tcl = -0.3)
    par(mfrow=c(1,1))
    total = apply(raw, 2, sum)
    cc.uniq = unique(cc);
    cols = match(cc, cc.uniq)
    #cols = (c(1:length(unique(cc)))-1)%/%3+1
    barplot(total/10^6, horiz = TRUE, names.arg = colnames(raw), las=1, col = cols, main='Total nb of reads quantified for features', xlab='number of reads (Million)')
    abline(v=c(1, 2, seq(5, 20, by=5)), col='red', lty=1, lwd=2.0);#abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
    #xx = raw
    #xx[which(xx==0)] = NA
    
    ## read distribution for each sample
    #par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
    #boxplot(log10(xx), las=3, col=cols, ylab='log10(nb of reads)', main='Read distribution for features')
    
    ### make DESeq object using read counts and design matrix
    countData = ceiling(raw)
    conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
    eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
    #dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
    dds <- dds[ rowSums(counts(dds)) > 0, ]
    dds <- estimateSizeFactors(dds)
    cpm = fpm(dds, robust = TRUE)
    
    ### boxplot of normalization read counts (cpm)
    cpm = as.matrix(cpm)
    par(mfrow=c(1,1))
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(16,3,2,0.2), tcl = -0.3)
    for(n in 1:ncol(cpm)){
      kk = which(cpm[,n]>0);
      if(n==1) boxplot(log2(cpm[kk, n]), horizontal = FALSE, las=2, at=(n), ylim=c(-6, 23), xlim=c(0, (ncol(cpm)+1)), names=as.character(colnames(cpm)[n]),
                       las=1, width=0.6, ylab='log2(fpm)', col=cols[n], main="Distribution of normalized signals (cpm)")
      else boxplot(log2(cpm[kk, n]), horizontal = FALSE, las=1, add=TRUE, at=(n), names=colnames(cpm)[n], width=0.6, col=cols[n])
      mtext(colnames(cpm)[n], side = 1, at=n, las=2)
    }
    
    ggs = sapply(rownames(cpm), function(x) unlist(strsplit(as.character(x), '_'))[1])
    
    sample.names = colnames(cpm)
    jj = order(sample.names)
    sample.names = sample.names[jj]
    for(n in 1:length(control.list))
    {
      # n = 1
      cat(control.list[n], "\n")
      kk = which(ggs==control.list[n])
      matplot(t(log2(cpm[kk, jj]+2^(-6))), type='b', lwd=0.4, cex=0.2, main=control.list[n], ylab="log2(cpm)", cex.main=1.5)
      mtext(sample.names, side = 1, at=c(1:length(sample.names)), las=2, cex=1.5)
      abline(h=-6, col='gray', lwd=1.5)
      
    }
    
  }
}
  
  

