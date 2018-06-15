#####################################
#### data version and analysis version
#####################################
version.data = '_R5213_v1';
version.analysis = paste0(version.data, '_analysis_v2')
nb.samples = 16;
resDir = "../Results/2017-10-06";
if(!dir.exists(resDir)) dir.create(resDir);
######################################
######################################
### import the raw data and clean the table
######################################
######################################
Import.DATA.table.exp.Design = FALSE
if(Import.DATA.table.exp.Design)
{
  library('xlsx')
  aa = read.table(paste('../DATA/out_5_Hagar_1MM_10X_2min_IncAbu_v3.txt', sep=''), header = TRUE, sep='\t')
  dim(aa)
  
  jj = grep('_In', colnames(aa))
  kk = grep('_Ab', colnames(aa))
  
  aa = data.frame(aa[, 1], aa[,kk], aa[,jj], stringsAsFactors = FALSE)
  colnames(aa)[1] = 'sgRNA'
  
  Change.Samples.Names = FALSE
  if(Change.Samples.Names)
  {
    sampleID.name = read.xlsx('../R4950/Infos_Samples_v2.xlsx', sheetIndex = 1)
    Changecolnames = function(x)
    {
      #x = colnames(aa)[2] 
      test = unlist(strsplit(as.character(x), '_'))
      ttest = paste0('F', test[2])
      nn = sampleID.name$Name[which(sampleID.name$Primer==ttest)]
      
      #return(paste(c(paste0('F', test[2]), test[3]), collapse = '_'))
      return(paste(c(test[3], as.character(nn), ttest), sep=' ',  collapse = '_'))
    }
    names = sapply(colnames(aa)[-1], Changecolnames)
    colnames(aa)[-1] = names
    kk = grep('Ab_', colnames(bb))
    jj = grep('In_', colnames(bb))
    #kk = kk[order(colnames(aa)[kk])]
    #samples = paste0('F', c(1, 2, 4, 5, 7, 8, 28, 29, 31, 10, 11, 13, 19, 20, 22, 14, 16, 17, 32, 34, 35, 23, 25, 26, 37, 38, 40))
    #kk = match(paste(samples, 'Ab', sep='_'), colnames(aa))
    #jj = match(paste(samples, 'In', sep='_'), colnames(aa))
    bb = data.frame(bb[,1], bb[,kk], bb[,jj], stringsAsFactors = FALSE)
    colnames(bb)[1] = 'sgRNA'
  }
  
  bb = aa;
  source("functions_Screen.R")
  gene = sapply(bb$sgRNA, find.gene.target)
  
  bb = data.frame(bb$sgRNA, gene, bb[, -1], stringsAsFactors = FALSE)
  colnames(bb)[1] = 'sgRNA'
  
  save(bb, file=paste0('../Rdata/Screen_countData', version.data, '.Rdata'))

  ################
  ##### clean sgRNA names
  Clean.sgRNA.Names = TRUE
  if(Clean.sgRNA.Names)
  {
    load(file=paste0('../Rdata/Screen_countData', version.data, '.Rdata'))
    sgrna.targets = NULL
    indexs = NULL
    
    for(n in 1:nrow(bb))
    #for(n in 1:100)
    {
      #n = 37;
      test = unlist(strsplit(as.character(bb$sgRNA[n]), ','))
      if(length(test)==1){
        targ = unlist(strsplit(as.character(test), '_'))[1]
        sgrna.targets = rbind(sgrna.targets, c(test, targ))
        indexs = c(indexs, n)
      }
      if(length(test)>1){
        for(ttest in test)
        {
          targ = unlist(strsplit(as.character(ttest), '_'))[1]
          sgrna.targets = rbind(sgrna.targets, c(ttest, targ))
          indexs = c(indexs, n)
        }
      }
      if(length(test)<1){
        print(n)
      }
    }
    
    xx = data.frame(sgrna.targets, bb[indexs, -c(1:2)], stringsAsFactors = FALSE)
    colnames(xx)[1:2] = c('sgRNA', 'Gene')
    
    bb = xx; 
    save(bb, file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
    
  }
  
  #write.table(bb, file='All_table_counts_barcodes.txxt', sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)
  #bb = read.table(file='All_table_counts_barcodes.txxt', sep='\t',header = TRUE, as.is = c(1, 2))
  #bb = data.frame(bb);
  #bb = as.matrix(bb) 
}

################
#### Control quality of screening data for sequencing and positive controls
################
Control.Quality.Screen.Data = FALSE
if(Control.Quality.Screen.Data)
{
  load(file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  
  cc = c("SFMIII_P5_P5_ABC", "SFMIII_NotP5_P5_ABC",  "SFMIII_P6_P7_ABC", "SFMIV_P6_P7_2_ABC", "Background", "Background", "Background",
         "SFMRIV", "SFMRIV", "SFMRIV",  "SFRIV_BFP", "SFRIV_BFP", "SFRIV_BFP", "SFRIV_GFP",  "SFRIV_GFP",  "SFRIV_GFP")
  
  Compare.Ramesh.background = TRUE ### add Ramesh's background to compare
  if(Compare.Ramesh.background)
  {
    bb0 = bb
    load(file=paste0('../../../RameshKarin/CrisprCas9_Screen_RK/Rdata/Screen_countData_sgRNA_Gene_clean_mapping_R4950_v2.Rdata'))
    rownames(bb) = bb$sgRNA
    ## keep both Ab and In
    xx = bb[match(bb0$sgRNA, bb$sgRNA), intersect(grep("Ab_|In_", colnames(bb)), grep("_bg", colnames(bb)))]
    #colnames(xx) = sapply(colnames(xx), function(x) paste0(unlist(strsplit(as.character(x), "_"))[-1], collapse = "_"))
    bb = data.frame(bb0, xx, stringsAsFactors = FALSE)
    cc = c(cc, rep(c("bg_s0", "bg_s2", "bg_s3"), each=3));
  }
  
  ### Check the read counts for samples
  source("functions_Screen.R")
  pdfname = paste0(resDir, "/SCREEN_Data_Qulity_Assessment_background_Comparison", version.analysis, ".pdf")
  pdf(pdfname, width = 16, height = 8)
  
  kk = grep('In', colnames(bb))
  raw = as.matrix(bb[, kk])
  rownames(raw) = bb$sgRNA;
  design.matrix = data.frame(samples = colnames(raw), conditions = as.factor(cc), stringsAsFactors = FALSE)
  
  #sels = grep("Background|bg", colnames(raw))
  sels = c(1:ncol(raw))
  
  Check.RNAseq.Quality(read.count = raw[, sels], design.matrix = design.matrix[sels,]);
  
  dev.off()
  
  
  #### chekcing individual examples across samples
  Checking.individual.examples.across.samples = FALSE
  if(Checking.individual.examples.across.samples)
  {
    control.list = c("Cbx7", "Rnf2", "Suz12","Rybp", "Eed", "Ezh2",  "Pcgf1", "Pcgf2", "Pcgf6", "Mtf2","Pcgf4", "Jarid2","Kdm2b",
                     "Dnmt1",  "Uhrf1", "Kdm3a", "Trp53")
                     #"Nf2", "Rarg", "Rxrb", "Morc2a", "Fam208a", "Trp53","Pds5a", "Zic3", "Tet2", "Pou5f1", "Mphosph8", "Mrpl28", 
                     #"Zfp11", "Foxp3", "Ep300", "Cacna1d", "Cdkn2a", "Cramp1l", "Jmjd1c", "Fam175a", "Prkce", "Nr3c2", "Foxa2","Ccr5"
                     #)
    control.list = unique(control.list)
    source("functions_Screen.R")
    
    Which2Compare = 'In';
    normalization = "total"
    kk = grep(Which2Compare, colnames(bb))
    
    pdfname = paste0(resDir, "/Controles_genes_", Which2Compare, "_", normalization,  version.analysis, ".pdf")
    pdf(pdfname, width = 20, height = 8)
    raw = as.matrix(bb[, kk])
    rownames(raw) = bb$sgRNA;
    design.matrix = data.frame(samples = colnames(raw), conditions = as.factor(cc), stringsAsFactors = FALSE)
    
    ## change the order of samples
    ii = c(grep("bg0|bg_s2|bg_s3", colnames(raw)), grep("Background", colnames(raw)))
    ii = c(ii, setdiff(c(1:ncol(raw)), ii))
    
    Check.positve.controls(read.count=raw[, ii], design.matrix=design.matrix[ii, ], control.list=control.list, norm = normalization)
    
    dev.off()
    
    
    source("functions_Screen.R")
    pdfname = paste0(resDir, "/SCREEN_Compare_backgrounds", version.analysis, ".pdf")
    pdf(pdfname, width = 20, height = 16)
    
    yy = as.matrix(fpm[, grep("Background|Ab_bg", colnames(fpm))])
    yy[which(yy==0)] = 2^-6;
    yy = log2(yy)
    pairs(yy, lower.panel=NULL, upper.panel=panel.fitting)
    
    dev.off()
    
  }
  
  #### check the correlation between read counts of UMI
  umi = as.matrix(bb[, grep("_In", colnames(bb))])
  pdfname = paste0(resDir, "/SCREEN_read_counts_vs_UMI", version.analysis, ".pdf")
  pdf(pdfname, width = 16, height = 14)
  par(mfrow=c(2,2))
  
  for(n in 1:ncol(raw))
  {
    #n = 1
    plot(log10(raw[,n]+0.1), log10(umi[, n]+0.1), xlab='log10(read counts)', ylab='log10(UMI)', main=colnames(raw)[n], cex=0.5);
    abline(0, 1, lwd=1.5, col='red')
  }
  
  dev.off()
  
  ### Check read counts per barcoded sgRNAs
  Check.readcount.for.barcode = FALSE
  if(Check.readcount.for.barcode)
  {
    datafile = "../DATA/out_6_Hagar_1MM_10X_2min_CrUMI_v3.txt"
    tab5rows <- read.table(datafile, header = TRUE, nrows = 5)
    classes <- sapply(tab5rows, class)
    xx <- read.table(datafile, header = TRUE, colClasses = classes)
    yy = as.matrix(xx[, -c(1:2)])
    kk = which(rowSums(yy)>0)
    
    raw = as.matrix(yy)
    rownames(raw) = paste0(xx$guide,"_", xx$barcode) 
    design.matrix = data.frame(samples = colnames(raw), conditions = as.factor(cc), stringsAsFactors = FALSE)
    
    control.list = c("Cbx7", "Rnf2", "Eed", "Suz12", "Ezh2", "Mtf2")
    
    source("functions_Screen.R")
    pdfname = paste0(resDir, "/SCREEN_Controles_list_readCounts_UMI", version.analysis, ".pdf")
    pdf(pdfname, width = 20, height = 10)
    
    Check.positve.controls.sgRNA.UMI(read.count=raw, design.matrix=design.matrix, control.list=control.list)
    
    dev.off() 
  }
 
}

############################
############################
#### Select the read count matrix and make bash script for mageck
############################
############################
Prepare.Run.mageck = FALSE
if(Prepare.Run.mageck)
{
  load(file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  DIR.mageck = paste0(resDir, "/mageck/");
  system(paste0('mkdir -p ', DIR.mageck))
  
  Using.Replicates = FALSE
  Which2Compare = 'In';
  Compare.Ramesh.background = FALSE ### add Ramesh's background to compare
  
  DIR.res = "Comparison_Incidence_eachSample_normalization_median_testRRA/" ## this is path relative to the script "run_mageck.sh"
  system(paste0('mkdir -p ', DIR.mageck, DIR.res))
  
  ## prepare table for mageck
  kk = grep(Which2Compare, colnames(bb))
  dds =  data.frame(bb[, c(1, 2, kk)], stringsAsFactors = FALSE)
  
  cc = c("SFMIII_P5_P5_ABC", "SFMIII_NotP5_P5_ABC",  "SFMIII_P6_P7_ABC", "SFMIV_P6_P7_2_ABC", "Background", "Background", "Background",
         "SFMRIV", "SFMRIV", "SFMRIV",  "SFRIV_BFP", "SFRIV_BFP", "SFRIV_BFP", "SFRIV_GFP",  "SFRIV_GFP",  "SFRIV_GFP")
  print(colnames(dds))
  
  if(Compare.Ramesh.background)
  {
    bb0 = bb
    load(file=paste0('../../../RameshKarin/CrisprCas9_Screen_RK/Rdata/Screen_countData_sgRNA_Gene_clean_mapping_R4950_v2.Rdata'))
    bb.ramesh = bb;
    bb = bb0;
    rownames(bb.ramesh) = bb.ramesh$sgRNA
    xx = bb.ramesh[match(bb$sgRNA, bb.ramesh$sgRNA), intersect(grep(Which2Compare, colnames(bb.ramesh)), grep("_bg", colnames(bb.ramesh)))];
    #colnames(xx) = sapply(colnames(xx), function(x) paste0(unlist(strsplit(as.character(x), "_"))[-1], collapse = "_"))
    
    dds = data.frame(dds, xx, stringsAsFactors = FALSE)
    
    cc = c(cc, rep(c("bg_s0", "bg_s2", "bg_s3"), each=3));
    
    print(colnames(dds))
  }
  
  ### make bash file for selected comparisons
  #DIR.CWD = getwd()
  #DIR.OUT = paste0(getwd(),'/mageck/'
 
  
  ### define comparisons for mageck
  find.conditions = function(x){
    return(unlist(strsplit(as.character(x), '_'))[1]);
  }
  #cc = unique(sapply(colnames(dds)[-c(1:2)], find.conditions)) 
  
  remove.samples = FALSE
  if(remove.samples){
    if(length(jj)>0) dds = dds[,-jj]
  }
  
  if(Using.Replicates)
  {
    conds = c('sgRNA', 'Gene', 
              "SFMIII_P5_P5_ABC", "SFMIII_NotP5_P5_ABC",  "SFMIII_P6_P7_ABC", "SFMIV_P6_P7_2_ABC", 
              paste0("Background_r", c(1:3)),
              paste0("SFMRIV_r", c(1:3)),  
              paste0("SFRIV_BFP_r", c(1:3)),  
              paste0("SFRIV_GFP_r", c(1:3))
    )
    
    #conds = colnames(dds)
    conds <- factor(conds)
    colnames(dds) = conds
    
    ## define pairwise comparisons
    comparison = list()
    n = 1
    for(ss in unique(cc)){
      if(ss != "Background"){
        comparison[[n]] = c(ss, "Background")
        n = n +1;
      }
    }
    
    ## examples
    #comparison = list( c('bg.s2', 'bg0'), c('bg.s3', 'bg0'),  c('bg.s3', 'bg.s2'),
    #                   c("maint.s2.R1", 'bg0'), c("maint.s2.R1", 'bg.s2'), 
    #                   c("maint.s2.R2", 'bg0'), c("maint.s2.R2", 'bg.s2'), 
    #                   c("maint.s3", 'bg0'), c("maint.s3", 'bg.s3') 
    #c("init.s2.R2", 'bg.s1'), c("init.s3.R2", 'bg.s1'), 
    #c("maint.s2", 'bg.f'), c("maint.s3", 'bg.f')
    #c("init.s2.R1", 'bg.f'), c("init.s3.R1", 'bg.f'), c("init.s2.R2", 'bg.f'), c("init.s3.R2", 'bg.f'),
    #c("maint.s2", 'bg0'), c("maint.s3", 'bg0'), 
    #c("init.s2.R1", 'bg0'), c("init.s3.R1", 'bg0'), c("init.s2.R2", 'bg0'), c("init.s3.R2", 'bg0'))
    #)
  }else{
    conds = c('sgRNA', 'Gene', 
              "SFMIII_P5_P5_ABC", "SFMIII_NotP5_P5_ABC",  "SFMIII_P6_P7_ABC", "SFMIV_P6_P7_2_ABC", 
              paste0("Background_r", c(1:3)),
              paste0("SFMRIV_", c("A", "B", "C")),  
              paste0("SFRIV_BFP_", c("A", "B", "C")),  
              paste0("SFRIV_GFP_", c("A", "B", "C"))
    )
    
    if(Compare.Ramesh.background) conds = c(conds, c(paste0("bg_s0_r", c(1:3)), paste0("bg_s2_r", c(1:3)), paste0("bg_s3_r", c(1:3))))
    
    colnames(dds) = factor(conds)
    
    ## define pairwise comparisons
    comparison = list()
    n = 1
    for(ss in unique(conds)){
      cat(ss, "\n")
      if(ss != "sgRNA" & ss !="Gene" & length(grep("Background|bg_s", ss))==0){
        if(Compare.Ramesh.background){
          comparison[[n]] = c(ss, "bg_s0")
        }else{
          comparison[[n]] = c(ss, "Background")
        }
       
        n = n +1;
      }
    }
  }
  
  ### save the read count table for the mageck
  sgrna_counnt = dds
  write.table(sgrna_counnt, file=paste0(DIR.mageck, 'sgrna_count_table.txt'), sep= '\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
  
  ### make bash script to run mageck
  for(n in 1:length(comparison))
  {
    #n = 1
    treatment = unlist(comparison[n])[1]; control = unlist(comparison[n])[2];
    samples.t = grep(treatment, colnames(sgrna_counnt));
    samples.c = grep(control, colnames(sgrna_counnt));
    if(length(samples.t)<1)  cat('ERROR for ', treatment, '\n');
    if(length(samples.c)<1)  cat('ERROR for ', control, '\n');
    if(length(samples.t)>=1 & length(samples.c)>=1)
    {
      treatment = gsub("\\.", "_", treatment);
      control = gsub("\\.", "_", control);
      logs = paste0(DIR.mageck, 'run_mageck.sh')
      
      output = paste0(DIR.res, Which2Compare, '_', treatment, '_VS_', control)
      
      samples.t = paste(colnames(sgrna_counnt)[samples.t], sep='', collapse = ',') 
      samples.c = paste(colnames(sgrna_counnt)[samples.c], sep='', collapse = ',') 
      if(n==1) {
        cmd = paste0('mageck test -k ', 'sgrna_count_table.txt', ' -t ', samples.t,  ' -c ',  samples.c, ' -n ',  output, 
                     ' --adjust-method pounds --sort-criteria pos --pdf-report --normcounts-to-file --norm-method median --remove-zero treatment')
        #cat(cmd, '\n')
        write(cmd, file=logs)
      }else{
        cmd = paste0('mageck test -k ', 'sgrna_count_table.txt', ' -t ', samples.t,  ' -c ',  samples.c, ' -n ',  output,  
                     ' --adjust-method pounds --sort-criteria pos --pdf-report --normcounts-to-file  --norm-method median --remove-zero treatment')
        #cat(cmd, '\n')
        write(cmd, file=logs, append = TRUE)
      }
    }  
  }
  
}

########################################################################
########################################################################
##### Test screen analysis using only RANK to avoid normalizatoin problems
########################################################################
########################################################################
Screen.Analysis.edgeR.DESeq2 = FALSE
if(Screen.Analysis.edgeR.DESeq2)
{
  load(file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  
  ## check the rank of guides
  Check.Rank.each.guide = FALSE
  if(Check.Rank.each.guide)
  {
    if(Compare.Ramesh.background)
    {
      bb0 = bb
      load(file=paste0('../../../RameshKarin/CrisprCas9_Screen_RK/Rdata/Screen_countData_sgRNA_Gene_clean_mapping_R4950_v2.Rdata'))
      bb.ramesh = bb;
      bb = bb0;
      rownames(bb.ramesh) = bb.ramesh$sgRNA
      xx = bb.ramesh[match(bb$sgRNA, bb.ramesh$sgRNA), grep("_bg", colnames(bb.ramesh))];
      #colnames(xx) = sapply(colnames(xx), function(x) paste0(unlist(strsplit(as.character(x), "_"))[-1], collapse = "_"))
      
      bb = data.frame(bb, xx, stringsAsFactors = FALSE)
      
      #dds = data.frame(dds, xx, stringsAsFactors = FALSE)
      #cc = c(cc, rep(c("bg_s0", "bg_s2", "bg_s3"), each=3));
      #print(colnames(dds))
    }
    
    sample.to.check = "SFMIII_P6_P7_ABC"
    kk = grep(sample.to.check, colnames(bb))
    
    yy = bb[, c(1, 2, kk, grep("Background|bg", colnames(bb)))]
    
    rownames(yy) = yy$sgRNA
    #Which2Compare = 'In';
    #o1 = order(-yy[, 4])
    #yy = yy[o1, ]
    
    index = which(yy$Gene=="Eed")
    print(index)
    yy[index, ]
    
    rr = cbind( (rank(yy$SFMIII_P6_P7_ABC_In)/nrow(yy))[index], 
    (rank(yy$Background_A_15_6_In)/nrow(yy))[index],
    (rank(yy$Background_B_15_6_In)/nrow(yy))[index],
    (rank(yy$Background_C_15_6_In)/nrow(yy))[index])
    
    print(rr)
    
    (rank(yy$SFMIII_P6_P7_ABC_Ab)/nrow(yy))[index]
    (rank(yy$Background_A_15_6_Ab)/nrow(yy))[index]
    (rank(yy$Background_B_15_6_Ab)/nrow(yy))[index]
    (rank(yy$Background_C_15_6_Ab)/nrow(yy))[index]
    
    require(gCrisprTools)
    data('fit')
    data('ann')
    output <- ct.generateResults(fit, ann, permutations = 100, RRAalphaCutoff = 1)
    head(output)
    p = seq(0, 1, length.out=20)
    fc = seq(-3, 3, length.out=20)
    fc[2] = NA
    fc[3] = -20
    stats = data.frame(
      Depletion.P=p,
      Enrichment.P=rev(p),
      fc=fc
    )
    ct.applyAlpha(stats,scoring="combined")
    
    
  }
  
  ### test gCrisprTools 
  Test.gCrisprTools.DESeq2 = FALSE
  if(Test.gCrisprTools.DESeq2)
  {
    require(gCrisprTools)
    
    data('fit')
    data('ann')
    genePvals <- ct.RRAaPvals(fit$p.value, ann, permute = 100)
    
    DIR.mageck = "../Results/2017-09-25/gCrisprTools/";
    DIR.res = "Comparison_usingRep_test/" ## this is path relative to the script "run_mageck.sh"
    #if(!dir.exists(resMageck)) dir.create(resMageck);
    system(paste0('mkdir -p ', DIR.mageck))
    system(paste0('mkdir -p ', DIR.mageck, DIR.res))
    ### make bash file for selected comparisons
    #DIR.CWD = getwd()
    #DIR.OUT = paste0(getwd(),'/mageck/')
    
    Which2Compare = 'In';
    
    kk = grep(Which2Compare, colnames(bb))
    dds =  data.frame(bb[, c(1, 2, kk)], stringsAsFactors = FALSE)
    
    cc = c("SFMIII_P5_P5_ABC", "SFMIII_NotP5_P5_ABC",  "SFMIII_P6_P7_ABC", "SFMIV_P6_P7_2_ABC", "Background", "Background", "Background",
           "SFMRIV", "SFMRIV", "SFMRIV",  "SFRIV_BFP", "SFRIV_BFP", "SFRIV_BFP", "SFRIV_GFP",  "SFRIV_GFP",  "SFRIV_GFP")
    
    cat(colnames(dds), "\n")
    
    rownames(dds) = dds$sgRNA;
    
    library(Biobase)
    library(limma)
    library(gCrisprTools)
    annot =  data.frame(sample.ID= colnames(dds)[-c(1:2)], condition=cc, stringsAsFactors = FALSE)
    rownames(annot) = annot[,1]
    annot = AnnotatedDataFrame(data=annot)
    es = ExpressionSet(as.matrix(dds[, -c(1,2)]), phenoData= annot)
    #library(edgeR)
    #data("es", package = "gCrisprTools")
    es
    ann = data.frame(dds[, c(1, 2, 2)], stringsAsFactors = FALSE)
    rownames(ann) = ann[, 1];
    colnames(ann) = c("ID", "geneID", "geneSymbol")
    
    sk <- relevel(as.factor(pData(es)$condition), "Background")
    names(sk) <- row.names(pData(es))
    sk
    
    es.floor <- ct.filterReads(es, read.floor = 1, sampleKey = sk)
    
    es <- ct.normalizeGuides(es, 'scale', annotation = ann, sampleKey = sk, plot.it = TRUE)
    
    ct.gRNARankByReplicate(es, sk)
    
    ct.gRNARankByReplicate(es, sk, annotation = ann, geneSymb = "Cbx7")
  }
  
  
 
  x5 = DGEList(counts=as.matrix(dds[, -c(1,2)]), group=cc)
  #y <- DGEList(counts=x,group=group)
  #x5$counts = 
  x5$samples = data.frame("SampleID"=colnames(x5$counts),
                          "group"=as.factor(cc), "lib.size"=colSums(x5$counts))
  y = x5
  y <- calcNormFactors(y, method = "TMM")
  y$samples
  y1 <- calcNormFactors(y, method = "RLE")
  y1$samples
  
  x5 = y;
  x5$genes = dds[,1:2]
  rownames(x5$genes) = dds[, 1]
  sel = rowSums(cpm(x5$counts)>5)>=2
  x5 = x5[sel,]
  
  plotMDS(x5, labels=x5$samples$group, xlim=c(-2,4),
          col=c(1,rep(c(2,3),each=4)), main="Shalem: MDS Plot")
  legend("topright", legend=c("Baseline", "DMSO", "PLX"), col=1:3, pch=15)
  
  cpm = cpm(x5)
  
  cpm[grep("Suz12", rownames(cpm)),]
  
}


########################################################################
#### Compare Hp screen an Cbx7 screen
########################################################################
xlist<-list.files(path=paste0('../../../RameshKarin/CrisprCas9_Screen_RK/mageck_Ab/R4950_v3_all'),
                  pattern = "*.gene_summary.txt", full.names = TRUE)
xlist = xlist[grep("Ab_maint_s2|Ab_maint_s3", xlist)]
bnames = basename(xlist)

#hits.rk = NULL;
for(n in 1:length(xlist)) 
{
  #n = 2
  ff = read.table(xlist[n], sep='\t', header = TRUE, as.is = 1);
  if(n==1){
    hits.rk = data.frame(cbind(as.character(ff$id), ff$pos.score), stringsAsFactors = FALSE);
    colnames(hits.rk)[1] = "gene";
  }else{
    mm = match(ff$id, hits.rk$gene)
    hits.rk = cbind(hits.rk[mm, ], ff$pos.score);
  }
  #summary[[n]] = data.frame(ff, stringsAsFactors = FALSE);
};

hits.rk = data.frame(hits.rk, stringsAsFactors = FALSE)
colnames(hits.rk)[-1] = bnames
rownames(hits.rk) = hits.rk$gene
hits.rk = hits.rk[, grep("bg_s2|bg_s3", colnames(hits.rk))]
xx = as.numeric((apply(as.matrix(hits.rk), 1, min)))
names(xx) = rownames(hits.rk)
hits.rk = xx;
hits.rk = hits.rk[order(hits.rk)]

ylist<-list.files(path=paste0('../Results/2017-09-26/mageck/Comparison_Incidences_eachSample_normalization_median'),
                  pattern = "*.gene_summary.txt", full.names = TRUE)
ylist = ylist[grep("SFM", ylist)]
bnames = basename(ylist)
bnames = sapply(bnames, function(x) unlist(strsplit(as.character(x), "[.]"))[1])

cutoff = 0.01
pdfname = paste0(resDir, "/Compare_HP_Cbx7_screen_candidates",  version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 10)
for(n in 1:length(ylist))
{
  #n = 1;
  out = read.table(ylist[n], sep='\t', header = TRUE)
  out = out[match(names(hits.rk), out$id), ]
  lims = -log10(min(c(out$pos.score, hits.rk)))
  jj = which(out$pos.score<cutoff | hits.rk<cutoff)
  plot(-log10(out$pos.score)[jj], -log10(hits.rk)[jj], cex=0.5, xlab="Cbx", ylab= "Hp1", 
       xlim = c(0, lims), ylim = c(0, lims), main = paste0('RRA score for ', bnames[n]), type='p')
  abline(0, 1, col='darkgray', lwd=2.0)
  abline(v=2, col='darkgray', lwd=1.0, lty=2);abline(h=2, col='darkgray', lwd=1.0, lty=2)
  text(-log10(out$pos.score)[jj], -log10(hits.rk)[jj], out$id[jj],  cex=0.6, pos = 4, offset = 0.3)
  #sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr<fdr.cutoff); points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col='darkblue', pch=16)
  #text(-log10(out1$pos.score)[sel[1:10]], -log10(out2$pos.score)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
  #sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff); points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col='blue', pch=0)
  #text(-log10(out1$pos.score)[sel[1:10]], -log10(out2$pos.score)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
  #sel = which(out1$pos.fdr> fdr.cutoff & out2$pos.fdr<fdr.cutoff); points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col='lightblue', pch=17)
  #text(-log10(out1$pos.score)[sel[1:10]], -log10(out2$pos.score)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
  #sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='lightblue', pch=0)
  #sel = which(out1$pos.fdr>= fdr.cutoff & out2$p
}
dev.off()

########################################################################
########################################################################
##### make a summary and display for selected pairwise comparisons
########################################################################
########################################################################
Which2Compare = 'Ab';
#paste0(Which2Compare, "_*.gene_summary.txt")
xlist<-list.files(path=paste0(getwd(), '/mageck_Ab/R4950_v3_all'), pattern = "*.gene_summary.txt", full.names = TRUE)
xlist = xlist[grep(Which2Compare, xlist)]
summary = list();
bnames = basename(xlist)
for(n in 1:length(xlist)) 
{
  ff = read.table(xlist[n], sep='\t', header = TRUE);
  summary[[n]] = data.frame(ff, stringsAsFactors = FALSE);
};
names(summary) = bnames

#ggs = unique(bb$Gene)
out.names = c('gene')
for(n in 1:length(comparison))
{
  treatment = unlist(comparison[n])[1]; control = unlist(comparison[n])[2];
  treatment = gsub("\\.", "_", treatment); control = gsub("\\.", "_", control);
  output = paste0(Which2Compare, '_', treatment, '_VS_', control);
  out.names = c(out.names, output)
}

Pairwise.Comparison = FALSE
if(Pairwise.Comparison)
{
  ##### volcanoplot for bg_s1, miain_s3 and init_s3 using bg0 as the reference
  Volcano.plot.single.comparion = FALSE
  if(Volcano.plot.single.comparion)
  {
    #TOplots = c("Ab_bg_s1_VS_bg0", "Ab_maint_s3_VS_bg0", "Ab_init_s3_R1_VS_bg0", "Ab_init_s2_R2_VS_bg0",
    #            "Ab_maint_s3_VS_bg_f", "Ab_init_s3_R1_VS_bg_f", "Ab_init_s2_R2_VS_bg_f")
    
    TOplots = c("Ab_bg_s1_VS_bg0", "Ab_maint_s3_VS_bg0")
    ggs = c("Ddx3", "Uba2", "Pou5f1", "Nanog", "Sox2", 
            "Setdb1", "Ehmt2", "Ehmt1", "Ago2", "Atf7ip", "Mbd1", "Pphln1",
            "Mphosph8", "Morc2a", "Phf8", "Fam208a", "Zscan10", "Dnmt1", "Uhrf1")           
    rra.cutoff = 0.01;cex=0.6;pch=16;
    col= 'gray'
    pdf("PLOTs/Vocanoplot_logFC_RRA_Negative_Positive_v6.pdf", width = 4, height = 3)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,1.5,3)+0.1, tcl = -0.3)
    for(n in 1:length(TOplots))
    {
      #n = 1;
      output = TOplots[n];
      out1 = basename(xlist[grep(output, xlist)])
      out1 = summary[[which(names(summary)==out1)]]
      lims = c(-7, 7);
      if(n>1) {lims = max(abs(out1$pos.lfc));lims = c(-lims, lims);}
      cex.text=0.3;cex.offset=0.2; cex=0.5;
      #ylims = c(0, 24)
      for(m in c(1))
      {
        plot(out1$pos.lfc, -log10(out1$pos.score), cex=0.5, type='n', main=output, xlab='log(fold change)', ylab='-log10(RRA)', xlim = lims)
        if(n ==1){
          kk = which(out1$pos.lfc>0);
          kk1 = kk[which(out1$pos.score[kk]<rra.cutoff)]; kk2 = kk[which(out1$pos.score[kk]>=rra.cutoff)]; 
          points(out1$pos.lfc[kk1], -log10(out1$pos.score[kk1]), cex=cex, col=col, pch=pch);
          points(out1$pos.lfc[kk2], -log10(out1$pos.score[kk2]), cex=0.2, col=col);
          jj = which(out1$pos.lfc<=0);
          jj1 = jj[which(out1$neg.score[jj]<rra.cutoff)]; jj2 = jj[which(out1$neg.score[jj]>=rra.cutoff)]; 
          points(out1$pos.lfc[jj1], -log10(out1$neg.score[jj1]), cex=cex, col=col, pch=pch)
          points(out1$pos.lfc[jj2], -log10(out1$neg.score[jj2]), cex=0.2, col=col)
          if(m==1){
            kk1 = kk1[which(!is.na(match(out1$id[kk1], ggs))==TRUE)];jj1 = jj1[which(!is.na(match(out1$id[jj1], ggs))==TRUE)];
            text(out1$pos.lfc[kk1], -log10(out1$pos.score[kk1]), out1$id[kk1], cex=cex.text, offset = cex.offset, pos = 2);
            text(out1$pos.lfc[jj1], -log10(out1$neg.score[jj1]), out1$id[jj1], cex=cex.text, offset = cex.offset, pos=2);
            points(out1$pos.lfc[kk1], -log10(out1$pos.score[kk1]), cex=cex, col='red', pch=pch);
            points(out1$pos.lfc[jj1], -log10(out1$neg.score[jj1]), cex=cex, col='red', pch=pch)
          }
        }else{
          kk1 = which(out1$pos.score[kk]<rra.cutoff); kk2 = which(out1$pos.score[kk]>=rra.cutoff); 
          points(out1$pos.lfc[kk1], -log10(out1$pos.score[kk1]), cex=cex, col=col, pch=pch);
          points(out1$pos.lfc[kk2], -log10(out1$pos.score[kk2]), cex=0.2, col=col);
          if(m==1) {
            jj = which(!is.na(match(out1$id, ggs))==TRUE);
            text(out1$pos.lfc[jj], -log10(out1$pos.score)[jj], out1$id[jj], cex=cex.text, offset = cex.offset, pos = 2);
            points(out1$pos.lfc[jj], -log10(out1$pos.score[jj]), cex=cex, col='red', pch=pch);
          }
        }
        abline(v=0, lwd=2.0, col='darkgray');  abline(h=2, lwd=2.0, col=col);
      }
      
    }
    dev.off()
  }
  
  ##########
  ##### Compare two samples using the same background bg0
  ##########
  ## two background bg_s1 and bg_f
  pdf("PLOTs/bg_s1_AND_bg_f_enrichement_depletion.pdf", width = 22, height = 10)
  par(mfrow=c(1,2))
  output = "Ab_bg_s1_VS_bg0";
  output = basename(xlist[grep(output, xlist)])
  out1 = summary[[which(names(summary)==output)]]
  plot(-log10(out1$pos.score), -log10(out1$neg.score), cex=0.5, main='bg_s1', xlab='-log10(RRA.pos)', ylab='-log10(RRA.neg)')
  abline(v=2, lwd=2.0, col='darkgray')
  abline(h=2, lwd=2.0, col='darkgray');abline(0, 1, col='black')
  jj = which(out1$neg.fdr<0.001); text(-log10(out1$pos.score)[jj], -log10(out1$neg.score)[jj], out1$id[jj], col='red', cex=0.6, pos=4, offset = 0.2)
  jj = which(out1$pos.fdr<0.001); text(-log10(out1$pos.score)[jj], -log10(out1$neg.score)[jj], out1$id[jj], col='blue', cex=0.6, pos=3, offset = 0.2)
  
  output = "Ab_bg_f_VS_bg0";
  output = basename(xlist[grep(output, xlist)])
  out2 = summary[[which(names(summary)==output)]]
  
  out2 = out2[match(out1$id, out2$id),]
  
  plot(-log10(out2$pos.score), -log10(out2$neg.score), cex=0.5, main='bg_f', xlab='-log10(RRA.pos)', ylab='-log10(RRA.neg)')
  abline(v=2, lwd=2.0, col='darkgray')
  abline(h=2, lwd=2.0, col='darkgray');abline(0, 1, col='black')
  jj = which(out2$neg.fdr<0.001); text(-log10(out2$pos.score)[jj], -log10(out2$neg.score)[jj], out2$id[jj], col='red', cex=0.6, pos=4, offset = 0.2)
  jj = which(out2$pos.fdr<0.001); text(-log10(out2$pos.score)[jj], -log10(out2$neg.score)[jj], out2$id[jj], col='blue', cex=0.6, pos=3, offset = 0.2)
  
  dev.off()
  
  pdf("PLOTs/bg_s1_AND_bg_f_logFC.pdf", width = 12, height = 12)
  #sel = which(out1$neg.fdr<0.01 & out2$neg.fdr<0.01)
  fdr.cutoff = 0.01; 
  plot((out1$neg.lfc), (out2$neg.lfc), cex=0.15, xlab="Ab_bg_s1_vs_bg0", ylab= "Ab_bg_f_vs_bg0", xlim = c(-9, 9), ylim = c(-9, 9), main = 'logFC (Ab) ')
  sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col='darkblue', pch=16)
  text((out1$neg.lfc)[sel[1:10]], (out2$neg.lfc)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
  sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='blue', pch=0)
  sel = which(out1$pos.fdr>= fdr.cutoff & out2$pos.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='lightblue', pch=17)
  sel = which(out1$neg.fdr< fdr.cutoff & out2$neg.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col='darkred', pch=16)
  text((out1$neg.lfc)[sel[1:10]], (out2$neg.lfc)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
  sel = which(out1$neg.fdr< fdr.cutoff & out2$neg.fdr>=fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='red', pch=0)
  sel = which(out1$neg.fdr>= fdr.cutoff & out2$neg.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='darksalmon', pch=17)
  abline(0, 1, col='red', lwd=2.0)
  abline(v=0, col='red', lwd=2.0, lty=3);abline(h=0, col='red', lwd=2.0, lty=3)
  legend('topright', col=c('darkblue', 'blue', 'lightblue', 'darkred', 'red', 'darksalmon'), 
         legend = c('signif in both', 'signif in bg_s1', 'signif in bg_f', 'signif in both', 'signif in bg_s1', 'signif in bg_f'), 
         pch = c(16, 0, 17, 16, 0, 17), bty = 'n')
  
  #plot(-log10(out1$neg.score[sel]), -log10(out2$neg.score[sel]), cex=0.6)
  #abline(v=2, lwd=2.0, col='darkgray', lty=1)
  #abline(h=2, lwd=2.0, col='darkgray', lty=2);abline(0, 1, col='black', lty=2)
 dev.off()
 
 pdf("PLOTs/bg_s1_AND_bg_f_RRA_scores.pdf", width = 12, height = 12)
 #sel = which(out1$neg.fdr<0.01 & out2$neg.fdr<0.01)
 fdr.cutoff = 0.01; 
 plot(-log10(out1$pos.score), -log10(out2$pos.score), cex=0.15, xlab="Ab_bg_s1_vs_bg0", ylab= "Ab_bg_f_vs_bg0", 
      xlim = c(0, 32), ylim = c(0, 32), main = 'RRA (Ab) ')
 abline(0, 1, col='red', lwd=2.0)
 abline(v=2, col='red', lwd=2.0, lty=3);abline(h=2, col='red', lwd=2.0, lty=3)
 sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr<fdr.cutoff); points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col='darkblue', pch=16)
 text(-log10(out1$pos.score)[sel[1:10]], -log10(out2$pos.score)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
 sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff); points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col='blue', pch=0)
 text(-log10(out1$pos.score)[sel[1:10]], -log10(out2$pos.score)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
 sel = which(out1$pos.fdr> fdr.cutoff & out2$pos.fdr<fdr.cutoff); points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col='lightblue', pch=17)
 text(-log10(out1$pos.score)[sel[1:10]], -log10(out2$pos.score)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
 #sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='lightblue', pch=0)
 #sel = which(out1$pos.fdr>= fdr.cutoff & out2$pos.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='lightblue', pch=17)
 #dev.off()
 
 sel = which(out1$neg.fdr< fdr.cutoff & out2$neg.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col='darkred', pch=16)
 text((out1$neg.lfc)[sel[1:10]], (out2$neg.lfc)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
 sel = which(out1$neg.fdr< fdr.cutoff & out2$neg.fdr>=fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='red', pch=0)
 sel = which(out1$neg.fdr>= fdr.cutoff & out2$neg.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='darksalmon', pch=17)
 abline(0, 1, col='red', lwd=2.0)
 abline(v=0, col='red', lwd=2.0, lty=3);abline(h=0, col='red', lwd=2.0, lty=3)
 legend('topright', col=c('darkblue', 'blue', 'lightblue', 'darkred', 'red', 'darksalmon'), 
        legend = c('signif in both', 'signif in bg_s1', 'signif in bg_f', 'signif in both', 'signif in bg_s1', 'signif in bg_f'), 
        pch = c(16, 0, 17, 16, 0, 17), bty = 'n')
 
 #plot(-log10(out1$neg.score[sel]), -log10(out2$neg.score[sel]), cex=0.6)
 #abline(v=2, lwd=2.0, col='darkgray', lty=1)
 #abline(h=2, lwd=2.0, col='darkgray', lty=2);abline(0, 1, col='black', lty=2)
 dev.off()
 
 ##########
 ##### The same sample based on two background bg0 and bg_f
 #########
 cc = c("maint.s2.R1", "maint.s2.R2",  "maint.s3")
 pdfname = paste0("PLOTs/Maint_RRA_scores_vs_bg0_bg_t_FDR_0.05", version.analysis, '.pdf')
 pdf(pdfname, width = 14, height = 14)
 fdr.cutoff = 0.05; 
 nb.text = 10;
 for(ss in cc)
 {
   #ss = cc[8];
   output = paste0(gsub("\\.", "_", ss), '_VS_bg0')
   output = basename(xlist[grep(output, xlist)])
   out2 = summary[[which(names(summary)==output)]]
   if(ss=="maint.s3"){
     output = paste0(gsub("\\.", "_", ss), '_VS_bg_s3')
   }else{
     output = paste0(gsub("\\.", "_", ss), '_VS_bg_s2')
   }
   output = basename(xlist[grep(output, xlist)])
   out1 = summary[[which(names(summary)==output)]]
   out2 = out2[match(out1$id, out2$id),]
   #sel = which(out1$neg.fdr<0.01 & out2$neg.fdr<0.01)
   cols = c('chartreuse3', 'blue', 'darksalmon')
   start = 30
   cex.text = 0.6
   for(ii in c(1:2))
   {
     if(ii==1){
       lims = c(0, 35)
     }else{
       lims = c(0, 6)
     }
     plot(-log10(out1$pos.score), -log10(out2$pos.score), cex=0.4, xlab= paste0(gsub("\\.", "_", ss), '_vs_bgf'), 
          ylab= paste0(gsub("\\.", "_", ss), '_vs_bg_0'), xlim = lims, ylim = lims, main = paste0(ss, '--RRA (Ab) '))
     abline(0, 1, col='darkgray', lwd=2.0)
     abline(v=2, col='darkgray', lwd=2.0, lty=3);abline(h=2, col='darkgray', lwd=2.0, lty=3)
     
     if(ii==1){
       sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr<fdr.cutoff); 
       if(length(grep('init', ss))>0){
         sel=which(out1$pos.score< 0.01 & out2$pos.score<0.01); 
       }
       sel0 = sel;nb.text= length(sel);
       if(nb.text>0){
         points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col=cols[1], pch=19)
         text(-log10(out1$pos.score)[sel[1:nb.text]], -log10(out2$pos.score)[sel[1:nb.text]], out1$id[sel[1:nb.text]],  cex=cex.text, pos = 4, offset = 0.2)
         #text(rep(22, nb.text), (start-(1*c(1:nb.text))),   out1$id[sel[order(out1$pos.rank[sel])]], cex = 1, col=cols[1])
         #text(rep(24, nb.text), (start-(1*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 1, col=cols[1])
       }
       sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff); 
       sel1 = sel;nb.text= length(sel);
       if(nb.text>0){
         points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col=cols[2], pch=15)
         #text(rep(26, nb.text), (start-(1*c(1:nb.text))),  out1$id[sel[order(out1$pos.rank[sel])]], cex = 1, col=cols[2])
         text(-log10(out1$pos.score)[sel[1:nb.text]], -log10(out2$pos.score)[sel[1:nb.text]], out1$id[sel[1:nb.text]],  cex=cex.text, pos = 4, offset = 0.3)
       }
       sel = which(out1$pos.fdr> fdr.cutoff & out2$pos.fdr<fdr.cutoff); 
       sel2 = sel; nb.text= length(sel);
       if(nb.text>0){
         points(-log10(out1$pos.score)[sel], -log10(out2$pos.score)[sel], cex=1.0, col=cols[3], pch=17)
         #text(rep(28, nb.text), (start-(1*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 1, col=cols[3])
         text(-log10(out1$pos.score)[sel[1:nb.text]], -log10(out2$pos.score)[sel[1:nb.text]], out1$id[sel[1:nb.text]],  cex=cex.text, pos = 3, offset = 0.3)
       }
     }else{
       kk = which((-log10(out1$pos.score)>2 | (-log10(out2$pos.score))>2))
       text(-log10(out1$pos.score)[kk], -log10(out2$pos.score)[kk], out1$id[kk], cex=0.6, pos = 1, offset = 0.2)
       
     }
     
     #sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='lightblue', pch=0)
     #sel = which(out1$pos.fdr>= fdr.cutoff & out2$pos.fdr<fdr.cutoff); points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='lightblue', pch=17)
     #dev.off()
     legend('topleft', col=cols, legend = c('signif in both', 'signif bg_f', 'signif bg0'), 
            pch = c(19, 15, 17), bty = 'n')
     
   }
   
   
 }
 dev.off()
 
}

### Compare two papulations with the same reference bg0
Compare.populations.using.same.reference = FALSE
if(Compare.populations.using.same.reference)
{
  output = "Ab_bg_s1_VS_bg0";
  output = basename(xlist[grep(output, xlist)])
  out1 = summary[[which(names(summary)==output)]]
  
  output = "Ab_bg_f_VS_bg0";
  output = basename(xlist[grep(output, xlist)])
  out2 = summary[[which(names(summary)==output)]]
  out2 = out2[match(out1$id, out2$id),]
  
  start = 10
  
  pdf("PLOTs/bg_s1_AND_bg_f_logFC_FDR_0.01.pdf", width = 12, height = 12)
  #sel = which(out1$neg.fdr<0.01 & out2$neg.fdr<0.01)
  fdr.cutoff = 0.01; 
  cols = c('darkblue', 'blue', 'lightblue', 'darkred', 'red', 'darksalmon');
  pchs = c(16, 0, 17, 16, 0, 17)
  lims = 10;
  plot((out1$neg.lfc), (out2$neg.lfc), cex=0.15, xlab="Ab_bg_s1_vs_bg0", ylab= "Ab_bg_f_vs_bg0", xlim = c(-lims, lims), ylim = c(-lims, lims), main = 'logFC (Ab) ')
  abline(0, 1, col='darkgray', lwd=2.0)
  abline(v=0, col='darkgray', lwd=2.0, lty=3);abline(h=0, col='darkgray', lwd=2.0, lty=3)
  
  sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr<fdr.cutoff); nb.text= length(sel);
  if(nb.text>0){
    points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col=cols[1], pch=pchs[1])
    text((out1$neg.lfc)[sel[1:nb.text]], (out2$neg.lfc)[sel[1:nb.text]], out1$id[sel[1:nb.text]], cex=0.8, pos = 4, offset = 0.3)
    text(rep(6, nb.text), (start-(0.4*c(1:nb.text))),   out1$id[sel[order(out1$pos.rank[sel])]], cex = 0.8, col=cols[1])
    text(rep(7, nb.text), (start-(0.4*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 0.8, col=cols[1])
  }
  sel = which(out1$pos.fdr< fdr.cutoff & out2$pos.fdr>=fdr.cutoff);
  nb.text= length(sel);
  if(nb.text>0){
    points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col=cols[2], pch=pchs[2])
    text((out1$neg.lfc)[sel[1:nb.text]], (out2$neg.lfc)[sel[1:nb.text]], out1$id[sel[1:nb.text]], cex=0.5, pos = 4, offset = 0.3)
    text(rep(8, nb.text), (start-(0.4*c(1:nb.text))),   out1$id[sel[order(out1$pos.rank[sel])]], cex = 0.6, col=cols[2])
    #text(rep(7, nb.text), (start-(0.5*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 0.7, col=cols[2])
  }
  sel = which(out1$pos.fdr>= fdr.cutoff & out2$pos.fdr<fdr.cutoff); 
  #points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='lightblue', pch=17)
  nb.text= length(sel);
  if(nb.text>0){
    points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col=cols[3], pch=pchs[3])
    text((out1$neg.lfc)[sel[1:nb.text]], (out2$neg.lfc)[sel[1:nb.text]], out1$id[sel[1:nb.text]], cex=0.5, pos = 4, offset = 0.3)
    text(rep(9, nb.text), (start-(0.4*c(1:nb.text))),   out1$id[sel[order(out1$pos.rank[sel])]], cex = 0.6, col=cols[3])
    #text(rep(7, nb.text), (start-(0.5*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 0.7, col=cols[2])
  }
  
  sel = which(out1$neg.fdr< fdr.cutoff & out2$neg.fdr<fdr.cutoff); #points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col='darkred', pch=16)
  nb.text= length(sel);
  if(nb.text>0){
    points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col=cols[4], pch=pchs[4])
    text((out1$neg.lfc)[sel[1:nb.text]], (out2$neg.lfc)[sel[1:nb.text]], out1$id[sel[1:nb.text]], cex=0.5, pos = 4, offset = 0.3)
    text(rep(-6, nb.text), (start-(0.5*c(1:nb.text))),   out1$id[sel[order(out1$pos.rank[sel])]], cex = 0.6, col=cols[4])
    text(rep(-7, nb.text), (start-(0.5*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 0.7, col=cols[4])
  }
  #text((out1$neg.lfc)[sel[1:10]], (out2$neg.lfc)[sel[1:10]], out1$id[sel[1:10]],  cex=0.8, pos = 4, offset = 0.3)
  sel = which(out1$neg.fdr< fdr.cutoff & out2$neg.fdr>=fdr.cutoff); #points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='red', pch=0)
  nb.text= length(sel);
  if(nb.text>0){
    points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col=cols[5], pch=pchs[5])
    text((out1$neg.lfc)[sel[1:nb.text]], (out2$neg.lfc)[sel[1:nb.text]], out1$id[sel[1:nb.text]], cex=0.5, pos = 4, offset = 0.3)
    text(rep(-8, nb.text), (start-(0.4*c(1:nb.text))),   out1$id[sel[order(out1$pos.rank[sel])]], cex = 0.6, col=cols[5])
    #text(rep(7, nb.text), (start-(0.5*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 0.7, col=cols[2])
  }
  sel = which(out1$neg.fdr>= fdr.cutoff & out2$neg.fdr<fdr.cutoff); #points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=0.8, col='darksalmon', pch=17)
  nb.text= length(sel);
  if(nb.text>0){
    points((out1$neg.lfc)[sel], (out2$neg.lfc)[sel], cex=1.0, col=cols[6], pch=pchs[6])
    text((out1$neg.lfc)[sel[1:nb.text]], (out2$neg.lfc)[sel[1:nb.text]], out1$id[sel[1:nb.text]], cex=0.5, pos = 4, offset = 0.3)
    text(rep(-9, nb.text), (start-(0.4*c(1:nb.text))),   out1$id[sel[order(out1$pos.rank[sel])]], cex = 0.6, col=cols[6])
    #text(rep(7, nb.text), (start-(0.5*c(1:nb.text))),   out2$id[sel[order(out2$pos.rank[sel])]], cex = 0.7, col=cols[2])
  }
  
  legend('topright', col=c('darkblue', 'blue', 'lightblue', 'darkred', 'red', 'darksalmon'), 
         legend = c('signif in both', 'signif in bg_s1', 'signif in bg_f', 'signif in both', 'signif in bg_s1', 'signif in bg_f'), 
         pch = c(16, 0, 17, 16, 0, 17), bty = 'n', cex=0.5, pt.cex = 0.5)
  
  #plot(-log10(out1$neg.score[sel]), -log10(out2$neg.score[sel]), cex=0.6)
  #abline(v=2, lwd=2.0, col='darkgray', lty=1)
  #abline(h=2, lwd=2.0, col='darkgray', lty=2);abline(0, 1, col='black', lty=2)
  dev.off()
}


##########################
#########  Compare time courses of samples
##########################
pdf("PLOTs/Kinetic_logFC_bg0_all_samples_FDR_0.25_bg_s1.pdf", width = 12, height = 12)
output = "_VS_bg0"; output = basename(xlist[grep(output, xlist)]);
outs = summary[match(output, names(summary))]
xx = outs[[2]];
kk = which(xx$pos.fdr<0.25)
plot(density(xx$pos.lfc[kk]), type='n')
#points(density(out2$pos.lfc[kk]), col='red', type='l')
#points(density(out1$pos.lfc[kk]), col='blue', type='l')
abline(v=0, col='darkgray', lwd=2.0)
for(n in 1:length(outs))
{
  xx = outs[[n]];
  points(density(xx$pos.lfc[kk]), col=n, type='l', lwd=2.0)
}
legend('topright', col=c(1:length(outs)), 
       legend = names(outs), 
       lty=1, bty = 'n', cex=1.0, pt.cex = 2)

dev.off()

###############################
######## make summary of all comparisons
###############################
res = NULL
for(gg in ggs)
{
  cat(gg, '\n')
  #treatment = unlist(comparison[n])[1]; control = unlist(comparison[n])[2];
  test = c()
  for(n in 2:length(out.names))
  {
    #n = 2
    output = basename(xlist[grep(out.names[n], xlist)])
    out = summary[[which(names(summary)==output)]]
    #out = read.table(file=paste0(getwd(), '/Tables/', output), sep = '\t', header = TRUE)
    kk = which(out$id == gg)
    test = c(test, c(out$pos.p.value[kk], out$pos.rank[kk]))
  }
  res = rbind(res, test)
}
kk = c(1:21)*2-1
xx = data.frame(ggs, res[,kk], res[,-kk],  stringsAsFactors = FALSE)
colnames(xx) = c(out.names[1], paste0(out.names[-1], '.pval'), paste0(out.names[-1], '.rank'))

res = xx;
o1 = order(res$Ab_maint_s2_VS_bg0.pval)
res = res[o1,]

write.table(res, file='Tables/Gene_Summary_pvals_ranks_all_comparisons.txt', sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)

pdf("PLOTs/Gene_summary_aross_all_comparisons.pdf", width = 10, height = 8)
par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,5,2,0.8)+0.1, tcl = -0.3)

jj = grep('pval', colnames(res))
pvals = -log10(as.matrix(res[ ,jj]))
samples = c("Ab_maint_s2", "Ab_maint_s3", "Ab_init_s2_R1", "Ab_init_s3", "Ab_init_s2_R2", "Ab_init_s3_R2")
for(n in 1:nrow(res))
{
  lims = c(0, 6.2)
  plot(c(1,9), lims, xlim = c(0.5, 10), ylim = lims, type='n', main=res$gene[n], xlab=NA, ylab='-log10(pvals)', axes = FALSE)
  mtext(at=c(1:9), side=1, text=c(out.names[c(2:4)], samples), las=3)
  axis(side=2)
  box();abline(h=2, col='darkgray', lwd=2.0)
  cex=1.0
  points(c(1:3), pvals[n, c(1:3)], cex=cex, col='black', type = 'b')
  points(c(4:9), pvals[n, c(4:9)], cex=cex, col='blue', type='b')
  points(c(4:9), pvals[n, c(10:15)], cex=cex, col='green', type='b')
  points(c(4:9), pvals[n, c(16:21)], cex=cex, col='red', type='b')
  legend('topright', col=c('blue', 'green', 'red'), legend = c('bg_s1', 'bg_f', 'bg0'), pch = 1, bty = 'n')
  
}

dev.off()

