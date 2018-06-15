#####################################
#### data version and analysis version
#####################################
version.data = '_R5213_v1';
version.analysis = paste0(version.data, '_2017_10_30')
nb.samples = 16;
resDir = "../Results/2017-10-30";
if(!dir.exists(resDir)) dir.create(resDir);

###############
### merge replicates for the same UMI
###############
Check.readcount.for.barcode = FALSE
if(Check.readcount.for.barcode)
{
  datafile = "../DATA/out_6_Hagar_1MM_10X_2min_CrUMI_v3.txt"
  tab5rows <- read.table(datafile, header = TRUE, nrows = 5)
  classes <- sapply(tab5rows, class)
  xx <- read.table(datafile, header = TRUE, colClasses = classes)
  
  ### merge replicates for each UMI
  yy = as.matrix(xx[, 3:6])
  
  #aa = yy[, ]
  Tomger = c("Background")
  kk = grep(Tomger, colnames(xx))
  yy = cbind(yy, apply(as.matrix(xx[, kk]), 1, sum))
  colnames(yy)[ncol(yy)] = paste0(Tomger, "_ABC")
  
  Tomger = c("SFMRIV")
  kk = grep(Tomger, colnames(xx))
  yy = cbind(yy, apply(as.matrix(xx[, kk]), 1, sum))
  colnames(yy)[ncol(yy)] = paste0(Tomger, "_ABC")
  
  Tomger = c("BFP")
  kk = grep(Tomger, colnames(xx))
  yy = cbind(yy, apply(as.matrix(xx[, kk]), 1, sum))
  colnames(yy)[ncol(yy)] = paste0("SFRIV_",Tomger, "_ABC")
  
  Tomger = c("GFP")
  kk = grep(Tomger, colnames(xx))
  yy = cbind(yy, apply(as.matrix(xx[, kk]), 1, sum))
  colnames(yy)[ncol(yy)] = paste0("SFRIV_",Tomger, "_ABC")
  
  Tomger = c("SFRIV")
  kk = grep(Tomger, colnames(xx))
  yy = cbind(yy, apply(as.matrix(xx[, kk]), 1, sum))
  colnames(yy)[ncol(yy)] = paste0(Tomger, "GFP_BFP_ABC")
  
  yy = data.frame(xx[, c(1, 2)], yy, stringsAsFactors = FALSE)
  colnames(yy)[ncol(yy)] = "SFRIV_GFP_BFP_ABC"
  
  kk = grep("Background", colnames(xx))
  yy = data.frame(yy, xx[, kk], stringsAsFactors = FALSE)
  colnames(yy)[12:14] = paste0("Background_", c("A", "B", "C"))
  
  #write.table(yy, file=paste0(resDir, "/out_CrUMI_pooled_replicates.txt"), sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
  
  
  ### check the distribution of read counts for control genes
  test = yy[grep("Rnf2", yy$guide), ]
  
  ### UMI counts for each guide 
  guides = unique(yy$guide)
  
  aa = matrix(NA, nrow = length(guides), ncol=(ncol(yy)-2)) ## UMI counts
  colnames(aa) = colnames(yy)[-c(1:2)]
  rownames(aa) = guides
  
  bb = aa; ## read counts
  colnames(bb) = paste0(colnames(bb), ".Ab")
  
  for(n in 1:nrow(aa)){
    # n = 1
    cat(n, "\n")
    kk = which(yy$guide==guides[n])
    test = as.matrix(yy[kk, -c(1:2)])
    # assign read counts
    bb[n, ] = apply(test, 2, sum)
    aa[n, ] = apply(test, 2, function(x) sum(x>0))
  }
  
  colnames(aa) = paste0(colnames(aa), ".In")
  raw = data.frame(aa, bb)
  
  write.table(raw, file=paste0(resDir, "/out_Ab_In_pooled_replicates.txt"), sep='\t', col.names = TRUE, row.names = TRUE, quote=FALSE)
  
  write.table(raw[, c(5:9)], file=paste0(resDir, "/out_In_pooled_replicates_Initiation.txt"), sep='\t', col.names = TRUE, row.names = TRUE, quote=FALSE)
  #rownames(raw) = paste0(xx$guide,"_", xx$barcode) 
  #design.matrix = data.frame(samples = colnames(raw), conditions = as.factor(cc), stringsAsFactors = FALSE)
  
  raw = read.table(file=paste0(resDir, "/out_Ab_In_pooled_replicates.txt"), sep="\t", header = TRUE, row.names = 1)
  
  Clean.sgRNA.Names = TRUE
  if(Clean.sgRNA.Names)
  {
    #load(file=paste0('../Rdata/Screen_countData', version.data, '.Rdata'))
    bb = raw;
    bb = data.frame(sgRNA=rownames(bb), bb, stringsAsFactors = FALSE)
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
    save(bb, file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping_pooled_initiation', version.data, '.Rdata'))
  }
  
}

###############
### Quality control and comapre the merged samples with separated samples
###############
load(file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping_pooled_initiation', version.data, '.Rdata'))

Compare.with.samples.before.pooled = FALSE
if(Compare.with.samples.before.pooled)
{
  aa = bb;
  load(file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  
  cc = "SFRIV"
  jj1 = intersect(grep(cc, colnames(bb)), grep("GFP_In", colnames(bb)))
  jj2 = intersect(grep(cc, colnames(aa)), grep("GFP_ABC.In", colnames(aa)))
  
  plot(bb[, jj1[1]], aa[match(bb$sgRNA, aa$sgRNA), jj2[1]], cex=0.7, log='xy');
  abline(0, 1, lwd=2.0, col='red')
  
}

## general quanlity controls for screen data
source("functions_Screen.R")

pdfname = paste0(resDir, "/SCREEN_Data_Qulity_Assessment_background_Comparison", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 8)

kk = grep('In', colnames(bb))
raw = as.matrix(bb[, kk])
rownames(raw) = bb$sgRNA;

#cc = c("SFMIII_P5_P5_ABC", "SFMIII_NotP5_P5_ABC",  "SFMIII_P6_P7_ABC", "SFMIV_P6_P7_2_ABC", "Background", "Background", "Background",
#       "SFMRIV", "SFMRIV", "SFMRIV",  "SFRIV_BFP", "SFRIV_BFP", "SFRIV_BFP", "SFRIV_GFP",  "SFRIV_GFP",  "SFRIV_GFP")

cc = sapply(colnames(raw), function(x) gsub("*.In", "", x))

design.matrix = data.frame(samples = colnames(raw), conditions = as.factor(cc), stringsAsFactors = FALSE)

#sels = grep("Background|bg", colnames(raw))
sels = c(1:ncol(raw))

Check.RNAseq.Quality(read.count = raw[, sels], design.matrix = design.matrix[sels,]);

dev.off()

###############
### test for hits using mageck
###############
Prepare.Run.mageck = FALSE
if(Prepare.Run.mageck)
{
  load(file=paste0('../Rdata/Screen_countData_sgRNA_Gene_clean_mapping_pooled_initiation', version.data, '.Rdata'))
  
  DIR.mageck = paste0(resDir, "/mageck/");
  system(paste0('mkdir -p ', DIR.mageck))
  DIR.res = "Comparison_Incidence_Initiation_pooled_RRA_0.5/" ## this is path relative to the script "run_mageck.sh"
  system(paste0('mkdir -p ', DIR.mageck, DIR.res))
  
  Which2Compare = 'In';
  
  ## prepare table for mageck
  kk = grep(Which2Compare, colnames(bb))[4:8]
  dds =  data.frame(bb[, c(1, 2, kk)], stringsAsFactors = FALSE)
  
  cc = sapply(colnames(dds)[-c(1:2)], function(x) gsub("*.In", "", x))
  
  print(colnames(dds))

  ### define comparisons for mageck
  find.conditions = function(x){
    return(unlist(strsplit(as.character(x), '_'))[1]);
  }
  #cc = unique(sapply(colnames(dds)[-c(1:2)], find.conditions)) 
  
  conds = c('sgRNA', 'Gene', cc)
  
  colnames(dds) = factor(conds)
  
  ## define pairwise comparisons
  comparison = list()
  n = 1
  for(ss in unique(conds)){
    cat(ss, "\n")
    if(ss != "sgRNA" & ss !="Gene" & length(grep("Background|bg_s", ss))==0)
    {
      comparison[[n]] = c(ss, "Background_ABC")
      n = n +1;
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
        cmd = paste0('mageck test -k ', 'sgrna_count_table.txt', ' -t ', samples.t,  ' -c ',  samples.c, ' -n ',  output, ' --adjust-method pounds --sort-criteria pos --pdf-report --normcounts-to-file --norm-method median --remove-zero treatment --additional-rra-parameters p=0.5')
        #cat(cmd, '\n')
        write(cmd, file=logs)
      }else{
        cmd = paste0('mageck test -k ', 'sgrna_count_table.txt', ' -t ', samples.t,  ' -c ',  samples.c, ' -n ',  output, ' --adjust-method pounds --sort-criteria pos --pdf-report --normcounts-to-file --norm-method median --remove-zero treatment --additional-rra-parameters p=0.5')
        #cat(cmd, '\n')
        write(cmd, file=logs, append = TRUE)
      }
    }  
  }
  
}

###############################
#### Volcano plots 
###############################
#Which2Compare = 'Ab';
#paste0(Which2Compare, "_*.gene_summary.txt")
ylist<-list.files(path=paste0(resDir, '/mageck/Comparison_Incidence_Initiation_pooled'),
                  pattern = "*.gene_summary.txt", full.names = TRUE)

#ylist = ylist[grep("SFM", ylist)]
bnames = basename(ylist)
bnames = sapply(bnames, function(x) unlist(strsplit(as.character(x), "[.]"))[1])

summary = list();
#bnames = basename(xlist)
for(n in 1:length(ylist))
{
  ff = read.table(ylist[n], sep='\t', header = TRUE);
  summary[[n]] = data.frame(ff, stringsAsFactors = FALSE);
};
names(summary) = bnames
ggs = c("Cbx7", "Rnf2", "Suz12", "Ezh2", "Eed", "Mtf2")

RRA.ranks.volcano.plot.single.comparion = FALSE
if(RRA.ranks.volcano.plot.single.comparion)
{
  ## RRA rank plot
  TOplots = names(summary)
  
  rra.cutoff = 0.01;cex=0.7;pch=16;
  col= 'gray'
  pdfname = paste0(resDir, "/RRAscores_ranks.pdf")
  pdf(pdfname, width = 8, height = 6)
  par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,1.5,3)+0.1, tcl = -0.3)
  for(n in 1:length(TOplots))
  {
    #n = 1;
    output = TOplots[n];
    #out1 = basename(xlist[grep(output, xlist)])
    out1 = summary[[which(names(summary)==output)]]
    
    cex.text=0.8;cex.offset=0.4; cex=0.7;
    ## RRA scores and ranks
    plot(out1$pos.rank, -log10(out1$pos.score), cex=cex, type='p', main=output, xlab='RRA ranks', ylab='-log10(RRA)')
    kk1 = which(out1$pos.rank<=10);
    text(out1$pos.rank[kk1], -log10(out1$pos.score[kk1]), out1$id[kk1], cex=cex.text, offset = cex.offset, pos = 4, col = 'blue');
    mm1 = match(ggs, out1$id)
    text(out1$pos.rank[mm1], -log10(out1$pos.score[mm1]), out1$id[mm1], cex=cex.text, offset = cex.offset, pos = 4, col = 'darkgreen');
  }
  dev.off()
  
  ## volcano plot
  Label.gene.names = FALSE
  
  pdfname = paste0(resDir, "/volcano_plots_without_geneNames_labeling.pdf")
  pdf(pdfname, width = 8, height = 6)
  par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,1.5,3)+0.1, tcl = -0.3)
  for(n in 1:length(TOplots))
  {
    #n = 1;
    output = TOplots[n];
    #out1 = basename(xlist[grep(output, xlist)])
    out1 = summary[[which(names(summary)==output)]]
    
    #cex.text=0.7;cex.offset=0.4; cex=0.8;
    rra.cutoff = 0.01;pch=16;
    col= 'gray'
    lims = c(-7, 7);
    if(n>1) {lims = max(abs(out1$pos.lfc));lims = c(-lims, lims);}
    cex.text=0.6;
    cex.offset=0.2; 
    cex=0.8;
    cex.bg = 0.2
    col = 'darkblue'; 
    pval.cutoff = 0.05;
    #ylims = c(0, 24)
    for(m in c(1))
    {
      plot(out1$pos.lfc, -log10(out1$pos.p.value), cex=cex.bg, type='p', main=output, xlab='log(fold change)', ylab='-log10(pvalue)', xlim = lims, col='gray')
      kk = which(out1$pos.lfc>0);
      kk1 = kk[which(out1$pos.p.value[kk]<pval.cutoff)]; 
      #kk2 = kk[which(out1$pos.p.value[kk]>=pval.cutoff)]; 
      points(out1$pos.lfc[kk1], -log10(out1$pos.p.value[kk1]), cex=cex, col=col, pch=pch);
      if(Label.gene.names){
        text(out1$pos.lfc[kk1], -log10(out1$pos.p.value[kk1]), out1$id[kk1], cex=cex.text, offset = cex.offset, pos = 4);
        mm1 = match(ggs, out1$id)
        text(out1$pos.lfc[mm1], -log10(out1$pos.p.value[mm1]), out1$id[mm1], cex=cex.text, offset = cex.offset, pos = 4, col = 'darkgreen'); 
      }
      abline(v=0, lwd=2.0, col='darkgray');  abline(h=-log10(pval.cutoff), lwd=2.0, col=col);
      
      plot(out1$pos.lfc, -log10(out1$pos.p.value), cex=cex.bg, type='p', main=output, xlab='log(fold change)', ylab='-log10(pvalue)',
           xlim = c(0, max(-log10(out1$pos.p.value))), col='gray')
      kk = which(out1$pos.lfc>0);
      kk1 = kk[which(out1$pos.p.value[kk]<pval.cutoff)]; 
      #kk2 = kk[which(out1$pos.p.value[kk]>=pval.cutoff)]; 
      points(out1$pos.lfc[kk1], -log10(out1$pos.p.value[kk1]), cex=cex, col=col, pch=pch);
      if(Label.gene.names){
        text(out1$pos.lfc[kk1], -log10(out1$pos.p.value[kk1]), out1$id[kk1], cex=cex.text, offset = cex.offset, pos = 4);
        mm1 = match(ggs, out1$id)
        text(out1$pos.lfc[mm1], -log10(out1$pos.p.value[mm1]), out1$id[mm1], cex=cex.text, offset = cex.offset, pos = 4, col = 'darkgreen');
      }
    
      abline(v=0, lwd=2.0, col='darkgray');  abline(h=-log10(pval.cutoff), lwd=2.0, col=col);
      
    }
  }
  dev.off()
  
}

###############################
#### Compare Hp screen an Cbx7 screen
###############################
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


ylist<-list.files(path=paste0(resDir, '/mageck/Comparison_Incidence_Initiation_pooled'),
                  pattern = "*.gene_summary.txt", full.names = TRUE)

#ylist = ylist[grep("SFM", ylist)]
bnames = basename(ylist)
bnames = sapply(bnames, function(x) unlist(strsplit(as.character(x), "[.]"))[1])

cutoff = 0.01
pdfname = paste0(resDir, "/Compare_HP_Cbx7_screen_candidates",  version.analysis, ".pdf")
pdf(pdfname, width = 8, height = 8)
for(n in 1:length(ylist))
{
  #n = 1;
  cex.text = 0.5;
  out = read.table(ylist[n], sep='\t', header = TRUE)
  out = out[match(names(hits.rk), out$id), ]
  lims = -log10(min(c(out$pos.score, hits.rk)))
  jj = which((out$pos.score<cutoff) | (hits.rk<cutoff))
  jj1 = which((out$pos.score<cutoff & out$pos.score<10^-5) | (hits.rk<cutoff & hits.rk<10^-5))
  plot(-log10(out$pos.score)[jj], -log10(hits.rk)[jj], cex=0.5, xlab="Cbx (RRA)", ylab= "Hp1 (RRA)", 
       xlim = c(0, lims), ylim = c(0, lims), main = paste0('RRA score for ', bnames[n]), type='p')
  abline(0, 1, col='darkgray', lwd=2.0)
  abline(v=2, col='darkgray', lwd=1.0, lty=2);abline(h=2, col='darkgray', lwd=1.0, lty=2)
  text(-log10(out$pos.score)[jj1], -log10(hits.rk)[jj1], out$id[jj1],  cex=cex.text, pos = 4, offset = 0.3)
  
  jj = which(out$pos.score<cutoff | hits.rk<cutoff)
  jj1 = which((out$pos.score<cutoff & out$pos.score>10^-5) | (hits.rk<cutoff & hits.rk>10^-5))
  plot(-log10(out$pos.score)[jj], -log10(hits.rk)[jj], cex=0.5, xlab="Cbx (RRA)", ylab= "Hp1 (RRA)", 
       xlim = c(0, 5), ylim = c(0, 5), main = paste0('RRA score for ', bnames[n]), type='p')
  abline(0, 1, col='darkgray', lwd=2.0)
  abline(v=2, col='darkgray', lwd=1.0, lty=2);abline(h=2, col='darkgray', lwd=1.0, lty=2)
  text(-log10(out$pos.score)[jj], -log10(hits.rk)[jj], out$id[jj],  cex=cex.text, pos = 4, offset = 0.3)
  
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



###############
### test for hits using gCrisprTools
###############
toCompare = c("SFMRIV_ABC.In", "SFRIV_BFP_ABC.In", "SFRIV_GFP_ABC.In","SFRIV_GFP_BFP_ABC.In")
require(gCrisprTools)

for(n in 1:length(toCompare))
{
  # n = 1
  jj1 = which(colnames(raw)==toCompare[n])
  jj2 = which(colnames(raw)=="Background_ABC.In")
  fcs = (raw[, jj1])/sum(raw[, jj1])/((raw[, jj2]+10^-2)/sum(raw[, jj2]))
  test = raw[, c(jj2, jj1)]
  test = test[order(-fcs), ]
  
}

### check control genes
control.list = c("Cbx7", "Rnf2", "Eed", "Suz12", "Ezh2", "Mtf2")

source("functions_Screen.R")
pdfname = paste0(resDir, "/SCREEN_Controles_list_readCounts_UMI", version.analysis, ".pdf")
pdf(pdfname, width = 20, height = 10)
Check.positve.controls.sgRNA.UMI(read.count=raw, design.matrix=design.matrix, control.list=control.list)
dev.off()




