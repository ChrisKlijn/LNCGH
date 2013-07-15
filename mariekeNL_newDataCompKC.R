library(KCsmart)
data(hsMirrorLocs)

projectDir <- '~/Projects/NKIProjects/LNCGH/newData/'
setwd(projectDir)

# Load data
load(file=file.path(projectDir, '../newDataImport.RData'))
sampleInfoTumor <- sampleInfo[sampleInfo$Type %in% "Tumor",]

# Exclude sample 592 and 456 as they failed qc

KC    <- segData[,!grepl('592|456', colnames(segData))]
KCLN  <- KC[,grepl('maploc|chrom|LN', colnames(KC))]
KCTum <- KC[,!grepl('LN', colnames(KC))]
colnames(KCLN)  <-  gsub('LN', '', colnames(KCLN))
colnames(KCTum) <-  gsub('Tumor', '', colnames(KCTum))
all.equal(colnames(KCLN),colnames(KCTum))

# Tumor vs. LN

KCcollTLN <- calcSpmCollection(data=KCTum, mirrorLocs=hsMirrorLocs, data2=KCLN)
KCcollTLNcomp <- compareSpmCollection(KCcollTLN, nperms=1000)
sigReg <- getSigRegionsCompKC(KCcollTLNcomp, fdr=.05)

# Tumor vs. LN - Triple Negative

KCcollTLN_TN <- calcSpmCollection(
  data=KCTum[,colnames(KCTum) %in% c('maploc',
                                     'chrom',
                                     sampleInfoTumor$NR[sampleInfoTumor$Mol_subtype %in% 'TN'])], 
  data2=KCLN[,colnames(KCTum) %in% c('maploc',
                                     'chrom',
                                     sampleInfoTumor$NR[sampleInfoTumor$Mol_subtype %in% 'TN'])], 
  mirrorLocs=hsMirrorLocs)
KCcollTLNcomp_TN <- compareSpmCollection(KCcollTLN_TN, nperms=1000)
sigReg_TN <- getSigRegionsCompKC(KCcollTLNcomp_TN, fdr=.05)

# Tumor vs. LN - Triple Negative

KCcollTLN_Lum <- calcSpmCollection(
  data=KCTum[,colnames(KCTum) %in% c('maploc',
                                     'chrom',
                                     sampleInfoTumor$NR[sampleInfoTumor$Mol_subtype %in% 'TN'])], 
  data2=KCLN[,colnames(KCTum) %in% c('maploc',
                                     'chrom',
                                     sampleInfoTumor$NR[sampleInfoTumor$Mol_subtype %in% 'TN'])], 
  mirrorLocs=hsMirrorLocs)
KCcollTLNcomp_Lum <- compareSpmCollection(KCcollTLN_Lum, nperms=1000)
sigReg_Lum <- getSigRegionsCompKC(KCcollTLNcomp_Lum, fdr=.05)


# Plot the results
pdf(file='Figures/compKC_all.pdf', width=15, height=5)
plot(KCcollTLNcomp, sigRegions=sigReg, 
     col1=colors()[122], col2=colors()[148],
     main='Difference Tumor - Lymph Node - All Samples')
abline(h=0)
dev.off()

pdf(file='Figures/compKC_TN.pdf', width=15, height=5)
plot(KCcollTLNcomp_TN, sigRegions=sigReg_TN, 
     col1=colors()[122], col2=colors()[148],
     main='Difference Tumor - Lymph Node - Triple Negative')
abline(h=0)
dev.off()

pdf(file='Figures/compKC_Lum.pdf', width=15, height=5)
plot(KCcollTLNcomp_Lum, sigRegions=sigReg_Lum, 
     col1=colors()[122], col2=colors()[148],
     main='Difference Tumor - Lymph Node - Luminal')
abline(h=0)
dev.off()

save.image('compKC_analysisImage.RData')
