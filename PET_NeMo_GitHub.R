# load the library  ####
library(stringr)
library(R.matlab)
library(viridis)
library("lattice")
library(robustbase)
library(matrixStats)
library(vioplot)
library(perm)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(ggsignif) 
library(RNifti)
library(oro.nifti)
library(ggpubr)

# Create the average PET images for controls at voxel level ####
PET_mean_HC<-array(NA,c(182,218,182))
PET_sd_HC<-array(NA,c(182,218,182))

for(i in 1:182){(cat(i))
  for(k in 1:218){
    for(m in 1:182){
      PET_n16_HC_v2<-NULL
      for(pat in 1:16){
        PET_n16_HC_v2<-rbind(PET_n16_HC_v2,PET_n16_HC[[pat]][i,k,m])
      }
      PET_mean_HC[i,k,m]<-mean(PET_n16_HC_v2)
      PET_sd_HC[i,k,m]<-sd(PET_n16_HC_v2)
      
    }
  }
}


# Create the average PET images for MS patients at voxel level ####

PET_mean_MS_PRL_n26<-array(NA,c(182,218,182))
PET_mean_MS_nonPRL_n18<-array(NA,c(182,218,182))

# z-scored PET
PET_mean_MS_PRL_n26_zscored<-array(NA,c(182,218,182))
PET_mean_MS_nonPRL_n18_zscored<-array(NA,c(182,218,182))

PET_n26_PRL_voxelwisezscoredHC <- vector("list", length(PET_n26_allsubjects))
PET_n18_nonPRL_voxelwisezscoredHC <- vector("list", length(PET_n18_allsubjects))

for (i in 1:length(PET_n26_allsubjects)) {PET_n26_PRL_voxelwisezscoredHC[[i]] <- array(NA, dim = c(182,218,182))}
for (i in 1:length(PET_n18_allsubjects)) {PET_n18_nonPRL_voxelwisezscoredHC[[i]] <- array(NA, dim = c(182,218,182))}

# empty matrices for the Student's t-stat results
PET_MSvsHC_tstat_estimate<-array(NA,c(182,218,182))
PET_MSvsHC_tstat_pvalue<-array(NA,c(182,218,182))

PET_PRLvsnonPRL_tstat_estimate<-array(NA,c(182,218,182))
PET_PRLvsnonPRL_tstat_pvalue<-array(NA,c(182,218,182))

# empty matrices for ANCOVA results
PET_MSvsHC_ancova_age_sex_estimate<-array(NA,c(182,218,182))
PET_MSvsHC_ancova_age_sex_pval<-array(NA,c(182,218,182))

PET_PRLvsnonPRL_ancova_estimate<-array(NA,c(182,218,182))
PET_PRLvsnonPRL_ancova_pval<-array(NA,c(182,218,182))

for(i in 1:182){(cat(i))
  for(k in 1:218){
    for(m in 1:182){
      
      # subjects with PRL
      PETsubjectswithPRL<-NULL
      PETsubjectswithPRL_zscored<-NULL
      for(pat in 1:length(PET_n26_allsubjects)){
        PETsubjectswithPRL<-rbind(PETsubjectswithPRL,PET_n26_allsubjects[[pat]][i,k,m])
        PET_n26_PRL_voxelwisezscoredHC[[pat]][i,k,m]<-(PET_n26_allsubjects[[pat]][i,k,m]-PET_mean_HC[i,k,m])/PET_sd_HC[i,k,m]
        PETsubjectswithPRL_zscored<-rbind(PETsubjectswithPRL_zscored,PET_n26_PRL_voxelwisezscoredHC[[pat]][i,k,m])
        
      }
      PET_mean_MS_PRL_n26[i,k,m]<-mean(PETsubjectswithPRL)
      PET_mean_MS_PRL_n26_zscored[i,k,m]<-mean(PETsubjectswithPRL_zscored)
      
      # subjects without PRL
      PETsubjectswithoutPRL<-NULL
      PETsubjectswithoutPRL_zscored<-NULL
      
      for(pat2 in 1:length(PET_n18_allsubjects)){
        PETsubjectswithoutPRL<-rbind(PETsubjectswithoutPRL,PET_n18_allsubjects[[pat2]][i,k,m])
        PET_n18_nonPRL_voxelwisezscoredHC[[pat2]][i,k,m]<-(PET_n18_allsubjects[[pat2]][i,k,m]-PET_mean_HC[i,k,m])/PET_sd_HC[i,k,m]
        PETsubjectswithoutPRL_zscored<-rbind(PETsubjectswithoutPRL_zscored,PET_n18_nonPRL_voxelwisezscoredHC[[pat2]][i,k,m])
        
      }
      PET_mean_MS_nonPRL_n18[i,k,m]<-mean(PETsubjectswithoutPRL)
      PET_mean_MS_nonPRL_n18_zscored[i,k,m]<-mean(PETsubjectswithoutPRL_zscored)
      
      # Student's t-test comparing HC vs MS
      PET_MS44<-NULL
      for(pat3 in 1:length(PET_n44_allsubjects)){
        PET_MS44<-rbind(PET_MS44,PET_n44_allsubjects[[pat3]][i,k,m])
      }
      
      PET_HC16<-NULL
      for(hcno in 1:length(PET_n16_HC)){
        PET_HC16<-rbind(PET_HC16,PET_n16_HC[[hcno]][i,k,m])
      }
      PET_MSvsHC_tstat_estimate[i,k,m]<-t.test(PET_MS44,PET_HC16)$statistic
      PET_MSvsHC_tstat_pvalue[i,k,m]<-t.test(PET_MS44,PET_HC16)$p.value
      
      # ANCOVA comparing HC vs MS 
      # Reorder MS_PET_demo based on the order of ID_n44
      MS_PET_demo_reordered <- MS_PET_demo[match(ID_n44, MS_PET_demo$ID), ]
      data1<-rbind(data.frame(metric=PET_HC16,age=HC_PET_demo2$Age,sex=HC_PET_demo2$Sex,group=rep("HC",16)),
                   data.frame(metric=PET_MS44,age=MS_PET_demo_reordered$Age.at.Base,sex=MS_PET_demo_reordered$Gender,group=rep("MS",44)))
      PET_MSvsHC_ancova_age_sex_estimate[i,k,m]<-summary(aov(metric~group+age+sex,data=data1))[[1]][, "F value"][1]
      PET_MSvsHC_ancova_age_sex_pval[i,k,m]<-summary(aov(metric~group+age,data=data1))[[1]][, "Pr(>F)"][1]
      
      # Student's t-test comparing MS patients with vs without PRL
      PET_PRL_n26<-NULL
      for(pat3 in 1:length(PET_n26_allsubjects)){
        PET_PRL_n26<-rbind(PET_PRL_n26,PET_n26_allsubjects[[pat3]][i,k,m])
      }
      
      PET_nonPRL_n18<-NULL
      for(hcno in 1:length(PET_n18_allsubjects)){
        PET_nonPRL_n18<-rbind(PET_nonPRL_n18,PET_n18_allsubjects[[hcno]][i,k,m])
      }
      PET_PRLvsnonPRL_tstat_estimate[i,k,m]<-t.test(PET_PRL_n26,PET_nonPRL_n18)$statistic
      PET_PRLvsnonPRL_tstat_pvalue[i,k,m]<-t.test(PET_PRL_n26,PET_nonPRL_n18)$p.value
      
      # ANCOVA comparing MS patients PRL vs nonPRL
      # Reorder MS_PET_demo based on the order of ID_n44
      ID_MS_combined<-c(ID_PET_n26_allsubjects_july2024,ID_PET_n18_allsubjects_july2024)
      MS_PET_demo_reordered2 <- MS_PET_demo[match(ID_MS_combined, as.numeric(MS_PET_demo$ID)), ];MS_PET_demo_reordered2$ID
      data1<-rbind(data.frame(metric=PET_PRL_n26,age=MS_PET_demo_reordered2$Age.at.Base[1:26],group=rep("PRL",26)),
                   data.frame(metric=PET_nonPRL_n18,age=MS_PET_demo_reordered2$Age.at.Base[27:44],group=rep("nonPRL",18)))
      PET_PRLvsnonPRL_ancova_estimate[i,k,m]<-summary(aov(metric~group+age,data=data1))[[1]][, "F value"][1]
      PET_PRLvsnonPRL_ancova_pval[i,k,m]<-summary(aov(metric~group+age,data=data1))[[1]][, "Pr(>F)"][1]
      
    }
  }
}



# load the WM mask which was created using 16 HC, the voxels within 1mm proximity from the CSF and GM were removed ####
wm_mask_aseg_0dot5thr_bin_16HC_ero1mm <- readNifti("wm_mask_aseg_0dot5thr_bin_16HC_ero1mm.nii.gz") # matrix with a dimension of 182x218x182

# z-scored PET at voxel level ####
# remove the voxels within 1mm proximity of the CSF and GM for the z-scored DVR
PET_mean_MS_PRL_n26_zscored_v2<-replace(PET_mean_MS_PRL_n26_zscored,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)
PET_mean_MS_nonPRL_n18_zscored_v2<-replace(PET_mean_MS_nonPRL_n18_zscored,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)

# remove the voxels within 1mm proximity of the CSF and GM for the t-statistics and p-value comparing HC vs MS
PET_MSvsHC_tstat_estimate_v2<-replace(PET_MSvsHC_tstat_estimate,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)
PET_MSvsHC_tstat_pvalue_v2<-replace(PET_MSvsHC_tstat_pvalue,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)

# remove the voxels within 1mm proximity of the CSF and GM for the t-statistics and p-value comparing patients with vs without PRL
PET_PRLvsnonPRL_tstat_estimate_v2<-replace(PET_PRLvsnonPRL_tstat_estimate,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)
PET_PRLvsnonPRL_tstat_pvalue_v2<-replace(PET_PRLvsnonPRL_tstat_pvalue,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)

# remove the voxels within 1mm proximity of the CSF and GM for the ANCOVA derived F-value and p-value comparing HC vs MS
PET_MSvsHC_ancova_age_sex_estimate_v2<-replace(PET_MSvsHC_ancova_age_sex_estimate,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)
PET_MSvsHC_ancova_age_sex_pval_v2<-replace(PET_MSvsHC_ancova_age_sex_pval,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)

# remove the voxels within 1mm proximity of the CSF and GM for the ANCOVA derived F-value and p-value comparing patients with vs without PRL
PET_PRLvsnonPRL_ancova_estimate_v2<-replace(PET_PRLvsnonPRL_ancova_estimate,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)
PET_PRLvsnonPRL_ancova_pvalue_v2<-replace(PET_PRLvsnonPRL_ancova_pval,which(wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==0),0)

# create nifti images
setwd("NEWDIRECTORY")

writeNIfTI(PET_mean_MS_PRL_n26_zscored_v2, "PET_mean_MS_PRL_n26_zscored_v2")
writeNIfTI(PET_mean_MS_nonPRL_n18_zscored_v2, "PET_mean_MS_nonPRL_n18_zscored_v2")

#t-stat
writeNIfTI(PET_MSvsHC_tstat_estimate_v2, "PET_MSvsHC_tstat_estimate_v2")
writeNIfTI(PET_MSvsHC_tstat_pvalue_v2, "PET_MSvsHC_tstat_pvalue_v2")

writeNIfTI(PET_PRLvsnonPRL_tstat_estimate_v2, "PET_PRLvsnonPRL_tstat_estimate_v2")
writeNIfTI(PET_PRLvsnonPRL_tstat_pvalue_v2, "PET_PRLvsnonPRL_tstat_pvalue_v2")

# ancova
writeNIfTI(PET_MSvsHC_ancova_age_sex_estimate_v2, "PET_MSvsHC_ancova_age_sex_estimate_v2")
writeNIfTI(PET_MSvsHC_ancova_age_sex_pval_v2, "PET_MSvsHC_ancova_age_sex_pval_v2")

writeNIfTI(PET_PRLvsnonPRL_ancova_estimate_v2, "PET_PRLvsnonPRL_ancova_estimate_v2")
writeNIfTI(PET_PRLvsnonPRL_ancova_pvalue_v2, "PET_PRLvsnonPRL_ancova_pvalue_v2")

# RUN this section on Terminal
cd NEWDIRECTORY

# create the PET images by removing the voxels with a value of zero
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_mean_MS_PRL_n26_zscored_v2.nii.gz PET_mean_MS_PRL_n26_zscored_v2.nii.gz
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_mean_MS_nonPRL_n18_zscored_v2.nii.gz PET_mean_MS_nonPRL_n18_zscored_v2.nii.gz
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_MSvsHC_tstat_estimate_v2.nii.gz PET_MSvsHC_tstat_estimate_v2.nii.gz
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_MSvsHC_tstat_pvalue_v2.nii.gz PET_MSvsHC_tstat_pvalue_v2.nii.gz
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_PRLvsnonPRL_tstat_estimate_v2.nii.gz PET_PRLvsnonPRL_tstat_estimate_v2.nii.gz
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_PRLvsnonPRL_tstat_pvalue_v2.nii.gz PET_PRLvsnonPRL_tstat_pvalue_v2.nii.gz

fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_MSvsHC_ancova_age_sex_estimate_v2 PET_MSvsHC_ancova_age_sex_estimate_v2
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_MSvsHC_ancova_age_sex_pval_v2 PET_MSvsHC_ancova_age_sex_pval_v2
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_PRLvsnonPRL_ancova_estimate_v2.nii.gz PET_PRLvsnonPRL_ancova_estimate_v2.nii.gz
fslmaths /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii.gz -mul 0 -add PET_PRLvsnonPRL_ancova_pvalue_v2.nii.gz PET_PRLvsnonPRL_ancova_pvalue_v2.nii.gz

# visualize - Run this section on Terminal
fsleyes /usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz PET_mean_MS_PRL_n26_zscored_v2.nii.gz PET_mean_MS_nonPRL_n18_zscored_v2.nii.gz PET_MSvsHC_tstat_estimate_v2.nii.gz PET_MSvsHC_tstat_pvalue_v2.nii.gz PET_PRLvsnonPRL_tstat_estimate_v2.nii.gz PET_PRLvsnonPRL_tstat_pvalue_v2.nii.gz
fsleyes /usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz PET_MSvsHC_ancova_age_sex_estimate_v2 PET_MSvsHC_ancova_age_sex_pval_v2 PET_PRLvsnonPRL_ancova_estimate_v2.nii.gz PET_PRLvsnonPRL_ancova_pvalue_v2.nii.gz

# FIGURE 1 - PANEL D ####


# load the matrices containing the PET metrics for HC and MS patienrs
load(file="PET_n16_HC_ontheFA1mmexcl.RData")
load(file="meanPET_n44_MS_ontheFA1mmexcl_july2024.RData")
load(file="meanPET_n26_MS_PRL_ontheFA1mmexcl_july2024.RData")
load(file="meanPET_n18_MS_woutPRL_ontheFA1mmexcl_july2024.RData")

# create a data frame containing the DVR metrics per subject and a Group column indicating the group assignment
PET_n16_HC_ontheFA1mmexcl_june2024_v2 <- data.frame(Value = c(PET_n16_HC_ontheFA1mmexcl), Group = 'HC')
meanPET_n44_MS_ontheFA1mmexcl_july2024_v2 <- data.frame(Value = c(meanPET_n44_MS_ontheFA1mmexcl_july2024), Group = 'MS')
meanPET_n26_MS_PRL_ontheFA1mmexcl_july2024_v2 <- data.frame(Value = c(meanPET_n26_MS_PRL_ontheFA1mmexcl_july2024), Group = 'PRLMS')
meanPET_n18_MS_woutPRL_ontheFA1mmexcl_july2024_v2 <- data.frame(Value = c(meanPET_n18_MS_woutPRL_ontheFA1mmexcl_july2024), Group = 'nonPRLMS')

# Combine the all data frames
combined_data_1mmexcl <- rbind(PET_n16_HC_ontheFA1mmexcl_june2024_v2, 
                               meanPET_n44_MS_ontheFA1mmexcl_july2024_v2,
                               meanPET_n26_MS_PRL_ontheFA1mmexcl_july2024_v2,
                               meanPET_n18_MS_woutPRL_ontheFA1mmexcl_july2024_v2)

# Count number of subjects in each group
hc_count <- nrow(PET_n16_HC_ontheFA1mmexcl_june2024_v2);hc_count
ms_count <- nrow(meanPET_n44_MS_ontheFA1mmexcl_july2024_v2);ms_count
prlms_count <- nrow(meanPET_n26_MS_PRL_ontheFA1mmexcl_july2024_v2);prlms_count
nonprlms_count <- nrow(meanPET_n18_MS_woutPRL_ontheFA1mmexcl_july2024_v2);nonprlms_count

# reorder the groups
combined_data_1mmexcl$Group <- factor(combined_data_1mmexcl$Group, levels = c("HC", "MS", "PRLMS","nonPRLMS"))

# combine PET images with demographics
combined_data_1mmexcl_withdemo<-cbind(combined_data_1mmexcl,rbind(HC_PET_demo3,MS_PET_demo4)[,2:3])

# perform ANCOVA to compare the groups and save the p-values
combined_data_1mmexcl_withdemo_ms_hc<-combined_data_1mmexcl_withdemo[which(combined_data_1mmexcl_withdemo$Group %in% c("HC","MS")),]
dim(combined_data_1mmexcl_withdemo_ms_hc)
p_hc_ms<-summary(aov(combined_data_1mmexcl_withdemo_ms_hc[,1]~Group+Age+Sex,data=combined_data_1mmexcl_withdemo_ms_hc))[[1]][, "Pr(>F)"][1]

combined_data_1mmexcl_withdemo_hc_PRLMS<-combined_data_1mmexcl_withdemo[which(combined_data_1mmexcl_withdemo$Group %in% c("HC","PRLMS")),]
dim(combined_data_1mmexcl_withdemo_hc_PRLMS)
p_hc_prlms<-summary(aov(combined_data_1mmexcl_withdemo_hc_PRLMS[,1]~Group+Age+Sex,data=combined_data_1mmexcl_withdemo_hc_PRLMS))[[1]][, "Pr(>F)"][1]

combined_data_1mmexcl_withdemo_hc_nonPRLMS<-combined_data_1mmexcl_withdemo[which(combined_data_1mmexcl_withdemo$Group %in% c("HC","nonPRLMS")),]
dim(combined_data_1mmexcl_withdemo_hc_nonPRLMS)
p_hc_nonprlms<-summary(aov(combined_data_1mmexcl_withdemo_hc_nonPRLMS[,1]~Group+Age+Sex,data=combined_data_1mmexcl_withdemo_hc_nonPRLMS))[[1]][, "Pr(>F)"][1]

combined_data_1mmexcl_withdemo_ms<-combined_data_1mmexcl_withdemo[which(combined_data_1mmexcl_withdemo$Group %in% c("PRLMS","nonPRLMS")),]
dim(combined_data_1mmexcl_withdemo_ms)
p_nonprl_prlms<-summary(aov(combined_data_1mmexcl_withdemo_ms[,1]~Group+Age+Sex,data=combined_data_1mmexcl_withdemo_ms))[[1]][, "Pr(>F)"][1]

# run the violinplot
textsize1<-6
ggplot(combined_data_1mmexcl, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.6) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Overlay boxplot
  theme(axis.title.x = element_text(size = 20),  # Adjust size for x-axis
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        legend.key.height = unit(2, "lines"),  # Adjust height as needed
        legend.key.width = unit(2, "lines"))+
  theme(legend.position = "none") +
  labs(title = "", 
       x = "", 
       y = "Average DVR in the WM") +
  scale_x_discrete(labels = c(paste0("HC \n (n=", hc_count, ")"), 
                              paste0("MS \n (n=", ms_count, ")"),
                              paste0("PRLMS \n (n=", prlms_count, ")"),
                              paste0("nonPRLMS \n (n=", nonprlms_count, ")"))) +
  geom_signif(comparisons = list(c("HC", "MS")), annotations=paste0("pBH=",round(p_hc_ms,3)), y_position = 1.4, tip_length = 0.03, textsize = textsize1)+
  geom_signif(comparisons = list(c("HC", "PRLMS")), annotations=paste0("pBH=",round(p_hc_prlms,3)), y_position = 1.45, tip_length = 0.03, textsize = textsize1)+
  geom_signif(comparisons = list(c("HC", "nonPRLMS")), annotations=paste0("pBH=",round(p_hc_nonprlms,3)), y_position = 1.5, tip_length = 0.03, textsize = textsize1)+
  geom_signif(comparisons = list(c("PRLMS", "nonPRLMS")), annotations=paste0("pBH=",round(p_nonprl_prlms,3)), y_position = 1.4, tip_length = 0.03, textsize = textsize1)


# FIGURE 3 - Average z-scored PET at subject level for different ChaCo level####

## Patients with PRL ####
# Create empty matrices
{

  PET_zscoredHC_outsidelesionmask_chacoZERO_n26<-NULL
  PET_zscoredHC_outsidelesionmask_chacoPRLZERO_n26<-NULL
  PET_zscoredHC_outsidelesionmask_chacoNONPRLZERO_n26<-NULL
  
  PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n26<-NULL
  PET_zscoredHC_outsidelesionmask_chacoPRL_0_0.1_n26<-NULL
  PET_zscoredHC_outsidelesionmask_chacoNONPRL_0_0.1_n26<-NULL
  
  PET_zscoredHC_outsidelesionmask_chaco09_1_n26<-NULL
  PET_zscoredHC_outsidelesionmask_chacoPRL09_1_n26<-NULL
  PET_zscoredHC_outsidelesionmask_chacoNONPRL09_1_n26<-NULL
  
}

k=1
for(i in ID_WM_Chaco_rimpos_n26){cat(k);k=k+1

# WM disruption due to rim+ and rim- lesions
ChaCoduetoRIMPOS<-WM_Chaco_rimpos_n26[[which(ID_WM_Chaco_rimpos_n26==i)]];dim(ChaCoduetoRIMPOS)
ChaCoduetoRIMNEG<-WM_Chaco_rimneg_n26[[which(ID_WM_Chaco_rimneg_n26==i)]]

# z-scored DVR for this subject
PET_voxelwise_zscoredHC<-PET_n26_PRL_voxelwisezscoredHC[[which(ID_PET_n26_allsubjects_july2024==i)]];dim(PET)

# identify the rim+ and rim- lesion masks for this patient
lesionmask1<-lesionmask_rimpos_n26_v2[[which(ID_lesionmask_rimpos_n26_v2==i)]];dim(lesionmask1)
lesionmask2<-lesionmask_rimneg_n26_v2[[which(ID_lesionmask_rimneg_n26_v2==i)]];dim(lesionmask2)

# WM disruption outside of the lesion mask and 1mm trimmed from CSF and GM
ChaCoduetoRIMPOS_outsidelesionmask<-ChaCoduetoRIMPOS[which(lesionmask1==0 & lesionmask2==0 & wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==1)];length(ChaCoduetoRIMPOS_outsidelesionmask)
ChaCoduetoRIMNEG_outsidelesionmask<-ChaCoduetoRIMNEG[which(lesionmask1==0 & lesionmask2==0 & wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==1)];length(ChaCoduetoRIMNEG_outsidelesionmask)

# z-scored DVR outside of the lesion mask and 1mm trimmed from CSF and GM
PET_outsidelesionmask<-PET[which(lesionmask1==0 & lesionmask2==0 & wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==1)];length(PET_outsidelesionmask)
PET_voxelwise_zscoredHC_outsidelesionmask<-PET_voxelwise_zscoredHC[which(lesionmask1==0 & lesionmask2==0 & wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==1)];length(PET_voxelwise_zscoredHC_outsidelesionmask)
summary(as.vector(PET_voxelwise_zscoredHC_outsidelesionmask))

# create the mean value for this subject
PET_voxelwise_zscoredHC_outsidelesionmask_1mmtrimmed_n26<-rbind(PET_voxelwise_zscoredHC_outsidelesionmask_1mmtrimmed_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask))

# calculate the mean for the voxels with WM disruption=0
PET_zscoredHC_outsidelesionmask_chacoZERO_n26<-rbind(PET_zscoredHC_outsidelesionmask_chacoZERO_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMPOS_outsidelesionmask==0 & ChaCoduetoRIMNEG_outsidelesionmask==0)]))
PET_zscoredHC_outsidelesionmask_chacoPRLZERO_n26<-rbind(PET_zscoredHC_outsidelesionmask_chacoPRLZERO_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMPOS_outsidelesionmask==0)]))
PET_zscoredHC_outsidelesionmask_chacoNONPRLZERO_n26<-rbind(PET_zscoredHC_outsidelesionmask_chacoNONPRLZERO_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMNEG_outsidelesionmask==0)]))

# calculate the mean for the voxels with WM disruption between 0 and 0.1
PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n26<-rbind(PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMPOS_outsidelesionmask> 0 & ChaCoduetoRIMPOS_outsidelesionmask<  0.1 | ChaCoduetoRIMNEG_outsidelesionmask> 0 & ChaCoduetoRIMNEG_outsidelesionmask< 0.1)]))
PET_zscoredHC_outsidelesionmask_chacoPRL_0_0.1_n26<-rbind(PET_zscoredHC_outsidelesionmask_chacoPRL_0_0.1_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMPOS_outsidelesionmask> 0 & ChaCoduetoRIMPOS_outsidelesionmask < 0.1)]))
PET_zscoredHC_outsidelesionmask_chacoNONPRL_0_0.1_n26<-rbind(PET_zscoredHC_outsidelesionmask_chacoNONPRL_0_0.1_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMNEG_outsidelesionmask> 0 & ChaCoduetoRIMNEG_outsidelesionmask<0.1)]))

# calculate the mean for the voxels with WM disruption between 0.9 and 1

PET_zscoredHC_outsidelesionmask_chaco09_1_n26<-rbind(PET_zscoredHC_outsidelesionmask_chaco09_1_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMPOS_outsidelesionmask>= 0.9 & ChaCoduetoRIMPOS_outsidelesionmask<  1 | ChaCoduetoRIMNEG_outsidelesionmask>= 0.9 & ChaCoduetoRIMNEG_outsidelesionmask< 1)]))
PET_zscoredHC_outsidelesionmask_chacoPRL09_1_n26<-rbind(PET_zscoredHC_outsidelesionmask_chacoPRL09_1_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMPOS_outsidelesionmask>= 0.9 & ChaCoduetoRIMPOS_outsidelesionmask <  1)]))
PET_zscoredHC_outsidelesionmask_chacoNONPRL09_1_n26<-rbind(PET_zscoredHC_outsidelesionmask_chacoNONPRL09_1_n26,mean(PET_voxelwise_zscoredHC_outsidelesionmask[which(ChaCoduetoRIMNEG_outsidelesionmask>= 0.9 & ChaCoduetoRIMNEG_outsidelesionmask< 1)]))

}

## Patients without PRL ####

#create empty matrices
{

  PET_zscoredHC_outsidelesionmask_chacoZERO_n18<-NULL
  PET_zscoredHC_outsidelesionmask_chaco0.9_1_n18<-NULL
  PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n18<-NULL
  
  
  PET_zscoredHC_outsidelesionmask_chaco_0_trimmed_n18<-NULL
  PET_zscoredHC_outsidelesionmask_chaco_0_trimmed_abovezscore0_n18<-NULL
  
  PET_zscoredHC_outsidelesionmask_chaco_0_0dot1_trimmed_n18<-NULL
  PET_zscoredHC_outsidelesionmask_chaco_0_0dot1_trimmed_abovezscore0_n18<-NULL
  
  PET_zscoredHC_outsidelesionmask_chaco0dot9_1_trimmed_n18<-NULL
  PET_zscoredHC_outsidelesionmask_chaco0dot9_1_trimmed_abovezscore0_n18<-NULL
}
k=1

for(i in ID_WM_Chaco_rimneg_n18){cat(k);k=k+1

# WM disruption for this patient
ChaCo<-WM_Chaco_rimneg_n18[[which(ID_WM_Chaco_rimneg_n18==i)]];dim(ChaCo)
# DVR metrics for this patient
PET_zscoredHC<-PET_n18_nonPRL_voxelwisezscoredHC[[which(ID_PET_n18_allsubjects_july2024==i)]];dim(PET)
# lesion mask for this patient
lesionmask<-lesionmask_rimneg_n18_v2[[which(ID_lesionmask_rimneg_n18_v2==i)]];dim(lesionmask)

# remove the lesion mask and the voxels within 1mm proximity of CSF and GM from the WM disruption mask
ChaCo_outsidelesionmask<-ChaCo[which(lesionmask==0 & wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==1)];length(ChaCo_outsidelesionmask)

# remove the lesion mask and the voxels within 1mm proximity of CSF and GM from the DVR image
PET_zscoredHC_outsidelesionmask<-PET_zscoredHC[which(lesionmask==0 & wm_mask_aseg_0dot5thr_bin_16HC_ero1mm==1)];length(PET_outsidelesionmask)

# DVR metrics where WM disruption is zero
PET_zscoredHC_outsidelesionmask_chacoZERO_n18<-rbind(PET_zscoredHC_outsidelesionmask_chacoZERO_n18,mean(PET_zscoredHC_outsidelesionmask[which(ChaCo_outsidelesionmask==0)]))
# DVR metrics where WM disruption is between 0-0.1
PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n18<-rbind(PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n18,mean(PET_zscoredHC_outsidelesionmask[which(ChaCo_outsidelesionmask>0 & ChaCo_outsidelesionmask < 0.1)])) 
# DVR metrics where WM disruption is between 0.9-1
PET_zscoredHC_outsidelesionmask_chaco0.9_1_n18<-rbind(PET_zscoredHC_outsidelesionmask_chaco0.9_1_n18,mean(PET_zscoredHC_outsidelesionmask[which(ChaCo_outsidelesionmask>=0.9 & ChaCo_outsidelesionmask < 1)]))

}


# subjects without and with PRL
nonPRLnumber<-18
PRLnumber<-26

# combine the datasets for subjects with and without PRL for the DVR metrics on the voxels with different level of disruption
PET_zscoredHC_1mmtrimmed_2024_10_02<-rbind(
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chacoZERO_n18,group=rep("T2_ChaCo0_n18",nonPRLnumber),paired=fornonPRL),    
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chacoZERO_n26,group=rep("T2_ChaCo0_n26",PRLnumber),paired=forPRL),    

  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n18,group=rep("T2_LowChaCo_n18",nonPRLnumber),paired=fornonPRL),    
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chaco_0_0.1_n26,group=rep("T2_LowChaCo",PRLnumber),paired=forPRL),    
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chacoPRL_0_0.1_n26,group=rep("PRL_LowChaCo",PRLnumber),paired=forPRL),
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chacoNONPRL_0_0.1_n26,group=rep("nonPRL_LowChaCo",PRLnumber),paired=forPRL),
  
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chaco0.9_1_n18,group=rep("T2_HighChaCo_n18",nonPRLnumber),paired=fornonPRL),    
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chaco09_1_n26,group=rep("T2_HighChaCo",PRLnumber),paired=forPRL),    
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chacoPRL09_1_n26,group=rep("PRL_HighChaCo",PRLnumber),paired=forPRL),
  data.frame(PETmetric=PET_zscoredHC_outsidelesionmask_chacoNONPRL09_1_n26,group=rep("nonPRL_HighChaCo",PRLnumber),paired=forPRL)

)

# reorder the groups
PET_zscoredHC_1mmtrimmed_2024_10_02$group <- factor(PET_zscoredHC_1mmtrimmed_2024_10_02$group,
                                                    levels = c("T2_ChaCo0_n18","T2_ChaCo0_n26", "T2_LowChaCo_n18","T2_LowChaCo","PRL_LowChaCo", "nonPRL_LowChaCo","T2_HighChaCo_n18","T2_HighChaCo","PRL_HighChaCo","nonPRL_HighChaCo"))


# compare subjects with vs without PRL where WM disruption is zero
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_ChaCo0_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_ChaCo0_n26"),],demo_ID_WM_Chaco_rimpos_n26_v2))


pval_chacozero_subjnoprl_subjwithprl_chacoduetoT2<-summary(aov(PETmetric~group+Age.at.Base+Gender+EDSS+newMStype,data=data1))[[1]][, "Pr(>F)"][1];pval_chacozero_subjnoprl_subjwithprl_chacoduetoT2

# compare subjects with vs without PRL where WM disruption due to all lesions is low
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_LowChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_LowChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))

pval_chacolow_subjnoprl_subjwithprl_chacoduetoT2<-summary(aov(PETmetric~group+Age.at.Base+Gender+EDSS+newMStype,data=data1))[[1]][, "Pr(>F)"][1];pval_chacolow_subjnoprl_subjwithprl_chacoduetoT2

# compare subjects with vs without PRL where WM disruption due to all lesions vs PRL, respectively is low
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_LowChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="PRL_LowChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))

pval_chacolow_subjnoprl_subjwithprl_chacoduetoprl<-summary(aov(PETmetric~group+Age.at.Base+Gender+EDSS+newMStype,data=data1))[[1]][, "Pr(>F)"][1];pval_chacolow_subjnoprl_subjwithprl_chacoduetoprl

# compare subjects with vs without PRL where WM disruption due to all lesions vs non-PRL, respectively is low
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_LowChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="nonPRL_LowChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))

pval_chacolow_subjnoprl_subjwithprl_chacoduetononprl<-summary(aov(PETmetric~group+Age.at.Base+Gender+EDSS+newMStype,data=data1))[[1]][, "Pr(>F)"][1];pval_chacolow_subjnoprl_subjwithprl_chacoduetononprl

# compare subjects with PRL where WM disruption due to PRL vs non-PRL, respectively is low
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="PRL_LowChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="nonPRL_LowChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))
pval_chacolow_subjwithprl_chacoduetoprl_vs_nonprl<-wilcox.test(data1[which(data1$group=="PRL_LowChaCo"),"PETmetric"],data1[which(data1$group=="nonPRL_LowChaCo"),"PETmetric"],paired=TRUE)$p.value

# compare subjects with vs without PRL where WM disruption due to all lesions is high
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_HighChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_HighChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))

pval_chacohigh_subjnoprl_subjwithprl_chacoduetoT2<-summary(aov(PETmetric~group+Age.at.Base+Gender+EDSS+newMStype,data=data1))[[1]][, "Pr(>F)"][1];pval_chacohigh_subjnoprl_subjwithprl_chacoduetoT2

# compare subjects with vs without PRL where WM disruption due to all lesions and PRL is high
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_HighChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="PRL_HighChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))

pval_chacohigh_subjnoprl_subjwithprl_chacoduetoPRL<-summary(aov(PETmetric~group+Age.at.Base+Gender+EDSS+newMStype,data=data1))[[1]][, "Pr(>F)"][1];pval_chacohigh_subjnoprl_subjwithprl_chacoduetoPRL

# compare subjects with vs without PRL where WM disruption due to all lesions and non-PRL is high
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_HighChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="nonPRL_HighChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))

pval_chacohigh_subjnoprl_subjwithprl_chacoduetononPRL<-summary(aov(PETmetric~group+Age.at.Base+Gender+EDSS+newMStype,data=data1))[[1]][, "Pr(>F)"][1];pval_chacohigh_subjnoprl_subjwithprl_chacoduetononPRL

# compare subjects with PRL where WM disruption due to PRL is high
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="PRL_HighChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="nonPRL_HighChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))
pval_chacohigh_subjwithprl_chacoduetoprl_vs_nonprl<-wilcox.test(data1[which(data1$group=="PRL_HighChaCo"),"PETmetric"],data1[which(data1$group=="nonPRL_HighChaCo"),"PETmetric"],paired=TRUE)$p.value

# start the figure

par(mar = c(5, 5, 5, 5) + 0.1)  # Adjust the margins as needed
scales::hue_pal()(4)
col3 = c("#FFD64C","#00616D",
         "#FFD64C","#00616D","#00AFBB", "#08EFFF",
         "#FFD64C","#00616D","#00AFBB", "#08EFFF")
levels1<-c("ChaCo = 0","ChaCo = 0", 
           "0<ChaCo<=0.1","0<ChaCo<=0.1","0<ChaCo<=0.1","0<ChaCo<=0.1",  
           "0.9<=ChaCo<=1","0.9<=ChaCo<=1","0.9<=ChaCo<=1","0.9<=ChaCo<=1")

# combine the p-values for the WM with low disruption and correct them using BH
p_lowdisruption<-c(pval_chacolow_subjnoprl_subjwithprl_chacoduetoT2,
                   pval_chacolow_subjnoprl_subjwithprl_chacoduetoprl,
                   pval_chacolow_subjnoprl_subjwithprl_chacoduetononprl,
                   pval_chacolow_subjwithprl_chacoduetoprl_vs_nonprl)
p_lowdisruption_adjusted<-p.adjust(p_lowdisruption,n=4,method="BH")

# combine the p-values for the highly disrupted WM and correct them using BH
p_highdisruption<-c(pval_chacohigh_subjnoprl_subjwithprl_chacoduetoT2,
                    pval_chacohigh_subjnoprl_subjwithprl_chacoduetoPRL,
                    pval_chacohigh_subjnoprl_subjwithprl_chacoduetononPRL,
                    pval_chacohigh_subjwithprl_chacoduetoprl_vs_nonprl);p_highdisruption
p_highdisruption_adjusted<-p.adjust(p_highdisruption,n=4,method="BH");p_highdisruption_adjusted

# fix the text size 
textsize1<-5

ggplot(PET_zscoredHC_1mmtrimmed_2024_10_02, aes(x=group, y=PETmetric,fill=group)) +   
  ylim(min(PET_zscoredHC_1mmtrimmed_2024_10_02$PETmetric),5)+
  # ylim(min(yalist3[,1]),max(yalist3[,1]))+
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white", show.legend = FALSE)+
  labs(title="",x="", y = "Subject Level Averaged z-scored DVR")+
  theme(axis.title.x = element_text(size = 20),  # Adjust size for x-axis
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.key.height = unit(2, "lines"),  # Adjust height as needed
        legend.key.width = unit(2, "lines"),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x=element_blank())+
  scale_fill_brewer(palette="RdBu")+ 
  scale_fill_manual(values = col3)+
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("T2_ChaCo0_n18", "T2_ChaCo0_n26")), annotations="No Disruption", y_position = 3.5, tip_length = 0.03, textsize = textsize1,fontface="bold")+
  geom_signif(comparisons = list(c("T2_LowChaCo_n18", "nonPRL_LowChaCo")), annotations="Low Level of Disruption", y_position = 3.5, tip_length = 0.03, textsize = textsize1,fontface="bold")+
  geom_signif(comparisons = list(c("T2_HighChaCo_n18", "nonPRL_HighChaCo")), annotations="High Level of Disruption", y_position = 3.5, tip_length = 0.03, textsize = textsize1,fontface="bold")+
  
  geom_signif(comparisons = list(c("T2_ChaCo0_n18", "T2_ChaCo0_n26")), annotations=paste("pBH=", round(pval_chacozero_subjnoprl_subjwithprl_chacoduetoT2,3)), y_position = 2.8, tip_length = 0.03, textsize = textsize1)+
  
  geom_signif(comparisons = list(c("T2_LowChaCo_n18", "T2_LowChaCo")), annotations=paste("pBH=", round(p_lowdisruption_adjusted[1],3)), y_position = 2.8, tip_length = 0.03, textsize = textsize1)+
  geom_signif(comparisons = list(c("T2_LowChaCo_n18", "PRL_LowChaCo")), annotations=paste("pBH=", round(p_lowdisruption_adjusted[2],3)), y_position = 3, tip_length = 0.03, textsize = textsize1)+
  geom_signif(comparisons = list(c("T2_LowChaCo_n18", "nonPRL_LowChaCo")), annotations=paste("pBH=", round(p_lowdisruption_adjusted[3],3)), y_position = 3.2, tip_length = 0.03, textsize = textsize1)+
  geom_signif(comparisons = list(c("PRL_LowChaCo", "nonPRL_LowChaCo")), annotations=paste("pBH=", round(p_lowdisruption_adjusted[4],3)), y_position = 2.8, tip_length = 0.03, textsize = textsize1)+
  
  geom_signif(comparisons = list(c("T2_HighChaCo_n18", "T2_HighChaCo")), annotations=paste("pBH=", round(p_highdisruption_adjusted[1],3)), y_position = 2.8, tip_length = 0.03, textsize = textsize1,col="black")+
  geom_signif(comparisons = list(c("T2_HighChaCo_n18", "PRL_HighChaCo")), annotations=paste("pBH=", round(p_highdisruption_adjusted[2],3)), y_position = 3, tip_length = 0.03, textsize = textsize1,col="black")+
  geom_signif(comparisons = list(c("T2_HighChaCo_n18", "nonPRL_HighChaCo")), annotations=paste("pBH=", round(p_highdisruption_adjusted[3],3)), y_position = 3.2, tip_length = 0.03, textsize = textsize1,col="black")+
  geom_signif(comparisons = list(c("PRL_HighChaCo", "nonPRL_HighChaCo")), annotations=paste("pBH=", round(p_highdisruption_adjusted[4],3)), y_position = 2.8, tip_length = 0.03, textsize = textsize1,col="black")
dev.off()


# FIGURE 4- DVR vs EDSS for various ChaCo level ####

# WM disruption=0 ####
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_ChaCo0_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_ChaCo0_n26"),],demo_ID_WM_Chaco_rimpos_n26_v2))

# Fit linear models for both groups
model_yes <- lm(PETmetric ~ EDSS + Gender + Age.at.Base+newMStype  , data = data1[which(data1$QSM.rim == "Yes"),])
model_no <- lm(PETmetric ~ EDSS + Gender + Age.at.Base+newMStype , data = data1[which(data1$QSM.rim == "No"),])
model_all <- lm(PETmetric ~ EDSS + Gender + Age.at.Base+newMStype , data = data1)

# Extract estimates and p-values for EDSS from both models
estimate_yes <- summary(model_yes)$coefficients["EDSS", "Estimate"]
p_value_yes <- summary(model_yes)$coefficients["EDSS", "Pr(>|t|)"];p_value_yes

estimate_no <- summary(model_no)$coefficients["EDSS", "Estimate"]
p_value_no <- summary(model_no)$coefficients["EDSS", "Pr(>|t|)"];p_value_no

estimate_all <- summary(model_all)$coefficients["EDSS", "Estimate"]
p_value_all <- summary(model_all)$coefficients["EDSS", "Pr(>|t|)"];p_value_all

adjustedP<-p.adjust(c(p_value_yes,p_value_no,p_value_all),method="BH",n=3);adjustedP

#start the figure for WM disruption=0
levelup<-1.5
data1$PseudoGroup <- ifelse(is.na(data1$QSM.rim), "Overall", as.character(data1$QSM.rim))

p1_edss<-ggplot(data1, aes(x = EDSS, y = PETmetric, color = QSM.rim)) +ylim(c(-1,4))+
  geom_point(size = 3) +  # Scatter points
  geom_smooth(method = "lm", se = TRUE, aes(fill = QSM.rim), alpha = 0.2) +  # Regression lines
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 1) +  # Overall regression line
  scale_fill_manual(values = c("No" = "#E69F00", "Yes" = "#56B4E9")) +  # Improved colors
  scale_color_manual(values = c("No" = "#E69F00", "Yes" = "#56B4E9"),
                     labels = c("Patients \n without PRL",
                                "Patients \n with PRL")) +
  labs(title = "",
       x = "EDSS", 
       y = "Average z-scored DVR in the NAWM \n on which there was no disruption") +
  guides(fill = FALSE) +  # Remove fill legend for smoother line
  theme(axis.title.x = element_text(size = 20),  # Adjust font size for x-axis title
        axis.title.y = element_text(size = 20),  # Adjust font size for y-axis title
        axis.text.x = element_text(size = 20),  # Adjust font size for x-axis text
        axis.text.y = element_text(size = 20),  # Adjust font size for y-axis text
        plot.title = element_text(size = 22, face = "bold"),  # Adjust plot title size and bold
        legend.key.height = unit(3, "lines"),  # Adjust legend key height
        legend.key.width = unit(3, "lines"),
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 15)) +  # Adjust legend text size
  annotate("text", x = 0, y = 3, 
           label = paste0("EDSS Estimate = ", round(estimate_yes, 3), 
                          ", pBH = ", format.pval(adjustedP[1], digits = 3)),
           color = "#56B4E9", size = 8, hjust = 0) +
  
  annotate("text", x = 0, y = 3.4 , 
           label = paste0("EDSS Estimate = ", round(estimate_no, 3), 
                          ", pBH = ", format.pval(adjustedP[2], digits = 3)),
           color = "#E69F00", size = 8, hjust = 0)+
  
  annotate("text", x = 0, y = 3.8 ,
           label = paste0("EDSS Estimate = ", round(estimate_all, 3),
                          ", pBH = ", format.pval(adjustedP[3], digits = 3)),
           color = "black", size = 8, hjust = 0)

# WM disruption is low ####
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_LowChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_LowChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))


# Fit linear models for both groups
model_yes <- lm(PETmetric ~ EDSS + Gender + Age.at.Base +newMStype  , data = data1[which(data1$QSM.rim == "Yes"),])
model_no <- lm(PETmetric ~ EDSS + Gender + Age.at.Base+newMStype, data = data1[which(data1$QSM.rim == "No"),])
model_all <- lm(PETmetric ~ EDSS + Gender + Age.at.Base+newMStype, data = data1)

# Extract estimates and p-values for EDSS from both models
estimate_yes <- summary(model_yes)$coefficients["EDSS", "Estimate"]
p_value_yes <- summary(model_yes)$coefficients["EDSS", "Pr(>|t|)"];p_value_yes

estimate_no <- summary(model_no)$coefficients["EDSS", "Estimate"]
p_value_no <- summary(model_no)$coefficients["EDSS", "Pr(>|t|)"];p_value_no

estimate_all <- summary(model_all)$coefficients["EDSS", "Estimate"]
p_value_all <- summary(model_all)$coefficients["EDSS", "Pr(>|t|)"];p_value_all

# adjust p-values
adjustedP<-p.adjust(c(p_value_yes,p_value_no,p_value_all),method="BH",n=3);adjustedP

# start the plot for WM with low disruption
levelup<-2
p2_edss<-ggplot(data1, aes(x = EDSS, y = PETmetric, color = QSM.rim)) +ylim(c(-1,4))+
  geom_point(size = 3) +  # Scatter points
  geom_smooth(method = "lm", se = TRUE, aes(fill = QSM.rim), alpha = 0.2) +  # Regression lines
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 1) +  # Overall regression line
  
  scale_fill_manual(values = c("No" = "#E69F00", "Yes" = "#56B4E9")) +  # Improved colors
  scale_color_manual(values = c("No" = "#E69F00", "Yes" = "#56B4E9"),
                     labels = c("Subjects \n without PRL", "Subjects \n with PRL")) +
  labs(title = "",
       x = "EDSS", 
       y = "Average z-scored DVR in the NAWM \n on which there was a low level of disruption") +
  guides(fill = FALSE) +  # Remove fill legend for smoother line
  theme(axis.title.x = element_text(size = 20),  # Adjust font size for x-axis title
        axis.title.y = element_text(size = 20),  # Adjust font size for y-axis title
        axis.text.x = element_text(size = 20),  # Adjust font size for x-axis text
        axis.text.y = element_text(size = 20),  # Adjust font size for y-axis text
        plot.title = element_text(size = 22, face = "bold"),  # Adjust plot title size and bold
        legend.key.height = unit(3, "lines"),  # Adjust legend key height
        legend.key.width = unit(3, "lines"),
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 15)) +  # Adjust legend text size
  # Add annotations for model estimates and p-values
  annotate("text", x = 0, y = 3, 
           label = paste0("EDSS Estimate = ", round(estimate_yes, 3), 
                          ", pBH = ", format.pval(adjustedP[1], digits = 3)),
           color = "#56B4E9", size = 8, hjust = 0) +
  
  annotate("text", x = 0, y = 3.4 , 
           label = paste0("EDSS Estimate = ", round(estimate_no, 3), 
                          ", pBH = ", format.pval(adjustedP[2], digits = 3)),
           color = "#E69F00", size = 8, hjust = 0)+
  annotate("text", x = 0, y = 3.8,
           label = paste0("EDSS Estimate = ", round(estimate_all, 3),
                          ", pBH = ", format.pval(adjustedP[3], digits = 3)),
           color = "black", size = 8, hjust = 0)


# High WM disruption ####
data1<-rbind(data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_HighChaCo_n18"),],demo_ID_WM_Chaco_rimneg_n18),
             data.frame(PET_zscoredHC_1mmtrimmed_2024_10_02[which(PET_zscoredHC_1mmtrimmed_2024_10_02$group=="T2_HighChaCo"),],demo_ID_WM_Chaco_rimpos_n26_v2))

# Fit linear models for both groups
model_yes <- lm(PETmetric ~ EDSS + Gender + Age.at.Base +newMStype , data = data1[which(data1$QSM.rim == "Yes"),])
model_no <- lm(PETmetric ~ EDSS + Gender + Age.at.Base+newMStype, data = data1[which(data1$QSM.rim == "No"),])
model_all <- lm(PETmetric ~ EDSS + Gender + Age.at.Base+newMStype, data = data1)

# Extract estimates and p-values for EDSS from both models
estimate_yes <- summary(model_yes)$coefficients["EDSS", "Estimate"]
p_value_yes <- summary(model_yes)$coefficients["EDSS", "Pr(>|t|)"];p_value_yes

estimate_no <- summary(model_no)$coefficients["EDSS", "Estimate"]
p_value_no <- summary(model_no)$coefficients["EDSS", "Pr(>|t|)"];p_value_no

estimate_all <- summary(model_all)$coefficients["EDSS", "Estimate"]
p_value_all <- summary(model_all)$coefficients["EDSS", "Pr(>|t|)"];p_value_all

# adjust the p values
adjustedP<-p.adjust(c(p_value_yes,p_value_no,p_value_all),method="BH",n=3);adjustedP

# start the figure for the WM with high disruption
levelup<-1.5
p3_edss<-ggplot(data1, aes(x = EDSS, y = PETmetric, color = QSM.rim)) +ylim(c(-1,4))+
  geom_point(size = 3) +  # Scatter points
  geom_smooth(method = "lm", se = TRUE, aes(fill = QSM.rim), alpha = 0.2) +  # Regression lines
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 1) +  # Overall regression line
  
  scale_fill_manual(values = c("No" = "#E69F00", "Yes" = "#56B4E9")) +  # Improved colors
  scale_color_manual(values = c("No" = "#E69F00", "Yes" = "#56B4E9"),
                     labels = c("Subjects \n without PRL", "Subjects \n with PRL")) +
  labs(title = "",
       x = "EDSS", 
       y = "Average z-scored DVR in the NAWM \n on which there was a high level of disruption") +
  guides(fill = FALSE) +  # Remove fill legend for smoother line
  # theme_minimal() +  # Cleaner theme
  theme(axis.title.x = element_text(size = 20),  # Adjust font size for x-axis title
        axis.title.y = element_text(size = 20),  # Adjust font size for y-axis title
        axis.text.x = element_text(size = 20),  # Adjust font size for x-axis text
        axis.text.y = element_text(size = 20),  # Adjust font size for y-axis text
        plot.title = element_text(size = 22, face = "bold"),  # Adjust plot title size and bold
        legend.key.height = unit(3, "lines"),  # Adjust legend key height
        legend.key.width = unit(3, "lines"),
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 15)) +  # Adjust legend text size
  # Add annotations for model estimates and p-values
  annotate("text", x = 0, y = 3, 
           label = paste0("EDSS Estimate = ", round(estimate_yes, 3), 
                          ", pBH = ", format.pval(adjustedP[1], digits = 3)),
           color = "#56B4E9", size = 8, hjust = 0) +
  
  annotate("text", x = 0, y = 3.4 , 
           label = paste0("EDSS Estimate = ", round(estimate_no, 3), 
                          ", pBH = ", format.pval(adjustedP[2], digits = 3)),
           color = "#E69F00", size = 8, hjust = 0)+
  
  annotate("text", x = 0, y = 3.8, 
           label = paste0("EDSS Estimate = ", round(estimate_all, 3), 
                          ", pBH = ", format.pval(adjustedP[3], digits = 3)),
           color = "black", size = 8, hjust = 0)


# combine the 3 figures (WM disruption is zero, low, and high)
grid.arrange(p1_edss,p2_edss,p3_edss,ncol=3,nrow=1)
