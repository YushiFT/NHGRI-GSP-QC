################################################################
# Topic: Compare GSP QC 
#  Date: 07/30/2019
#    By: Yushi F.T.
################################################################

library(ggplot2)
library(dplyr)

# Step 1: Combine All Chromosomes ##############################
# Change the files as postqc or postqc_revisegt
# postqc: variants list after QC, with raw GT
# postqc_revisegt: variants list after QC, with revised GT
################################################################

metadata22 <- read.table('../temp_new/qced_22_sample_qc_info_postqc_revisegt.txt',header=T)
metadata21 <- read.table('../temp_new/qced_21_sample_qc_info_postqc_revisegt.txt',header=T)
metadata20 <- read.table('../temp_new/qced_20_sample_qc_info_postqc_revisegt.txt',header=T)
metadata19 <- read.table('../temp_new/qced_19_sample_qc_info_postqc_revisegt.txt',header=T)
metadata18 <- read.table('../temp_new/qced_18_sample_qc_info_postqc_revisegt.txt',header=T)
metadata17 <- read.table('../temp_new/qced_17_sample_qc_info_postqc_revisegt.txt',header=T)
metadata16 <- read.table('../temp_new/qced_16_sample_qc_info_postqc_revisegt.txt',header=T)
metadata15 <- read.table('../temp_new/qced_15_sample_qc_info_postqc_revisegt.txt',header=T)
metadata14 <- read.table('../temp_new/qced_14_sample_qc_info_postqc_revisegt.txt',header=T)
metadata13 <- read.table('../temp_new/qced_13_sample_qc_info_postqc_revisegt.txt',header=T)
metadata12 <- read.table('../temp_new/qced_12_sample_qc_info_postqc_revisegt.txt',header=T)
metadata11 <- read.table('../temp_new/qced_11_sample_qc_info_postqc_revisegt.txt',header=T)
metadata10 <- read.table('../temp_new/qced_10_sample_qc_info_postqc_revisegt.txt',header=T)
metadata9 <- read.table('../temp_new/qced_9_sample_qc_info_postqc_revisegt.txt',header=T)
metadata8 <- read.table('../temp_new/qced_8_sample_qc_info_postqc_revisegt.txt',header=T)
metadata7 <- read.table('../temp_new/qced_7_sample_qc_info_postqc_revisegt.txt',header=T)
metadata6 <- read.table('../temp_new/qced_6_sample_qc_info_postqc_revisegt.txt',header=T)
metadata5 <- read.table('../temp_new/qced_5_sample_qc_info_postqc_revisegt.txt',header=T)
metadata4 <- read.table('../temp_new/qced_4_sample_qc_info_postqc_revisegt.txt',header=T)
metadata3 <- read.table('../temp_new/qced_3_sample_qc_info_postqc_revisegt.txt',header=T)
metadata2 <- read.table('../temp_new/qced_2_sample_qc_info_postqc_revisegt.txt',header=T)
metadata1 <- read.table('../temp_new/qced_1_sample_qc_info_postqc_revisegt.txt',header=T)


# Before QC-ed variants
metadata22 <- read.table('../work/qced_22_sample_qc_info_preqc.txt',header=T)
metadata21 <- read.table('../work/qced_21_sample_qc_info_preqc.txt',header=T)
metadata20 <- read.table('../work/qced_20_sample_qc_info_preqc.txt',header=T)
metadata19 <- read.table('../work/qced_19_sample_qc_info_preqc.txt',header=T)
metadata18 <- read.table('../work/qced_18_sample_qc_info_preqc.txt',header=T)
metadata17 <- read.table('../work/qced_17_sample_qc_info_preqc.txt',header=T)
metadata16 <- read.table('../work/qced_16_sample_qc_info_preqc.txt',header=T)
metadata15 <- read.table('../work/qced_15_sample_qc_info_preqc.txt',header=T)
metadata14 <- read.table('../work/qced_14_sample_qc_info_preqc.txt',header=T)
metadata13 <- read.table('../work/qced_13_sample_qc_info_preqc.txt',header=T)
metadata12 <- read.table('../work/qced_12_sample_qc_info_preqc.txt',header=T)
metadata11 <- read.table('../work/qced_11_sample_qc_info_preqc.txt',header=T)
metadata10 <- read.table('../work/qced_10_sample_qc_info_preqc.txt',header=T)
metadata9 <- read.table('../work/qced_9_sample_qc_info_preqc.txt',header=T)
metadata8 <- read.table('../work/qced_8_sample_qc_info_preqc.txt',header=T)
metadata7 <- read.table('../work/qced_7_sample_qc_info_preqc.txt',header=T)
metadata6 <- read.table('../work/qced_6_sample_qc_info_preqc.txt',header=T)
metadata5 <- read.table('../work/qced_5_sample_qc_info_preqc.txt',header=T)
metadata4 <- read.table('../work/qced_4_sample_qc_info_preqc.txt',header=T)
metadata3 <- read.table('../work/qced_3_sample_qc_info_preqc.txt',header=T)
metadata2 <- read.table('../work/qced_2_sample_qc_info_preqc.txt',header=T)
metadata1 <- read.table('../work/qced_1_sample_qc_info_preqc.txt',header=T)

# Merge chromosome 22, 21
collist <- c(1,11,12,13,14,15,16,20,21,22,23)
temp <- merge(metadata22[,collist], metadata21[,collist], by = 's' )
metadata <- as.data.frame(temp$s)
colnames(metadata) <- 's'

metadata$n_called <- temp$sample_qc.n_called.x + temp$sample_qc.n_called.y
metadata$n_not_called <- temp$sample_qc.n_not_called.x + temp$sample_qc.n_not_called.y
metadata$n_filtered <- temp$sample_qc.n_filtered.x + temp$sample_qc.n_filtered.y
metadata$n_hom_ref <- temp$sample_qc.n_hom_ref.x + temp$sample_qc.n_hom_ref.y
metadata$n_hom_var <- temp$sample_qc.n_hom_var.x + temp$sample_qc.n_hom_var.y
metadata$n_het <- temp$sample_qc.n_het.x + temp$sample_qc.n_het.y
metadata$n_insertion <- temp$sample_qc.n_insertion.x + temp$sample_qc.n_insertion.y
metadata$n_deletion <- temp$sample_qc.n_deletion.x + temp$sample_qc.n_deletion.y
metadata$n_transition <- temp$sample_qc.n_transition.x + temp$sample_qc.n_transition.y
metadata$n_transversion <- temp$sample_qc.n_transversion.x + temp$sample_qc.n_transversion.y

# Merge chromosome 20
temp <- merge(metadata, metadata20[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 19
temp <- merge(metadata, metadata19[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 18
temp <- merge(metadata, metadata18[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 17
temp <- merge(metadata, metadata17[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 16
temp <- merge(metadata, metadata16[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion

# Merge chromosome 15
temp <- merge(metadata, metadata15[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 14
temp <- merge(metadata, metadata14[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 13
temp <- merge(metadata, metadata13[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 12
temp <- merge(metadata, metadata12[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 11
temp <- merge(metadata, metadata11[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 10
temp <- merge(metadata, metadata10[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 9
temp <- merge(metadata, metadata9[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 8
temp <- merge(metadata, metadata8[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 7
temp <- merge(metadata, metadata7[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 6
temp <- merge(metadata, metadata6[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 5
temp <- merge(metadata, metadata5[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion



# Merge chromosome 4
temp <- merge(metadata, metadata4[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 3
temp <- merge(metadata, metadata3[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 2
temp <- merge(metadata, metadata2[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion


# Merge chromosome 1
temp <- merge(metadata, metadata1[,collist], by = 's' )
metadata$n_called <- temp$n_called + temp$sample_qc.n_called
metadata$n_not_called <- temp$n_not_called + temp$sample_qc.n_not_called
metadata$n_filtered <- temp$n_filtered + temp$sample_qc.n_filtered
metadata$n_hom_ref <- temp$n_hom_ref + temp$sample_qc.n_hom_ref
metadata$n_hom_var <- temp$n_hom_var + temp$sample_qc.n_hom_var
metadata$n_het <- temp$n_het + temp$sample_qc.n_het
metadata$n_insertion <- temp$n_insertion + temp$sample_qc.n_insertion
metadata$n_deletion <- temp$n_deletion + temp$sample_qc.n_deletion
metadata$n_transition <- temp$n_transition + temp$sample_qc.n_transition
metadata$n_transversion <- temp$n_transversion + temp$sample_qc.n_transversion



# Step 2: Calculate ratios ################################################
metadata$sample_qc.r_ti_tv <- metadata$n_transition / metadata$n_transversion
metadata$sample_qc.r_het_hom_var <- metadata$n_het / metadata$n_hom_var
metadata$sample_qc.r_insertion_deletion <- metadata$n_insertion / metadata$n_deletion

# Step 3: Annotate center and population #############################################
# Annotate center and population information
pop_meta <- read.table('../data/ccdgf2_predicted_ethnicity_PC1-15.txt',header=T)
center_data <- read.csv('../data/sample_center_platform.csv', header=T)
metadata <- merge(metadata, pop_meta[,1:2], by.x = "s", by.y = "Sample")
metadata <- merge(metadata, center_data, by.x = "s", by.y = "sample_id",
                  all.x=T, all.y=F)

# Re-order by population with center as subgroup
rownames(metadata) <- metadata$s
meta_AFR_NYGC <- subset(metadata,metadata$Pop=='AFR' & metadata$center_platform=='NYGC' )
meta_AFR_WashU <- subset(metadata,metadata$Pop=='AFR' & metadata$center_platform=='WASHU')
meta_AFR_Baylor <- subset(metadata,metadata$Pop=='AFR' & metadata$center_platform=='BAYLOR' )
meta_AFR_Broad <- subset(metadata,metadata$Pop=='AFR' & metadata$center_platform=='BROAD' )
meta_AFR_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='AFR' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_AFR_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='AFR' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_AFR <- rbind(meta_AFR_Baylor, meta_AFR_Baylor_HiSeqX)
meta_AFR <- rbind(meta_AFR, meta_AFR_Baylor_NovaSeq)
meta_AFR <- rbind(meta_AFR, meta_AFR_Broad)
meta_AFR <- rbind(meta_AFR, meta_AFR_NYGC)
meta_AFR <- rbind(meta_AFR, meta_AFR_WashU)


meta_AMR_NYGC <- subset(metadata,metadata$Pop=='AMR' & metadata$center_platform=='NYGC' )
meta_AMR_WashU <- subset(metadata,metadata$Pop=='AMR' & metadata$center_platform=='WASHU')
meta_AMR_Baylor <- subset(metadata,metadata$Pop=='AMR' & metadata$center_platform=='BAYLOR' )
meta_AMR_Broad <- subset(metadata,metadata$Pop=='AMR' & metadata$center_platform=='BROAD' )
meta_AMR_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='AMR' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_AMR_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='AMR' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_AMR <- rbind(meta_AMR_Baylor, meta_AMR_Baylor_HiSeqX)
meta_AMR <- rbind(meta_AMR, meta_AMR_Baylor_NovaSeq)
meta_AMR <- rbind(meta_AMR, meta_AMR_Broad)
meta_AMR <- rbind(meta_AMR, meta_AMR_NYGC)
meta_AMR <- rbind(meta_AMR, meta_AMR_WashU)


meta_EAS_NYGC <- subset(metadata,metadata$Pop=='EAS' & metadata$center_platform=='NYGC' )
meta_EAS_WashU <- subset(metadata,metadata$Pop=='EAS' & metadata$center_platform=='WASHU')
meta_EAS_Baylor <- subset(metadata,metadata$Pop=='EAS' & metadata$center_platform=='BAYLOR' )
meta_EAS_Broad <- subset(metadata,metadata$Pop=='EAS' & metadata$center_platform=='BROAD' )
meta_EAS_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='EAS' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_EAS_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='EAS' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_EAS <- rbind(meta_EAS_Baylor, meta_EAS_Baylor_HiSeqX)
meta_EAS <- rbind(meta_EAS, meta_EAS_Baylor_NovaSeq)
meta_EAS <- rbind(meta_EAS, meta_EAS_Broad)
meta_EAS <- rbind(meta_EAS, meta_EAS_NYGC)
meta_EAS <- rbind(meta_EAS, meta_EAS_WashU)


meta_EUR_NYGC <- subset(metadata,metadata$Pop=='EUR' & metadata$center_platform=='NYGC' )
meta_EUR_WashU <- subset(metadata,metadata$Pop=='EUR' & metadata$center_platform=='WASHU')
meta_EUR_Baylor <- subset(metadata,metadata$Pop=='EUR' & metadata$center_platform=='BAYLOR' )
meta_EUR_Broad <- subset(metadata,metadata$Pop=='EUR' & metadata$center_platform=='BROAD' )
meta_EUR_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='EUR' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_EUR_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='EUR' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_EUR <- rbind(meta_EUR_Baylor, meta_EUR_Baylor_HiSeqX)
meta_EUR <- rbind(meta_EUR, meta_EUR_Baylor_NovaSeq)
meta_EUR <- rbind(meta_EUR, meta_EUR_Broad)
meta_EUR <- rbind(meta_EUR, meta_EUR_NYGC)
meta_EUR <- rbind(meta_EUR, meta_EUR_WashU)


meta_FIN_NYGC <- subset(metadata,metadata$Pop=='FIN' & metadata$center_platform=='NYGC' )
meta_FIN_WashU <- subset(metadata,metadata$Pop=='FIN' & metadata$center_platform=='WASHU')
meta_FIN_Baylor <- subset(metadata,metadata$Pop=='FIN' & metadata$center_platform=='BAYLOR' )
meta_FIN_Broad <- subset(metadata,metadata$Pop=='FIN' & metadata$center_platform=='BROAD' )
meta_FIN_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='FIN' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_FIN_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='FIN' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_FIN <- rbind(meta_FIN_Baylor, meta_FIN_Baylor_HiSeqX)
meta_FIN <- rbind(meta_FIN, meta_FIN_Baylor_NovaSeq)
meta_FIN <- rbind(meta_FIN, meta_FIN_Broad)
meta_FIN <- rbind(meta_FIN, meta_FIN_NYGC)
meta_FIN <- rbind(meta_FIN, meta_FIN_WashU)

meta_PUR_NYGC <- subset(metadata,metadata$Pop=='PUR' & metadata$center_platform=='NYGC' )
meta_PUR_WashU <- subset(metadata,metadata$Pop=='PUR' & metadata$center_platform=='WASHU')
meta_PUR_Baylor <- subset(metadata,metadata$Pop=='PUR' & metadata$center_platform=='BAYLOR' )
meta_PUR_Broad <- subset(metadata,metadata$Pop=='PUR' & metadata$center_platform=='BROAD' )
meta_PUR_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='PUR' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_PUR_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='PUR' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_PUR <- rbind(meta_PUR_Baylor, meta_PUR_Baylor_HiSeqX)
meta_PUR <- rbind(meta_PUR, meta_PUR_Baylor_NovaSeq)
meta_PUR <- rbind(meta_PUR, meta_PUR_Broad)
meta_PUR <- rbind(meta_PUR, meta_PUR_NYGC)
meta_PUR <- rbind(meta_PUR, meta_PUR_WashU)


meta_SAS_NYGC <- subset(metadata,metadata$Pop=='SAS' & metadata$center_platform=='NYGC' )
meta_SAS_WashU <- subset(metadata,metadata$Pop=='SAS' & metadata$center_platform=='WASHU')
meta_SAS_Baylor <- subset(metadata,metadata$Pop=='SAS' & metadata$center_platform=='BAYLOR' )
meta_SAS_Broad <- subset(metadata,metadata$Pop=='SAS' & metadata$center_platform=='BROAD' )
meta_SAS_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='SAS' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_SAS_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='SAS' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_SAS <- rbind(meta_SAS_Baylor, meta_SAS_Baylor_HiSeqX)
meta_SAS <- rbind(meta_SAS, meta_SAS_Baylor_NovaSeq)
meta_SAS <- rbind(meta_SAS, meta_SAS_Broad)
meta_SAS <- rbind(meta_SAS, meta_SAS_NYGC)
meta_SAS <- rbind(meta_SAS, meta_SAS_WashU)


meta_Other_NYGC <- subset(metadata,metadata$Pop=='Other' & metadata$center_platform=='NYGC' )
meta_Other_WashU <- subset(metadata,metadata$Pop=='Other' & metadata$center_platform=='WASHU')
meta_Other_Baylor <- subset(metadata,metadata$Pop=='Other' & metadata$center_platform=='BAYLOR' )
meta_Other_Broad <- subset(metadata,metadata$Pop=='Other' & metadata$center_platform=='BROAD' )
meta_Other_Baylor_HiSeqX <- subset(metadata,metadata$Pop=='Other' & metadata$center_platform=='BAYLOR_HiSeqX' )
meta_Other_Baylor_NovaSeq <- subset(metadata,metadata$Pop=='Other' & metadata$center_platform=='BAYLOR_NovaSeq' )
meta_Other <- rbind(meta_Other_Baylor, meta_Other_Baylor_HiSeqX)
meta_Other <- rbind(meta_Other, meta_Other_Baylor_NovaSeq)
meta_Other <- rbind(meta_Other, meta_Other_Broad)
meta_Other <- rbind(meta_Other, meta_Other_NYGC)
meta_Other <- rbind(meta_Other, meta_Other_WashU)


meta_new <- rbind(meta_AFR, meta_AMR)
meta_new <- rbind(meta_new, meta_EAS)
meta_new <- rbind(meta_new, meta_EUR)
meta_new <- rbind(meta_new, meta_FIN)
meta_new <- rbind(meta_new, meta_PUR)
meta_new <- rbind(meta_new, meta_SAS)
meta_new <- rbind(meta_new, meta_Other)



# Re-order by center with population as subgroup

meta_Baylor <- rbind(meta_AFR_Baylor, meta_AMR_Baylor)
meta_Baylor <- rbind(meta_Baylor, meta_EAS_Baylor)
meta_Baylor <- rbind(meta_Baylor, meta_EUR_Baylor)
meta_Baylor <- rbind(meta_Baylor, meta_FIN_Baylor)
meta_Baylor <- rbind(meta_Baylor, meta_PUR_Baylor)
meta_Baylor <- rbind(meta_Baylor, meta_SAS_Baylor)
meta_Baylor <- rbind(meta_Baylor, meta_Other_Baylor)

meta_Broad <- rbind(meta_AFR_Broad, meta_AMR_Broad)
meta_Broad <- rbind(meta_Broad, meta_EAS_Broad)
meta_Broad <- rbind(meta_Broad, meta_EUR_Broad)
meta_Broad <- rbind(meta_Broad, meta_FIN_Broad)
meta_Broad <- rbind(meta_Broad, meta_PUR_Broad)
meta_Broad <- rbind(meta_Broad, meta_SAS_Broad)
meta_Broad <- rbind(meta_Broad, meta_Other_Broad)

meta_NYGC <- rbind(meta_AFR_NYGC, meta_AMR_NYGC)
meta_NYGC <- rbind(meta_NYGC, meta_EAS_NYGC)
meta_NYGC <- rbind(meta_NYGC, meta_EUR_NYGC)
meta_NYGC <- rbind(meta_NYGC, meta_FIN_NYGC)
meta_NYGC <- rbind(meta_NYGC, meta_PUR_NYGC)
meta_NYGC <- rbind(meta_NYGC, meta_SAS_NYGC)
meta_NYGC <- rbind(meta_NYGC, meta_Other_NYGC)

meta_Baylor_HiSeqX <- rbind(meta_AFR_Baylor_HiSeqX, meta_AMR_Baylor_HiSeqX)
meta_Baylor_HiSeqX <- rbind(meta_Baylor_HiSeqX, meta_EAS_Baylor_HiSeqX)
meta_Baylor_HiSeqX <- rbind(meta_Baylor_HiSeqX, meta_EUR_Baylor_HiSeqX)
meta_Baylor_HiSeqX <- rbind(meta_Baylor_HiSeqX, meta_FIN_Baylor_HiSeqX)
meta_Baylor_HiSeqX <- rbind(meta_Baylor_HiSeqX, meta_PUR_Baylor_HiSeqX)
meta_Baylor_HiSeqX <- rbind(meta_Baylor_HiSeqX, meta_SAS_Baylor_HiSeqX)
meta_Baylor_HiSeqX <- rbind(meta_Baylor_HiSeqX, meta_Other_Baylor_HiSeqX)

meta_Baylor_NovaSeq <- rbind(meta_AFR_Baylor_NovaSeq, meta_AMR_Baylor_NovaSeq)
meta_Baylor_NovaSeq <- rbind(meta_Baylor_NovaSeq, meta_EAS_Baylor_NovaSeq)
meta_Baylor_NovaSeq <- rbind(meta_Baylor_NovaSeq, meta_EUR_Baylor_NovaSeq)
meta_Baylor_NovaSeq <- rbind(meta_Baylor_NovaSeq, meta_FIN_Baylor_NovaSeq)
meta_Baylor_NovaSeq <- rbind(meta_Baylor_NovaSeq, meta_PUR_Baylor_NovaSeq)
meta_Baylor_NovaSeq <- rbind(meta_Baylor_NovaSeq, meta_SAS_Baylor_NovaSeq)
meta_Baylor_NovaSeq <- rbind(meta_Baylor_NovaSeq, meta_Other_Baylor_NovaSeq)

meta_WashU <- rbind(meta_AFR_WashU, meta_AMR_WashU)
meta_WashU <- rbind(meta_WashU, meta_EAS_WashU)
meta_WashU <- rbind(meta_WashU, meta_EUR_WashU)
meta_WashU <- rbind(meta_WashU, meta_FIN_WashU)
meta_WashU <- rbind(meta_WashU, meta_PUR_WashU)
meta_WashU <- rbind(meta_WashU, meta_SAS_WashU)
meta_WashU <- rbind(meta_WashU, meta_Other_WashU)


meta_new_center <- rbind(meta_Baylor, meta_Baylor_HiSeqX)
meta_new_center <- rbind(meta_new_center, meta_Baylor_NovaSeq)
meta_new_center <- rbind(meta_new_center, meta_Broad)
meta_new_center <- rbind(meta_new_center, meta_NYGC)
meta_new_center <- rbind(meta_new_center, meta_WashU)


colnames(meta_new)
meta_new$id = c(1:nrow(meta_new))

colnames(meta_new_center)
meta_new_center$id = c(1:nrow(meta_new_center))


c(-1,1)*4*mad(meta_AFR$sample_qc.r_ti_tv) + median(meta_AFR$sample_qc.r_ti_tv)
c(-1,1)*4*mad(meta_AMR$sample_qc.r_ti_tv) + median(meta_AMR$sample_qc.r_ti_tv)
c(-1,1)*4*mad(meta_EAS$sample_qc.r_ti_tv) + median(meta_EAS$sample_qc.r_ti_tv)
c(-1,1)*4*mad(meta_EUR$sample_qc.r_ti_tv) + median(meta_EUR$sample_qc.r_ti_tv)
c(-1,1)*4*mad(meta_FIN$sample_qc.r_ti_tv) + median(meta_FIN$sample_qc.r_ti_tv)
c(-1,1)*4*mad(meta_PUR$sample_qc.r_ti_tv) + median(meta_PUR$sample_qc.r_ti_tv)
c(-1,1)*4*mad(meta_SAS$sample_qc.r_ti_tv) + median(meta_SAS$sample_qc.r_ti_tv)
c(-1,1)*4*mad(meta_Other$sample_qc.r_ti_tv) + median(meta_Other$sample_qc.r_ti_tv)

nrow(subset(meta_AFR, (meta_AFR$sample_qc.r_ti_tv < median(meta_AFR$sample_qc.r_ti_tv) - 4*mad(meta_AFR$sample_qc.r_ti_tv)) | (meta_AFR$sample_qc.r_ti_tv > median(meta_AFR$sample_qc.r_ti_tv) + 4*mad(meta_AFR$sample_qc.r_ti_tv)) ))
nrow(subset(meta_AMR, (meta_AMR$sample_qc.r_ti_tv < median(meta_AMR$sample_qc.r_ti_tv) - 4*mad(meta_AMR$sample_qc.r_ti_tv)) | (meta_AMR$sample_qc.r_ti_tv > median(meta_AMR$sample_qc.r_ti_tv) + 4*mad(meta_AMR$sample_qc.r_ti_tv)) ))
nrow(subset(meta_EAS, (meta_EAS$sample_qc.r_ti_tv < median(meta_EAS$sample_qc.r_ti_tv) - 4*mad(meta_EAS$sample_qc.r_ti_tv)) | (meta_EAS$sample_qc.r_ti_tv > median(meta_EAS$sample_qc.r_ti_tv) + 4*mad(meta_EAS$sample_qc.r_ti_tv)) ))
nrow(subset(meta_EUR, (meta_EUR$sample_qc.r_ti_tv < median(meta_EUR$sample_qc.r_ti_tv) - 4*mad(meta_EUR$sample_qc.r_ti_tv)) | (meta_EUR$sample_qc.r_ti_tv > median(meta_EUR$sample_qc.r_ti_tv) + 4*mad(meta_EUR$sample_qc.r_ti_tv)) ))
nrow(subset(meta_FIN, (meta_FIN$sample_qc.r_ti_tv < median(meta_FIN$sample_qc.r_ti_tv) - 4*mad(meta_FIN$sample_qc.r_ti_tv)) | (meta_FIN$sample_qc.r_ti_tv > median(meta_FIN$sample_qc.r_ti_tv) + 4*mad(meta_FIN$sample_qc.r_ti_tv)) ))
nrow(subset(meta_PUR, (meta_PUR$sample_qc.r_ti_tv < median(meta_PUR$sample_qc.r_ti_tv) - 4*mad(meta_PUR$sample_qc.r_ti_tv)) | (meta_PUR$sample_qc.r_ti_tv > median(meta_PUR$sample_qc.r_ti_tv) + 4*mad(meta_PUR$sample_qc.r_ti_tv)) ))
nrow(subset(meta_SAS, (meta_SAS$sample_qc.r_ti_tv < median(meta_SAS$sample_qc.r_ti_tv) - 4*mad(meta_SAS$sample_qc.r_ti_tv)) | (meta_SAS$sample_qc.r_ti_tv > median(meta_SAS$sample_qc.r_ti_tv) + 4*mad(meta_SAS$sample_qc.r_ti_tv)) ))
nrow(subset(meta_Other, (meta_Other$sample_qc.r_ti_tv < median(meta_Other$sample_qc.r_ti_tv) - 4*mad(meta_Other$sample_qc.r_ti_tv)) | (meta_Other$sample_qc.r_ti_tv > median(meta_Other$sample_qc.r_ti_tv) + 4*mad(meta_Other$sample_qc.r_ti_tv)) ))


c(-1,1)*4*mad(meta_AFR$sample_qc.r_het_hom_var) + median(meta_AFR$sample_qc.r_het_hom_var)
c(-1,1)*4*mad(meta_AMR$sample_qc.r_het_hom_var) + median(meta_AMR$sample_qc.r_het_hom_var)
c(-1,1)*4*mad(meta_EAS$sample_qc.r_het_hom_var) + median(meta_EAS$sample_qc.r_het_hom_var)
c(-1,1)*4*mad(meta_EUR$sample_qc.r_het_hom_var) + median(meta_EUR$sample_qc.r_het_hom_var)
c(-1,1)*4*mad(meta_FIN$sample_qc.r_het_hom_var) + median(meta_FIN$sample_qc.r_het_hom_var)
c(-1,1)*4*mad(meta_PUR$sample_qc.r_het_hom_var) + median(meta_PUR$sample_qc.r_het_hom_var)
c(-1,1)*4*mad(meta_SAS$sample_qc.r_het_hom_var) + median(meta_SAS$sample_qc.r_het_hom_var)
c(-1,1)*4*mad(meta_Other$sample_qc.r_het_hom_var) + median(meta_Other$sample_qc.r_het_hom_var)

nrow(subset(meta_AFR, (meta_AFR$sample_qc.r_het_hom_var < median(meta_AFR$sample_qc.r_het_hom_var) - 4*mad(meta_AFR$sample_qc.r_het_hom_var)) | (meta_AFR$sample_qc.r_het_hom_var > median(meta_AFR$sample_qc.r_het_hom_var) + 4*mad(meta_AFR$sample_qc.r_het_hom_var)) ))
nrow(subset(meta_AMR, (meta_AMR$sample_qc.r_het_hom_var < median(meta_AMR$sample_qc.r_het_hom_var) - 4*mad(meta_AMR$sample_qc.r_het_hom_var)) | (meta_AMR$sample_qc.r_het_hom_var > median(meta_AMR$sample_qc.r_het_hom_var) + 4*mad(meta_AMR$sample_qc.r_het_hom_var)) ))
nrow(subset(meta_EAS, (meta_EAS$sample_qc.r_het_hom_var < median(meta_EAS$sample_qc.r_het_hom_var) - 4*mad(meta_EAS$sample_qc.r_het_hom_var)) | (meta_EAS$sample_qc.r_het_hom_var > median(meta_EAS$sample_qc.r_het_hom_var) + 4*mad(meta_EAS$sample_qc.r_het_hom_var)) ))
nrow(subset(meta_EUR, (meta_EUR$sample_qc.r_het_hom_var < median(meta_EUR$sample_qc.r_het_hom_var) - 4*mad(meta_EUR$sample_qc.r_het_hom_var)) | (meta_EUR$sample_qc.r_het_hom_var > median(meta_EUR$sample_qc.r_het_hom_var) + 4*mad(meta_EUR$sample_qc.r_het_hom_var)) ))
nrow(subset(meta_FIN, (meta_FIN$sample_qc.r_het_hom_var < median(meta_FIN$sample_qc.r_het_hom_var) - 4*mad(meta_FIN$sample_qc.r_het_hom_var)) | (meta_FIN$sample_qc.r_het_hom_var > median(meta_FIN$sample_qc.r_het_hom_var) + 4*mad(meta_FIN$sample_qc.r_het_hom_var)) ))
nrow(subset(meta_PUR, (meta_PUR$sample_qc.r_het_hom_var < median(meta_PUR$sample_qc.r_het_hom_var) - 4*mad(meta_PUR$sample_qc.r_het_hom_var)) | (meta_PUR$sample_qc.r_het_hom_var > median(meta_PUR$sample_qc.r_het_hom_var) + 4*mad(meta_PUR$sample_qc.r_het_hom_var)) ))
nrow(subset(meta_SAS, (meta_SAS$sample_qc.r_het_hom_var < median(meta_SAS$sample_qc.r_het_hom_var) - 4*mad(meta_SAS$sample_qc.r_het_hom_var)) | (meta_SAS$sample_qc.r_het_hom_var > median(meta_SAS$sample_qc.r_het_hom_var) + 4*mad(meta_SAS$sample_qc.r_het_hom_var)) ))
nrow(subset(meta_Other, (meta_Other$sample_qc.r_het_hom_var < median(meta_Other$sample_qc.r_het_hom_var) - 4*mad(meta_Other$sample_qc.r_het_hom_var)) | (meta_Other$sample_qc.r_het_hom_var > median(meta_Other$sample_qc.r_het_hom_var) + 4*mad(meta_Other$sample_qc.r_het_hom_var)) ))


c(-1,1)*4*mad(meta_AFR$sample_qc.r_insertion_deletion) + median(meta_AFR$sample_qc.r_insertion_deletion)
c(-1,1)*4*mad(meta_AMR$sample_qc.r_insertion_deletion) + median(meta_AMR$sample_qc.r_insertion_deletion)
c(-1,1)*4*mad(meta_EAS$sample_qc.r_insertion_deletion) + median(meta_EAS$sample_qc.r_insertion_deletion)
c(-1,1)*4*mad(meta_EUR$sample_qc.r_insertion_deletion) + median(meta_EUR$sample_qc.r_insertion_deletion)
c(-1,1)*4*mad(meta_FIN$sample_qc.r_insertion_deletion) + median(meta_FIN$sample_qc.r_insertion_deletion)
c(-1,1)*4*mad(meta_PUR$sample_qc.r_insertion_deletion) + median(meta_PUR$sample_qc.r_insertion_deletion)
c(-1,1)*4*mad(meta_SAS$sample_qc.r_insertion_deletion) + median(meta_SAS$sample_qc.r_insertion_deletion)
c(-1,1)*4*mad(meta_Other$sample_qc.r_insertion_deletion) + median(meta_Other$sample_qc.r_insertion_deletion)

nrow(subset(meta_AFR, (meta_AFR$sample_qc.r_insertion_deletion < median(meta_AFR$sample_qc.r_insertion_deletion) - 4*mad(meta_AFR$sample_qc.r_insertion_deletion)) | (meta_AFR$sample_qc.r_insertion_deletion > median(meta_AFR$sample_qc.r_insertion_deletion) + 4*mad(meta_AFR$sample_qc.r_insertion_deletion)) ))
nrow(subset(meta_AMR, (meta_AMR$sample_qc.r_insertion_deletion < median(meta_AMR$sample_qc.r_insertion_deletion) - 4*mad(meta_AMR$sample_qc.r_insertion_deletion)) | (meta_AMR$sample_qc.r_insertion_deletion > median(meta_AMR$sample_qc.r_insertion_deletion) + 4*mad(meta_AMR$sample_qc.r_insertion_deletion)) ))
nrow(subset(meta_EAS, (meta_EAS$sample_qc.r_insertion_deletion < median(meta_EAS$sample_qc.r_insertion_deletion) - 4*mad(meta_EAS$sample_qc.r_insertion_deletion)) | (meta_EAS$sample_qc.r_insertion_deletion > median(meta_EAS$sample_qc.r_insertion_deletion) + 4*mad(meta_EAS$sample_qc.r_insertion_deletion)) ))
nrow(subset(meta_EUR, (meta_EUR$sample_qc.r_insertion_deletion < median(meta_EUR$sample_qc.r_insertion_deletion) - 4*mad(meta_EUR$sample_qc.r_insertion_deletion)) | (meta_EUR$sample_qc.r_insertion_deletion > median(meta_EUR$sample_qc.r_insertion_deletion) + 4*mad(meta_EUR$sample_qc.r_insertion_deletion)) ))
nrow(subset(meta_FIN, (meta_FIN$sample_qc.r_insertion_deletion < median(meta_FIN$sample_qc.r_insertion_deletion) - 4*mad(meta_FIN$sample_qc.r_insertion_deletion)) | (meta_FIN$sample_qc.r_insertion_deletion > median(meta_FIN$sample_qc.r_insertion_deletion) + 4*mad(meta_FIN$sample_qc.r_insertion_deletion)) ))
nrow(subset(meta_PUR, (meta_PUR$sample_qc.r_insertion_deletion < median(meta_PUR$sample_qc.r_insertion_deletion) - 4*mad(meta_PUR$sample_qc.r_insertion_deletion)) | (meta_PUR$sample_qc.r_insertion_deletion > median(meta_PUR$sample_qc.r_insertion_deletion) + 4*mad(meta_PUR$sample_qc.r_insertion_deletion)) ))
nrow(subset(meta_SAS, (meta_SAS$sample_qc.r_insertion_deletion < median(meta_SAS$sample_qc.r_insertion_deletion) - 4*mad(meta_SAS$sample_qc.r_insertion_deletion)) | (meta_SAS$sample_qc.r_insertion_deletion > median(meta_SAS$sample_qc.r_insertion_deletion) + 4*mad(meta_SAS$sample_qc.r_insertion_deletion)) ))
nrow(subset(meta_Other, (meta_Other$sample_qc.r_insertion_deletion < median(meta_Other$sample_qc.r_insertion_deletion) - 4*mad(meta_Other$sample_qc.r_insertion_deletion)) | (meta_Other$sample_qc.r_insertion_deletion > median(meta_Other$sample_qc.r_insertion_deletion) + 4*mad(meta_Other$sample_qc.r_insertion_deletion)) ))

# Total outliers
nrow(subset(meta_AFR, (meta_AFR$sample_qc.r_ti_tv < median(meta_AFR$sample_qc.r_ti_tv) - 4*mad(meta_AFR$sample_qc.r_ti_tv)) | (meta_AFR$sample_qc.r_ti_tv > median(meta_AFR$sample_qc.r_ti_tv) + 4*mad(meta_AFR$sample_qc.r_ti_tv)) |
                      (meta_AFR$sample_qc.r_het_hom_var < median(meta_AFR$sample_qc.r_het_hom_var) - 4*mad(meta_AFR$sample_qc.r_het_hom_var)) | (meta_AFR$sample_qc.r_het_hom_var > median(meta_AFR$sample_qc.r_het_hom_var) + 4*mad(meta_AFR$sample_qc.r_het_hom_var)) | 
                      (meta_AFR$sample_qc.r_insertion_deletion < median(meta_AFR$sample_qc.r_insertion_deletion) - 4*mad(meta_AFR$sample_qc.r_insertion_deletion)) | (meta_AFR$sample_qc.r_insertion_deletion > median(meta_AFR$sample_qc.r_insertion_deletion) + 4*mad(meta_AFR$sample_qc.r_insertion_deletion)) ))

nrow(subset(meta_AMR, (meta_AMR$sample_qc.r_ti_tv < median(meta_AMR$sample_qc.r_ti_tv) - 4*mad(meta_AMR$sample_qc.r_ti_tv)) | (meta_AMR$sample_qc.r_ti_tv > median(meta_AMR$sample_qc.r_ti_tv) + 4*mad(meta_AMR$sample_qc.r_ti_tv)) |
              (meta_AMR$sample_qc.r_het_hom_var < median(meta_AMR$sample_qc.r_het_hom_var) - 4*mad(meta_AMR$sample_qc.r_het_hom_var)) | (meta_AMR$sample_qc.r_het_hom_var > median(meta_AMR$sample_qc.r_het_hom_var) + 4*mad(meta_AMR$sample_qc.r_het_hom_var)) | 
              (meta_AMR$sample_qc.r_insertion_deletion < median(meta_AMR$sample_qc.r_insertion_deletion) - 4*mad(meta_AMR$sample_qc.r_insertion_deletion)) | (meta_AMR$sample_qc.r_insertion_deletion > median(meta_AMR$sample_qc.r_insertion_deletion) + 4*mad(meta_AMR$sample_qc.r_insertion_deletion)) ))

nrow(subset(meta_EAS, (meta_EAS$sample_qc.r_ti_tv < median(meta_EAS$sample_qc.r_ti_tv) - 4*mad(meta_EAS$sample_qc.r_ti_tv)) | (meta_EAS$sample_qc.r_ti_tv > median(meta_EAS$sample_qc.r_ti_tv) + 4*mad(meta_EAS$sample_qc.r_ti_tv)) |
              (meta_EAS$sample_qc.r_het_hom_var < median(meta_EAS$sample_qc.r_het_hom_var) - 4*mad(meta_EAS$sample_qc.r_het_hom_var)) | (meta_EAS$sample_qc.r_het_hom_var > median(meta_EAS$sample_qc.r_het_hom_var) + 4*mad(meta_EAS$sample_qc.r_het_hom_var)) | 
              (meta_EAS$sample_qc.r_insertion_deletion < median(meta_EAS$sample_qc.r_insertion_deletion) - 4*mad(meta_EAS$sample_qc.r_insertion_deletion)) | (meta_EAS$sample_qc.r_insertion_deletion > median(meta_EAS$sample_qc.r_insertion_deletion) + 4*mad(meta_EAS$sample_qc.r_insertion_deletion)) ))

nrow(subset(meta_EUR, (meta_EUR$sample_qc.r_ti_tv < median(meta_EUR$sample_qc.r_ti_tv) - 4*mad(meta_EUR$sample_qc.r_ti_tv)) | (meta_EUR$sample_qc.r_ti_tv > median(meta_EUR$sample_qc.r_ti_tv) + 4*mad(meta_EUR$sample_qc.r_ti_tv)) |
              (meta_EUR$sample_qc.r_het_hom_var < median(meta_EUR$sample_qc.r_het_hom_var) - 4*mad(meta_EUR$sample_qc.r_het_hom_var)) | (meta_EUR$sample_qc.r_het_hom_var > median(meta_EUR$sample_qc.r_het_hom_var) + 4*mad(meta_EUR$sample_qc.r_het_hom_var)) | 
              (meta_EUR$sample_qc.r_insertion_deletion < median(meta_EUR$sample_qc.r_insertion_deletion) - 4*mad(meta_EUR$sample_qc.r_insertion_deletion)) | (meta_EUR$sample_qc.r_insertion_deletion > median(meta_EUR$sample_qc.r_insertion_deletion) + 4*mad(meta_EUR$sample_qc.r_insertion_deletion)) ))

nrow(subset(meta_FIN, (meta_FIN$sample_qc.r_ti_tv < median(meta_FIN$sample_qc.r_ti_tv) - 4*mad(meta_FIN$sample_qc.r_ti_tv)) | (meta_FIN$sample_qc.r_ti_tv > median(meta_FIN$sample_qc.r_ti_tv) + 4*mad(meta_FIN$sample_qc.r_ti_tv)) |
              (meta_FIN$sample_qc.r_het_hom_var < median(meta_FIN$sample_qc.r_het_hom_var) - 4*mad(meta_FIN$sample_qc.r_het_hom_var)) | (meta_FIN$sample_qc.r_het_hom_var > median(meta_FIN$sample_qc.r_het_hom_var) + 4*mad(meta_FIN$sample_qc.r_het_hom_var)) | 
              (meta_FIN$sample_qc.r_insertion_deletion < median(meta_FIN$sample_qc.r_insertion_deletion) - 4*mad(meta_FIN$sample_qc.r_insertion_deletion)) | (meta_FIN$sample_qc.r_insertion_deletion > median(meta_FIN$sample_qc.r_insertion_deletion) + 4*mad(meta_FIN$sample_qc.r_insertion_deletion)) ))

nrow(subset(meta_PUR, (meta_PUR$sample_qc.r_ti_tv < median(meta_PUR$sample_qc.r_ti_tv) - 4*mad(meta_PUR$sample_qc.r_ti_tv)) | (meta_PUR$sample_qc.r_ti_tv > median(meta_PUR$sample_qc.r_ti_tv) + 4*mad(meta_PUR$sample_qc.r_ti_tv)) |
              (meta_PUR$sample_qc.r_het_hom_var < median(meta_PUR$sample_qc.r_het_hom_var) - 4*mad(meta_PUR$sample_qc.r_het_hom_var)) | (meta_PUR$sample_qc.r_het_hom_var > median(meta_PUR$sample_qc.r_het_hom_var) + 4*mad(meta_PUR$sample_qc.r_het_hom_var)) | 
              (meta_PUR$sample_qc.r_insertion_deletion < median(meta_PUR$sample_qc.r_insertion_deletion) - 4*mad(meta_PUR$sample_qc.r_insertion_deletion)) | (meta_PUR$sample_qc.r_insertion_deletion > median(meta_PUR$sample_qc.r_insertion_deletion) + 4*mad(meta_PUR$sample_qc.r_insertion_deletion)) ))

nrow(subset(meta_SAS, (meta_SAS$sample_qc.r_ti_tv < median(meta_SAS$sample_qc.r_ti_tv) - 4*mad(meta_SAS$sample_qc.r_ti_tv)) | (meta_SAS$sample_qc.r_ti_tv > median(meta_SAS$sample_qc.r_ti_tv) + 4*mad(meta_SAS$sample_qc.r_ti_tv)) |
              (meta_SAS$sample_qc.r_het_hom_var < median(meta_SAS$sample_qc.r_het_hom_var) - 4*mad(meta_SAS$sample_qc.r_het_hom_var)) | (meta_SAS$sample_qc.r_het_hom_var > median(meta_SAS$sample_qc.r_het_hom_var) + 4*mad(meta_SAS$sample_qc.r_het_hom_var)) | 
              (meta_SAS$sample_qc.r_insertion_deletion < median(meta_SAS$sample_qc.r_insertion_deletion) - 4*mad(meta_SAS$sample_qc.r_insertion_deletion)) | (meta_SAS$sample_qc.r_insertion_deletion > median(meta_SAS$sample_qc.r_insertion_deletion) + 4*mad(meta_SAS$sample_qc.r_insertion_deletion)) ))

nrow(subset(meta_Other, (meta_Other$sample_qc.r_ti_tv < median(meta_Other$sample_qc.r_ti_tv) - 4*mad(meta_Other$sample_qc.r_ti_tv)) | (meta_Other$sample_qc.r_ti_tv > median(meta_Other$sample_qc.r_ti_tv) + 4*mad(meta_Other$sample_qc.r_ti_tv)) |
              (meta_Other$sample_qc.r_het_hom_var < median(meta_Other$sample_qc.r_het_hom_var) - 4*mad(meta_Other$sample_qc.r_het_hom_var)) | (meta_Other$sample_qc.r_het_hom_var > median(meta_Other$sample_qc.r_het_hom_var) + 4*mad(meta_Other$sample_qc.r_het_hom_var)) | 
              (meta_Other$sample_qc.r_insertion_deletion < median(meta_Other$sample_qc.r_insertion_deletion) - 4*mad(meta_Other$sample_qc.r_insertion_deletion)) | (meta_Other$sample_qc.r_insertion_deletion > median(meta_Other$sample_qc.r_insertion_deletion) + 4*mad(meta_Other$sample_qc.r_insertion_deletion)) ))



# Filter outliers
meta_afr <- subset(meta_AFR, (meta_AFR$sample_qc.r_ti_tv >= median(meta_AFR$sample_qc.r_ti_tv) - 4*mad(meta_AFR$sample_qc.r_ti_tv)) & (meta_AFR$sample_qc.r_ti_tv <= median(meta_AFR$sample_qc.r_ti_tv) + 4*mad(meta_AFR$sample_qc.r_ti_tv)) &
                  (meta_AFR$sample_qc.r_het_hom_var >= median(meta_AFR$sample_qc.r_het_hom_var) - 4*mad(meta_AFR$sample_qc.r_het_hom_var)) & (meta_AFR$sample_qc.r_het_hom_var <= median(meta_AFR$sample_qc.r_het_hom_var) + 4*mad(meta_AFR$sample_qc.r_het_hom_var)) & 
                  (meta_AFR$sample_qc.r_insertion_deletion >= median(meta_AFR$sample_qc.r_insertion_deletion) - 4*mad(meta_AFR$sample_qc.r_insertion_deletion)) & (meta_AFR$sample_qc.r_insertion_deletion <= median(meta_AFR$sample_qc.r_insertion_deletion) + 4*mad(meta_AFR$sample_qc.r_insertion_deletion)) )

meta_amr <- subset(meta_AMR, (meta_AMR$sample_qc.r_ti_tv >= median(meta_AMR$sample_qc.r_ti_tv) - 4*mad(meta_AMR$sample_qc.r_ti_tv)) & (meta_AMR$sample_qc.r_ti_tv <= median(meta_AMR$sample_qc.r_ti_tv) + 4*mad(meta_AMR$sample_qc.r_ti_tv)) &
                  (meta_AMR$sample_qc.r_het_hom_var >= median(meta_AMR$sample_qc.r_het_hom_var) - 4*mad(meta_AMR$sample_qc.r_het_hom_var)) & (meta_AMR$sample_qc.r_het_hom_var <= median(meta_AMR$sample_qc.r_het_hom_var) + 4*mad(meta_AMR$sample_qc.r_het_hom_var)) & 
                  (meta_AMR$sample_qc.r_insertion_deletion >= median(meta_AMR$sample_qc.r_insertion_deletion) - 4*mad(meta_AMR$sample_qc.r_insertion_deletion)) & (meta_AMR$sample_qc.r_insertion_deletion <= median(meta_AMR$sample_qc.r_insertion_deletion) + 4*mad(meta_AMR$sample_qc.r_insertion_deletion)) )

meta_eas <- subset(meta_EAS, (meta_EAS$sample_qc.r_ti_tv >= median(meta_EAS$sample_qc.r_ti_tv) - 4*mad(meta_EAS$sample_qc.r_ti_tv)) & (meta_EAS$sample_qc.r_ti_tv <= median(meta_EAS$sample_qc.r_ti_tv) + 4*mad(meta_EAS$sample_qc.r_ti_tv)) &
                  (meta_EAS$sample_qc.r_het_hom_var >= median(meta_EAS$sample_qc.r_het_hom_var) - 4*mad(meta_EAS$sample_qc.r_het_hom_var)) & (meta_EAS$sample_qc.r_het_hom_var <= median(meta_EAS$sample_qc.r_het_hom_var) + 4*mad(meta_EAS$sample_qc.r_het_hom_var)) & 
                  (meta_EAS$sample_qc.r_insertion_deletion >= median(meta_EAS$sample_qc.r_insertion_deletion) - 4*mad(meta_EAS$sample_qc.r_insertion_deletion)) & (meta_EAS$sample_qc.r_insertion_deletion <= median(meta_EAS$sample_qc.r_insertion_deletion) + 4*mad(meta_EAS$sample_qc.r_insertion_deletion)) )

meta_eur <- subset(meta_EUR, (meta_EUR$sample_qc.r_ti_tv >= median(meta_EUR$sample_qc.r_ti_tv) - 4*mad(meta_EUR$sample_qc.r_ti_tv)) & (meta_EUR$sample_qc.r_ti_tv <= median(meta_EUR$sample_qc.r_ti_tv) + 4*mad(meta_EUR$sample_qc.r_ti_tv)) &
                  (meta_EUR$sample_qc.r_het_hom_var >= median(meta_EUR$sample_qc.r_het_hom_var) - 4*mad(meta_EUR$sample_qc.r_het_hom_var)) & (meta_EUR$sample_qc.r_het_hom_var <= median(meta_EUR$sample_qc.r_het_hom_var) + 4*mad(meta_EUR$sample_qc.r_het_hom_var)) & 
                  (meta_EUR$sample_qc.r_insertion_deletion >= median(meta_EUR$sample_qc.r_insertion_deletion) - 4*mad(meta_EUR$sample_qc.r_insertion_deletion)) & (meta_EUR$sample_qc.r_insertion_deletion <= median(meta_EUR$sample_qc.r_insertion_deletion) + 4*mad(meta_EUR$sample_qc.r_insertion_deletion)) )

meta_fin <- subset(meta_FIN, (meta_FIN$sample_qc.r_ti_tv >= median(meta_FIN$sample_qc.r_ti_tv) - 4*mad(meta_FIN$sample_qc.r_ti_tv)) & (meta_FIN$sample_qc.r_ti_tv <= median(meta_FIN$sample_qc.r_ti_tv) + 4*mad(meta_FIN$sample_qc.r_ti_tv)) &
                  (meta_FIN$sample_qc.r_het_hom_var >= median(meta_FIN$sample_qc.r_het_hom_var) - 4*mad(meta_FIN$sample_qc.r_het_hom_var)) & (meta_FIN$sample_qc.r_het_hom_var <= median(meta_FIN$sample_qc.r_het_hom_var) + 4*mad(meta_FIN$sample_qc.r_het_hom_var)) & 
                  (meta_FIN$sample_qc.r_insertion_deletion >= median(meta_FIN$sample_qc.r_insertion_deletion) - 4*mad(meta_FIN$sample_qc.r_insertion_deletion)) & (meta_FIN$sample_qc.r_insertion_deletion <= median(meta_FIN$sample_qc.r_insertion_deletion) + 4*mad(meta_FIN$sample_qc.r_insertion_deletion)) )

meta_pur <- subset(meta_PUR, (meta_PUR$sample_qc.r_ti_tv >= median(meta_PUR$sample_qc.r_ti_tv) - 4*mad(meta_PUR$sample_qc.r_ti_tv)) & (meta_PUR$sample_qc.r_ti_tv <= median(meta_PUR$sample_qc.r_ti_tv) + 4*mad(meta_PUR$sample_qc.r_ti_tv)) &
                  (meta_PUR$sample_qc.r_het_hom_var >= median(meta_PUR$sample_qc.r_het_hom_var) - 4*mad(meta_PUR$sample_qc.r_het_hom_var)) & (meta_PUR$sample_qc.r_het_hom_var <= median(meta_PUR$sample_qc.r_het_hom_var) + 4*mad(meta_PUR$sample_qc.r_het_hom_var)) & 
                  (meta_PUR$sample_qc.r_insertion_deletion >= median(meta_PUR$sample_qc.r_insertion_deletion) - 4*mad(meta_PUR$sample_qc.r_insertion_deletion)) & (meta_PUR$sample_qc.r_insertion_deletion <= median(meta_PUR$sample_qc.r_insertion_deletion) + 4*mad(meta_PUR$sample_qc.r_insertion_deletion)) )

meta_sas <- subset(meta_SAS, (meta_SAS$sample_qc.r_ti_tv >= median(meta_SAS$sample_qc.r_ti_tv) - 4*mad(meta_SAS$sample_qc.r_ti_tv)) & (meta_SAS$sample_qc.r_ti_tv <= median(meta_SAS$sample_qc.r_ti_tv) + 4*mad(meta_SAS$sample_qc.r_ti_tv)) &
                  (meta_SAS$sample_qc.r_het_hom_var >= median(meta_SAS$sample_qc.r_het_hom_var) - 4*mad(meta_SAS$sample_qc.r_het_hom_var)) & (meta_SAS$sample_qc.r_het_hom_var <= median(meta_SAS$sample_qc.r_het_hom_var) + 4*mad(meta_SAS$sample_qc.r_het_hom_var)) & 
                  (meta_SAS$sample_qc.r_insertion_deletion >= median(meta_SAS$sample_qc.r_insertion_deletion) - 4*mad(meta_SAS$sample_qc.r_insertion_deletion)) & (meta_SAS$sample_qc.r_insertion_deletion <= median(meta_SAS$sample_qc.r_insertion_deletion) + 4*mad(meta_SAS$sample_qc.r_insertion_deletion)) )

#meta_other <- subset(meta_Other, (meta_Other$sample_qc.r_ti_tv >= median(meta_Other$sample_qc.r_ti_tv) - 4*mad(meta_Other$sample_qc.r_ti_tv)) & (meta_Other$sample_qc.r_ti_tv <= median(meta_Other$sample_qc.r_ti_tv) + 4*mad(meta_Other$sample_qc.r_ti_tv)) &
#                    (meta_Other$sample_qc.r_het_hom_var >= median(meta_Other$sample_qc.r_het_hom_var) - 4*mad(meta_Other$sample_qc.r_het_hom_var)) & (meta_Other$sample_qc.r_het_hom_var <= median(meta_Other$sample_qc.r_het_hom_var) + 4*mad(meta_Other$sample_qc.r_het_hom_var)) & 
#                    (meta_Other$sample_qc.r_insertion_deletion >= median(meta_Other$sample_qc.r_insertion_deletion) - 4*mad(meta_Other$sample_qc.r_insertion_deletion)) & (meta_Other$sample_qc.r_insertion_deletion <= median(meta_Other$sample_qc.r_insertion_deletion) + 4*mad(meta_Other$sample_qc.r_insertion_deletion)) )

meta_clean <- rbind(meta_afr, meta_amr)
meta_clean <- rbind(meta_clean, meta_eas)
meta_clean <- rbind(meta_clean, meta_eur)
meta_clean <- rbind(meta_clean, meta_fin)
meta_clean <- rbind(meta_clean, meta_pur)
meta_clean <- rbind(meta_clean, meta_sas)
meta_clean <- rbind(meta_clean, meta_Other)

colnames(meta_clean)
meta_clean$id = c(1:nrow(meta_clean))


meta_washu <- subset(meta_WashU, (meta_WashU$sample_qc.r_ti_tv >= median(meta_WashU$sample_qc.r_ti_tv) - 4*mad(meta_WashU$sample_qc.r_ti_tv)) & (meta_WashU$sample_qc.r_ti_tv <= median(meta_WashU$sample_qc.r_ti_tv) + 4*mad(meta_WashU$sample_qc.r_ti_tv)) &
                     (meta_WashU$sample_qc.r_het_hom_var >= median(meta_WashU$sample_qc.r_het_hom_var) - 4*mad(meta_WashU$sample_qc.r_het_hom_var)) & (meta_WashU$sample_qc.r_het_hom_var <= median(meta_WashU$sample_qc.r_het_hom_var) + 4*mad(meta_WashU$sample_qc.r_het_hom_var)) & 
                     (meta_WashU$sample_qc.r_insertion_deletion >= median(meta_WashU$sample_qc.r_insertion_deletion) - 4*mad(meta_WashU$sample_qc.r_insertion_deletion)) & (meta_WashU$sample_qc.r_insertion_deletion <= median(meta_WashU$sample_qc.r_insertion_deletion) + 4*mad(meta_WashU$sample_qc.r_insertion_deletion)) )
meta_nygc <- subset(meta_NYGC, (meta_NYGC$sample_qc.r_ti_tv >= median(meta_NYGC$sample_qc.r_ti_tv) - 4*mad(meta_NYGC$sample_qc.r_ti_tv)) & (meta_NYGC$sample_qc.r_ti_tv <= median(meta_NYGC$sample_qc.r_ti_tv) + 4*mad(meta_NYGC$sample_qc.r_ti_tv)) &
                     (meta_NYGC$sample_qc.r_het_hom_var >= median(meta_NYGC$sample_qc.r_het_hom_var) - 4*mad(meta_NYGC$sample_qc.r_het_hom_var)) & (meta_NYGC$sample_qc.r_het_hom_var <= median(meta_NYGC$sample_qc.r_het_hom_var) + 4*mad(meta_NYGC$sample_qc.r_het_hom_var)) & 
                     (meta_NYGC$sample_qc.r_insertion_deletion >= median(meta_NYGC$sample_qc.r_insertion_deletion) - 4*mad(meta_NYGC$sample_qc.r_insertion_deletion)) & (meta_NYGC$sample_qc.r_insertion_deletion <= median(meta_NYGC$sample_qc.r_insertion_deletion) + 4*mad(meta_NYGC$sample_qc.r_insertion_deletion)) )
meta_baylor <- subset(meta_Baylor, (meta_Baylor$sample_qc.r_ti_tv >= median(meta_Baylor$sample_qc.r_ti_tv) - 4*mad(meta_Baylor$sample_qc.r_ti_tv)) & (meta_Baylor$sample_qc.r_ti_tv <= median(meta_Baylor$sample_qc.r_ti_tv) + 4*mad(meta_Baylor$sample_qc.r_ti_tv)) &
                     (meta_Baylor$sample_qc.r_het_hom_var >= median(meta_Baylor$sample_qc.r_het_hom_var) - 4*mad(meta_Baylor$sample_qc.r_het_hom_var)) & (meta_Baylor$sample_qc.r_het_hom_var <= median(meta_Baylor$sample_qc.r_het_hom_var) + 4*mad(meta_Baylor$sample_qc.r_het_hom_var)) & 
                     (meta_Baylor$sample_qc.r_insertion_deletion >= median(meta_Baylor$sample_qc.r_insertion_deletion) - 4*mad(meta_Baylor$sample_qc.r_insertion_deletion)) & (meta_Baylor$sample_qc.r_insertion_deletion <= median(meta_Baylor$sample_qc.r_insertion_deletion) + 4*mad(meta_Baylor$sample_qc.r_insertion_deletion)) )
meta_baylor_hiseqx <- subset(meta_Baylor_HiSeqX, (meta_Baylor_HiSeqX$sample_qc.r_ti_tv >= median(meta_Baylor_HiSeqX$sample_qc.r_ti_tv) - 4*mad(meta_Baylor_HiSeqX$sample_qc.r_ti_tv)) & (meta_Baylor_HiSeqX$sample_qc.r_ti_tv <= median(meta_Baylor_HiSeqX$sample_qc.r_ti_tv) + 4*mad(meta_Baylor_HiSeqX$sample_qc.r_ti_tv)) &
                     (meta_Baylor_HiSeqX$sample_qc.r_het_hom_var >= median(meta_Baylor_HiSeqX$sample_qc.r_het_hom_var) - 4*mad(meta_Baylor_HiSeqX$sample_qc.r_het_hom_var)) & (meta_Baylor_HiSeqX$sample_qc.r_het_hom_var <= median(meta_Baylor_HiSeqX$sample_qc.r_het_hom_var) + 4*mad(meta_Baylor_HiSeqX$sample_qc.r_het_hom_var)) & 
                     (meta_Baylor_HiSeqX$sample_qc.r_insertion_deletion >= median(meta_Baylor_HiSeqX$sample_qc.r_insertion_deletion) - 4*mad(meta_Baylor_HiSeqX$sample_qc.r_insertion_deletion)) & (meta_Baylor_HiSeqX$sample_qc.r_insertion_deletion <= median(meta_Baylor_HiSeqX$sample_qc.r_insertion_deletion) + 4*mad(meta_Baylor_HiSeqX$sample_qc.r_insertion_deletion)) )
meta_baylor_novaseq <- subset(meta_Baylor_NovaSeq, (meta_Baylor_NovaSeq$sample_qc.r_ti_tv >= median(meta_Baylor_NovaSeq$sample_qc.r_ti_tv) - 4*mad(meta_Baylor_NovaSeq$sample_qc.r_ti_tv)) & (meta_Baylor_NovaSeq$sample_qc.r_ti_tv <= median(meta_Baylor_NovaSeq$sample_qc.r_ti_tv) + 4*mad(meta_Baylor_NovaSeq$sample_qc.r_ti_tv)) &
                     (meta_Baylor_NovaSeq$sample_qc.r_het_hom_var >= median(meta_Baylor_NovaSeq$sample_qc.r_het_hom_var) - 4*mad(meta_Baylor_NovaSeq$sample_qc.r_het_hom_var)) & (meta_Baylor_NovaSeq$sample_qc.r_het_hom_var <= median(meta_Baylor_NovaSeq$sample_qc.r_het_hom_var) + 4*mad(meta_Baylor_NovaSeq$sample_qc.r_het_hom_var)) & 
                     (meta_Baylor_NovaSeq$sample_qc.r_insertion_deletion >= median(meta_Baylor_NovaSeq$sample_qc.r_insertion_deletion) - 4*mad(meta_Baylor_NovaSeq$sample_qc.r_insertion_deletion)) & (meta_Baylor_NovaSeq$sample_qc.r_insertion_deletion <= median(meta_Baylor_NovaSeq$sample_qc.r_insertion_deletion) + 4*mad(meta_Baylor_NovaSeq$sample_qc.r_insertion_deletion)) )
meta_broad <- subset(meta_Broad, (meta_Broad$sample_qc.r_ti_tv >= median(meta_Broad$sample_qc.r_ti_tv) - 4*mad(meta_Broad$sample_qc.r_ti_tv)) & (meta_Broad$sample_qc.r_ti_tv <= median(meta_Broad$sample_qc.r_ti_tv) + 4*mad(meta_Broad$sample_qc.r_ti_tv)) &
                     (meta_Broad$sample_qc.r_het_hom_var >= median(meta_Broad$sample_qc.r_het_hom_var) - 4*mad(meta_Broad$sample_qc.r_het_hom_var)) & (meta_Broad$sample_qc.r_het_hom_var <= median(meta_Broad$sample_qc.r_het_hom_var) + 4*mad(meta_Broad$sample_qc.r_het_hom_var)) & 
                     (meta_Broad$sample_qc.r_insertion_deletion >= median(meta_Broad$sample_qc.r_insertion_deletion) - 4*mad(meta_Broad$sample_qc.r_insertion_deletion)) & (meta_Broad$sample_qc.r_insertion_deletion <= median(meta_Broad$sample_qc.r_insertion_deletion) + 4*mad(meta_Broad$sample_qc.r_insertion_deletion)) )


meta_clean_center <- rbind(meta_baylor, meta_baylor_hiseqx)
meta_clean_center <- rbind(meta_clean_center, meta_baylor_novaseq)
meta_clean_center <- rbind(meta_clean_center, meta_broad)
meta_clean_center <- rbind(meta_clean_center, meta_nygc)
meta_clean_center <- rbind(meta_clean_center, meta_washu)

colnames(meta_clean_center)
meta_clean_center$id = c(1:nrow(meta_clean_center))




# Step 4: Generate plots ########################################
# Plots by population with center as subgroup
ggplot(meta_new, 
       aes(x=id, y=sample_qc.r_ti_tv, col=Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ti / Tv") +
  ylim(1.840,2.045) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_new, 
       aes(x=id, y=sample_qc.r_het_hom_var, col=Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Het. / Hom. Var.") +
  ylim(0.00,3.15) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_new, 
       aes(x=id, y=sample_qc.r_insertion_deletion, col=Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ins. / Del.") +
  ylim(0.80,1.27) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 


# Plots by center with population as subgroup
ggplot(meta_new_center, 
       aes(x=id, y=sample_qc.r_ti_tv, col=center_platform)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ti / Tv") +
  ylim(1.840,2.045) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_new_center, 
       aes(x=id, y=sample_qc.r_het_hom_var, col=center_platform)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Het. / Hom. Var.") +
  ylim(0.00,3.15) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_new_center, 
       aes(x=id, y=sample_qc.r_insertion_deletion, col=center_platform)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ins. / Del.") +
  ylim(0.80,1.27) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 


# Plots by center and population together
bp1 <- ggplot(meta_new_center, 
       aes(x = id, y=sample_qc.r_ti_tv, col=Pop )) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ti / Tv") +
  ylim(1.840,2.045) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(0.7)),
        axis.title = element_text(size=rel(1.5))) 
bp1 + facet_grid(. ~ center_platform, scales='free_x' )

bp2 <- ggplot(meta_new_center, 
       aes(x=id, y=sample_qc.r_het_hom_var, col = Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Het. / Hom. Var.") +
  ylim(0.00,3.15) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(0.7)),
        axis.title = element_text(size=rel(1.5))) 
bp2 + facet_grid(. ~ center_platform, scales='free_x' )

bp3 <- ggplot(meta_new_center, 
       aes(x=id, y=sample_qc.r_insertion_deletion, col = Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ins. / Del.") +
  ylim(0.80,1.27) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(0.7)),
        axis.title = element_text(size=rel(1.5))) 
bp3 + facet_grid(. ~ center_platform, scales='free_x' )


# Step 5: Remove outliers ########################################
# Plots by population with center as subgroup
ggplot(meta_clean, 
       aes(x=id, y=sample_qc.r_ti_tv, col=Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ti / Tv") +
  ylim(1.70,2.045) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_clean, 
       aes(x=id, y=sample_qc.r_het_hom_var, col=Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Het. / Hom. Var.") +
  ylim(0.00,3.15) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_clean, 
       aes(x=id, y=sample_qc.r_insertion_deletion, col=Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ins. / Del.") +
  ylim(0.80,1.45) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 


# Plots by center with population as subgroup
ggplot(meta_clean_center, 
       aes(x=id, y=sample_qc.r_ti_tv, col=center_platform)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ti / Tv") +
  ylim(1.840,2.045) +
  ggtitle("Before QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_clean_center, 
       aes(x=id, y=sample_qc.r_het_hom_var, col=center_platform)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Het. / Hom. Var.") +
  ylim(0.00,3.15) +
  ggtitle("After QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 

ggplot(meta_clean_center, 
       aes(x=id, y=sample_qc.r_insertion_deletion, col=center_platform)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ins. / Del.") +
  ylim(0.80,1.27) +
  ggtitle("After QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(1.5))) 


# Plots by center and population together
bp1 <- ggplot(meta_clean_center, 
       aes(x=id, y=sample_qc.r_ti_tv, col=Pop )) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ti / Tv") +
  ylim(1.840,2.045) +
  ggtitle("After QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(0.7)),
        axis.title = element_text(size=rel(1.5))) 
bp1 + facet_grid(. ~ center_platform, scales='free_x' )

bp2 <- ggplot(meta_clean_center, 
       aes(x=id, y=sample_qc.r_het_hom_var, col = Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Het. / Hom. Var.") +
  ylim(0.00,3.15) +
  ggtitle("After QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(0.7)),
        axis.title = element_text(size=rel(1.5))) 
bp2 + facet_grid(. ~ center_platform, scales='free_x' )

bp3 <- ggplot(meta_clean_center, 
       aes(x=id, y=sample_qc.r_insertion_deletion, col = Pop)) + 
  geom_point(size=3.0) + 
  xlab("Sample ID") + ylab("Ratio Ins. / Del.") +
  ylim(0.80,1.27) +
  ggtitle("After QC (Chromosome 1~22)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(0.7)),
        axis.title = element_text(size=rel(1.5))) 
bp3 + facet_grid(. ~ center_platform, scales='free_x' )






