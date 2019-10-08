##########################################
# Topic: Tidy cohort and center data
#  Date: 08/29/2019
#    By: Yushi F.T.
##########################################

##########################################
raw_data <- read.table("../data/CCDG_Freeze_2_Manifest_Expanded_2019-08-07.txt",
                       fill = T, header=T)
temp <- raw_data[,c(1,2,3,5,6,55)]

raw_data <- read.csv('../data/CCDG_Freeze2_manifest_simple.csv')
temp <- raw_data
write.csv(temp, '../data/CCDG_Freeze_2_Manifest_Expanded.csv',
          row.names = F)
##########################################



#################################################
library(dplyr)

center <- read.csv('../data/CCDG_Freeze_2_Manifest_Expanded.csv', header=T)
sample <- read.table('../work/qced_1_sample_qc_info_postqc.txt', header=T)
new_data_1 <- as.data.frame(sample$s)
new_data_2 <- as.data.frame(sample$s)
new_data_3 <- as.data.frame(sample$s)
colnames(new_data_1) <- c('sample_id')
colnames(new_data_2) <- c('sample_id')
colnames(new_data_3) <- c('sample_id')


new_data_1 <- merge(new_data_1, center, by.x='sample_id', by.y = 'sample_id')
new_data_2 <- merge(new_data_2, center, by.x='sample_id', by.y = 'alternate_sample_id_1')
new_data_3 <- merge(new_data_3, center, by.x='sample_id', by.y = 'subject_id')

new_temp_1 <- new_data_1[,c(1,2,3,6)]
new_temp_2 <- new_data_2[,c(1,2,3,6)]
new_temp_3 <- new_data_3[,c(1,2,3,6)]

unique <- rbind(new_temp_1, new_temp_2)
unique <- rbind(unique, new_temp_3)

unique <- distinct(unique)
unique <- unique[!duplicated(unique$sample_id),]



# Checking missing data from center2

center2 <- read.csv('../data/phenotype_ccdgf2_tidy.csv',header=T)

missing <- sample$s[!(sample$s %in% unique$sample_id) & (sample$s %in% center2$sample_id)]

missing_temp <- center2[center2$sample_id %in% missing,c(1,3,4)]

missing_temp$center <- "WashU"
missing_temp$platform <- NA

colnames(unique) <- colnames(missing_temp)
out_data <- rbind(unique, missing_temp)

## Summary platform data by sequencing centers ####################
with(out_data,
     table(platform, center, useNA = "ifany"))
summary(out_data$center)
# Baylor: 12800
# Broad:  11762
# NYGC:   20059
# WashU:  15924
# For samples in Baylor center, 7206 sequenced by HiSeqX,
# 1169 sequenced by NovaSeq, 4425 marked with NA


## Create mark for center with sequencing platform ################
out_data$center_platform <- c(1:nrow(out_data))
out_data$center_platform <- ifelse(out_data$center=='NYGC', 'NYGC', out_data$center_platform)
out_data$center_platform <- ifelse(out_data$center=='WashU', 'WASHU', out_data$center_platform)
out_data$center_platform <- ifelse(out_data$center=='Broad', 'BROAD', out_data$center_platform)
out_data$center_platform <- ifelse((out_data$center=='Baylor'), 'BAYLOR', out_data$center_platform)
out_data$center_platform <- ifelse((out_data$center=='Baylor') & (out_data$platform=='HiSeqX'), 'BAYLOR_HiSeqX', out_data$center_platform)
out_data$center_platform <- ifelse((out_data$center=='Baylor') & (out_data$platform=='NovaSeq'), 'BAYLOR_NovaSeq', out_data$center_platform)
out_data$center_platform <- ifelse((out_data$center=='Baylor') & is.na(out_data$center_platform), 'BAYLOR', out_data$center_platform)

table(out_data$center_platform)

write.csv(out_data, '../data/sample_center_platform.csv', 
          row.names=F)



## Create mark for seperate cohorts  ####################################
out_data$mark <- c(1:nrow(out_data))
out_data$mark <- ifelse(out_data$center=='NYGC', 'NYGC', out_data$mark)
out_data$mark <- ifelse(out_data$center=='WashU', 'WASHU', out_data$mark)
out_data$mark <- ifelse(out_data$center=='Baylor', 'BAYLOR', out_data$mark)
out_data$mark <- ifelse(out_data$center=='Broad', 'BROAD', out_data$mark)


out_data$mark <- ifelse(out_data$study=='aric', 'BAYLOR_ARIC', out_data$mark)
out_data$mark <- ifelse(out_data$study=='sol', 'BAYLOR_SOL', out_data$mark)
out_data$mark <- ifelse(out_data$study=='ssc', 'NYGC_SSC', out_data$mark)
out_data$mark <- ifelse(out_data$study=='gala2', 'NYGC_GALAII', out_data$mark)
out_data$mark <- ifelse(out_data$study=='af-lmu', 'BROAD_AFLMU', out_data$mark)
out_data$mark <- ifelse(out_data$study=='taichi', 'BROAD_TAICHI', out_data$mark)
out_data$mark <- ifelse(out_data$study=='virgo', 'BROAD_VIRGO', out_data$mark)
out_data$mark <- ifelse(out_data$study=='cleveland', 'WASHU_CLEVELAND', out_data$mark)

out_data$mark <- ifelse(out_data$sample_id=='SSC06708', 'withdraw', out_data$mark)
out_data$mark <- ifelse(out_data$sample_id=='SSC06703', 'withdraw', out_data$mark)
out_data$mark <- ifelse(out_data$sample_id=='SSC06699', 'withdraw', out_data$mark)
out_data$mark <- ifelse(out_data$sample_id=='SSC06709', 'withdraw', out_data$mark)

write.table(out_data[,c(1,4)], '../data/sample_marker.txt', row.names=F, quote=FALSE)


## Check the data ###################################

marker <- read.table('../data/sample_marker.txt', header=T)
sample <- read.table('../work/qced_1_sample_qc_info_postqc.txt', header=T)
new_data <- as.data.frame(sample$s)
colnames(new_data) <- c('sample_id')

new_data <- merge(new_data, marker, by='sample_id')




