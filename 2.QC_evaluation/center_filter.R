##########################################
# Topic: Tidy center data
#  Date: 08/19/2019
#    By: Yushi F.T.
##########################################

center_data <- read.table('../data/phenotype_ccdgf2.txt', header=T, fill=T)
center_nygc <- subset(center_data, center_data$center == 'NYGC')
center_washu <- subset(center_data, center_data$center == 'WashU')
center_baylor <- subset(center_data, center_data$center == 'Baylor')
center_broad <- subset(center_data, center_data$center == 'Broad')
center_page <- subset(center_data, center_data$center == 'PAGE')
center_na <- subset(center_data, center_data$center == 'NA')

center_new <- rbind(center_nygc, center_washu)
center_new <- rbind(center_new, center_baylor)
center_new <- rbind(center_new, center_broad)
center_new <- rbind(center_new, center_page)

write.csv(center_new, '../data/phenotype_ccdgf2_tidy.csv', row.names=F)


