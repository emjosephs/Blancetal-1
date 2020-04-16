# This script takes the raw RNAseq expression data and pulls out the lines for which we have complete data (no zero, no missing values) for each tissue type
# This script was used to generate the files i "../data/Raw_expression/" and after mean-centering "../data/Mean_centered_expression/"

# Read in libraries
library(data.table)
library(dplyr)
library(stringr)

## Read in Expression Data
expression_data <- fread("data/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded.txt")

##Subset data to include only data from one tissue type
Gene <- expression_data
dv <- str_split_fixed(Gene$`<Trait>`, "_", 5)
dv <- as.data.frame(dv)
dv <- dv[,2:5]
dv$V6 <- seq(1:1960)
dv[dv==""] <- NA
dv$V5 <- as.character(dv$V5)
s1 <- subset(dv, V5 != 'NA')
s1 <- s1[,2:5]
colnames(s1) <- c("V1", "V2", "V3", "V4")
s2 <- dv[rowSums(is.na(dv)) > 0,]
sj <- s2[,5]
s2 <- s2[,1:3]
s2 <- cbind(s2,sj)
colnames(s2) <- c("V1", "V2", "V3", "V4")
dv1 <- rbind(s1, s2)
dv1 <- dv1[order(dv1$V4),]
Gene <- Gene[ ,2:37128]
dv <- cbind(dv1, Gene)

# Subset to desired tissue
subs <- subset(dv, V1 == "Kern")

## Read in file that has line names in order
names <- read.table("data/LineNames.txt")

## Get the names of all of the genes
col_names <- colnames(subs[,5:ncol(subs)])

#remove duplicates (picking the first one)
subs1 = subs[!duplicated(subs$V2),]

## Order Values -- We need to pull of the expression value for each of the names in the jen263 list
df1 <- dplyr::left_join(names, subs1, by = c('V1'='V2'))

## Replace the old names with the ones that match from the Kinship matrix
conversion <- read.csv("data/name_conversions.csv", header = F)
new_names <- conversion$V2
row.names(df1) <- new_names

#remove NAs
df1 <- df1[!is.na(df1$V1.y),]
df1 <- df1[,5:ncol(df1)]

# Remove Columns that have any zeros
df1 <- df1[,which(as.numeric(colSums(df1 != 0)) == nrow(df1))]

# Save raw expression data
write.table(df1, file = "data/Raw_expression/Kern.txt")

# Mean center the data
for (i in 1:ncol(df1)){
  df1[,i] <- scale(df1[,i], scale = FALSE)
}

# Remove last row
df1 <- df1[-nrow(df1),]

# Save mean centered expression
write.table(df1, file = "data/Mean_centered_expression/Kern.txt")
