library(tidyverse)
library(ggplot2)
library(patchwork)
library(edgeR)

## Loading in data
adata_cardiac <- read.csv('adata_cardiac_obs.csv') %>% rename(cell = X)

## Calculating cell type proportions for CM and CF cells
samples <- unique(adata_cardiac$individual)
days <- unique(adata_cardiac$diffday)

df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df) <- c("day", "sample", "proportion", "cell_type")

for (i in days){
  for (j in samples){
    total_rows <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j,])
    CM_rows <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j & adata_cardiac$type=='CM',])
    CF_rows <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j & adata_cardiac$type=='CF',])
    fraction_cm <- CM_rows/total_rows
    fraction_cf <- CF_rows/total_rows
    df[nrow(df) + 1,] <-  list(as.numeric(gsub('day', '', i)), as.character(j), fraction_cm, 'cm')
    df[nrow(df) + 1,] <-  list(as.numeric(gsub('day', '', i)), as.character(j), fraction_cf, 'cf')
  }
}

## Plotting proportion of CM/CF cells
ggplot(data=df[df$cell_type=='cm' & !is.na(df$proportion),], aes(x=day, y=proportion, color=sample)) + 
  geom_point()+
  geom_line(alpha=0.8) +
  theme_classic()+
  scale_color_viridis_d(option="C") +
  xlab('Day') + 
  ylab('Proportion') + 
  theme(legend.position = "none") +
  ggtitle('Proportion of CM cells') +
  geom_hline(yintercept=0.5) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7,9,11,15))

ggplot(data=df[df$cell_type=='cf' & !is.na(df$proportion),], aes(x=day, y=proportion, color=sample)) + 
  geom_point() +
  geom_line(alpha=0.8) +
  theme_classic()+
  scale_color_viridis_d(option="C") +
  xlab('Day') + 
  ylab('Proportion') + 
  theme(legend.position = "none") + 
  ggtitle('Proportion of CF cells') +
  geom_hline(yintercept=0.5) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7,9,11,15))


## Calculating combined CM and CF Proportion
samples <- unique(adata_cardiac$individual)
days <- unique(adata_cardiac$diffday)

df_total <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df_total) <- c("day", "sample", "proportion")

for (i in days){
  for (j in samples){
    total_rows <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j,])
    CM_rows <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j & adata_cardiac$type=='CM',])
    CF_rows <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j & adata_cardiac$type=='CF',])
    fraction_cm <- CM_rows/total_rows
    fraction_cf <- CF_rows/total_rows
    df_total[nrow(df_total) + 1,] <-  list(as.numeric(gsub('day', '', i)), as.character(j), fraction_cf + fraction_cm)
  }
}

## Plotting Combined CM and CF Proportion
ggplot(data=df_total[!is.na(df_total$proportion),], aes(x=day, y=proportion, color=sample)) + 
  geom_point() +
  geom_line(alpha=0.8) +
  scale_color_viridis_d(option="C") +
  theme_classic()+
  theme(legend.position = "none") +
  xlab('Day')+
  ylab('Proportion CM + CF') + 
  scale_x_continuous(breaks = c(0, 1, 3, 5, 7,9,11,15))


## Classifying cells 

# Classifying based on D15 proportion 
# Cell line considered biased towards a cell fate if proportion of cells adopting that fate > 0.5

# CM-biased lines
df[df$day==15 & df$cell_type=='cm' & df$proportion>0.5, 2]
length(df[df$day==15 & df$cell_type=='cm' & df$proportion>0.5, 2])
# CF-biased lines
df[df$day==15 & df$cell_type=='cf' & df$proportion>0.5, 2]
length(df[df$day==15 & df$cell_type=='cf' & df$proportion>0.5, 2])


# Classifying based on D15 proportions
# Cell line considered biased towards cm or cf fate if proportion cm > proportion cf at d15 and vice versa

# CM-biased lines
idx = df[df$day==15 & df$cell_type=='cm' , 3] > df[df$day==15 & df$cell_type=='cf' , 3]
df[df$day==15 & df$cell_type=='cf', 2][idx]
length(df[df$day==15 & df$cell_type=='cf', 2][idx])
# CF-biased lines
idx = df[df$day==15 & df$cell_type=='cm' , 3] > df[df$day==15 & df$cell_type=='cf' , 3]
df[df$day==15 & df$cell_type=='cf', 2][!idx]
length(df[df$day==15 & df$cell_type=='cf', 2][!idx])


# Classifying based on D11 proportion 
# Cell line considered biased towards a cell fate if proportion of cells adopting that fate > 0.5

#CM-biased lines
df[df$day==11 & df$cell_type=='cm' & df$proportion>0.5, 2]
length(df[df$day==11 & df$cell_type=='cm' & df$proportion>0.5, 2])
# CF-biased lines
df[df$day==11 & df$cell_type=='cf' & df$proportion>0.5, 2]
length(df[df$day==11 & df$cell_type=='cf' & df$proportion>0.5, 2])


# Classifying based on D11 proportion 
# Cell line considered biased towards cm or cf fate if proportion cm > proportion cf at d15 and vice versa

# CM-biased lines
idx = df[df$day==11 & df$cell_type=='cm' , 3] > df[df$day==15 & df$cell_type=='cf' , 3]
df[df$day==11 & df$cell_type=='cf', 2][idx]
length(df[df$day==11 & df$cell_type=='cf', 2][idx])
# CF-biased lines
idx = df[df$day==11 & df$cell_type=='cm' , 3] > df[df$day==15 & df$cell_type=='cf' , 3]
df[df$day==11 & df$cell_type=='cf', 2][!idx]
length(df[df$day==11 & df$cell_type=='cf', 2][!idx])

## Working with later time points

# Getting data for cell type proportions across different time points
samples <- unique(adata_cardiac$individual)
days <- unique(adata_cardiac$diffday)
types <- unique(adata_cardiac$type)

df_types <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_types) <- c("day", "sample", "proportion", "cell_type")

for (i in days){
  for (j in samples){
    total_rows <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j,])
    for (k in types){
      type_row <- nrow(adata_cardiac[adata_cardiac$diffday==i & adata_cardiac$individual==j & adata_cardiac$type==k,])
      fraction <- type_row/total_rows
      df_types[nrow(df_types) + 1,] <-  list(as.numeric(gsub('day', '', i)), as.character(j), fraction, k)
    }
  }
}

# Plotting
ggplot(data=df_types[df_types$sample=='18520' ,], aes(x=day, y=proportion)) + geom_col(aes(fill=cell_type))

do.call(patchwork::wrap_plots, lapply(samples, function(x) {
  ggplot(df_types[df_types$sample==x ,], aes(x=day, y=proportion)) + 
    geom_col(aes(fill=cell_type)) +
    theme_minimal() + 
    theme(legend.position = "none") + 
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
})) 



### Aggregating time points and classifying bias

# Classifying based on D7-15 proportion 
# Cell line considered biased towards cm or cf fate if proportion cm > proportion cf at d15 and vice versa
df_to_agg <- df[df$day >= 7 & df$cell_type=='cm',] 
df_agg_cm= aggregate(df_to_agg,
                     by = list(df_to_agg$sample),
                     FUN = mean)

df_to_agg <- df[df$day >= 7 & df$cell_type=='cf',] 
df_agg_cf= aggregate(df_to_agg,
                     by = list(df_to_agg$sample),
                     FUN = mean)

# CM-biased lines
idx = df_agg_cm[, 4] > df_agg_cf[, 4]
df_agg_cm[, 1][idx]
length(df_agg_cm[, 1][idx])

# CF-biased lines
idx = df_agg_cm[, 4] < df_agg_cf[, 4]
df_agg_cm[, 1][idx]
length(df_agg_cm[, 1][idx])



#### Working with Bulk and Pseudobulk data

## edgeR Normalization

# Loading in data set
bulk <- read.csv('bulk_counts.tsv', header=TRUE, sep='\t', check.names = FALSE)
pseudobulk <- read.csv('pseduobulk_counts.tsv', header=TRUE, sep='\t', check.names = FALSE)

# Renaming column names so that they can be merged together
colnames(bulk)[-(1)] <- paste("bulk", colnames(bulk)[-(1)] , sep = "_")
colnames(pseudobulk)[-(1)] <- paste("pseudobulk", colnames(pseudobulk)[-(1)] , sep = "_")

# Merging data
combined <- merge(bulk, pseudobulk, by='gene') %>% 
  column_to_rownames(var="gene") 

# Merge adds 'x' and 'y' characters to tail of colnames. Remove and convert to matrix for edger input.
colnames(combined) <- gsub('.x','',names(combined))
colnames(combined) <- gsub('.y','',names(combined))
combined <- as.matrix(combined)

# Read into EdgeR and normalize using TMM
y <- DGEList(counts=combined)
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.sizes = TRUE, method='TMM') 

# Get all D0 time points for bulk and pseudobulk
idx <- grep("_0", colnames(z))
d0_data <- z[,idx]

## Calculating Cell Type Bias from Days 7-15

# Read in bulk and pseudobulk cell type data
bulk_cell_types <- read.csv('~/rotations/battle/bulk.inferred.tsv', sep='\t', header=TRUE)
pseudobulk_cell_types <- read.csv('~/rotations/battle/pseudobulk_indday.full.true.tsv', sep='\t', header=TRUE)

# get day 7-11 data
bulk_terminal <- bulk_cell_types[grep("_7|_11|_15", bulk_cell_types$Mixture), c('Mixture', 'CM', 'CF')]
pseudobulk_terminal <- pseudobulk_cell_types[grep("_7|_11|_15", pseudobulk_cell_types$sample), c('sample', 'CM', 'CF')]

# get day 0 data
bulk_d0 <- bulk_cell_types[grep("_0", bulk_cell_types$Mixture), c('Mixture', 'IPSC')]
pseudobulk_d0 <- pseudobulk_cell_types[grep("_0", pseudobulk_cell_types$sample), c('sample', 'IPSC')]

## Bulk CM/CF Bias
# Make a new column with just the cell_line (ie. 18499_11 convereted to 18499) 
bulk_terminal <- bulk_terminal %>% 
  separate(Mixture, c("Line", "Day"))

# Aggregate by cell line
bulk_terminal_aggregated = aggregate(bulk_terminal,
                                     by = list(bulk_terminal$Line),
                                     FUN = mean)

# Determine if CM or CF biased
bulk_terminal_aggregated <- bulk_terminal_aggregated %>% 
  mutate(bias = ifelse(CM>CF, 'CM', 'CF'))

## Pseudobulk CM/CF Bias
# Make a new column with just the cell_line (ie. 18499_11 convereted to 18499) 
pseudobulk_terminal <- pseudobulk_terminal %>% 
  separate(sample, c("Line", "Day"))

# Aggregate by cell line
pseduobulk_terminal_aggregated = aggregate(pseudobulk_terminal,
                                           by = list(pseudobulk_terminal$Line),
                                           FUN = mean)

# Determine if CM or CF biased
pseduobulk_terminal_aggregated <- pseduobulk_terminal_aggregated %>% 
  mutate(bias = ifelse(CM>CF, 'CM', 'CF'))

## Get list of CM-biased bulk and pseudobulk lines
bulk_terminal_aggregated[bulk_terminal_aggregated$bias == 'CM', 'Group.1']
pseduobulk_terminal_aggregated[pseduobulk_terminal_aggregated$bias == 'CM', 'Group.1']


