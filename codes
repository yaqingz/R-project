---
  title: Study of clearance and genetic mutation
author: Yaqing Zhang
email: yaqing.zhang@outlook.com
date: '`r format(Sys.time(), "%d %B, %Y")`'
mainfont: 
  sansfont: 
  mathfont: 
  monofont: 
  fontsize: 
  output: 
  html_document: 
  code_folding: show
number_sections: yes
toc: yes
toc_depth: 3
toc_float:
  collapsed: false
smooth_scroll: true
---
# Environment setup & Load the packages
```{r}
library(tidyr)    
library(knitr)
library(dplyr)
library(ggplot2)  
```

# Data management
## Load the Data Sets
```{r PK data}
# Define PK data path and import PK data
dir_pk <-'PK_data/'
file_pk <-'data_200mg.csv'
input_path_pk <- paste(dir_pk, file_pk, sep = "")
data_pk <- read.csv(file = input_path_pk, header = T, as.is = T)

# Define genetic data path and import genetic data
dir_gn <-'SNP_data/'
file_gn <-'genetics_data.csv'
input_path_gn <- paste(dir_gn, file_gn, sep = "")
data_gn <- read.csv(file = input_path_gn, header = T, as.is = T)
```

## Input Check and table shaping
```{r, eval = TRUE}

str(data_pk) # Check the PK data
knitr::kable(x       = data_pk,
             align   = 'c',
             caption = "data of PK")

str(data_gn) # Check the genetic data
knitr::kable(x       = data_gn,
             align   = 'c',
             caption = "data of gene")
```

## Data set cleaning
```{r}
# Change the column names
colnames(data_pk)[1]<-"ID"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[2]<-"TIME_0.5h"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[3]<-"TIME_1h"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[4]<-"TIME_2h"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[5]<-"TIME_3h"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[6]<-"TIME_6h"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[7]<-"TIME_12h"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[8]<-"TIME_24h"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[10]<-"WT"
r<-paste(colnames(data_pk), collapse = ', ')
colnames(data_pk)[11]<-"HT"
r<-paste(colnames(data_pk), collapse = ', ')

# Change column classes
data_pk$TIME_24h <-as.numeric(data_pk$TIME_24h)
data_pk$TIME_0.5h <-as.numeric(data_pk$TIME_0.5h)
class(data_pk$TIME_24h)
class(data_pk$TIME_0.5h)

# Define 'Sex' factor
data_pk$Sex <- as.factor(data_pk$Sex)
f <- factor(data_pk$Sex, levels=c("M","F"), labels=c("male","female"))
class(data_pk$Sex)

# The code below will be used in the calculation of AUC
AUC_pk <- data_pk
```
## Reshape data sets
```{r}
# Gather the data
data_pk <- gather(data_pk, key = TIME,value = Conc, -ID, -Sex, -HT, -WT)
data_pk <- data_pk[order(data_pk$ID), ]

# Gather the genetic data and define mutation types
data_gn <- gather(data_gn, key = MUTATION,value = GENE, -ID)
data_gn <- data_gn[order(data_gn$ID), ]
data_gn$GENE <- as.character(data_gn$GENE)
data_gn$GENOTYPE <- ifelse(data_gn$GENE == '0', 'No mutation', NA)
data_gn$GENOTYPE <- ifelse(data_gn$GENE == '1', 'Heterozygous', data_gn$GENOTYPE)
data_gn$GENOTYPE <- ifelse(data_gn$GENE == '2', 'Homozygous', data_gn$GENOTYPE)
data_gn <- mutate(data_gn, GENOTYPE)
data_gn = data_gn[,-3]
```

## Merge data sets
```{r}
data_all <- merge(x   = data_pk,
                  y   = data_gn,
                  all = TRUE,
                  by  = 'ID')

str(data_all)

knitr::kable(x       = head(data_all, n=5),
             align   = 'c',
             caption = 'All Data')
data_all$TIME <- gsub('TIME_', '', data_all$TIME)
data_all$TIME <- gsub('h', '', data_all$TIME)
data_all$TIME <- as.numeric(data_all$TIME)
```

# Data summary
```{r}
# Time-changing continuous
summary(data_all[, c('TIME', 'Conc')])

# Non-time changing continuous
summary(data_all[, c('Sex', 'HT', 'WT')])
```

# Process calculation
```{r}
# Delete redundant data
data_unique <- data_all[!duplicated(data_all$ID),]
```

## Drug concentration-Time Profiles
```{r}
# Plotting a graph of drug concentation and time
p0 <- ggplot() + geom_line(aes(y = Conc, TIME, colour = ID),
                           data = data_all, stat="identity")
p0
```

## Calculation of AUC
```{r}
# The calculation of area under curve(AUC) will be based on time-concentration curve. AUC will be divided to several ladders, and the estimated AUC is the addtion results of these ladders
AUC_pk$TIME_24h <- as.numeric(AUC_pk$TIME_24h)

# In the calculation, the missing concentration data (listed as NA) in 24h will be replace by an average concentration of other patients
AUC_pk[is.na(AUC_pk)] <- mean(AUC_pk$TIME_24h,na.rm=TRUE)

# The ladder area after 24h will be assumed to be the triangle area comstituted by the 24h concentration and the extension of concentration decrease line 
# Equation:y=kx+b
k <- (AUC_pk$TIME_24h-AUC_pk$TIME_12h)/12
b <- AUC_pk$TIME_24h-k*24

# The time when concentration decrease to 0 will be calculated by the equeation below
t <- (-b)/k

# Define the function of calculation AUC
time1<-c(2:7)
time2<-c(3:8)
for (i in AUC_pk){
  area1<-(AUC_pk[,time1]+AUC_pk[,time2])
}
area1
area2<-matrix(c(0.5,1,1,3,6,12),byrow=TRUE,nrow=15,ncol=6)
area2
area3<-area1*area2*0.5
area3
area4<-rowSums (area3)+0.5*(t-24)*AUC_pk$TIME_24h
area4

# Output AUC results
AUC_pk$AUC<-area4

# Merge data sets
AUC_pk = AUC_pk[,-(2:11)]
data_all <- merge(x   = AUC_pk,
                  y   = data_all,
                  all = TRUE,
                  by  = 'ID')
```

##Calculation of clearance(CL) and weight	adjusted	clearance(CLw)
```{r}
# CL=DOSE/AUC
data_all$CL <- 200/data_all$AUC

# CLw=CL*(WT/70)^0.75
data_all$CLw <- data_all$CL*(data_all$WT/70)^0.75
```

##Compare the standard deviation of CL and CLw
```{r}
# Calculate the standard deviation and standard error
sd_CL <- sd(data_all$CL)
sd_CLw <- sd(data_all$CLw)
se_CL <- sd_CL/sqrt(sd_CL)
se_CLw <- sd_CLw/sqrt(sd_CLw)
# Compare
cat('SD of CL:', sd_CL, '\n')
cat('SD of CLw:', sd_CLw, '\n')
cat('SE of CL:', se_CL, '\n')
cat('SE of CLw:', se_CLw, '\n')
```

# Covariate testing of PK parameters
##The classfication of the genetic	variations
```{r}
# In previous chunk, mutation types are classified as Heterozygous, Homozygous and No mutatioin. 
# Data of heterozygous is collected below, it will be used in comparison to no mutation samples
data_hete <- filter(data_all, GENOTYPE !='Homozygous')
data_hete$MUTATION <- as.factor(data_hete$MUTATION)

# Data of homozygous is collected below
data_homo <- filter(data_all, GENOTYPE !='Heterozygous')
data_homo$MUTATION <- as.factor(data_homo$MUTATION)

# Data of all mutation samples is collected below
data_muta <- filter(data_all, GENOTYPE !='No mutation')
data_muta$MUTATION <- as.factor(data_muta$MUTATION)
```

##T-test of the CLw by different genotype
### T-test of heterozygous genotype
```{r Ttest of heterozygous data}
# t-test of heterozygous data
ttest_hete <- data_hete %>%
  filter(GENOTYPE != 'No mutation') %>% 
  group_by(MUTATION) %>% 
  do(p_value = t.test(x = data_hete$CLw[data_hete$GENOTYPE == 'No mutation'],
                      y = .$CLw,
                      paired = FALSE)$p.value) %>% 
  transmute(MUTATION        = MUTATION,
            p_value     = unlist(p_value),
            Significant = ifelse(p_value <= 0.05, 'yes', 'no'),
            p_value     = as.character(signif(p_value, digits = 3)))

ttest_hete %>% 
  kable(caption   = 't-tests between heterozygous and no mutation', 
        row.names = FALSE,
        align     = 'c')

```

### Plot of the data with significant different
```{r}
# A boxplot will be created here to compare the clearance of samples with G105A heterozygous mutation and no mutation samples
hete_G1051A <- subset(data_hete, data_hete$MUTATION == "G1051A")
hete_G1051A$MUTATION <- as.factor(hete_G1051A$MUTATION)

p1 <- ggplot(data =hete_G1051A,aes(x = GENOTYPE, y = CL,fill=GENOTYPE))+geom_boxplot(outlier.shape=1,color = "darkred",alpha = 0.5,width=0.25)+ggtitle("G1051A - heterozygous and no mutation")
p1
```

# T-test of homozygous genotype
```{r}
# T-test of the homozygous data
ttest_homo <- data_homo %>%
  filter(GENOTYPE != 'No mutation') %>% 
  group_by(MUTATION) %>% 
  do(p_value = t.test(x = data_homo$CLw[data_homo$GENOTYPE == 'No mutation'],
                      y = .$CLw,
                      paired = FALSE)$p.value) %>% 
  transmute(MUTATION        = MUTATION,
            p_value     = unlist(p_value),
            Significant = ifelse(p_value <= 0.05, 'yes', 'no'),
            p_value     = as.character(signif(p_value, digits = 3)))

ttest_homo %>% 
  kable(caption   = 't-tests between homozygous and no mutation', 
        row.names = FALSE,
        align     = 'c')

```

### Plot of the data with significant different
```{r}
# Boxplots will be created here to compare the clearance of samples with T134A and G1051A homozygous mutation and no mutation samples
homo_T134A <- subset(data_homo, data_homo$MUTATION == "T134A")
homo_T134A$MUTATION <- as.factor(homo_T134A$MUTATION)

p2 <- ggplot(data =homo_T134A,aes(x = GENOTYPE, y = CL,fill=GENOTYPE))+geom_boxplot(outlier.shape=1,color = "darkred",alpha = 0.5,width=0.25)+ ggtitle("T134A - homozygous and no mutation") 
p2

# G1051A
homo_G1051A <- subset(data_homo, data_homo$MUTATION == "G1051A")
homo_G1051A$MUTATION <- as.factor(homo_G1051A$MUTATION)

p3 <- ggplot(data =homo_G1051A,aes(x = GENOTYPE, y = CL,fill=GENOTYPE))+geom_boxplot(outlier.shape=1,color = "darkred",alpha = 0.5,width=0.25) + ggtitle("G1051A - homozygous and no mutation")
p3
```

# T-test of both two types of mutation
```{r Ttest of mutation data}
#T-test of the mutation data
ttest_muta <- data_muta %>%
  filter(GENOTYPE != 'Homozygous') %>% 
  group_by(MUTATION) %>% 
  do(p_value = t.test(x = data_muta$CLw[data_muta$GENOTYPE == 'Homozygous'],
                      y = .$CLw,
                      paired = FALSE)$p.value) %>% 
  transmute(MUTATION        = MUTATION,
            p_value     = unlist(p_value),
            Significant = ifelse(p_value <= 0.05, 'yes', 'no'),
            p_value     = as.character(signif(p_value, digits = 3)))

ttest_muta %>% 
  kable(caption   = 't-tests between heterozygous and homozygous', 
        row.names = FALSE,
        align     = 'c')

```
