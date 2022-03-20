# BiocManager::install("RTCGA.clinical")
# BiocManager::install("RTCGA.mRNA")
library(RTCGA.clinical)
library(RTCGA.mRNA)
library(tibble)

# use data from RTCGA.clinical package
dim(BRCA.clinical)
names(BRCA.clinical)
# use the helper function from RTCGA.clinical package to extract the column of interest
clin <- survivalTCGA(BRCA.clinical, OV.clinical, GBM.clinical, 
                     extract.cols="admin.disease_code")
# preview the dataset
head(clin)
# count of each disease
table(clin$admin.disease_code)
# contingency table
xtabs(~admin.disease_code+patient.vital_status, data=clin)

# fit the survival model of disease
coxph(Surv(times, patient.vital_status)~admin.disease_code, data=clin)
sfit <- survfit(Surv(times, patient.vital_status)~admin.disease_code, data=clin)
summary(sfit, times=seq(0,365*5,365))
# plot
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE)


# example of looking at multiple features
clinkid <- survivalTCGA(KIPAN.clinical, 
                     extract.cols=c("admin.disease_code", "patient.gender"))
head(clinkid)
xtabs(~admin.disease_code+patient.gender, data=clinkid)
coxph(Surv(times, patient.vital_status)~admin.disease_code+patient.gender, data=clinkid)

sfit <- survfit(Surv(times, patient.vital_status)~admin.disease_code+patient.gender, data=clinkid)
summary(sfit, times=seq(0,365*5,365))
# plot
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE)

# kmeans
data <- rbind(BRCA.mRNA, OV.mRNA)[, 2:17815]
std <- apply(data, 2, sd, na.rm =TRUE)
top200 <- order(std, decreasing = TRUE)[1:200]
gene_names <- colnames(data)
selected_genes <- gene_names[top200]
expr <- expressionsTCGA(BRCA.mRNA, OV.mRNA,
                        extract.cols = selected_genes)
expr <- na.omit(expr)
km.res <- kmeans(expr[3:202], centers = 2)
str(km.res)
table(true=expr$dataset,cluster=km.res$cluster)


