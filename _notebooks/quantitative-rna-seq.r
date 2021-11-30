library('edgeR')
library('readr')
library('magrittr')





##### Using edgeR from an ExpressionSet object #####

# STEP 1 LOAD THE COUNT DATA

load(file.path(getwd(), "rdata/modencodefly_eset.RData"))

# STEP 2 SPECIFY EXPERIMENTS OF INTEREST

experiments_of_interest <- c("L1Larvae", "L2Larvae")

columns_of_interest <- which(phenoData(modencodefly.eset)[['stage']] %in% experiments_of_interest )

# STEP 3 FORM THE GROUPING FACTOR

grouping <- droplevels(phenoData(modencodefly.eset)[['stage']][columns_of_interest] )

# STEP 4 FORM THE SUBSET OF COUNT DATA

counts_of_interest <- exprs(modencodefly.eset)[,columns_of_interest]

# STEP 5 CREATE THE DGE OBJECT

eset_dge <- edgeR::DGEList(counts = counts_of_interest, group = grouping)

# STEP 6 PERFORM DIFFERENTIAL EXRESSION ANALYSIS 

design <- model.matrix(~ grouping)
eset_dge <- edgeR::estimateDisp(eset_dge, design)
fit <- edgeR::glmQLFit(eset_dge, design)
result <- edgeR::glmQLFTest(fit, coef=2)
topTags(result)


