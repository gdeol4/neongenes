# Introduction

### What is smoc2?

### A look at the dataset

Introduction to dataset: Smoc2 Smoc2, or Secreted modular
calcium-binding protein 2, has been shown to have increased expression
in kidney fibrosis, which is characterized by an excess of extracellular
matrix in the space between tubules and capillaries within the kidney.
However, it is unknown how Smoc2 functions in the induction and
progression of fibrosis.

Introduction to dataset: Experimental design There are four sample
groups being tested: normal, control mice, referred to as wild type
mice, with and without fibrosis and Smoc2 over-expressing mice with and
without fibrosis. There are three biological replicates for all normal
samples and four replicates for all fibrosis samples. Initially, we will
explore the effect of fibrosis on gene expression using ‘Wild type’
samples during lectures and ‘Smoc2 over-expression’ data during
exercises.

### Tools used

Differential expression analysis: tools While there are a large number
of statistical packages developed for DE analysis, DESeq2 and EdgeR are
two of the most popular tools. Both tools use the negative binomial
model to model the raw count data, employ similar methods, and
typically, yield similar results. Both are pretty stringent and have a
good balance between sensitivity and specificity. Limma-Voom is another
set of tools often used together for DE analysis, while also a good
method, it can be a bit less sensitive for smaller sample sizes. We will
be using DESeq2 for the DE analysis, which has an extensive vignette
available to help with questions.

``` r
#Load library for DESeq2
library(DESeq2)

#Load library for RColorBrewer
library(RColorBrewer)

#Load library for pheatmap
library(pheatmap)

#Load library for tidyverse
library(tidyverse)
```

### Overview of workflow

Differential expression analysis: DESeq2 workflow We will be using the
DESeq2 tool to perform the differential expression analysis, but what
are the steps we will need to perform? Displayed in the workflow are the
steps in the differential expression analysis, separated into quality
control and DE analysis steps. To start, we will take the count matrix
containing the raw counts for each sample and perform quality control
steps. First, the counts will be normalized to account for differences
in library depth. Then, principal component analysis and hierarchical
clustering using heatmaps will be performed to identify potential sample
outliers and major sources of variation in the data. After QC, the DE
analysis will be performed, including the modeling of counts, shrinking
the log2 fold changes, and testing for differential expression. We will
cover each of these steps in more detail as we progress through the
workflow.

# Importing read counts associated with genes

### Importing the data

``` r
data <- read_csv("C:/Users/gurka/Downloads/fibrosis_smoc2_rawcounts/fibrosis_smoc2_rawcounts.csv")
```

    ## New names:
    ## Rows: 47729 Columns: 8
    ## ── Column specification
    ## ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (7): smoc2_fibrosis1, smoc2_fibrosis4, smoc2_normal1, smoc2_normal3, smoc2_fibrosis3, smoc2_normal4, smoc2_fibrosis2
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
head(data)
```

    ## # A tibble: 6 × 8
    ##   ...1               smoc2_fibrosis1 smoc2_fibrosis4 smoc2_normal1 smoc2_normal3 smoc2_fibrosis3 smoc2_normal4 smoc2_fibrosis2
    ##   <chr>                        <dbl>           <dbl>         <dbl>         <dbl>           <dbl>         <dbl>           <dbl>
    ## 1 ENSMUSG00000102693               0               0             0             0               0             0               0
    ## 2 ENSMUSG00000064842               0               0             0             0               0             0               0
    ## 3 ENSMUSG00000051951              72              30             0             3              36             1              51
    ## 4 ENSMUSG00000102851               0               0             0             0               0             0               0
    ## 5 ENSMUSG00000103377               0               0             1             0               0             0               0
    ## 6 ENSMUSG00000104017               0               0             0             0               0             0               0

``` r
data <- subset (data[-c(1)])
```

### Putting together the metadata

Preparation for differential expression analysis: metadata In addition
to our raw counts, we require sample metadata. At the very least, we
need to know which of our samples correspond to each condition. To
generate our metadata, we create a vector for each column and combine
the vectors into a data frame. The sample names are added as the row
names.

Preparation for differential expression analysis: metadata After we have
the counts and metadata files, we can start the differential expression
analysis workflow.

``` r
#create genotype vector
genotype = c('smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe','smoc2_oe')

#create condition vector
condition = c('fibrosis','fibrosis','fibrosis','fibrosis','normal','normal','normal')

#create dataframe

smoc2_metadata = data.frame(genotype, condition)

#Assign the row names of the data frame
rownames(smoc2_metadata) = c('smoc2_fibrosis1','smoc2_fibrosis2','smoc2_fibrosis3','smoc2_fibrosis4', 'smoc2_normal1','smoc2_normal3', 'smoc2_normal4' )
```

# Normalization

### Preparing the data for DESeq2

After we have our data loaded, we need to make sure it’s in a specific
format so that it will be accepted as input to DESeq2.

Bringing in data for DESeq2: sample order DESeq2 requires the sample
names in the metadata and counts datasets to be in the same order.
Therefore, the row names of our metadata need to be in the same order as
the column names of our counts data. As we can see in these images, the
sample names for the counts are in a nice order while the names in the
metadata are not.

Bringing in data for DESeq2: sample order Alternatively, we can explore
the row names and column names with the rownames() and colnames()
functions.

Bringing in data for DESeq2: sample order By looking at our sample names
in both datasets, we can see that the order is not the same, but it’s
not always clear, so using the all() function with the “double equal to”
sign can check if all of the row names of the metadata are in the same
order as the column names of the raw counts data. The all() function
returns a FALSE value. Now we know that our samples are not in the same
order, so we need to reorder the data to use it with DESeq2.

Matching order between vectors To easily reorder the rows of the
metadata to match the order of columns in the counts data, we can use
the match() function. The match() function takes two vectors as input.
The first is a vector of values in the order we want, and the second is
a vector of values we would like to reorder. In our example, the column
names of the raw counts data is the vector with the order we want, so it
will be in the first position of the match() function. The row names of
the metadata is the vector to be reordered, so it will be in the second
position. The output shows how we would need to reorder the rows of the
metadata to be in the same order as the columns in the count data. For
instance the 6th row would need to come first, followed by the 9th row,
then the 1st row, so on and so forth.

Reordering with the match() function Now, we can use the output of the
match() function to reorder the rows of the metadata to be in the same
order as the columns in the count data. To do this we can save the
indices output by match() to a variable, in this case called idx. Then,
we can rearrange the metadata by using the square brackets and adding
idx to the rows position. The samples should now be in the same order
for both datasets.

Checking the order To check the order we can use the all() function
again. Since they match, we can now use these datasets to create the
DESeq2 object needed to start the DESeq2 workflow.

8.  Creating the DESeq2 object To create the DESeq2 object, use the
    DESeqDataSetFromMatrix() function. This function takes as input the
    raw counts, associated metadata, and a design formula detailing
    which conditions in the metadata we want to use for differential
    expression analysis. We will talk in more detail about the design
    formula later. This function will create a DESeq2 object, of the
    class Ranged Summarized Experiment. This is a list-like object with
    slots available for the data it will generate throughout the
    analysis. Currently, it only has a few of the slots filled with the
    count data, metadata, and design information.

``` r
#check if the column orders are the same with the all function
all(rownames(smoc2_metadata) ==  colnames(data))
```

    ## [1] FALSE

``` r
#determine the matching order between the vectors using match()
match(colnames(data), rownames(smoc2_metadata))
```

    ## [1] 1 4 5 6 3 7 2

``` r
#use the match() function to reorder the column of the raw counts
reorder_idx = match(colnames(data), rownames(smoc2_metadata))

#reorder the metadata table
reordered_smoc2_meta = smoc2_metadata[reorder_idx, ]

#view the new table
reordered_smoc2_meta
```

    ##                 genotype condition
    ## smoc2_fibrosis1 smoc2_oe  fibrosis
    ## smoc2_fibrosis4 smoc2_oe  fibrosis
    ## smoc2_normal1   smoc2_oe    normal
    ## smoc2_normal3   smoc2_oe    normal
    ## smoc2_fibrosis3 smoc2_oe  fibrosis
    ## smoc2_normal4   smoc2_oe    normal
    ## smoc2_fibrosis2 smoc2_oe  fibrosis

``` r
#check to see if the column and row names are in the same order
all(rownames(reordered_smoc2_meta) ==  colnames(data))
```

    ## [1] TRUE

Good work, we now have a DESeq2 object storing our raw counts and
metadata that we can use to explore the data with DESeq2 functions and
to use for performing the differential expression analysis.

### Creating the DESeq2 DataSet (DDS)

``` r
# Create a DESeq2 object
dds_smoc2 = DESeqDataSetFromMatrix(countData = data,
                                   colData = reordered_smoc2_meta,
                                   design = ~ condition)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors

### Estimating size factors

Got It! 1. Count normalization Now that we have our DESeq2 object
created with the raw counts and metadata stored inside, we can start the
DESeq2 workflow.

2.  DESeq workflow - normalization The first step in the workflow is to
    normalize the raw counts to assess sample-level quality control
    metrics.

3.  Count normalization But what does it mean to normalize the raw
    counts? The raw counts represent the number of reads aligning to
    each gene and should be proportional to the expression of the RNA in
    the sample; however, there are factors other than RNA expression
    that can influence the number of reads aligning to each gene. We can
    adjust the count data to remove the influence of these factors on
    the overall counts using normalization methods. The main factors
    often considered during normalization of count data are library
    depth, gene length, and RNA composition.

4.  Library depth normalization Differences in library size between
    samples can lead to many more reads being aligned to genes in one
    sample versus another sample. In this example, sample A has nearly
    twice the reads, represented as small rectangles, aligning to each
    of the genes as sample B only because sample A has nearly twice the
    number of reads sequenced. Therefore, we need to adjust the counts
    assigned to each gene based on the size of the library prior to
    doing differential expression analysis.

5.  Gene length normalization Another normalization factor often
    adjusted for is gene length. A longer gene generates a longer
    transcript, which generates more fragments for sequencing.
    Therefore, a longer gene will often have more counts than a shorter
    gene with the same level of expression. In this example, gene X is
    twice as long as gene Y and, due to the difference in length, is
    assigned twice as many reads. Since DE analysis compares expression
    levels of the same genes between conditions, we do not need to
    normalize for gene length. However, if you were to compare the
    expression levels of different genes, you would need to account for
    lengths of the genes.

6.  Library composition effect When adjusting for library size, the
    composition of the library is also important. A few highly
    differentially expressed genes can skew many normalization methods
    that are not resistant to these outliers. In this image, we can see
    that the green DE gene takes up a large proportion of reads for
    Sample A. If we just divided our counts by the total number of
    reads, normalization for the majority of genes would be skewed by
    the highly expressed DE gene. For this reason, when performing a DE
    analysis, we need to use a method that is resistant to these outlier
    genes.

7.  Normalized counts: calculation To calculate the normalized counts
    with DESeq2, we can use the function estimateSizeFactors() on the
    DESeq2 object, dds_wt, and assign the output to a slot in the DESeq2
    object, by re-assigning to dds_wt. DESeq2 will use these size
    factors to normalize the raw counts. The raw counts for each sample
    are divided by the associated sample-specific size factor for
    normalization. To view the size factors used for normalization, we
    can use the sizeFactors() function.

``` r
#estimating the size factors and feeding them back to the dds object by reassigning to the dds object
dds_smoc2 = estimateSizeFactors(dds_smoc2)

#viewing the size factors
sizeFactors(dds_smoc2)
```

    ## smoc2_fibrosis1 smoc2_fibrosis4   smoc2_normal1   smoc2_normal3 smoc2_fibrosis3   smoc2_normal4 smoc2_fibrosis2 
    ##       1.4319832       1.0826799       0.7106482       0.7989734       1.2480024       0.8482426       1.1189642

## Extracting normalized counts

7.  DESeq2 normalization DESeq2 uses a ‘median of ratios’ method of
    normalization. This method adjusts the raw counts for library size
    and is resistant to large numbers of differentially expressed genes.

8.  Normalized counts: extraction Once the size factors have been
    calculated and added to the DESeq2 object, the normalized counts can
    be extracted from it. We can extract the normalized counts from the
    DESeq2 object using the counts() function while specifying that we
    would like the normalized counts with the normalized = TRUE
    argument. If the default was left as normalized = FALSE, then we
    would extract the raw counts from the object.

``` r
#once factors are calculated and reassigned then the normalized counts can be extracted using the counts function with normalized=TRUE argument
normalized_smoc2_counts = counts(dds_smoc2, normalized=TRUE)
```

``` r
head(normalized_smoc2_counts)
```

    ##      smoc2_fibrosis1 smoc2_fibrosis4 smoc2_normal1 smoc2_normal3 smoc2_fibrosis3 smoc2_normal4 smoc2_fibrosis2
    ## [1,]         0.00000         0.00000      0.000000      0.000000          0.0000      0.000000         0.00000
    ## [2,]         0.00000         0.00000      0.000000      0.000000          0.0000      0.000000         0.00000
    ## [3,]        50.27992        27.70902      0.000000      3.754818         28.8461      1.178908        45.57786
    ## [4,]         0.00000         0.00000      0.000000      0.000000          0.0000      0.000000         0.00000
    ## [5,]         0.00000         0.00000      1.407166      0.000000          0.0000      0.000000         0.00000
    ## [6,]         0.00000         0.00000      0.000000      0.000000          0.0000      0.000000         0.00000

### Transform the normalized counts

Congratulations, we now have our normalized counts, which we can use to
accurately compare gene expression between samples. We will be using the
normalized counts to explore similarities in gene expression between
each of our samples, with the expection that our biological replicates
are more similar to each other and the different conditions (wild type
and fibrosis) are more different.

# Unsupervised clustering analysis

1.  Unsupervised clustering analyses Now that we have our normalized
    counts, we can continue on in the differential expression analysis
    workflow.

2.  Unsupervised clustering analyses With our counts normalized for
    library size, we can now compare the counts between the different
    samples. We can explore how similar the samples are to each other
    with regards to gene expression to assess the quality of our
    experiment. To do this we use visualization methods for unsupervised
    clustering analyses, including hierarchical clustering heatmaps and
    principal component analysis or PCA. We perform these QC methods to
    get an idea of how similar the biological replicates are to each
    other and to identify outlier samples and major sources of variation
    present in the dataset.

### log transformations

3.  Unsupervised clustering analyses: log transformation When using
    these visualization methods, we should first log transform the
    normalized counts to improve the visualization of the clustering.
    For RNA-Seq data, DESeq2 uses a variance stabilizing transformation
    (VST), which is a logarithmic transformation that moderates the
    variance across the mean. We can transform the normalized counts by
    using the DESeq2 vst() function on the DESeq2 object. The blind=TRUE
    argument specifies that the transformation should be blind to the
    sample information given in the design formula; this argument should
    be specified when performing quality assessment.

``` r
#Unsupervised clustering analyses: log transformation

#transform the normalized counts
vsd_smoc2 = vst(dds_smoc2, blind=TRUE)
```

### Hierarchical clustering with a correlation heatmap

4.  Hierarchical clustering with correlation heatmaps Hierarchical
    clustering with heatmaps is used to assess the similarity in gene
    expression between the different samples in a dataset. This
    technique is used to explore how similar replicates are to each
    other and whether the samples belonging to different sample groups
    cluster separately. The heatmap is created by using the gene
    expression correlation values for all pairwise combinations of
    samples in the dataset, with the value 1 being perfect correlation.
    The hierarchical tree shows which samples are more similar to each
    other and the colors in the heatmap depict the correlation values.
    We expect the biological replicates to cluster together and sample
    conditions to cluster apart. Since the majority of genes should not
    be differentially expressed, samples should generally have high
    correlations with each other. Samples with correlation values below
    0-point-8 may require further investigation to determine whether
    these samples are outliers or have contamination.

5.  Hierarchical clustering with correlation heatmaps To create a
    hierarchical heatmap, we are going to extract the VST-transformed
    normalized counts as a matrix from the vsd object using the assay()
    function. Then, we can compute the pairwise correlation values
    between each pair of samples using the cor() function. Using View()
    we can view the correlation values between each of the sample pairs.

6.  Hierarchical clustering with correlation heatmaps After generating
    the correlation values, we can use the pheatmap package to create
    the heatmap. The annotation argument selects which factors in the
    metadata to include as annotation bars. We use the select() function
    from the dplyr package to select the condition column in the
    wildtype metadata.

``` r
#Hierarchical clustering with correlation heatmaps

#extract the vst matrix (of transformed counts) from the object
vsd_mat_smoc2 = assay(vsd_smoc2)

#compute the pairwise correlation values between samples
vsd_cor_smoc2 = cor(vsd_mat_smoc2)

#view correlation statistics
vsd_cor_smoc2
```

    ##                 smoc2_fibrosis1 smoc2_fibrosis4 smoc2_normal1 smoc2_normal3 smoc2_fibrosis3 smoc2_normal4 smoc2_fibrosis2
    ## smoc2_fibrosis1       1.0000000       0.9948457     0.9631782     0.9637385       0.9949621     0.9597660       0.9958622
    ## smoc2_fibrosis4       0.9948457       1.0000000     0.9668704     0.9671808       0.9951846     0.9635529       0.9946819
    ## smoc2_normal1         0.9631782       0.9668704     1.0000000     0.9935479       0.9678895     0.9940744       0.9651063
    ## smoc2_normal3         0.9637385       0.9671808     0.9935479     1.0000000       0.9691965     0.9942610       0.9656578
    ## smoc2_fibrosis3       0.9949621       0.9951846     0.9678895     0.9691965       1.0000000     0.9650676       0.9952111
    ## smoc2_normal4         0.9597660       0.9635529     0.9940744     0.9942610       0.9650676     1.0000000       0.9618977
    ## smoc2_fibrosis2       0.9958622       0.9946819     0.9651063     0.9656578       0.9952111     0.9618977       1.0000000

7.  Hierarchical clustering with correlation heatmaps The output from
    the heatmap shows that the biological replicates cluster together
    and the conditions separate. This is encouraging since our
    differentially expressed genes between our conditions are likely to
    be driving this separation. Also, all correlation values are
    expectedly high without any outlier samples. If our replicates did
    not cluster as expected, we could plot the heatmap with all of the
    metadata and see whether any other factor corresponds to the
    separation of the samples. If so, you would want to see if you get
    similar results with the Principal Component Analysis, which we will
    be covering in the next lesson. If identified by both methods, we
    can account for it in the DESeq2 model. If you have an outlier
    identified, you would also want to check the PCA. If you see the
    outlier with both methods, you could investigate that sample more
    and decide whether to remove it from the analysis.

``` r
#plot the heatmap
pheatmap(vsd_cor_smoc2, annotation = select(smoc2_metadata, condition))
```

![](pre-post_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

### Principal component analysis

To continue with the quality assessment of our samples, in the first
part of this exercise, we will perform PCA to look how our samples
cluster and whether our condition of interest corresponds with the
principal components explaining the most variation in the data.

When performing quality assessment of our count data, we need to
transform the normalized counts for better visualization of the variance
for unsupervised clustering analyses. To assess the similarity of the
smoc2 samples using hierarchical heatmaps, transform the normalized
counts and perform hierarchical clustering analysis. Assume all
libraries have been loaded, the DESeq2 object created, and the size
factors have been stored in the DESeq2 object, dds_smoc2.

Well done, now let’s assess the heatmap to determine whether our samples
cluster as expected and if we have any outliers. If we had additional
metadata for possible sources of variation, then we would want to
include these factors as annotation bars as well.

Got It! 1. Principal component analysis The next type of unsupervised
clustering method is principal component analysis, or PCA.

2.  Principal component analysis (PCA): Overview PCA is a technique used
    to emphasize the variation present in a dataset. PCA finds the
    principal components of a dataset, with the first principal
    component, or PC1, representing the greatest amount of variance in
    the data.

3.  Principal component analysis (PCA): Theory To understand this a bit
    better, we can think of a dataset with two samples. We could plot
    the normalized counts of every gene for one sample on the x-axis and
    the other sample on the y-axis. In this example, Gene A has four
    counts for sample 1 plotted on the x-axis and 5 counts for sample 2
    on the y-axis. We can plot the other genes similarly.

4.  Principal component analysis (PCA): Theory We can draw a line
    through the dataset where there exists the most variation, or where
    there is the largest spread. In this example, the line with largest
    spread is between genes B and C. This line represents the first
    principal component. The second most variation in the dataset,
    represented as PC2, must be perpendicular to PC1, in order to best
    describe the variance in the dataset not included in PC1. In this
    example, PC2 is drawn between genes A and D. The spread is much
    smaller for PC2. In reality, your dataset will have more samples and
    many more genes. The number of principal components is equal to the
    number of samples, n, in the dataset, so finding the largest amount
    of variation, PC1, means plotting a line through n-dimensional
    space.

5.  Principal component analysis (PCA): Theory The most variant genes
    for a principal component have the most influence on that principal
    component’s direction. In our example, the most variant genes for
    PC1, genes B and C, would affect the direction of the line more than
    genes A and D.

6.  Principal component analysis (PCA): Theory We give quantitative
    scores to genes based on how much they influence the different PCs.

7.  Principal component analysis (PCA): Theory A ‘per sample’ PC value
    is computed by taking the product of the influence and the
    normalized read count for each gene and summing across all genes.

8.  Principal component analysis (PCA): Theory For PCA we generally plot
    these per sample PC values. Samples that cluster together have more
    similar gene expression profiles than samples that cluster apart,
    especially for the most variant genes.

9.  Principal component analysis (PCA): Theory This is a good method to
    explore the quality of the data as we hope to see replicates cluster
    together and conditions to separate on PC1. Sample outliers and
    major sources of variation can also be identified with this method.

10. Principal component analysis (PCA): Theory We can perform PCA using
    DESeq2’s plotPCA() function to plot the first two PCs. This function
    takes as input the transformed vsd object, and we can use the
    intgroup argument to specify what factor in the metadata to use to
    color the plot. We can see that the sample groups, normal and
    fibrosis, separate well on PC1. This means that our condition
    corresponds to PC1, which represents 88% of the variance in our
    data, while 4% is explained by PC2. This is great since it seems
    that a lot of the variation in gene expression in the dataset can
    likely be explained by the differences between sample groups.
    However, if the samples do not separate by PC1 or PC2, the effect of
    the condition could be small or there are sometimes other larger
    sources of variation present. You can color the PCA by other
    factors, such as age, sex, batch, etcetera, to identify other major
    sources of variation that could correspond to one of the large
    principal components. We’ll talk later about how we can account for
    these sources of variation in the model. Just to note, if you would
    like to explore PCs other than PC1 or PC2, the prcomp() base
    function allows a more thorough analysis of PCs.

PCA analysis To continue with the quality assessment of our samples, in
the first part of this exercise, we will perform PCA to look how our
samples cluster and whether our condition of interest corresponds with
the principal components explaining the most variation in the data. In
the second part, we will answer questions about the PCA plot.

To assess the similarity of the smoc2 samples using PCA, we need to
transform the normalized counts then perform the PCA analysis. Assume
all libraries have been loaded, the DESeq2 object created, and the size
factors have been stored in the DESeq2 object, dds_smoc2.

``` r
#plot pca
plotPCA(vsd_smoc2, intgroup="condition")
```

![](pre-post_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Well done! If there is an outlier on the heatmap, then we would want to
see it with PCA as well. We similarly hope our biological replicates
cluster together and conditions separate by PC1 and/or PC2. If we don’t
see this, there may be sources of variation present in our data and if
these sources are present in our metadata, then we can explore these
sources of variation by coloring the PCA by these factors.
Alternatively, we might not see sample groups separate if the condition
of interest does not cause a big change in gene expression.

### Running the DE analysis

Got It! 1. Overview of the DE analysis Now that we have explored the
quality of our samples, removed any outlier samples, and assessed the
major sources of variation present, we can start the differential
expression analysis with DESeq2.

2.  Review the dataset/question Remember, by performing the differential
    expression analysis with DESeq2, we are trying to identify which
    genes have significant differences in expression between the normal
    and fibrosis sample groups.

3.  Overview of the DE analysis The differential expression analysis
    with DESeq2 consists of roughly three steps: fitting the raw counts
    for each gene to the DESeq2 negative binomial model and testing for
    differential expression, shrinking the log2 fold changes, and
    extracting and visualizing the results.

4.  DESeq2 workflow: Model We will start by fitting the raw counts to
    the DESeq2 model. To do this, DESeq2 will first need to estimate the
    size factors, if they haven’t already been estimated, and estimate
    the variation in expression across replicates for each gene. After
    these calculations, the data can be fit to the negative binomial
    model.

5.  DESeq2 workflow: Model To perform these calculations and fit the
    negative binomial model requires only two functions. The first
    function creates the DESeq2 object, which is the same function we
    used previously during count normalization. We wouldn’t need to
    re-create the object unless we removed samples or found additional
    sources of variation during QC using PCA and the correlation
    heatmap. When we created the DESeq2 object, we provided the raw
    counts, metadata, and design formula. We have explored the raw
    counts and metadata, but what exactly was specified with the design
    formula? The design formula is an important part of modeling,
    telling DESeq2 the known major sources of variation to control for,
    or regress out, as well as, the condition of interest to use for
    differential expression testing.

6.  DESeq2 workflow: Design formula For example, if this is your
    metadata and you know that strain, sex, and treatment are major
    sources of variation in the data as determined by the PCA and
    heatmap, then all of these factors should be included in the design
    formula. If my condition of interest is treatment, then it would
    come last in the formula with the other factors preceding it in any
    order. Therefore, the design formula would be: tilde, strain plus
    sex plus treatment. The tilde tells DESeq2 to model the counts using
    the following formula, so it should always proceed your factors.
    Also, the factor names in the design formula need to exactly match
    the column names in the metadata.

7.  DESeq2 workflow: Design formula DESeq2 also allows for complex
    designs. For instance, using the same metadata, if we wanted to know
    the effect of sex on the effect of treatment, we could use an
    interaction term. In this case, we could regress out the variation
    due to strain, sex and treatment and test for genes that
    significantly differ in their treatment effect due to sex using the
    interaction term, sex colon treatment, as the last term in the
    formula. For more information about specifying complex designs, I
    recommend reading through the DESeq2 vignette or the Bioconductor
    support site.

8.  DESeq2 workflow: Running Once you have the DESeq2 object with the
    raw counts, metadata, and the appropriate design formula, you can
    perform model fitting with the function DESeq(). As it runs it will
    output the completed steps as it fills in the different slots in the
    DESeq2 object. The final DESeq2 object will contain all of the
    information needed for performing the differential expression
    testing between specific sample groups.

``` r
# Run DESeq2 analysis 
dds_smoc2 <- DESeq(dds_smoc2)
```

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

### Dispersion

After fitting the model in the previous exercise, let’s explore the fit
of our smoc2 data to the negative binomial model by plotting the
dispersion estimates using the plotDispEsts() function. Remember that
the dispersion estimates are used to model the raw counts; if the
dispersions don’t follow the assumptions made by DESeq2, then the
variation in the data could be poorly estimated and the DE results could
be less accurate.

The assumptions DESeq2 makes are that the dispersions should generally
decrease with increasing mean and that they should more or less follow
the fitted line

Got It! 1. DESeq2 model Now that we have run the DE analysis, we could
explore our results. However, before proceeding, we should explore how
well our data fit the model.

2.  DESeq2 model The goal of the differential expression analysis is to
    determine whether a gene’s mean expression between sample groups is
    different given the variation within groups. This is determined by
    testing the probability of the log2 fold changes between groups
    being significantly different from zero. The log2 fold changes are
    found by the log of the one sample group mean, shown here as the
    treatment group, divided by the mean of the other sample group,
    shown here as the control group. Therefore, to model the counts
    requires information about the mean and variation in the data. To
    explore the variation in our data, we will observe the variance in
    gene expression relative to the mean. Variance is the square of the
    standard deviation, representing how far away the expression of the
    individual samples, as shown by the dark red and blue circles, are
    from the means, shown in pink and light blue.

DESeq2 model - mean-variance relationship For RNA-Seq data, the variance
is generally expected to increase with the gene’s mean expression. To
observe this relationship, we can calculate the means and variances for
every gene of the normal samples using the apply() function.

4.  DESeq2 model - dispersion Then we can create a data frame for
    plotting with ggplot2. We plot the mean and variance values for each
    gene using log10 scales. Each black dot represents a gene.

5.  DESeq2 model - dispersion We see the variance in gene expression
    increases with the mean. This is expected for RNA-Seq data. Also,
    note how the range in values for variance is greater for lower mean
    counts than higher mean counts. This is also expected for RNA-Seq
    count data. A measure of the variance for a given mean is described
    by a metric called dispersion in the DESeq2 model. The DESeq2 model
    uses dispersion to assess the variability in expression when
    modeling the counts.

6.  DESeq2 model - dispersion The DESeq2 model calculates dispersion as
    being indirectly related to the mean and directly related to the
    variance of the data using the formula displayed. So, an increase in
    variance will increase dispersion, while an increase in mean will
    decrease dispersion. For any two genes with the same mean
    expression, the only difference in dispersion will be based on
    differences in variance. To check the fit of our data to the DESeq2
    model, it can be useful to look at the dispersion estimates.

7.  DESeq2 model - dispersion To plot the dispersions relative to the
    means for each gene, we can use the plotDispEsts() function on the
    DESeq2 object. Each black dot is a gene with associated mean and
    dispersion values. We expect dispersion values to decrease with
    increasing mean, which is what we see. With only a few replicates
    for RNA-Seq experiments, gene-wise estimates of dispersion are often
    inaccurate, so DESeq2 uses information across all genes to determine
    the most likely estimates of dispersion for a given mean expression
    value, shown with the red line in the figure. Genes with
    inaccurately small estimates of variation could yield many false
    positives, or genes that are identified as DE, when they are really
    not. Therefore, the original gene-wise dispersion estimates, shown
    as the black dots in the figure, are shrunken towards the curve to
    yield more accurate estimates of dispersion, shown as blue dots. The
    more accurate, shrunken dispersion estimates are used to model the
    counts for determining the differentially-expressed genes. Extremely
    high dispersion values, shown surrounded by blue circles, are not
    shrunken, due to the likelihood that the gene may have higher
    variability than others for biological or technical reasons and
    reducing the variation could result in false positives. The strength
    of the shrinkage is dependent on the distance from the curve and
    sample size. Larger numbers of replicates can estimate the mean and
    variation more accurately, so yield less shrinkage.

8.  DESeq2 model - dispersion Worrisome plots would include a cloud of
    data that doesn’t follow the curve or dispersions that don’t
    decrease with increasing mean. These problems can often be explained
    by sample outliers or contamination. Examples of worrisome
    dispersion plots are shown in the figures.

``` r
# Calculating mean for each gene (each row)
mean_counts = apply(data[, 1:3], 1, mean)
head(mean_counts)
```

    ## [1]  0.0000000  0.0000000 34.0000000  0.0000000  0.3333333  0.0000000

``` r
# Calculating variance for each gene( each row)
variance_counts = apply(data[, 1:3], 1, var)
head(variance_counts)
```

    ## [1]    0.0000000    0.0000000 1308.0000000    0.0000000    0.3333333    0.0000000

``` r
#Plotting relationship between mean and variance:

# Creating data frame with mean and vairance for every gene
df = data.frame(mean_counts, variance_counts)
head(df)
```

    ##   mean_counts variance_counts
    ## 1   0.0000000       0.0000000
    ## 2   0.0000000       0.0000000
    ## 3  34.0000000    1308.0000000
    ## 4   0.0000000       0.0000000
    ## 5   0.3333333       0.3333333
    ## 6   0.0000000       0.0000000

### Relationship between mean and variance

``` r
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene")
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous x-axis

![](pre-post_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### Plotting dispersion estimates

``` r
#plot dispersion estimates
plotDispEsts(dds_smoc2)
```

![](pre-post_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### extracting the results of DE analysis

DESeq2 model - exploring dispersions NOTE: It may take a bit longer to
load this exercise.

After fitting the model in the previous exercise, let’s explore the fit
of our smoc2 data to the negative binomial model by plotting the
dispersion estimates using the plotDispEsts() function. Remember that
the dispersion estimates are used to model the raw counts; if the
dispersions don’t follow the assumptions made by DESeq2, then the
variation in the data could be poorly estimated and the DE results could
be less accurate.

The assumptions DESeq2 makes are that the dispersions should generally
decrease with increasing mean and that they should more or less follow
the fitted line. Nice work! This is always a good plot to check because
it gives us an idea about the fit of our data to the model. Now, let’s
observe the plot and answer questions about whether these dispersion
estimates meet expectations for RNA-seq data.

Got It! 1. DESeq2 model - contrasts Now that we have explored the fit of
our data to the model, we can extract the DE testing results.

2.  DESEq2 workflow We can also make more accurate estimates of the
    foldchanges, which represent the expression of one samplegroup
    relative to another.

3.  DESeq2 workflow When we ran the DESeq() function, the size factors
    were calculated to generate the normalized counts, the shrunken
    dispersions were estimated prior to model fitting, then,
    differential expression testing was performed. However, before we
    extract our results, let’s explore the model a bit more to
    understand the results.

4.  DESeq2 Negative Binomial Model We know that RNA-Seq data can be
    represented well using the negative binomial model, as it accounts
    for the additional variation in the data added by the small number
    of biological replicates. The size factors, indicated by Sij,
    normalized counts, indicated by Qij, and shrunken dispersions,
    denoted by alpha-i, are used as input to the negative binomial model
    to fit the raw count data.

5.  DESeq2 Negative Binomial Model For each gene, the model uses the
    log2 normalized counts, denoted on the left side of the equation to
    determine the log2 foldchange estimates, indicated by the beta-ir,
    for the samples of the condition of interest, represented by xjr,
    for each gene. In addition, to determining the log2 foldchanges, the
    associated standard error is also output.

``` r
# The syntax for DESeq2 contrasts is

# results(dds,
#         contrast = c("condition_factor", "level_to_compare", "base_level"),
#         alpha = 0.05)

smoc2_results = results(dds_smoc2,
                        contrast = c("condition", "fibrosis", "normal"),
                        alpha = 0.05)
```

``` r
head(smoc2_results)
```

    ## log2 fold change (MLE): condition fibrosis vs normal 
    ## Wald test p-value: condition fibrosis vs normal 
    ## DataFrame with 6 rows and 6 columns
    ##    baseMean log2FoldChange     lfcSE      stat      pvalue        padj
    ##   <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
    ## 1  0.000000             NA        NA        NA          NA          NA
    ## 2  0.000000             NA        NA        NA          NA          NA
    ## 3 22.478090        4.49814  0.829291  5.424085 5.82520e-08 2.54201e-07
    ## 4  0.000000             NA        NA        NA          NA          NA
    ## 5  0.201024       -1.59170  3.816946 -0.417009 6.76672e-01          NA
    ## 6  0.000000             NA        NA        NA          NA          NA

### MA plot

``` r
plotMA(smoc2_results, ylim=c(-8,8))
```

![](pre-post_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

### LFC shrinkage

6.  DESeq2 contrasts By default, DESeq2 will perform the Wald test for
    pairwise comparisons to test for differences in expression between
    two sample groups for the condition of interest, in this case,
    condition. The sample groups for condition are fibrosis and normal.
    The results of the testing can be extracted using the results()
    function and specifying a significance level, or alpha value, using
    the alpha argument. You can choose the alpha based on how stringent
    you want to be with your analysis. Lower alpha values indicate less
    probability of identifying a gene as DE when it is actually not. We
    will use a standard alpha of 0-point-05. The top of the output shows
    “log2 foldchange condition normal versus fibrosis”, indicating that
    fibrosis is the base level of comparison. This means that all log2
    fold changes represent the normal group relative to the fibrosis
    group, which doesn’t seem the most intuitive. Instead of using the
    default results we can perform any pairwise comparison between the
    sample groups by supplying our own contrast.

7.  DESeq2 contrasts Contrasts, which specify the sample groups to
    compare, can be given directly to the results() function using the
    contrast argument. Within the combine function, we need to specify
    the condition of interest, level to compare, and base level. To
    specify the normal sample group as the base level for condition, we
    could create the contrast defining normal as the base level and
    fibrosis as the level to compare. Note that the condition of
    interest and sample groups for the condition need to match the names
    in the metadata.

8.  DESeq2 contrasts Now the results give log2 fold changes of the
    fibrosis group relative to normal.

9.  DESeq2 LFC shrinkage To explore our results a bit, the MA plot can
    be helpful. The MA plot shows the mean of the normalized counts
    versus the log2 fold changes for all genes tested. We can use the
    DESeq2 function plotMA() to create the plot and the genes that are
    significantly DE are colored red. Note the large log2 foldchanges,
    particularly for genes with lower mean count values. These fold
    changes are unlikely to be as accurate for genes that have little
    information associated with them, such as genes with low numbers of
    counts or high dispersion values.

10. LFC shrinkage To improve the estimated fold changes we can use log2
    foldchange shrinkage. For genes with low amounts of information
    available, shrinkage uses information from all genes to generate
    more likely, lower, log2 fold change estimates, similar to what we
    did with dispersions. DESeq2 has the lfcShrink() function to
    generate the shrunken log2 foldchanges. We need to specify the
    DESeq2 object, the contrast, and our results object. We can then
    create the MA plot again.

11. LFC shrinkage Now we see more restricted log2 foldchange values,
    especially for lowly expressed genes. These shrunken log2
    foldchanges should be more accurate; however, shrinking the log2
    foldchanges will not affect the number of differentially expressed
    genes returned, only the log2 fold change values. Now that we have
    accurate fold changes, we can extract the significant DE genes and
    perform further visualizations of results.

``` r
smoc2_results = lfcShrink(dds_smoc2,
                          contrast = c("condition", "fibrosis", "normal"),
                          type = 'ashr',
                          res = smoc2_results)
```

    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041

Now that we have extracted our results, we can get a nice overview of
the number of differentially expressed genes there are for our
designated alpha level using the summary() function. It will output the
numbers/percentages of up- and down-regulated genes, as well as, give
information about independent filtering and outliers removed.

### Create a results dataframe

``` r
summary(smoc2_results)
```

    ## 
    ## out of 29556 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 5776, 20%
    ## LFC < 0 (down)     : 5332, 18%
    ## outliers [1]       : 15, 0.051%
    ## low counts [2]     : 7207, 24%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

Well done! Now we implemented a log2 fold change threshold of 1.25 fold
(log2 0.32) when testing for significant genes. Now we can use these
results to subset only the significant genes with p-adjusted values less
than 0.05. Now that we have extracted our results, we can get a nice
overview of the number of differentially expressed genes there are for
our designated alpha level using the summary() function. It will output
the numbers/percentages of up- and down-regulated genes, as well as,
give information about independent filtering and outliers removed.

``` r
smoc2_results = results(dds_smoc2, contrast = c("condition", "fibrosis", "normal"), lfcThreshold = 0.32, alpha = 0.05)
```

``` r
# Save results as a data frame
smoc2_res_all <- data.frame(smoc2_results)
```

Got It! 1. DESeq2 results In this lesson, we will explore our
differential expression results.

2.  DESeq2 results table To get descriptions for the columns in the
    results table, we can use the mcols() function. The first column is
    the mean value across all samples, followed by the shrunken log2
    foldchanges, standard error of the fold change estimates, the Wald
    statistics output from the Wald test for differential expression,
    the Wald test p-value, and the Benjamini-Hochberg adjusted p-value.

3.  DESeq2 results table Now let’s look at the values in the results
    table and identify the differentially expressed genes. To determine
    significant DE genes, we will be using the p-values adjusted for
    multiple test correction in the last column. The reason for this is
    that for every gene tested with an alpha of 0-point-05, there is a
    5% chance that the gene is called as DE when it is not, yielding
    false positives. If we were to test the roughly 47,000 genes in the
    raw counts file, we would have about 5% or over 2,000 genes as false
    positives. It would be difficult to identify the true positives, or
    genes that are called DE when they truly are, from the false.
    Therefore, multiple test correction is performed by DESeq2 using the
    Benjamini-Hochberg, or BH-method, to adjust p-values for multiple
    testing and control the proportion of false positives relative to
    true. Using the BH-method and an alpha value of 0-point-05, if we
    had 1,000 genes identified as DE, we would expect 5% of the DE genes
    to be false positives, or 50 genes. To reduce the number of genes
    tested, DESeq2 automatically filters out genes unlikely to be truly
    differentially expressed prior to testing, such as genes with zero
    counts across all samples, genes with low mean values across all
    samples, and genes with extreme count outliers. We can see the
    filtered genes in the results tables represented by an NA in the
    p-adjusted column.

4.  Significant DE genes - summary DESeq2’s summary() function provides
    the number of differentially expressed genes for our alpha level and
    information about the number of genes filtered. Our results give
    over 10,000 genes as DE, which is the sum of the DE genes with log2
    fold changes less than 0 and those with fold changes greater than 0.
    This is a lot of genes to sift through. If we wanted to return the
    genes most likely to be biologically relevant, we could also include
    a log2 fold change threshold. Oftentimes, a log2 fold change
    threshold isn’t preferred. However, it can be helpful when dealing
    with such large numbers of DE genes.

5.  Significant DE genes - fold-change threshold To test for significant
    genes using both an alpha value threshold and a log2 foldchange
    threshold different from 0, we need to re-run the results function.
    Let’s use a small 1-point-25-foldchange threshold, which equals
    0-point-32 on the log2 scale, by adding the lfcThreshold argument to
    our results() function. While using any log2 fold change cut-off
    increases the risk of losing biologically relevant genes, by using a
    very small log2 foldchange threshold, we are hoping to reduce the
    risk that the genes more biologically meaningful.

6.  Significant DE genes - summary Now that we have the results, we need
    to re-shrink the foldchanges, then run the summary() function again.
    Now, we have returned just over 6,000 DE genes.

7.  Results - annotate To better understand which genes the results
    pertain to, we can use the annotables package to quickly obtain gene
    names for the Ensembl gene IDs using the table of gene annotations
    for the Grch38 mouse genome build.

8.  Results - extract To annotate the genes with gene names and
    descriptions, we need to first turn our results table into a data
    frame using the data-dot-frame() function. Then, after changing the
    row names to a column, we can merge the gene names and descriptions
    with our results using the left_join() function and merging by
    Ensembl gene IDs. Now we have our entire results table.

9.  Significant DE genes - arrange To extract the significant DE genes,
    we’ll subset the results table, using the subset() function, for
    genes with p-adjusted values less than the alpha value of
    0-point-05. We should see all p-adjusted values are less than
    0-point-05 and log2 foldchanges are greater than the absolute value
    of 0-point-32. We will use the arrange() function to order the genes
    by p-adjusted values to generate the final table of significant
    results. We can explore this table for interesting or expected genes
    with a high probability of being related to kidney fibrosis.

10. Significant DE genes If we look up many of these top genes, we will
    find known roles associated with fibrosis, which is an encouraging
    and exciting result for us.

To reduce the number of DE genes that we are returning and to reduce the
likelihood of the DE genes being biologically meaningful, we are going
to use a small log2 fold change threshold to determine the DE genes.

``` r
# Subset the results to only return the significant genes with p-adjusted values less than 0.05
smoc2_res_sig <- subset(smoc2_res_all, padj < 0.05)
```

``` r
head(smoc2_res_sig)
```

    ##      baseMean log2FoldChange      lfcSE       stat       pvalue         padj
    ## 3    22.47809      4.4981432 0.82929064   5.038213 4.698974e-07 2.829520e-06
    ## 17   12.06950     -2.3959881 0.60066414  -3.456155 5.479409e-04 2.278479e-03
    ## 33 1380.35712     -0.8942696 0.09748260  -5.890996 3.838750e-09 2.815588e-08
    ## 35 2522.97515     -1.9163511 0.14944995 -10.681510 1.242343e-26 2.749899e-25
    ## 40   11.55182      2.3983021 0.70397545   2.952237 3.154811e-03 1.160401e-02
    ## 45 1921.19192     -0.9062709 0.08444206  -6.942877 3.841942e-12 3.594718e-11

4.  Visualizing results - Volcano plot In addition to the MA plot
    explored previously, another useful plot providing a global view of
    the results is the volcano plot, which shows the fold changes
    relative to the adjusted p-values for all genes. First, using all
    results, wt_res_all, convert the row names to a column called
    ensgene, then create a column of logical values indicating if the
    gene is DE using the mutate() function, with p-adjusted value
    threshold less than 0-point-05. Then, use ggplot2 to plot the log2
    foldchange values versus the -log10 adjusted p-value. The points for
    the genes should then be colored by whether they are significant
    using the threshold column.

5.  Visualizing results - Volcano plot We can zoom in on the volcano
    plot to visualize better the significance cut-off using the ylim()
    function within ggplot2.

To explore the results, visualizations can be helpful to see a global
view of the data, as well as, characteristics of the significant genes.
Usually, we expect to see significant genes identified across the range
of mean values, which we can plot using the MA plot. If we only see
significant genes with high mean values, then it could indicate an issue
with our data. The volcano plot helps us to get an idea of the range of
fold changes needed to identify significance in our data.

``` r
# Generate logical column 
smoc2_res_all <- data.frame(smoc2_results) %>% mutate(threshold = padj < 0.05)
```

``` r
# Create the volcano plot
ggplot(smoc2_res_all) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))
```

    ## Warning: Removed 25395 rows containing missing values (geom_point).

![](pre-post_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->
