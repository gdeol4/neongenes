Got It! 1. Differential gene expression overview Differential gene
expression analysis is a powerful technique to determine whether genes
are expressed at significantly different levels between two or more
samplegroups. We will use the DESeq2 package to model the gene counts
and identify differentially expressed genes.

2.  Differential expression analysis: Goal In this image, we see a heat
    map of genes as rows colored by number of counts. These genes
    represent the genes with large expression differences or fold
    changes between sample groups. To determine which genes are
    differentially expressed, one might ask ‘why not just identify the
    genes with the largest fold changes in expression between sample
    groups?’.

3.  Differential expression analysis: Goal To get at the answer, let’s
    observe the plot of normalized counts for gene A. The points
    represent the gene A expression levels for five biological
    replicates for ‘untreated’ and ‘treated’ conditions. The mean
    expression for the ‘treated’ condition is over twice that of the
    untreated. However, there appears to be greater variation in the
    ‘treated’ condition and the difference in expression may not be
    significant. We need to account for variation in the data when we
    determine whether genes are differentially expressed. Therefore, the
    goal of differential expression analysis is to determine for each
    gene whether the differences in expression between groups is
    significant given the amount of variation within groups, or between
    the biological replicates.

4.  Introduction to dataset To explore the workflow, we will be using a
    publicly available RNA-Seq dataset from Gerarduzzi et al from the
    journal JCI Insight. In this paper, the goal of the RNA-Seq
    experiment was to explore why mice over-expressing the Smoc2 gene,
    or producing more Smoc2 mRNA than normal, are more likely to develop
    kidney fibrosis.

5.  Introduction to dataset: Smoc2 Smoc2, or Secreted modular
    calcium-binding protein 2, has been shown to have increased
    expression in kidney fibrosis, which is characterized by an excess
    of extracellular matrix in the space between tubules and capillaries
    within the kidney. However, it is unknown how Smoc2 functions in the
    induction and progression of fibrosis.

6.  Introduction to dataset: Experimental design There are four sample
    groups being tested: normal, control mice, referred to as wild type
    mice, with and without fibrosis and Smoc2 over-expressing mice with
    and without fibrosis. There are three biological replicates for all
    normal samples and four replicates for all fibrosis samples.
    Initially, we will explore the effect of fibrosis on gene expression
    using ‘Wild type’ samples during lectures and ‘Smoc2
    over-expression’ data during exercises.

7.  RNA-Seq count distribution To test whether the expression of genes
    between two or more groups is significantly different, we need an
    appropriate statistical model. An appropriate statistical model is
    determined by the count distribution. When we plot the distribution
    of counts for a single sample, we can visualize key features of
    RNA-Seq count data, including a large proportion of genes with low
    counts and many genes with zero counts. Also note the long right
    tail, which is due to there being no limit for maximum expression in
    RNA-Seq data. If there was no expression variation between
    biological replicates, a frequently used count distribution known as
    the Poisson distribution, would be an appropriate model. But, there
    is always biological variation, and this additional variation
    present in RNA-Seq data can be modeled well using the negative
    binomial model, which we will be using as part of DESeq2.

8.  Preparation for differential expression analysis: DESeq2 object To
    start the differential expression analysis we use the
    `DESeqDataSetFromMatrix()` function, which takes a raw count matrix
    as input, along with the metadata and a design formula to create the
    DESeq2 object. The design formula given should contain major
    expected sources of variation to control for and the condition of
    interest as the last term in the formula. If the raw count data is a
    Summarized Experiment from the htseq-count tool, or generated by
    pseudo-alignment tools, DESeq2 has other functions to use to create
    the DESeq2 object as detailed in the vignette.

9.  Preparation for differential expression analysis: metadata In
    addition to our raw counts, we require sample metadata. At the very
    least, we need to know which of our samples correspond to each
    condition. To generate our metadata, we create a vector for each
    column and combine the vectors into a data frame. The sample names
    are added as the row names.

10. Preparation for differential expression analysis: metadata After we
    have the counts and metadata files, we can start the differential
    expression analysis workflow.

11. Let’s practice! Let’s practice exploring counts and getting our
    files ready for analysis.

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](2022-05-11-performing-rna-seq-with-r_files/figure-gfm/pressure-1.png)<!-- -->

``` r
install.packages('readr')
```

    ## Installing package into 'C:/Users/gurka/AppData/Local/R/win-library/4.2'
    ## (as 'lib' is unspecified)

    ## package 'readr' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\gurka\AppData\Local\Temp\RtmpimOLRT\downloaded_packages

``` r
library(readr)
```

``` r
data <- read_csv("C:/Users/gurka/Downloads/fibrosis_smoc2_rawcounts/fibrosis_smoc2_rawcounts.csv")
```

    ## New names:
    ## Rows: 47729 Columns: 8
    ## ── Column specification
    ## ────────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (7): smoc2_fibrosis1, smoc2_fibrosis4, smoc2_normal1,
    ## smoc2_normal3, smoc2...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ Specify
    ## the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
head(data)
```

    ## # A tibble: 6 × 8
    ##   ...1               smoc2_fibrosis1 smoc2_fibrosis4 smoc2_normal1 smoc2_normal3
    ##   <chr>                        <dbl>           <dbl>         <dbl>         <dbl>
    ## 1 ENSMUSG00000102693               0               0             0             0
    ## 2 ENSMUSG00000064842               0               0             0             0
    ## 3 ENSMUSG00000051951              72              30             0             3
    ## 4 ENSMUSG00000102851               0               0             0             0
    ## 5 ENSMUSG00000103377               0               0             1             0
    ## 6 ENSMUSG00000104017               0               0             0             0
    ## # … with 3 more variables: smoc2_fibrosis3 <dbl>, smoc2_normal4 <dbl>,
    ## #   smoc2_fibrosis2 <dbl>

``` r
data <- subset (data[-c(1)])
```

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

1.  Differential expression analysis In chapter 1 we learned about the
    goals of differential expression analysis, features of RNA-Seq count
    data, and the types of data required for performing differential
    expression or DE analysis. In this chapter, we will explore the
    workflow using DESeq2 and prepare our fibrosis experimental data for
    differential expression analysis.

2.  Differential expression analysis: tools While there are a large
    number of statistical packages developed for DE analysis, DESeq2 and
    EdgeR are two of the most popular tools. Both tools use the negative
    binomial model to model the raw count data, employ similar methods,
    and typically, yield similar results. Both are pretty stringent and
    have a good balance between sensitivity and specificity. Limma-Voom
    is another set of tools often used together for DE analysis, while
    also a good method, it can be a bit less sensitive for smaller
    sample sizes. We will be using DESeq2 for the DE analysis, which has
    an extensive vignette available to help with questions.

3.  Differential expression analysis: DESeq2 vignette We can use the
    `vignette()` function to open the DESeq2 vignette, which contains
    detailed information about the workflow. This should often be the
    first place to look when you have questions regarding the tool or
    workflow. For example, if we want to know how to get help for
    specific questions regarding DESeq2, we can click on the second
    bullet point under the Standard workflow.

4.  Differential expression analysis: DESeq2 vignette This link will
    take us to the help section, where the documentation suggests
    posting questions to the Bioconductor support site if answers can’t
    be found in the vignette.

5.  Differential expression analysis: DESeq2 workflow We will be using
    the DESeq2 tool to perform the differential expression analysis, but
    what are the steps we will need to perform? Displayed in the
    workflow are the steps in the differential expression analysis,
    separated into quality control and DE analysis steps. To start, we
    will take the count matrix containing the raw counts for each sample
    and perform quality control steps. First, the counts will be
    normalized to account for differences in library depth. Then,
    principal component analysis and hierarchical clustering using
    heatmaps will be performed to identify potential sample outliers and
    major sources of variation in the data. After QC, the DE analysis
    will be performed, including the modeling of counts, shrinking the
    log2 fold changes, and testing for differential expression. We will
    cover each of these steps in more detail as we progress through the
    workflow.

6.  Bringing in data for DESeq2 To prepare our data for the DESeq2
    workflow, we need the raw counts of the number of reads aligning to
    each gene and the associated sample metadata brought into R. Let’s
    start with the raw counts. We can use the read-dot-csv() function
    like we did previously to bring it in. Then take a quick peek at the
    data frame with the View() function. Just to note, in addition to
    output from standard alignment and counting tools, DESeq2 will also
    take counts from pseudo-alignment tools like Kallisto and Salmon.
    However, these abundance estimates need to be formatted properly for
    DESeq2 by using the tximport package as discussed in the vignette.

7.  Bringing in data for DESeq2: metadata We can also load in the
    metadata file we created earlier using the read-dot-csv() function.
    Now we have both of the files needed for performing the differential
    expression analysis.

8.  Let’s practice! Let’s explore the DESeq2 vignette a bit more on our
    own and practice bringing in our data.
