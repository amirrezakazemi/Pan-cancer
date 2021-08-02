# pan-cancer
this research project is about introducing a machine learning technique to cluster cancer samples based on their somatic mutation profiles. the problem is unsupervised clustering. we performed statistical biomarker analyses to evaluate our subtyping approach.
## preprocessing
gene annotation and obtaining mutation profile of samples. corresponding codes are in "preprocessing" folder<br/>
codes:<br/> 
add_gene.r & add_genemotif.r: in thses two scripts we use ensembl data to annotate dnas and idenify start and end of each gene. input: icgc data, output: a dataframe containing new columns.<br/>
count_genes.r & count_genemotifs.r: these two scripts are used to make a matrix that rows are samples and columns are genes/genemotifs. input: icgc data , output: the matrix<br/>
health samples.r: to filter only patients<br/>
intgr_donor.r: preprocessing for clinical analysis. <br/>
summarize_data.r: to filter only somatic mutations in icgc dataset <br/>
integrate.r: preprocess ensembl dataset and filter coding genes.<br/>

## extract meaningful features.
the mutation profile of cancer samples is highly dispersing, so we need to identify genes that are significantly mutated and obtained genes are our candidate features. since our project is on all cancer types we identified genes in each cancer separately. to do this, we performed a background statistical method (negative binomial distribution) to fit mutation data of each cancer. then based on a fixed p-value we obtained significantly mutated genes.corresponding codes are in "sgnf_genes" folder.<br/>
codes:<br/>
decide_distr.r: to choose the best fitting distribution on mutational data for each cancertype<br/>
genes_bg.r: fit the negative binomial distribution and use the threshold to filter and extract the most significant genes in each cancertype.<br/>
genemotiffs_bg.r: fit the beta binomial distribution and use the threshold to filter and extract the most significant genemotifs. after these two experiments we only use genes as our features. so the input for clustering process is a matrix containing number of mutations in each significant gene for each sample.<br/>


## clustering
after obtaining features, we used a classical machine learning method to cluster samples. a hierarchical mixture model is performed in this step. we also examined other clustering methods such as k-means and db-scan. corresponding codes are in "clustering" folder.<br/>
codes:<br/>
clustring1.r: this code is performing the process of clustring with kmeans and mclust<br/>
subclusters.r: this code is for setting the threshold of 95 percent. we examined 2 threshold and identified two sets of subtypes and then we chose the best one.<br/>
cluster_validation.r: performing statistical analysis on our subtypes to validate the clustering performance by mclust package.<br/>
cluster_analysis.r: cluster samples using mclust and plotting our subtypes and save the result in bio_res folder.<br/>

## statistical biomarker analysis
in this step, we attempted to validate our clustering using biological interpretations. different biomarker analyses such as enriched genes, enriched gene motifs, enriched gene ontologies, signature and survival analysis are performed and clusters had distinguished characteristics. corresponding codes are in "biomarker_analysis" folder.<br/>
codes:<br/>
statistical.R:input is our preprocessed data with labels and output is the plot of distribution of cancertypes in each subtype.<br/>
survival_analysis.R: we used survminer and ggsurvplot to plot and evaluate the survival time of patients. the result will be saved in "bio_res" folder. <br/>
top_motifs.R: this code is to compute the top 100 gene-motifs in each subtype and then compute the intersection between each two subtypes. the result will be save as a figure and input is sample_data(summarized_data computing in preprocess step).<br/>
visualization.R: this code is to plot a bar figure in order to show the number of samples in each subtype, the cancertype distribution in each subtype. each cancertype will be shown by a color. <br/>
webgetstalt.R: using the most significant genes in each subtype(saved in top100_exp_codings), we used the R tool to analyze the subtypes in pathway, ontology and disease.
sig_motif.R: applying fisher's exact test to extract top 100 motifs in each subtype.<br/>
pairwise_motif.R: this analysis is to show that subtypes mutated in a gene, have been mutated in different motifs and thus their samples are diferent. we evaluate this analysis for all pairs and filter the best results.<br/>
mutation_rate.R: we used sample_data as inputs and plot three types of mutational rate analysis, for feature genes, for allcoding genes and for all non-coding genes. mutational rate is the fraction of subtypes' samples that are mutated in a specific gene.<br/>
motif_analysis.R: in this analysis, we first choose some genes that are highly mutated in many subtypes.(using top100 genes analysis) and then compare the mutation motifs for that gene in all subtypes to indicate that the mutation are occuring in different context.<br/>
gene_rate.R: two types of gene rate analysis in the script, 1. considering only samples that have 3/5 mutations in a specific gene (importance of number of samples).  2.differential analysis(importance of different rate of mutation).<br/>
fisher_genes.R: the fisher method is written in this script and is applicable in other scripts.<br/>
consequence_type.R: the plot of consequence type will be generated by this script. the summarized data has cosequnce type for each sample and we plotted a bar figure for each subtype using this script<br/>
clinical_analysis.R: two plots in this script: 1. gender distribution for each subtype(male/female), 2. region distribution for each subtype. the input is not all samples and we filter those samples whom their clinical data were available.






