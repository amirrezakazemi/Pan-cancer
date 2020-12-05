# Pan-cancer
This research project is about introducing a machine learning technique to cluster cancer samples based on their somatic mutation profiles. The problem is unsupervised clustering. We performed statistical biomarker analyses to evaluate our subtyping approach.
## Extract meaningful features.
The mutation profile of cancer samples is highly dispersing, so we need to identify genes that are significantly mutated and obtained genes are our candidate features. Since our project is on all cancer types we identified genes in each cancer separately. To do this, we performed a background statistical method (negative binomial distribution) to fit mutation data of each cancer. then based on a fixed P-value we obtained significantly mutated genes.
## Clustering
After obtaining features, we used a classical machine learning method to cluster samples. A hierarchical mixture model is performed in this step. Corresponding codes are in the analysis folder. We also examined other clustering methods such as k-means and db-scan.
## Statistical biomarker analysis
In this step, we attempted to validate our clustering using biological interpretations. Different biomarker analyses such as enriched genes, enriched gene motifs, enriched gene ontologies, signature and survival analysis are performed and clusters had distinguished characteristics.
