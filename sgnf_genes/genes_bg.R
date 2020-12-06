{
  library(readr)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(VGAM)
  library(MASS)
  library(fitdistrplus)
  library(highcharter)
}
#using data in "data/counts/gene_count.tsv.gz" we computed the most significant genes
#for each cancer type which is obtained by using fitdistplus. the result is stored in "sgnf_genes/allcancer".

gene_count = read_tsv("data/counts/gene_count.tsv.gz")
cancers = gene_count$cancer_type %>% unique()

find_sigenes <- function(cancer) {
  print(cancer)
  gene_count %>% dplyr::filter(cancer_type == cancer) -> df
  # descdist(df$count, boot = 100, discrete = TRUE)
  # descdist(df$count, boot = 100, discrete = FALSE)
  
  # head and neck -> lognormal
  # esophagus -> poisson
  
  if (cancer == "Skin") {
    ## fit gamma
    fit <- fitdist(df$count, "gamma") 
    shape = coef(fit)[1]
    rate = coef(fit)[2]
    # plot(ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 1000) +
    #        stat_function(fun = dgamma, args = list(shape = unname(shape), rate = unname(rate)), color = "blue", alpha = 0.8) +
    #        ggtitle(cancer))
    
    pval = 1 - pgamma(df$count, shape = shape, rate = rate)
  }
  else {
    tb = as.data.frame(table(df$count))
    colnames(tb) = c("y", "w")
    tb$y = as.numeric(tb$y)
    tb$w = as.numeric(tb$w)
    fit <- vglm(y ~ 1, negbinomial(deviance = TRUE), data = tb,
                weights = w, crit = "coef") 
    mu = Coef(fit)[1]
    size = Coef(fit)[2]
    
    plot(ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 500) +
      stat_function(fun = dnbinom, args = list(size = unname(size), mu = unname(mu)), color = "blue", alpha = 0.8) +
      ggtitle(cancer))
    
    pval = 1 - pnbinom(df$count, size = size, mu = mu)
  }
  df = cbind(df, pval)
  sgnf_genes = df %>% dplyr::filter(pval < treshold) %>% dplyr::select(gene_id, cancer_type, count)
  print(cat("number of sgnf in ", cancer, " : ", nrow(sgnf_genes), "of ", length(unique(df$gene_id)), "  ---------------------------------------"))
  return(sgnf_genes)
}

treshold = 0.001
sigenes_mat = data.frame(matrix(ncol = 3, nrow = 0))
colnames(sigenes_mat) <- c("gene_id", "cancer_type", "count")

for (cancer in cancers) {
  sgnf_genes = find_sigenes(cancer)
  sigenes_mat = rbind(sigenes_mat, sgnf_genes)
}

p1<- ggplot(sigenes_mat) + geom_tile(aes(y = gene_id, x = factor(cancer_type), fill = count)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_minimal()

plot(p1)
ggsave(file.path(paste0("../plots/", mat_name), paste0("sgnf_genes", ".png")), height = 2000, units = "mm")


sigenes_mat %>% dplyr::select(gene_id) %>% unique() -> sgnf_genes
write_tsv(sgnf_genes, file.path("sgnf_genes", "allcancer_0.05.tsv"))



# analysis ----------------------------------------------------------------

threshold = 0.005
sfng_genes = read_tsv("sgnf_genes/allcancer_0.005.tsv")
sgnf_genes = sfng_genes$gene_id

cosmic_genes = read_tsv("sgnf_genes/from_dataset/census_genes.tsv")
cosmic_genes = cosmic_genes$gene_id

sum(sgnf_genes %in% cosmic_genes)
length(sgnf_genes)


