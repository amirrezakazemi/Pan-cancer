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

genemotif_count = read_tsv("data/counts/genemotif_count.tsv.gz")
cancers = genemotif_count$cancer_type %>% unique()

find_sgnf_genemotifs <- function(cancer) {
  print(cancer)
  genemotif_count %>% dplyr::filter(cancer_type == cancer) -> df
  # descdist(df$count, boot = 100, discrete = TRUE)
  # descdist(df$count, boot = 100, discrete = FALSE)
  
  # beta binomial for now
  n = max(df$count)
  fit <- vglm(cbind(df$count, n - df$count) ~ 1, betabinomial, data = df, trace = TRUE) 
  mu = Coef(fit)[1]
  rho = Coef(fit)[2]
  
  x = seq(0, 1, length.out = length(df$count))
  ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 100) +
    stat_function(fun = VGAM::dbetabinom, args = list(n, mu, rho), color = "blue", alpha = 0.8) + 
    ggtitle(cancer)
  # geom_line(aes(x, VGAM::dbetabinom(x, n, mu, rho)), color = "blue")
  
  # ks.test(df$count,'pbetabinom',n, mu, rho)
  # goftest::cvm.test(df$count, 'pbetabinom', n, mu, rho)
  
  pval = 1 - pbetabinom(df$count, n, rho = rho, prob = mu)
  
  df = cbind(df, pval)
  sgnf_genes = df %>% dplyr::filter(pval < treshold) %>% dplyr::select(genemotif, cancer_type, count)
  print(cat("number of sgnf in ", cancer, " : ", nrow(sgnf_genes), "of ", length(unique(df$gene_id)), "  ---------------------------------------"))
  return(sgnf_genes)
}

treshold = 0.001
sigenes_mat = data.frame(matrix(ncol = 3, nrow = 0))
colnames(sigenes_mat) <- c("genemotif", "cancer_type", "count")

for (cancer in cancers) {
  sgnf_genes = find_sgnf_genemotifs(cancer)
  sigenes_mat = rbind(sigenes_mat, sgnf_genes)
}

ggplot(sigenes_mat) + geom_tile(aes(y = gene_id, x = factor(cancer_type), fill = count)) +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()

sigenes_mat %>% dplyr::select(gene_id) %>% unique() -> sgnf_genes_final
write_tsv(sgnf_genes_final, file.path("sgnf_genes", "allcancer_0.05.tsv"))



