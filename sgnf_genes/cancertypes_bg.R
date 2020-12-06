library(readr)
library(dplyr)
library(MASS)
library(VGAM)
library(fitdistrplus)
library(ggplot2)

gene_count = read_tsv("data/counts/gene_count.tsv.gz")
cancers = gene_count$cancer_type %>% unique()

find_sigenes <- function(cancer) {
  print(cancer)
  gene_count %>% dplyr::filter(cancer_type == cancer) -> df
  # descdist(df$count, boot = 100, discrete = TRUE)
  # fit <- fitdistr(df$count, "negative binomial")
  # size = fit$estimate[1]
  # mu = fit$estimate[2]
  tb = as.data.frame(table(df$count))
  colnames(tb) = c("y", "w")
  tb$y = as.numeric(tb$y)
  fit <- vglm(y ~ 1, negbinomial(deviance = TRUE), data = tb,
              weights = w, crit = "coef") 
  print(Coef(fit))
  mu = Coef(fit)[1]
  size = Coef(fit)[2]
  ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 500) +
      stat_function(fun = dnbinom, args = list(size = unname(size), mu = unname(mu)), color = "blue", alpha = 0.8)
  pval = 1 - pnbinom(df$count, size = size, mu = mu)
  
  # if (cancer == "Skin") {
    # fit <- fitdistr(df$count, "poisson")
    # lambda = fit$estimate
    # ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 1000) +
    #   stat_function(fun = dpois, args = list(lambda = unname(lambda)), color = "blue", alpha = 0.8)
  # }
  
  df = cbind(df, pval)
  sgnf_genes = df %>% dplyr::filter(pval < treshold) %>% dplyr::select(gene_id, cancer_type, pval)
  print(cat("number of sgnf in ", cancer, " : ", nrow(sgnf_genes), "  ---------------------------------------"))
  return(sgnf_genes)
}

treshold = 0.001
sigenes_mat = data.frame(matrix(ncol = 2, nrow = 0))
colnames(sigenes_mat) <- c("gene_id", "cancer_type")

for (cancer in cancers) {
  sgnf_genes = find_sigenes(cancer)
  sigenes_mat = rbind(sigenes_mat, sgnf_genes)
}

sigenes_mat %>% group_by(gene_id) %>% dplyr::summarise(n = dplyr::n()) %>% filter(n == max(n)) -> sgnf

sigenes_mat %>% dplyr::select(gene_id) %>% unique() -> sgnf_genes
write_tsv(sgnf_genes, file.path("sgnf_genes", "allcancer_0.002.tsv"))

write_csv(sigenes_mat, file.path("sgnf_genes", "sgnf_genes_pvals.csv"))

# Genemotif_col -----------------------------------------------------------

genemotif = read_tsv("data/genemotifs/Colorectal_genemotif.tsv.gz")
genemotif %>% 
  dplyr::select(icgc_sample_id, genemotif_3mer) %>% 
  unique() %>% 
  group_by(genemotif_3mer) %>% 
  summarise(count = n()) -> genemotif_count

df = genemotif_count
descdist(df$count, boot = 100, discrete = TRUE) # gamma

fit <- fitdistr(df$count, "beta binomial")
shape = fit$estimate[1]
rate = fit$estimate[2]
ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 25) +
  stat_function(fun = dnbinom, args = list(shape = unname(shape), rate = unname(rate)), color = "blue", alpha = 0.8)
pval = 1 - pgamma(df$count, shape = shape, rate = rate)
df2 = cbind(df, pval = pval)
sgnf_genemotifs = df2 %>% dplyr::filter(pval < 0.001) %>% dplyr::select(genemotif_3mer) # 44k


# gene --------------------------------------------------------------------

gene_count %>% dplyr::filter(cancer_type == "Colorectal") -> df
descdist(df$count, boot = 100, discrete = TRUE)
fit <- fitdistr(df$count, "negative binomial")
size = fit$estimate[1]
mu = fit$estimate[2]
ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 200) +
  stat_function(fun = dnbinom, args = list(size = unname(size), mu = unname(mu)), color = "blue", alpha = 0.8)
ks.test(df$count, "pnbinom", size = size, mu = mu)
goftest::cvm.test(df$count, "pnbinom", size = size, mu = mu)
goftest::ad.test(df$count, "pnbinom", size = size, mu = mu)

pval = 1 - pnbinom(df$count, size = size, mu = mu)

df = cbind(df, pval)
sgnf_genes = df %>% dplyr::filter(pval < 0.01) %>% dplyr::select(gene_id, cancer_type, pval) # 523


# both --------------------------------------------------------------------

sgnf_genemotifs %>% 
  mutate(gene = sub("_*.", "", genemotif_3mer)) %>% 
  filter(gene %in% sgnf_genes$gene_id) %>% 
  View()

