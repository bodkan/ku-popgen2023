library(slendr)
init_env()


# Exercise #3 ------------------------------------------------------------------

chimp <- population("CHIMP", time = 7e6, N = 5000)
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# this is a thing we added to the base model without gene flow
gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 50000)

model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,  # <-------------------------- we need to compile this too
  generation_time = 30
)

plot_model(model, proportions = TRUE)
plot_model(model, proportions = TRUE, log = TRUE)

ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)

# for comparison of gene flow and non-gene flow models:
# ts_save(ts, "/tmp/ex3.trees")
# ts <- ts_load("ex3.trees", model)
# ts <- ts_load("ex1.trees", model)


# Exercise #4 ------------------------------------------------------------------

library(dplyr)
library(ggplot2)

# first get a (named!) list of individuals in each population

samples <- ts_names(ts, split = "pop")
samples

# compute diversity

pi <- ts_diversity(ts, sample_sets = samples)

arrange(pi, diversity)

# compute divergence

div <- ts_divergence(ts, sample_sets = samples)

arrange(div, divergence)


# outgroup f3

ts_f3(ts, B = "AFR_1", C = "EUR_1", A = "CHIMP_1")

samples[["AFR"]] %>% head(5)
samples[["EUR"]] %>% head(5)

ts_f3(ts, B = samples["AFR"], C = samples["EUR"], A = "CHIMP_1")

ts_f3(ts, B = samples["AFR"], C = samples["NEA"], A = "CHIMP_1")

ts_f3(ts, B = samples["EUR"], C = samples["NEA"], A = "CHIMP_1")


# admixture detection -----------------------------------------------------

# D(AFR, EUR; NEA, CHIMP)

ts_f4(ts, W = "AFR_1", X = "AFR_2", Y = "NEA_1", Z = "CHIMP_1")

ts_f4(ts, W = "AFR_1", X = "EUR_1", Y = "NEA_1", Z = "CHIMP_1")

set.seed(42)
afr_samples <- samples$AFR %>% sample(25)
eur_samples <- samples$EUR %>% sample(25)

f4_afr <- lapply(afr_samples, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()
f4_eur <- lapply(eur_samples, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()

f4_afr$pop <- "AFR"
f4_eur$pop <- "EUR"

rbind(f4_afr, f4_eur) %>%
  ggplot(aes(pop, f4, color = pop)) +
  geom_boxplot() +
  geom_jitter() +
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("f4(AFR, EUR; NEA, CHIMP)")
