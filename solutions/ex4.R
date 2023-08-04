library(dplyr)

library(slendr)
init_env()

chimp <- population("CHIMP", time = 7e6, N = 5000)
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

gf <- gene_flow(from = nea, to = eur, rate = 0.05, start = 50000, end = 40000)

model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,
  generation_time = 30
)

plot_model(model, sizes = FALSE, proportions = TRUE)
plot_model(model, log = TRUE, proportions = TRUE)

ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)
ts


# extract samples as a list of names --------------------------------------

samples <- ts_samples(ts) %>%
  split(., .$pop) %>%
  lapply(pull, "name")
samples

# compute diversity -------------------------------------------------------

pi <- ts_diversity(ts, sample_sets = samples)

arrange(pi, diversity)

# compute divergence ------------------------------------------------------

div <- ts_divergence(ts, sample_sets = samples)

arrange(div, desc(divergence))


# outgroup f3 -------------------------------------------------------------

ts_f3(ts, A = "AFR_1", B = "EUR_1", C = "CHIMP_1")

samples$AFR
samples$EUR
ts_f3(ts, A = samples$AFR, B = samples$EUR, C = "CHIMP_1")


ts_f3(ts, A = list(afr = samples$AFR), B = list(eur = samples$EUR), C = "CHIMP_1")

ts_f3(ts, A = list(afr = samples$AFR), B = list(nea = samples$NEA), C = "CHIMP_1")


# admixture detection -----------------------------------------------------

# D(AFR, EUR; NEA, CHIMP)

ts_f4(ts, W = "AFR_1", X = "AFR_2", Y = "NEA_1", Z = "CHIMP_1")

ts_f4(ts, W = "AFR_1", X = "EUR_1", Y = "NEA_1", Z = "CHIMP_1")

afr_samples <- samples$AFR %>% sample(25)
eur_samples <- samples$EUR %>% sample(25)

f4_afr <- lapply(afr_samples, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()
f4_eur <- lapply(eur_samples, function(x) ts_f4(ts, W = "AFR_1", X = x, Y = "NEA_1", Z = "CHIMP_1")) %>% bind_rows()

f4_afr$pop <- "AFR"
f4_eur$pop <- "EUR"

library(ggplot2)

rbind(f4_afr, f4_eur) %>%
  ggplot(aes(pop, f4, color = pop)) +
  geom_boxplot() +
  geom_jitter() +
  geom_hline(yintercept = 0, linetype = 2)
