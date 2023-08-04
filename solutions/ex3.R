library(dplyr)

library(slendr)
init_env()

chimp <- population("CHIMP", time = 7e6, N = 5000)
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# this is a thing we added to the base model without gene flow
gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 50000, end = 40000)

model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,                            # <--- we need to compile this too
  generation_time = 30
)

ts <- msprime(model, sequence_length = 10e6, recombination_rate = 1e-8) %>%
  ts_mutate(1e-8)

samples <- ts_samples(ts)

samples <- ts_samples(ts) %>%
  split(., .$pop) %>%
  lapply(pull, "name")


ts_diversity(ts, sample_sets = samples)

?ts_diversity

# verify visually
plot_model(model)
plot_model(model, sizes = FALSE, proportions = TRUE)
plot_model(model, sizes = FALSE, log = TRUE, proportions = TRUE)
