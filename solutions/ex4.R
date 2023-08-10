library(dplyr)
library(ggplot2)

library(slendr)
init_env()

chimp <- population("CHIMP", time = 7e6, N = 5000)
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

gf <- gene_flow(from = nea, to = eur, rate = 0.05, start = 55000, end = 50000)

model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,
  generation_time = 30
)


# define sampling events
nea_samples <- schedule_sampling(model, times = c(70000, 40000), list(nea, 1))
present_samples <- schedule_sampling(model, times = 0, list(chimp, 1), list(afr, 5), list(eur, 10))
emh_samples <- schedule_sampling(model, times = runif(n = 40, min = 10000, max = 40000), list(eur, 1))

ts <- msprime(model, samples = rbind(nea_samples, present_samples, emh_samples),
              sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)


# extract table with names and times of sampled Europeans (ancient and present day)
eur_inds <- ts_samples(ts) %>% filter(pop == "EUR")
eur_inds

# compute f4-ration statistic (this will take ~30s)
nea_ancestry <- ts_f4ratio(ts, X = eur_inds$name, "NEA_1", "NEA_2", "AFR_1", "CHIMP_1")
nea_ancestry

# add the computed proportions to the data frame
eur_inds$ancestry <- nea_ancestry$alpha

eur_inds %>%
  ggplot(aes(time, ancestry)) +
    geom_point() +
    geom_smooth(method = "lm", linetype = 2, color = "red", linewidth = 0.5) +
    xlim(40000, 0) + coord_cartesian(ylim = c(0, 0.1)) +
    labs(x = "time [years ago]", y = "Neanderthal ancestry proportion")
