library(slendr)
init_env()


# Exercise #1 ------------------------------------------------------------------

# chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 5000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)

# Neanderthal population splitting at 600 ky ago from modern humans
# (becomes extinct by 40 ky ago)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# Neanderthal introgression event (3% admixture between 55-50 kya)
gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 50000)

# compile the entire model into a single object
model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  gene_flow = gf,
  generation_time = 30
)

plot_model(model, proportions = TRUE)
plot_model(model, proportions = TRUE, log = TRUE)

ts <-
  msprime(model, sequence_length = 100e6, recombination_rate = 1e-8) %>%
  ts_mutate(mutation_rate = 1e-8)

# ts_save(ts, "ex3.trees")
# ts <- ts_load("ex3.trees", model)

# Exercise #4 ------------------------------------------------------------------

library(dplyr)
library(ggplot2)

# first get a (named!) list of individuals in each population

samples <- ts_names(ts, split = "pop")
samples

lapply(samples, head, 5) # get just 5 samples from each population

# compute nucleotide diversity (pi) in each population

pi <- ts_diversity(ts, sample_sets = samples)

arrange(pi, diversity)

# heterozygosity in a single sample
ts_diversity(ts, "NEA_1")

# compute pairwise population divergence

div <- ts_divergence(ts, sample_sets = samples)

arrange(div, divergence)


# outgroup f3

ts_f3(ts, B = "AFR_1", C = "EUR_1", A = "CHIMP_1")

samples[["AFR"]] %>% head(5)
samples[["EUR"]] %>% head(5)

ts_f3(ts, B = samples["AFR"], C = samples["EUR"], A = "CHIMP_1")

ts_f3(ts, B = samples["AFR"], C = samples["NEA"], A = "CHIMP_1")

ts_f3(ts, B = samples["EUR"], C = samples["NEA"], A = "CHIMP_1")


# outgroup f3 as a linear combination of f2 statistics

# branch-based version
ts_f3(ts, B = "AFR_1", C = "AFR_2", A = "CHIMP_1", mode = "branch")

my_f3 <- (
  ts_f2(ts, A = "AFR_1", B = "CHIMP_1", mode = "branch")$f2 +
  ts_f2(ts, A = "AFR_2", B = "CHIMP_1", mode = "branch")$f2 -
  ts_f2(ts, A = "AFR_1", B = "AFR_2", mode = "branch")$f2
) / 2
my_f3


# admixture detection -----------------------------------------------------

# D(AFR, EUR; NEA, CHIMP)      ~  (BABA - ABBA) / (BABA + ABBA)


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



# unique quartets

# # install a combinatorics R package
# install.packages("combinat")

quartet <- c("AFR_1", "EUR_1", "NEA_1", "CHIMP_1")

quartets <- combinat::permn(quartet)

# how many permutations there are in total?
#   4! = 4 * 3 * 2 * 1 = 24

length(quartets)

# loop across all quartets, computing the corresponding f4 statistic
all_f4s <- lapply(quartets, function(q) ts_f4(ts, q[1], q[2], q[3], q[4], mode = "branch"))

# bind the list of f4 results into a single data frame
all_f4s <- bind_rows(all_f4s) %>% arrange(abs(f4))
print(all_f4s, n = Inf)

distinct(all_f4s, f4, .keep_all = TRUE)
distinct(all_f4s, abs(f4), .keep_all = TRUE)
