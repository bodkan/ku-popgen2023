library(slendr)
init_env()

# chimpanzee outgroup
chimp <- population("CHIMP", time = 7e6, N = 5000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 15000)
eur <- population("EUR", parent = afr, time = 70e3, N = 3000)

# Neanderthal population splitting at 600 ky ago from modern humans
# (becomes extinct by 40 ky ago)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# compile the entire model into a single object
model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  generation_time = 30
)

# verify visually
plot_model(model)
plot_model(model, sizes = FALSE)
plot_model(model, sizes = FALSE, log = TRUE)


# notes -------------------------------------------------------------------

# models can also be serialized to a permanent location

model_dir <- "/tmp/ex1_model"

model <- compile_model(
  populations = list(chimp, nea, afr, eur),
  generation_time = 30,
  path = model_dir
)

list.files(model_dir)




# you can ignore the rest -------------------------------------------------

# ts <-
#   msprime(model, sequence_length = 100e6, recombination_rate = 1e-8) %>%
#   ts_mutate(mutation_rate = 1e-8)
# 
# ts_save(ts, "/tmp/ex1.trees")
