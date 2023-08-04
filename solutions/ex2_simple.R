library(slendr)
init_env()

simulate_afs <- function(Ne) {
  # create a slendr model with a given N size
  pop <- population("pop", N = Ne, time = 1)
  model <- compile_model(pop, generation_time = 1,
                         simulation_length = 100000)
  
  # simulate a tree sequence
  ts <-
    msprime(model, sequence_length = 10e6, recombination_rate = 1e-8) %>%
    ts_mutate(mutation_rate = 1e-8)
  
  # get a random sample of names of 10 individuals
  samples <- ts_samples(ts) %>%
    dplyr::sample_n(10) %>%
    dplyr::pull(name)
  
  # compute the AFS vector
  afs <- ts_afs(ts, list(samples), polarised = TRUE)
  
  afs
}

afs_observed <- c(2520, 1449, 855, 622, 530, 446, 365, 334, 349, 244,
                  264, 218,  133, 173, 159, 142, 167, 129, 125, 143)

# we know that the Ne is between 1000 and 30000
afs_Ne1k <- simulate_afs(1000)
afs_Ne5k <- simulate_afs(5000)
afs_Ne6k <- simulate_afs(6000)
afs_Ne10k <- simulate_afs(10000)
afs_Ne20k <- simulate_afs(20000)
afs_Ne30k <- simulate_afs(30000)

plot(afs_observed, type = "b", col = "black", lwd = 3,
     xlab = "allele count bin", ylab = "count", ylim = c(0, 13000))
lines(afs_Ne1k, lwd = 2, col = "blue")
lines(afs_Ne5k, lwd = 2, col = "green")
lines(afs_Ne6k, lwd = 2, col = "pink")
lines(afs_Ne10k, lwd = 2, col = "purple")
lines(afs_Ne20k, lwd = 2, col = "orange")
lines(afs_Ne30k, lwd = 2, col = "cyan")
legend("topright",
       legend = c("observed AFS", "Ne = 1000", "Ne = 5000",
                  "Ne = 6000", "Ne = 10000", "Ne = 20000", "Ne = 30000"),
       fill = c("black", "blue", "green", "pink", "purple", "orange", "cyan"))

