library(slendr)
init_env()

simulate_afs <- function(Ne) {
  # create a slendr model with a given N size
  pop <- population("pop", N = Ne, time = 100000)
  model <- compile_model(pop, generation_time = 1, direction = "backward")
  
  # simulate a tree sequence
  ts <-
    msprime(model, sequence_length = 10e6, recombination_rate = 1e-8) %>%
    ts_mutate(mutation_rate = 1e-8)
  
  # get a random sample of names of 10 individuals
  samples <- ts_samples(ts) %>% dplyr::sample_n(10) %>% dplyr::pull(name)
  
  # compute the AFS vector
  afs <- ts_afs(ts, list(samples), polarised = TRUE)
  
  afs
}

afs_observed <- c(2520, 1449, 855, 622, 530, 446, 365, 334, 349, 244,
                  264, 218,  133, 173, 159, 142, 167, 129, 125, 143)

# regularly spaced values of potential Ne
Ne_grid <- seq(from = 1000, to = 30000, by = 1000)
Ne_grid

afs_grid <- lapply(Ne_grid, simulate_afs) # saveRDS(afs_grid, "data/ex2_grid.rds")
afs_grid <- readRDS("data/ex2_grid.rds")
names(afs_grid) <- Ne_grid

# plot the observed AFS and overlay the simulated AFS vectors on top of it
plot(afs_observed, type = "b", col = "red", lwd = 3,
     xlab = "allele count bin", ylab = "count")
for (i in seq_along(Ne_grid)) {
  lines(afs_grid[[i]], lwd = 0.5)
}
legend("topright", legend = c("observed AFS", "simulated AFS"),
       fill = c("red", "black"))

# compute mean-squared error of the AFS produced by each Ne value across the grid
errors <- sapply(afs_grid, function(sim_afs) {
  sum((sim_afs - afs_observed)^2) / length(sim_afs)
})

# plot the errors, highlight the most likely value
plot(Ne_grid, errors, ylab = "error")
abline(v = Ne_grid[which.min(errors)], col = "green")
abline(v = TRUE_NE, col = "purple")
legend("topright", legend = paste("closest inferred Ne =",
                                  Ne_grid[which.min(errors)]), fill = "green")

# Plot the AFS again, highlighting the most likely spectrum
plot(afs_observed, type = "b", col = "red", lwd = 1,
     xlab = "allele count bin", ylab = "count")
for (i in seq_along(Ne_grid)) {
  color <- if (i == which.min(errors)) "green" else "gray"
  width <- if (i == which.min(errors)) 2 else 0.5
  lines(afs_grid[[i]], lwd = width, col = color)
}
legend("topright", legend = c("observed AFS", paste("Ne =", Ne_grid[which.min(errors)])),
       fill = c("red", "green"))
