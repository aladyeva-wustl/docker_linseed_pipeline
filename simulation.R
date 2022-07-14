#' Simulation builder initialization: simulation of basis, proportions and gene expression data without noise.
#' The motivation behind this builder approach is that every add... method creates a copy of the input simulation,
#' so that it is possible to create a single base simulation with many different child simulations.
#' For example, with different noise applied to the same data.
#' The add... functions inherently support magrittr's pipe operator (%>%).
#' @examples
#' sim_initial <- createSimulation(12000, 100, 3)
#' sim_noise_3 <- addNoise(sim_initial, 3)
#' sim_noise_5 <- sim_initial %>% addNoise(5) %>% addBasisSamples()
createSimulation <- function(n_genes, n_samples, n_cell_types, with_marker_genes = F) {
  if (with_marker_genes) {
    n_genes <- n_genes - n_cell_types
  }
  basis <- generateBasis(n_genes, n_cell_types)
  if (with_marker_genes) {
    tmp <- appendMarkerGenesToBasis(basis)
    basis <- tmp[["basis"]]
    marker_gene_names <- tmp[["marker_gene_names"]]
  }

  proportions <- generateProportions(n_samples, n_cell_types)
  data <- basis %*% proportions

  data[data < 0] <- 0

  res <- list(
    basis = basis,
    proportions = proportions,
    data = data,
    mixed_sample_names = colnames(data)
  )
  if (with_marker_genes) {
    res$marker_gene_names <- marker_gene_names
  }

  res
}

#' Add N pure samples (basis + noise) per cell type to the simulation
addPureSamples <- function(simulation, samples_per_cell_type, noise_deviation = 0.05) {
  pure_samples <- generatePureSamples(
    simulation$basis,
    samples_per_cell_type,
    noise_deviation = noise_deviation
  )
  simulation$pure_sample_names <- colnames(pure_samples$data)
  simulation$data <- cbind(simulation$data, pure_samples$data)
  simulation$proportions <- cbind(simulation$proportions, pure_samples$proportions)
  simulation$pure_samples_noise_deviation <- noise_deviation
  simulation
}

#' Add samples, which are pure basis (without noise) to the simulation
addBasisSamples <- function(simulation) {
  basis_samples <- generateBasisSamples(simulation$basis)
  simulation$basis_sample_names <- colnames(basis_samples$data)
  simulation$data <- cbind(simulation$data, basis_samples$data)
  simulation$proportions <- cbind(simulation$proportions, basis_samples$proportions)
  simulation
}

#' Add noise to the simulation's data
addNoise <- function(simulation, noise_deviation, protect_genes = c(), protect_samples = c()) {
  simulation$data <- applyNoiseToData(simulation$data, noise_deviation, protect_genes, protect_samples)
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  simulation$data_noise_deviation <- noise_deviation
  simulation
}


#' Create a simulation with noise on data, backwards compatibility wrapper
generateData <- function(genes, samples, cellTypes, noiseDeviation = 0) {
  simulation <- createSimulation(genes, samples, cellTypes)
  if (noiseDeviation > 0) {
    simulation <- addNoise(simulation, noiseDeviation)
  }
  simulation
}


# All underlying logic is below
generateBasis <- function(n_genes, n_cell_types, sd_ = 0.2) {
  data_1 <- c(
    rnorm(n_genes * 2, mean = 4, sd = 0.75),
    rnorm(n_genes * 3, mean = 10, sd = 1.5)
  )
  basis <- matrix(0, nrow = n_genes, ncol = n_cell_types)
  sds_ <- rnorm(n_cell_types, mean = sd_, sd = 0.02)

  for (i in 1:n_genes) {
    basis[i, 1] <- data_1[sample(seq_len(length(data_1)), 1)] * rnorm(1, mean = 1, sd = sds_[1])
    for (j in 2:n_cell_types) {
      basis[i, j] <- basis[i, 1] * rnorm(1, mean = 1, sd = sds_[j])
    }
  }

  basis <- limma::normalizeBetweenArrays(basis, method = "quantile")
  basis <- 2^basis
  rownames(basis) <- paste0("gene_", 1:n_genes)
  colnames(basis) <- paste0("cell_type_", 1:n_cell_types)

  basis
}

appendMarkerGenesToBasis <- function(basis) {
  n_cell_types <- ncol(basis)
  marker_expressions <- diag(1, n_cell_types, n_cell_types)
  rownames(marker_expressions) <- paste0("marker_gene_", 1:n_cell_types)
  list(
    basis = rbind(basis, marker_expressions),
    marker_gene_names = rownames(marker_expressions)
  )
}

#' Generate n points uniformly distributed on a k-dimensional standard simplex
sampleFromSimplexUniformly <- function(n, k, M = 100000) {
  X <- matrix(0, nrow = k + 1, ncol = n)
  X[k + 1,] <- M

  X[2:k,] <- replicate(n, sample(1:(M - 1), k - 1))
  X <- apply(X, 2, sort)
  Y <- (X - X[c(k + 1, 1:k),])[2:(k + 1),]
  return(Y / M)
}

generateProportions <- function(n_samples, n_cell_types) {
  proportions <- sampleFromSimplexUniformly(n_samples, n_cell_types)
  colnames(proportions) <- paste0("sample_", 1:n_samples)
  rownames(proportions) <- paste0("cell_type_", 1:n_cell_types)

  proportions
}

applyNoiseToData <- function(data, noise_deviation, protect_genes = c(), protect_samples = c()) {
  noise_mask <- matrix(rnorm(length(data), sd = noise_deviation), nrow = nrow(data), ncol = ncol(data))
  rownames(noise_mask) <- rownames(data)
  colnames(noise_mask) <- colnames(data)
  noise_mask <- 2^noise_mask
  noise_mask[protect_genes,] <- 0
  noise_mask[, protect_samples] <- 0
  noisy_data <- data + noise_mask
  noisy_data[noisy_data < 0] <- 0
  noisy_data
}

generatePureSamples <- function(basis, samples_per_cell_type, noise_deviation = 0.05) {
  data <- basis[, rep(seq_len(ncol(basis)), each = samples_per_cell_type)] * matrix(
    rnorm(nrow(basis) * ncol(basis) * samples_per_cell_type, mean = 1, sd = noise_deviation),
    nrow = nrow(basis),
    ncol = ncol(basis) * samples_per_cell_type
  )
  proportions <- diag(1, ncol(basis), ncol(basis))[rep(seq_len(ncol(basis)), each = samples_per_cell_type),]
  colnames(proportions) <- colnames(basis)
  rownames(proportions) <- paste0("pure_sample_", seq_len(nrow(proportions)))
  proportions <- t(proportions)
  colnames(data) <- colnames(proportions)
  list(proportions = proportions, data = data)
}

generateBasisSamples <- function(basis) {
  proportions <- diag(1, ncol(basis), ncol(basis))
  rownames(proportions) <- colnames(basis)
  colnames(proportions) <- paste0("basis_sample_", seq_len(ncol(proportions)))
  data <- basis
  colnames(data) <- colnames(proportions)
  list(proportions = proportions, data = data)
}
