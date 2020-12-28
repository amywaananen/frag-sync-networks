#### MODEL RUN ####
library(igraph)
# library(mateable)
library(dplyr)
library(ggplot2)
library(sn)
# library(devtools)
# install_github("stuartWagenius/mateable@v0.3.1")
library(mateable)
# library(mobsim)
library(spatstat)

## CONSTANTS ####
N <- 100 # initial population size
SDiv <- 12 # S-allele diversity
MEAN_START_DATE <- 0
SD_START_DATE <- 15
SKEW_START_DATE <- 3
MEAN_DURATION <- 10
SD_DURATION <- 3
X_RANGE <- c(0, 1000)
Y_RANGE <- c(0, 1000)

## FUNCTIONS ####

# todo: update simulation to 
# (1) create populations at different levels of aggregation (3 levels - random, regular, poisson aggregated), 
# then (2) apply different levels of habitat amount and fragmentation

# 1a. Aggregated: We used the
# function sim_thomas_community from the R package mobsim to
# simulate species distributions following the Thomas process (May,
# Gerstner, McGlinn, Xiao, & Chase, 2018). We simulated a full factorial
# design of the parameter combinations nP in {10; 100; 1,000},
# nC in {1; 2; 5; 10}, and σ in {0.01, 0.02, 0.05, 0.1}.

# 1b. Random: For random distributions, we use
# the Poisson process, which assumes complete spatial randomness
# (CSR) without any interactions between individuals in the simulated
# arena. We simulated random distributions with a density of 10, 100,
# and 1,000 points in a square arena of 1 × 1 units.

# 1c. We simulated realizations
# of the Strauss process using the Metropolis–Hastings algorithm
# as implemented in the function rmh in the R package spatstat
# (Baddeley et al., 2015). Again, we conducted simulations for a full
# factorial design of the parameter values nP in {10; 100; 1,000} and
# γ in {0, 0.01, 0.1, 1}.
# install.packages('spatstat')

#2. We simulated fractal raster maps with defined habitat amount
# and fragmentation per se (Campos, Rosas, de Oliveira, & Gomes,
# 2013; Körner & Jeltsch, 2008; With, 1997; With, Gardner, &
# Turner, 1997). For this purpose, we used the midpoint displacement
# algorithm (Saupe, 1988) as implemented in the R package
# FieldSim.

# consider which levels of correlation are relevant - is it just space and time / space and genetics?


# Load functions

# function to create population coordinates by Poisson Thomas Process
makeThomasPop <- function(kappa, sigma, mu, N){
  pop <- rThomas(kappa, sigma, mu)
  dplyr::slice_sample(.data = coords(pop),n = N)
}
# makeThomasPop(15, 0.02, 5, 50)

source("makePop.R") # returns function runMakePop

sampleSAlleles <- function(ind, probdf = probsByX, varb = "x") {
  if (probdf == "uniform") {
    S1 <- sample(1:12, 1)
    S2 <- S1
    while (S2 == S1) {
      S2 <- sample(1:12, 1)
    }
  } else {
    S1 <- sample(1:12, 1, prob = probdf[, ind[varb]])
    S2 <- S1
    while (S2 == S1) {
      S2 <- sample(1:12, 1, prob = probdf[, ind[varb]])
    }
  }
  return(c(S1, S2))
}

netFun <- function(df, main = "", plot = F) {
  finallist <- list()

  nullPop <- df[, c("x", "y", "start1", "end1", "S1Null", "S2Null")]
  stCorPop <- df[, c("x", "y", "start2", "end2", "S1Null", "S2Null")]
  sgCorPop <- df[, c("x", "y", "start1", "end1", "S1SG", "S2SG")]
  tgCorPop <- df[, c("x", "y", "start1", "end1", "S1TG", "S2TG")]
  stgCorPop <- df[, c("x", "y", "start2", "end2", "S1STG", "S2STG")]

  pops <- list(nullPop, stCorPop, sgCorPop, tgCorPop, stgCorPop)
  popNames <- c("nullPop", "stCorPop", "sgCorPop", "tgCorPop", "stgCorPop")
  templist <- list(popType = character(), mods = numeric(), edge_density = numeric(), transitivity = numeric())

  for (i in 1:length(pops)) {
    dfi <- pops[[i]]
    colnames(dfi) <- c("x", "y", "s", "e", "s1", "s2")

    dfi$start <- as.Date("2019-07-04") + dfi$s
    dfi$end <- as.Date("2019-07-04") + dfi$e
    scene <- makeScene(dfi)
    sync <- synchrony(scene, method = "augspurger", subject = "pairwise")
    comp <- compatibility(scene, subject = "pair", method = "si_echinacea")$pair
    dist <- as.matrix(exp(dist(df[, c("x", "y")], upper = TRUE) * -.13)) # should omp decay with distance?
    # dist <- proximity(scene, method = 'maxProp', subject = 'pair')$pair
    omp <- sync * comp * dist
    gg <- graph_from_adjacency_matrix(omp, weighted = T, mode = "undirected")

    # Identify community clusters
    # modulos<-cluster_walktrap(gg)
    modulos <- fastgreedy.community(gg)

    if (plot) {
      dfi <- dfi %>%
        group_by(s) %>%
        mutate(e = 1:n())

      plot.igraph(gg,
        layout = as.matrix(dfi[, c("s", "e")]),
        vertex.label = NA, vertex.size = 0.1 * igraph::degree(gg),
        edge.width = (edge_attr(gg)$weight),
        arrow.mode = "-",
        edge.curved = 0.3, edge.color = "gray50",
        mark.groups = membership(modulos), mark.border = NA
      )
      mtext(main, line = -4, asp = F)
      mtext("Number of modules:", side = 1, line = -3, cex = 0.75)
      mtext(length(modulos), side = 1, line = -2, cex = 0.75)

      mtext("Edge density:", side = 1, line = -5.5, cex = 0.75)
      mtext(round(edge_density(gg), digits = 2), side = 1, line = -4.5, cex = 0.75)

      mtext("Transitivity:", side = 1, line = -8, cex = 0.75)
      mtext(round(transitivity(gg), digits = 2), side = 1, line = -7, cex = 0.75)
    }

    templist$popType[[i]] <- popNames[i]
    templist$mods[[i]] <- modularity(modulos)
    templist$edge_density[[i]] <- edge_density(gg)
    templist$transitivity[[i]] <- transitivity(gg)
  }
  templist
}


# pop <- runMakePop(15, 15, 15, 15, xCoords = ech2012$Ecoord, yCoords = ech2012$Ncoord)
# pop1 <- runMakePop(1, 1, 1, 1, xCoords = ech2012$Ecoord, yCoords = ech2012$Ncoord)
# pop15 <- runMakePop(100, 100, 100, 100, xCoords = ech2012$Ecoord, yCoords = ech2012$Ncoord)
# 
# ggplot(pop1, aes(x = x, y = start2)) +
#   geom_point() +
#   theme_classic()
# ggplot(pop15, aes(x = x, y = start2)) +
#   geom_point() +
#   theme_classic()
# ggplot(pop1, aes(x = x, y = y, fill = start2, col = start2)) +
#   geom_point()
# ggplot(pop15, aes(x = x, y = y, fill = start2, col = start2)) +
#   geom_point()

# out <- netFun(pop, plot = T)

## SIMULATE POPULATIONS ####


# We simulated a full factorial
# design of the parameter combinations nP in {10; 100; 1,000},
# nC in {1; 2; 5; 10}, and σ in {0.01, 0.02, 0.05, 0.1}.

# cluster parameters
clusterParams <- expand.grid(kappa = c(10, 100, 1000), mu = c(1,2,5,10), sigma = c(.01, .02, .05, .1))
head(clusterParams)
nrow(clusterParams) # 48 populations

allSd <- expand.grid(seq(0, 100, by = 50), seq(0, 100, by = 50), seq(0, 100, by = 50), seq(0, 100, by = 50))
# allSims <- expand.grid(seq(0,100, by = 50), seq(0,100, by = 50), seq(0,50, length.out = 3))


allParameters <- expand.grid(kappa = c(10, 100, 1000), mu = c(1,2,5,10), sigma = c(.01, .02, .05, .1), 
                             SD_ST = seq(0, 100, by = 50), SD_SG = seq(0, 100, by = 50), SD_TG = seq(0, 100, by = 50), SD_STG = seq(0, 100, by = 50))

syncLevels <- seq(0, 30, length.out = 3)
nboot <- 3

datalist <- list()

for (j in 1:nrow(allSd)) {
  SDs <- allSd[j, ] # loop through all combinations of desired SD values (0:5 for SG, ST, TG, and STG)
  for (k in 1:length(syncLevels)) {
    SD_START_DATE <- syncLevels[k]
    for (m in 1:nboot) { # this will make nboot populations at the given SD values
      pop <- runMakePop(SDs[[1]], SDs[[2]], SDs[[3]], SDs[[4]]) # make a population with SD values
      dat <- netFun(pop, plot = F)
      dat[["i"]] <- rep(m, 5)
      datalist[[length(datalist) + 1]] <- dat
    }
  }
}


dc <- function(list) {
  (
    do.call(cbind, list)
  )
}

# dc(datalist[[1]])

test <- lapply(datalist, dc)
test2 <- do.call(rbind, test)
test2 <- as.data.frame(test2)
sdReps <- allSd %>% slice(rep(1:n(), each = nboot * 5 * length(syncLevels)))
sdReps$sd_start <- rep(rep(syncLevels, each = 5 * nboot), times = nrow(allSd))

colnames(sdReps) <- c("SD_ST", "SD_SG", "SD_TG", "SD_STG", "sd_start")
test2$i <- as.numeric(test2$i)
test3 <- cbind(test2, sdReps)
test3$mods <- as.numeric(test3$mods)
test3$edge_density <- as.numeric(test3$edge_density)
test3$transitivity <- as.numeric(test3$transitivity)

big_data <- test3

ggplot(big_data, aes(x = mods)) +
  geom_histogram(color = "black", fill = "gray", bins = 35) +
  theme_bw()

ggplot(big_data, aes(x = edge_density)) +
  geom_histogram(color = "black", fill = "gray", bins = 35) +
  theme_bw() +
  xlab("Edge Density")

ggplot(big_data, aes(x = transitivity)) +
  geom_histogram(color = "black", fill = "gray", bins = 35) +
  theme_bw()

# table(big_data$mods)
# table(big_data$edge_density)
# table(big_data$transitivity)

nrow(big_data)


ggplot(big_data, aes(x = transitivity, group = popType, fill = popType)) +
  geom_density(alpha = 0.5)

ggplot(big_data, aes(x = transitivity, group = popType, fill = popType)) +
  geom_histogram(alpha = 0.5)

ggplot(big_data[big_data$popType %in% "stgCorPop", ], aes(y = transitivity, x = SD_STG)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(big_data[big_data$popType %in% "stCorPop", ], aes(y = transitivity, x = SD_ST)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(big_data[big_data$popType %in% "stCorPop", ], aes(y = mods, x = SD_ST)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(big_data, aes(y = transitivity, x = SD_STG)) +
  facet_wrap(vars(popType)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(big_data, aes(y = mods, x = SD_STG)) +
  facet_wrap(vars(popType)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(big_data, aes(y = edge_density, x = SD_STG)) +
  # facet_wrap(vars(popType))+
  geom_point() +
  stat_smooth(aes(color = popType), method = "lm") +
  theme_bw()

ggplot(big_data, aes(y = edge_density, x = SD_SG)) +
  # facet_wrap(vars(popType))+
  geom_point() +
  stat_smooth(aes(color = popType), method = "lm") +
  theme_bw()

ggplot(big_data, aes(y = edge_density, x = SD_ST)) +
  # facet_wrap(vars(popType))+
  geom_point() +
  stat_smooth(aes(color = popType), method = "lm") +
  theme_bw()

big_data$SD_ST <- factor(big_data$SD_ST, levels = c(100, 50, 0))

# Modularity
ggplot(big_data, aes(y = mods, x = sd_start)) +
  facet_grid(rows = vars(SD_ST), cols = vars(SD_SG)) +
  geom_jitter(size = 0.5, width = 0.2, height = 0) +
  stat_smooth(aes(color = popType), method = "lm", se = FALSE) +
  theme_bw() +
  xlab("Standard Deviation of Start Date") +
  ylab("Modularity")

# Edge density
ggplot(big_data, aes(y = edge_density, x = sd_start)) +
  facet_grid(rows = vars(SD_ST), cols = vars(SD_SG)) +
  geom_jitter(size = 0.5, width = 0.2, height = 0) +
  stat_smooth(aes(color = popType), method = "lm", se = FALSE) +
  theme_bw() +
  xlab("Standard Deviation of Start Date") +
  ylab("Edge Density")

# Transitivity
ggplot(big_data, aes(y = transitivity, x = sd_start)) +
  facet_grid(rows = vars(SD_ST), cols = vars(SD_SG)) +
  geom_jitter(size = 0.5, width = 0.2, height = 0) +
  stat_smooth(aes(color = popType), method = "lm", se = FALSE) +
  theme_bw() +
  xlab("Standard Deviation of Start Date") +
  ylab("Transitivity")
