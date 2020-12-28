runMakePop <- function(SD_ST, # standard deviation of distribution for sampling timing by space
                       SD_SG, # standard deviation of distribution for sampling S-alleles by space
                       SD_TG, # standard deviation of distribution for sampling S-alleles by time
                       SD_STG, # standard deviation of distribution for sampling S-alleles by time,
                       kappa, # parameter for Thomas poisson process
                       sigma, # parameter for Thomas poisson process
                       mu, # parameter for Thomas poisson process
                       xCoords, # pre-supplied spatial coordinates
                       yCoords # pre-supplied spatial coordinates
                       ) {
  
  #### STEP 0: Assign spatial coordinates ####
  if (!any(is.na(xCoords))) {
    N <- length(xCoords)
    pop <- data.frame(id = 1:N, x = xCoords, y = yCoords) # if coordinates are supplied in xCoord and yCoord, use these
  } else if(!any(!is.na(kappa))){
    pop <- makeThomasPop(kappa, sigma, mu, N)
    pop$id <- 1:N
  } else { # otherwise, assign coordinates randomly
    pop <- data.frame(id = 1:N, x = sample(X_RANGE[1]:X_RANGE[2], N, replace = TRUE), y = sample(Y_RANGE[1]:Y_RANGE[2], N, replace = TRUE)) # Step 0
  }

  #### STEP 1: Assign start date and duration randomly ####
  pop$start1 <- as.integer(round(sn::rsn(n = N, 0, omega = SD_START_DATE, alpha = SKEW_START_DATE), 0)) # assign start dates randomly from skew normal distribution with mean of 0
  pop$start1 <- pop$start1 + abs(min(pop$start1)) + 1 # shift start dates by absolute value of the minimum plus one so that all start dates are positive
  pop$duration <- round(rnorm(n = N, mean = MEAN_DURATION, sd = SD_DURATION)) # sample duration from normal distribution with given mean and standard deviation
  pop$end1 <- pop$start1 + pop$duration # add duration to start date to get end date

  #### STEP 2: Assign start date and duration weighted by location ####

  # This loop reassigns the start dates assigned randomly, though individuals maintain their duration.

  startRanked <- sort(pop$start1) # make an ordered list of the randomly assigned start dates
  pop[order(pop$x), "xRank"] <- 1:N # add a column to show the rank of individuals in terms of their x-position
  pop$start2ind <- NA # add a  column for indexing from the startRanked vector

  for (a in 1:nrow(pop)) {
    sample <- round(rnorm(1, mean = pop[a, "xRank"], sd = SD_ST)) # draw a number from a normal distribution centered around the individuals xRank with given standard deviation
    availableIndices <- seq(1, N)[!seq(1, N) %in% unique(pop$start2ind)] # identify indices that haven't been assigned yet
    diffFromAvailableIndices <- abs(sample - availableIndices) # calculate the difference between the sampled number and the available indices
    pop[a, "start2ind"] <- availableIndices[which(diffFromAvailableIndices %in% min(diffFromAvailableIndices))[1]] # assign the individual the index that is the least different from the sampled number
    pop[a, "start2"] <- startRanked[pop[a, "start2ind"]] # get the start date corresponding to the start index and assign to column 'start2'
  }

  pop$end2 <- pop$start2 + pop$duration # assign end date 'end2' by adding duration to start date

  # # visualize spatiotemporal gradient
  # par(mfrow = c(2,1))
  # plot(pop$x, pop$start1)
  # plot(pop$x, pop$start2)

  #### STEP 3: Assign S-alleles randomly ####

  # Use function sampleSAlleles to randomly assign S alleles from pool of S alleles
  # The number of S-alleles is determined by SDiv under 'Model Initializations'
  # This (should) also influence the mean compatibility rate of the population,
  # which the population size will also affect.

  pop[, c("S1Null", "S2Null")] <- t(apply(pop, 1, sampleSAlleles, probdf = "uniform")) # assign S1 and S2 alleles, labeled S1/S2 Null

  #### STEP 4: Assign S-alleles weighted by x-coordinate ####

  # General approach here is  to rank s-alleles and x-coordinates, then match them up,
  # similar to how flowering start dates were assigned when spatially autocorrelated
  # and skipping to the next available allele when S1 and S2 are the same.

  S1alleles <- sort(pop$S1Null) # make a sorted list of S1 alleles
  S2alleles <- sort(pop$S2Null) # maek a sorted list of S2 alleles
  pop$S1SG <- NA
  pop$S2SG <- NA

  for (b in 1:nrow(pop)) {
    sample <- round(rnorm(1, mean = pop[b, "xRank"], sd = SD_SG)) # draw a sample from normal distribution centered at xRank
    availableIndices <- seq(1, N)[!seq(1, N) %in% unique(pop$S1ind)] # identify available indices by excluding indices that are already assigned
    diffFromAvailableIndices <- abs(sample - availableIndices)
    pop[b, "S1ind"] <- availableIndices[which(diffFromAvailableIndices %in% min(diffFromAvailableIndices))[1]] # see above, find minimally different from sampled xRank
    pop[b, "S1SG"] <- S1alleles[pop[b, "S1ind"]] # assign S1 allele based on index

    while (is.na(pop[b, "S2SG"]) | pop[b, "S2SG"] == pop[b, "S1SG"]) {
      sample <- round(rnorm(1, mean = pop[b, "xRank"], sd = SD_SG)) # Draw sample based on xRank for S2
      if (length(intersect(seq(1, N)[!seq(1, N) %in% unique(pop$S2ind)], seq(1, N)[!S2alleles %in% pop[b, "S1SG"]])) == 0) { # this if clause is to deal with issues where the all remaining S2 alleles are the same as the assigned S1 allele, leaving no options for S2
        options <- pop[!is.na(pop$S1SG) & !pop$S1SG %in% pop[b, "S1SG"] & !pop$S2SG %in% pop[b, "S1SG"], "xRank"] # find options to reassign S1 allele based on xRanks that have been assigned an S1 allele that is not the same as the S1 allele currently being dealt with; also need to be sure that the S2 allele doesn't match the current S1 allele so that swapping doesn't create new issues
        swapRank <- options[which.min(abs(options - pop[b, "xRank"]))] # pick the allele to swap the S1 allele with based on minimizing the difference in xRank (i.e., swap with the nearest individual that has been assigned with a different S1 allele)
        swap1 <- pop[b, "S1SG"] # save the original S1 allele in index i for the swap
        pop[b, "S1SG"] <- pop[pop$xRank %in% swapRank, "S1SG"] # do the swap, give the swap rank allele to index i
        pop[pop$xRank %in% swapRank, "S1SG"] <- swap1 # and give the swapped individual the original S1 allele that was in index i
      }
      availableIndices <- intersect(seq(1, N)[!seq(1, N) %in% unique(pop$S2ind)], seq(1, N)[!S2alleles %in% pop[b, "S1SG"]])
      diffFromAvailableIndices <- abs(sample - availableIndices)
      pop[b, "S2ind"] <- availableIndices[which(diffFromAvailableIndices %in% min(diffFromAvailableIndices))[1]]
      pop[b, "S2SG"] <- S2alleles[pop[b, "S2ind"]]
    }
  }

  # Do some tests to make sure this works:
  # # 1. check that all alleles are the same in null & spatial assignment
  # all(sort(pop$S1Null) == sort(pop$S1SG))
  # all(sort(pop$S2Null) == sort(pop$S2SG))
  #
  # # 2. check that there aren't any cases where s1 = s2
  # nrow(pop[pop$S1Null == pop$S2Null,])
  # nrow(pop[pop$S1SG == pop$S2SG,])
  #
  # # 3. check that there's a spatial pattern of allele assignment
  # par(mfrow = c(2,2))
  # plot(pop$x, pop$S1Null)
  # plot(pop$xRank, pop$S2Null)
  # plot(pop$x, pop$S1SG)
  # plot(pop$x, pop$S2SG)

  #### STEP 4: Assign S-alleles weighted by timing of flowering where timing of flowering is not spatially autocorrelated ####

  # Same loop as before, but varying: SD_XX variable, column to subset for rank, column to assign S1 and S2 alleles

  pop$S1ind <- NA
  pop$S2ind <- NA
  pop$S1TG <- NA
  pop$S2TG <- NA
  pop[order(pop$start1), "start1Rank"] <- 1:N

  for (c in 1:nrow(pop)) {
    sample <- round(rnorm(1, mean = pop[c, "start1Rank"], sd = SD_TG))
    availableIndices <- seq(1, N)[!seq(1, N) %in% unique(pop$S1ind)]
    diffFromAvailableIndices <- abs(sample - availableIndices)
    pop[c, "S1ind"] <- availableIndices[which(diffFromAvailableIndices %in% min(diffFromAvailableIndices))[1]]
    pop[c, "S1TG"] <- S1alleles[pop[c, "S1ind"]]

    while (is.na(pop[c, "S2TG"]) | pop[c, "S2TG"] == pop[c, "S1TG"]) {
      sample <- round(rnorm(1, mean = pop[c, "start1Rank"], sd = SD_TG))

      if (length(intersect(
        seq(1, N)[!seq(1, N) %in% unique(pop$S2ind)], # available indices
        seq(1, N)[!S2alleles %in% pop[c, "S1TG"]]
      )) == 0) { # compared to indices that match the current S1 allele

        # this if clause is to deal with issues where the all remaining S2 alleles
        # are the same as the assigned S1 allele, leaving no options for S2

        options <- pop[!is.na(pop$S1TG) & !pop$S1TG %in% pop[c, "S1TG"] & !pop$S2TG %in% pop[c, "S1TG"], "start1Rank"] # find options to reassign S1 allele based on start1Ranks that have been assigned an S1 allele that is not the same as the S1 allele currently being dealt with; also need to be sure that the S2 allele doesn't match the current S1 allele so that swapping doesn't create new issues
        swapRank <- options[which.min(abs(options - pop[c, "start1Rank"]))] # pick the allele to swap the S1 allele with based on minimizing the difference in start1Rank (i.e., swap with the nearest individual that has been assigned with a different S1 allele)
        swap1 <- pop[c, "S1TG"] # save the original S1 allele in index i for the swap
        pop[c, "S1TG"] <- pop[pop$start1Rank %in% swapRank, "S1TG"] # do the swap, give the swap rank allele to index i
        pop[pop$start1Rank %in% swapRank, "S1TG"] <- swap1 # and give the swapped individual the original S1 allele that was in index i
      }
      availableIndices <- intersect(seq(1, N)[!seq(1, N) %in% unique(pop$S2ind)], seq(1, N)[!S2alleles %in% pop[c, "S1TG"]])
      diffFromAvailableIndices <- abs(sample - availableIndices)
      pop[c, "S2ind"] <- availableIndices[which(diffFromAvailableIndices %in% min(diffFromAvailableIndices))[1]]
      pop[c, "S2TG"] <- S2alleles[pop[c, "S2ind"]]
    }
  }

  # # 1. check that all alleles are the same in null & spatial assignment
  # all(sort(pop$S1Null) == sort(pop$S1TG))
  # all(sort(pop$S2Null) == sort(pop$S2TG))
  #
  # # 2. check that there aren't any cases where s1 = s2
  # nrow(pop[pop$S1Null == pop$S2Null,])
  # nrow(pop[pop$S1TG == pop$S2TG,])
  #
  # # 3. check that there's a spatial pattern of allele assignment
  # par(mfrow = c(2,2))
  # plot(pop$start1, pop$S1Null)
  # plot(pop$start1Rank, pop$S2Null)
  # plot(pop$start1, pop$S1TG)
  # plot(pop$start1Rank, pop$S2TG)

  #### STEP 5: Assign S-Alleles weighted by time of flowering where timing is spatially autocorrelated ####
  pop$S1ind <- NA
  pop$S2ind <- NA
  pop$S1STG <- NA
  pop$S2STG <- NA
  pop[order(pop$start2), "start2Rank"] <- 1:N

  # Same loop as before, but varying: SD_XX variable, column to subset for rank, column to assign S1 and S2 alleles

  for (d in 1:nrow(pop)) {
    sample <- round(rnorm(1, mean = pop[d, "start2Rank"], sd = SD_STG))
    availableIndices <- seq(1, N)[!seq(1, N) %in% unique(pop$S1ind)]
    diffFromAvailableIndices <- abs(sample - availableIndices)
    pop[d, "S1ind"] <- availableIndices[which(diffFromAvailableIndices %in% min(diffFromAvailableIndices))[1]]
    pop[d, "S1STG"] <- S1alleles[pop[d, "S1ind"]]
    while (is.na(pop[d, "S2STG"]) | pop[d, "S2STG"] == pop[d, "S1STG"]) {
      sample <- round(rnorm(1, mean = pop[d, "start2Rank"], sd = SD_STG))
      if (length(intersect(seq(1, N)[!seq(1, N) %in% unique(pop$S2ind)], seq(1, N)[!S2alleles %in% pop[d, "S1STG"]])) == 0) { # this if clause is to deal with issues where the all remaining S2 alleles are the same as the assigned S1 allele, leaving no options for S2
        options <- pop[!is.na(pop$S1STG) & !pop$S1STG %in% pop[d, "S1STG"] & !pop$S2STG %in% pop[d, "S1STG"], "start2Rank"] # find options to reassign S1 allele based on start2Ranks that have been assigned an S1 allele that is not the same as the S1 allele currently being dealt with; also need to be sure that the S2 allele doesn't match the current S1 allele so that swapping doesn't create new issues
        # print(options)
        swapRank <- options[which.min(abs(options - pop[d, "start2Rank"]))] # pick the allele to swap the S1 allele with based on minimizing the difference in start2Rank (i.e., swap with the nearest individual that has been assigned with a different S1 allele)
        swap1 <- pop[d, "S1STG"] # save the original S1 allele in index i for the swap
        pop[d, "S1STG"] <- pop[pop$start2Rank %in% swapRank, "S1STG"] # do the swap, give the swap rank allele to index i
        pop[pop$start2Rank %in% swapRank, "S1STG"] <- swap1 # and give the swapped individual the original S1 allele that was in index i
      }
      availableIndices <- intersect(seq(1, N)[!seq(1, N) %in% unique(pop$S2ind)], seq(1, N)[!S2alleles %in% pop[d, "S1STG"]])
      diffFromAvailableIndices <- abs(sample - availableIndices)
      # print(diffFromAvailableIndices)
      pop[d, "S2ind"] <- availableIndices[which(diffFromAvailableIndices %in% min(diffFromAvailableIndices))[1]]
      pop[d, "S2STG"] <- S2alleles[pop[d, "S2ind"]]
    }
  }

  # output from this script is pop

  # colnames(pop)
  pop <- pop[, c(
    "id", "x", "y",
    "duration",
    "start1", "end1",
    "start2", "end2",
    "S1Null", "S2Null",
    "S1SG", "S2SG",
    "S1TG", "S2TG",
    "S1STG", "S2STG"
  )]
}
