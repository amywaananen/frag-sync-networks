simulateScene <- function(size = 30, meanSD = "2012-07-12", sdSD = 6, meanDur = 11,
                          sdDur = 3, skSD = 0, xRange = c(0, 100), yRange = c(0, 100),
                          distro = "unif", sAlleles = 10) {
  md <- as.integer(as.Date(meanSD, "%Y-%m-%d")) # mean start date of flowering
  sd <- as.integer(md + round(sn::rsn(
    n = size, 0, omega = sdSD, # sample start dates from skew normal distribution; omega = standard deviation of start date
    alpha = skSD
  ), 0)) # alpha is a the skew of the start date (defined by function of mean and standard deviation)
  ed <- as.integer(sd + abs(round(
    rnorm(size, meanDur, sdDur), # select end date by choosing a duration from a normal disribution with given mean and standard deviation
    0
  )))
  if (distro != "unif") {
    warning("distro must be unif")
  } # spatial distribution is uniform
  xv <- runif(size, min = xRange[1], max = xRange[2]) # adjust here to create clumped distribution?
  yv <- runif(size, min = yRange[1], max = yRange[2])
  sM <- sample(x = 1:sAlleles, size = size, replace = TRUE) # sample allele 1 for all individuals from all alleles with replacement
  if (sAlleles == 2) {
    sP <- 3 - sM # if there are onlt two alleles, the next allele is whatever the first one wasn't
  }
  else {
    sP <- sapply(sM, FUN = function(x) {
      sample(
        (1:sAlleles)[-x], # otherwise, sample the second allele from all alleles, excluding the one chosen for the first allele
        1
      )
    })
  }
  df <- data.frame(
    id = 1:size, start = sd, end = ed, x = xv, # put information together in a dataframe
    y = yv, s1 = sM, s2 = sP
  )
  makeScene(df,
    startCol = "start", endCol = "end", xCol = "x", # make scene using mateable makeScene function
    yCol = "y", idCol = "pla", dateFormat = "1970-01-01"
  )
}


# R help creating clumped vs random: https://stat.ethz.ch/pipermail/r-help/2006-May/106099.html
# install.packages('mobsim')

compatibility <- function(scene, method, subject = "all",
                          averageType = "mean") {
  method <- match.arg(method, c("si_echinacea", "dioecious"))
  subject <- match.arg(subject, c(
    "population", "pairwise",
    "individual", "all"
  ), several.ok = T)
  averageType <- match.arg(averageType, c("mean", "median"))

  if (averageType == "mean") {
    average <- mean
  } else if (averageType == "median") {
    average <- median
  }

  if (is.list(scene) & !is.data.frame(scene)) {
    potential <- lapply(scene, compatibility, method, subject, averageType)
  } else {
    if (method == "si_echinacea") {
      pairCompat <- pair_si_ech(scene$s1, scene$s2)
      attr(pairCompat, "idOrder") <- scene$id

      indCompat <- data.frame(id = scene$id, compatibility = -1)
      if (averageType == "mean") {
        indCompat$compatibility <- rowMeans(pairCompat)
      } else if (averageType == "median") {
        indCompat$compatibility <- row_medians(pairCompat)
      }

      popCompat <- average(indCompat$compatibility)
    } else if (method == "dioecious") {
      pairCompat <- pair_dioecious(scene$s1)
      attr(pairCompat, "idOrder") <- scene$id

      indCompat <- data.frame(id = scene$id, compatibility = -1)
      if (averageType == "mean") {
        indCompat$compatibility <- rowMeans(pairCompat)
      } else if (averageType == "median") {
        indCompat$compatibility <- row_medians(pairCompat)
      }

      popCompat <- average(indCompat$compatibility)
    }

    # return
    potential <- list()
    if ("population" %in% subject) {
      potential$pop <- popCompat
    }
    if ("individual" %in% subject) {
      potential$ind <- indCompat
    }
    if ("pairwise" %in% subject) {
      potential$pair <- pairCompat
    }
    if ("all" %in% subject) {
      potential$pop <- popCompat
      potential$ind <- indCompat
      potential$pair <- pairCompat
    }
    attr(potential, "t") <- FALSE
    attr(potential, "s") <- FALSE
    attr(potential, "c") <- TRUE
    potential
  }
}



makeScene <- function(df, multiYear = FALSE, startCol = "start", endCol = "end",
                      xCol = "x", yCol = "y", s1Col = "s1", s2Col = "s2",
                      idCol = "id", otherCols = NULL, dateFormat = "%Y-%m-%d",
                      split = NULL) {
  if (multiYear) {
    if (dateFormat == "%Y") {
      dates <- as.Date(as.character(df[, startCol]), dateFormat)
    } else {
      dates <- as.Date(df[, startCol], dateFormat)
    }
    df$year <- as.numeric(format(dates, "%Y"))
    years <- levels(as.factor(df$year))
    newScene <- list()
    for (i in 1:length(years)) {
      newScene[[as.character(years[i])]] <-
        makeScene(
          df[df$year %in% years[i], ], F, startCol, endCol, xCol, yCol,
          s1Col, s2Col, idCol, otherCols, dateFormat, split
        )
    }
  } else if (!is.null(split)) {
    splitTo <- levels(as.factor(df[, split]))
    newScene <- list()
    for (i in 1:length(splitTo)) {
      newScene[[as.character(splitTo[i])]] <-
        makeScene(
          df[df[, split] %in% splitTo[i], ], F, startCol, endCol, xCol, yCol,
          s1Col, s2Col, idCol, otherCols, dateFormat
        )
    }
  } else {
    newScene <- data.frame(id = character(nrow(df)))

    if (idCol %in% names(df)) {
      newScene$id <- df[, idCol]
    } else {
      newScene$id <- 1:nrow(df)
    }

    attr(newScene, "t") <- FALSE
    attr(newScene, "s") <- FALSE
    attr(newScene, "mt") <- FALSE
    attr(newScene, "originalNames") <- names(df)

    if (all(c(startCol, endCol) %in% names(df))) {
      attr(newScene, "t") <- TRUE
      newScene$start <- as.integer(as.Date(df[, startCol], dateFormat))
      firstDay <- min(newScene$start)
      newScene$start <- newScene$start - firstDay + 1
      newScene$end <- as.integer(as.Date(df[, endCol], dateFormat)) - firstDay + 1
      newScene$duration <- newScene$end - newScene$start + 1
      origin <- as.Date(firstDay - 1, "1970-01-01")

      attr(newScene, "origin") <- origin
    }

    if (all(c(xCol, yCol) %in% names(df))) {
      attr(newScene, "s") <- TRUE
      newScene$x <- df[, xCol]
      newScene$y <- df[, yCol]
    }

    if (all(c(s1Col, s2Col) %in% names(df))) {
      attr(newScene, "mt") <- TRUE
      newScene$s1 <- as.factor(df[, s1Col])
      newScene$s2 <- as.factor(df[, s2Col])
    }

    if (!is.null(otherCols)) {
      newScene[, otherCols] <- df[, otherCols]
    }
    # not going to add this for now because it's unlikely we'll make our
    # own generics or use oop
    # class(newScene) <- "matingScene"
  }
  newScene
}
