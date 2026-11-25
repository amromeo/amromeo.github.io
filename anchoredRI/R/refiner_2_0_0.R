# SPDX-License-Identifier: GPL-3.0-or-later
# refineR 2.0.0 by Matthias Beck, Tatjana Ammer, Christopher M Rank and
# Andre Schuetzenmeister. See inst/NOTICE and inst/licenses/refineR-GPL-3.txt.
# Source: anchorRefine's Stage 0 vendor/refineR at revision
# 55b9d09ba7e17b13e341e67b0aa03855351bc25b.
# The search and endpoint functions below retain the upstream calculations.
# A private environment isolates the original names and their calls.
.refineR_200 <- local({
  # From utils.R (refineR 2.0.0).

  # Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
  BoxCox <- function(x, lambda) {

      if(abs(lambda) < 1e-20)
          x <- log(x)
      else
          x <- (x^lambda-1) / lambda
      x
  }

  # Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
  invBoxCox <- function(x, lambda) {

      if(abs(lambda) < 1e-20)
          x <- exp(x)
      else
          x <- (lambda*x+1)^(1/lambda) # if(lambda*x+1) is negative & lambda < 1-> result is NaN
      x
  }

  # Christopher Rank \email{christopher.rank@@roche.com}
  pnormApprox <- function(q, pNormVal, mean = 0, oneOverSd = 1, oneOverH = 10) {

      index <- ((q-mean)*oneOverSd + 5.2)*oneOverH
      index[index < 0] <- 0
      index[index >= (length(pNormVal)-2)] <- length(pNormVal)-1.5

      indexInt <- as.integer(index)

      w <- index-indexInt

      res <- (1-w)*pNormVal[indexInt+1] + w*pNormVal[indexInt+2]

      return(res)
  }

  # Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}
  ashDensity <- function(x, ab, nbin, m, kopt = c(2, 1), normToAB = FALSE) {

      invisible(capture.output(d <- ash1(bins = bin1(x = x, ab = ab, nbin = nbin), m = m, kopt = kopt)))

      if(normToAB)
          PCorr <- 1.0
      else
          PCorr <- sum(x >= ab[1] & x < ab[2])/length(x)

      list(x = d$x, y = d$y*PCorr, m = m, ab = ab, nbin = nbin, PCorr = PCorr)
  }

  # Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}
  findRoundingBase <- function(x) {

      ab <- as.numeric(quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

  # generate table with difference between neighboring values
      diffVal  <- table(round(diff(sort(unique(x[x>=ab[1] & x<=ab[2]]))), digits=10))
  # remove number of those that have a difference smaller than e-10
      diffVal <- diffVal[names(diffVal) != "0"]

  # select differences that occur in >= 1% of cases
      diffVal  <- diffVal[diffVal>=0.01*sum(diffVal)]

  # find multiples of rounding base
      if(length(diffVal) > 1)
      {
          for(i in 1:length(diffVal))
          {
              if(diffVal[i] > 0)
              {
                  ratio <- as.numeric(names(diffVal))/as.numeric(names(diffVal[i]))

                  selection <- which(ratio>1.5 & abs(round(ratio, digits=9)-round(ratio))<1e-20)

                  diffVal[i] <- diffVal[i] + sum(diffVal[selection])

                  diffVal[selection] <- 0
              }
          }
      }

  # select differences that occur in >= 10% of cases
      diffVal  <- diffVal[diffVal>=0.1*sum(diffVal)]

      roundingBase <- NA

  # if dataset has only finite number of unique values (e.g. rounded data)
      if (length(diffVal) > 0) {

  # determine rounding base (maximum step size of discrete values)
          roundingBase <- max(as.numeric(names(diffVal)))
      }

      return(roundingBase)
  }

  # Matthias Beck \email{matthias.beck.mb1@@roche.com}, Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}
  cdfTruncatedBoxCox <- function(mu, sigma, lambda, shift, RIperc){

      stopifnot(is.numeric(RIperc), min(RIperc)>=0, max(RIperc)<=1)
      stopifnot(is.numeric(mu), is.numeric(sigma), is.numeric(lambda), is.numeric(shift))
      if(any(!is.finite(c(mu, sigma, lambda, shift))) || sigma <= 0)
      {
          return(rep(NA, length(RIperc)))
      }

      RIperc <- sort(RIperc)

  # formula for truncated normal distribution
      if(lambda >= 0) {
          RI <- pnorm(-1/lambda, mean=mu, sd=sigma) + RIperc*(1 - pnorm(-1/lambda, mean=mu, sd=sigma))
          RI <- qnorm(RI, mean=mu, sd=sigma)
      } else {
  # for lambda < 0
          RI <- qnorm(RIperc, mean=mu, sd=sigma)
      }

      RI <- invBoxCox(RI, lambda = lambda)
      RI <- RI + shift
      RI[RI<0 | is.na(RI)] <- 0
      return(RI) # remove negative values
  }
  # From algoInvHist.R (refineR 2.0.0).

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  generateHistData <- function(x, ab, roundingBase) {

      HistData <- list()
      HistData$counts <- HistData$breakL <- HistData$breakR <- NULL

  # calculate number of bins (without overlapping) depending on size of data
      NBins <- min(40, max(9, 11*(sum(x >=ab[1] & x <= ab[2])/5000)^(1/5)))

  # calculate number of bins and number of total bins (with overlapping)
      overlapFactor <- round(256/NBins)

      NBinsTotal    <- round(NBins)*overlapFactor

  # if dataset has only finite number of unique values (e.g. rounded data)
      if(!is.na(roundingBase))
      {
  # determine bin size, which is a multiple of rounding base
          abStepSize <- ceiling(diff(ab)/256/roundingBase)*roundingBase

  # determine overlap factor
          overlapFactor <- max(1, round(diff(ab)/abStepSize/NBins))

  # calculate number of unique values to adjust NBins
          xAB <- unique(x[x>=ab[1] & x<=ab[2]])

          if(length(xAB) < NBins){
              NBins <- ceiling(length(xAB)/overlapFactor)
          }

  # calculate number of bins and number of total bins (with overlapping)
          NBinsTotal <- round(NBins)*overlapFactor

  # adapt ab to discrete grid
          ab[1] <- max(0.5*roundingBase, round(ab[1]/roundingBase)*roundingBase - 0.5*roundingBase)
          ab[2] <- ab[1] + (NBinsTotal+overlapFactor-1)*abStepSize
      }

  # define vector of breaks and perform histogramming
      histBreaks <- seq(ab[1], ab[2], length.out = NBinsTotal+overlapFactor)
      fullHist   <- hist(x[x > histBreaks[1] & x <= histBreaks[NBinsTotal+overlapFactor]], breaks = histBreaks, plot = FALSE)

  # combine bins...
      for (i in 1:NBinsTotal) {
          HistData$counts <- c(HistData$counts, sum(fullHist$counts[i:(i+overlapFactor-1)]))
          HistData$breakL <- c(HistData$breakL, fullHist$breaks[i])
          HistData$breakR <- c(HistData$breakR, fullHist$breaks[i+overlapFactor])
      }

  # add bin for lower border
      if(min(x) <= ab[1])
      {
  # -1e-20 is used to include zeros although breakL is >= 1e-20
          HistData$counts <- c(HistData$counts, sum(x > -1e-20 & x <= ab[1]))
          HistData$breakL <- c(HistData$breakL, max(min(x), 1e-20))
          HistData$breakR <- c(HistData$breakR, ab[1])
          fullHist$counts <- c(sum(x > -1e-20 & x <= ab[1]), fullHist$counts)
          fullHist$breaks <- c(max(min(x), 1e-20), fullHist$breaks)

      } else
      {
          HistData$counts <- c(HistData$counts, 0)
          HistData$breakL <- c(HistData$breakL, ab[1])
          HistData$breakR <- c(HistData$breakR, ab[1])
          fullHist$counts <- c(0, fullHist$counts)
          fullHist$breaks <- c(ab[1], fullHist$breaks)
      }

  # add bin for upper border
      if(max(x) > ab[2])
      {
          HistData$counts <- c(HistData$counts, sum(x > ab[2] & x <= 1e20))
          HistData$breakL <- c(HistData$breakL, ab[2])
          HistData$breakR <- c(HistData$breakR, min(max(x), 1e20))
          fullHist$counts <- c(fullHist$counts, sum(x > ab[2] & x <= 1e20))
          fullHist$breaks <- c(fullHist$breaks, min(max(x), 1e20))

      } else
      {
          HistData$counts <- c(HistData$counts, 0)
          HistData$breakL <- c(HistData$breakL, ab[2])
          HistData$breakR <- c(HistData$breakR, ab[2])
          fullHist$counts <- c(fullHist$counts, 0)
          fullHist$breaks <- c(fullHist$breaks, ab[2])
      }

      fullHist$mids   <- 0.5*(fullHist$breaks[2:length(fullHist$breaks)] + fullHist$breaks[1:(length(fullHist$breaks)-1)])

  # calculate mids of histogram bins
      HistData$mids   <- 0.5*(HistData$breakL + HistData$breakR)

  # sort bins in order of midpoints
      sortIndex <- sort(HistData$mids, index.return = TRUE)$ix
      HistData$counts <- HistData$counts[sortIndex]
      HistData$mids   <- HistData$mids[sortIndex]
      HistData$breakL <- HistData$breakL[sortIndex]
      HistData$breakR <- HistData$breakR[sortIndex]

  # add further parameters to list
      HistData$overlapFactor <- overlapFactor
      HistData$NBins             <- round(NBins)
      HistData$abOr           <- ab
      HistData$NData            <- length(x)

      HistData$fullHist <- fullHist

      return(HistData)
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  defineSearchRegions <- function(x, lambdaVec, roundingBase, abEst = NULL) {

      if (length(x) > 100000) {
  # add jitter to the data in case of rounded data
          if (!is.na(roundingBase))
              x <- x + runif(n=length(x), min =-0.5*roundingBase, max=0.5*roundingBase-1e-20)

  # calculate 50000 quantiles type=7 for data (to get better performance)
          qData <- as.numeric(quantile(x = x, probs = seq(0, 1, length.out = 50000), type = 7))

      } else {
  # duplicate dataset and add jitter to the data in case of rounded data
          if (!is.na(roundingBase)) {

              x <- rep(x, times=floor(100000/length(x)))
              x <- x + runif(n=length(x), min =-0.5*roundingBase, max=0.5*roundingBase-1e-20)
              qData <- x

          }else {
                qData <- as.numeric(quantile(x = x, probs = seq(0, 1, length.out = 50000), type = 7))
          }
      }

  # only call estimate ab in first definition of search regions
      if(is.null(abEst))
          abEst <- estimateAB(x=qData)# restrict data for estimation of start values to region around main peak

      qData <- qData[qData >= abEst$abSearchReg[1] & qData <= abEst$abSearchReg[2]]

      muEst      <- matrix(nrow = length(lambdaVec), ncol = 3)
      sigmaEst <- matrix(nrow = length(lambdaVec), ncol = 3)

  #determine m for density estimation
      densM <- min(48, max(24, ceiling(50 - (1/10)*length(x)/1000)))

      c <- 1

  # loop over lambda
      for (lambda in lambdaVec) {

  # BoxCox transformation and density estimation...
  # apply BoxCox transformation to quantiles/data and remove NAs an infinity values
          qDataTr <- suppressWarnings(BoxCox(x = qData, lambda = lambda))
          qDataTr <- qDataTr[!is.na(qDataTr) & is.finite(qDataTr)]

  #get range of transformed data for density estimation
          abTr     <- range(qDataTr)
          abTrDiff <- diff(abTr)
          abTr[1]  <- abTr[1]-0.1*abTrDiff
          abTr[2]  <- abTr[2]+0.1*abTrDiff

  # find main peak of distribution
          mainPeak <- findMainPeak(x=qDataTr, ab=abTr, mStart=densM)
          peakInd  <- mainPeak$peakInd
          modEst   <- mainPeak$modEst
          qDensTr  <- mainPeak$Dens

  # find width of peak for densities 0.5 to 0.9 of max(density) in 0.01 steps
          steps       <- seq(0.5, 0.95, by = 0.05)
          DensRange <- steps*qDensTr$y[peakInd]

          widthsL <- NULL
          widthsR <- NULL

          finalSteps <- NULL
          for (i in 1:length(DensRange)) {
              d <- DensRange[i]
              if(length(which(qDensTr$x < modEst & qDensTr$y <= d)) != 0 & (length(which(qDensTr$x > modEst & qDensTr$y <= d)) != 0)){
  # get width (distance) from peak to the left to the first points that is smaller or equal to d
                  widthsL <- c(widthsL, modEst - max(qDensTr$x[which(qDensTr$x < modEst & qDensTr$y <= d)]))
  # get width (distance) from peak to the right to the first point that is smaller or equal to d
                  widthsR <- c(widthsR, (min(qDensTr$x[which(qDensTr$x > modEst & qDensTr$y <= d)]) - modEst))
                  finalSteps <- c(finalSteps, steps[i])
              }
          }
  # get points left and right of the peak
          left     <- modEst - widthsL
          right     <- modEst + widthsR

  # leftRight Matrix: cols: 1=density range, 2=left value, 3=right value, 4=widthsL, 5=widthsR, 6=mean of left and right --> mu
          leftRight     <- matrix(c(finalSteps, left, right, widthsL, widthsR), nrow = length(left), ncol = 5, byrow = FALSE)
          leftRight     <- cbind(leftRight, apply(leftRight[, 2:3], 1, mean))

          muTmp     <- c(leftRight[, 6], modEst)

  # estimate range of muEst as range of centers of c(widths, modEst)
          muRange     <- range(muTmp)
  #add margins to muRange
          diffMuRange    <- diff(muRange)
          muEstRange    <- c(muRange[1] - 0.5*diffMuRange, muRange[2] + 0.5*diffMuRange)
          muEst[c,]     <- c(lambda, muEstRange)

  # estimate range of sigmaEst as range of sds calculated from widths
  # minimal width (as width is highly influenced by pathological samples), denominator is derived from HWHM formula
          leftRight <- cbind(leftRight, apply(leftRight[, c(1, 4:5)], 1, function(x) {
                              return(min(x[c(2, 3)])/(sqrt((-2)*log(x[1]))))
                          }))
          sigmaTmp  <- leftRight[, 7]

  # add margins to sigmaRange on both sides
          sigmaRange         <- range(sigmaTmp)
          diffSigmaRange     <- diff(sigmaRange)
          sigmaEstRange    <- c(max(0, sigmaRange[1] - 2.5*diffSigmaRange), sigmaRange[2]+1*diffSigmaRange)
          sigmaEst[c,]     <- c(lambda, sigmaEstRange)

          c = c + 1
      }

      return(list("abEst" = abEst, "muEst" = muEst, "sigmaEst" = sigmaEst))
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  testParam <- function(lambdaVec, bestParam, Data, HistData, startValues, NIter, alpha = 0.01, alphaMcb = 0.1) {

      h <- 0.002
      pNormLookup <- list(h = h, val = pnorm(q = seq(from = -5.2, to = 5.2+h, by = h)))

      c <- 1
  # loop over all lambda values
      for (lambda in lambdaVec) {
  # generate mu and sigma vectors
          muVec <- startValues$muEst[c, 2:3]
          muVec <- seq(muVec[1], muVec[2], length.out = 9)

          sigmaVec <- startValues$sigmaEst[c, 2:3]
          sigmaVec <- seq(sigmaVec[1], sigmaVec[2], length.out = 13)

          muUnique         <- sigmaUnique <- NULL
          currentBestMu   <- currentBestSigma <- NA
          currentBestCost <- 1e9

  # loop for optimization using adapted grid search
          iter = 1
          continueIteration = TRUE
          while (continueIteration & iter <= NIter) {

              testedParam <- c(NA, NA, NA, NA, 1e9, NA)
              muUnique    <- sort(unique(c(muUnique,    muVec)),    method="quick")
              sigmaUnique <- sort(unique(c(sigmaUnique, sigmaVec)), method="quick")

              if (iter > 1) {
  # if in second iteration mu and sigma are still NA, no value in the grid will fit -> stop iteration
  # else get optimized grid
                  if (is.na(currentBestMu) & is.na(currentBestSigma) & iter == 2) {
                      continueIteration = FALSE
                  } else {
                      muVec      <- optimizeGrid(currentBestMu, muUnique, iter, sigmLimit = FALSE)
                      sigmaVec <- optimizeGrid(currentBestSigma, sigmaUnique, iter, sigmLimit = TRUE)
                  }
              }

              if ((iter == 1) | (iter > 1 & !is.na(currentBestMu) & !is.na(currentBestSigma))) {

                  testedParam <- calculateCostHist(lambda = lambda, muVec = muVec, sigmaVec = sigmaVec, HistData = HistData, alpha = alpha, alphaMcb = alphaMcb, pNormLookup = pNormLookup)

              }

  # update currentBest Parameters is better parameters (smaller costs) are found
              if (testedParam[5] < currentBestCost) {
                  currentBestMu    <- testedParam[2]
                  currentBestSigma <- testedParam[3]
                  currentBestCost  <- testedParam[5]
              }

  # update overall bestParam if tested parameters are superior
              if (testedParam[5] < bestParam[5]) {
                  bestParam    <- testedParam

              }
              iter = iter + 1
          }
          c = c + 1
      }

      return(bestParam)
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  calculateCostHist <- function(lambda, muVec, sigmaVec, HistData, alpha = 0.01, alphaMcb = 0.1, pNormLookup) {

      countsData    <- HistData$counts

  # Box Cox transformation of histogram breaks and histogram range
      breakL <- suppressWarnings(BoxCox(HistData$breakL, lambda))
      breakR <- suppressWarnings(BoxCox(HistData$breakR, lambda))

      breakMinMax <- c(min(breakL), max(breakR))

  # parameters for pnormApprox function
      oneOverH <- 1/pNormLookup$h
      pNormVal <- pNormLookup$val

  # limits for surrounding peak areas needed fo calculation of P
      pLimitMin <- c(0.5, 0.55, 0.60, 0.65, 0.70)
      pLimitMax <- c(0.9, 0.94, 0.95, 0.96, 0.97)

  # precalculated factor for confidence region
      qNormFactor <- -qnorm(p = alpha/2)

  # maximum number of bins expected to be below confidence region
      maxCountsBelow <- as.integer(ceiling(length(countsData)*alphaMcb/2))
  # initial bestParam
      bestParam <- c(NA, NA, NA, NA, 1e9, NA)

      for (mu in muVec) {
          for (sigma in sigmaVec)    {

  # theoretical prediction of bin counts
              pCorr <- pnorm(q=breakMinMax, mean=mu, sd=sigma)
              pCorr <- 1/(pCorr[2]-pCorr[1])

              countsPred <- pCorr*HistData$NData*(pnormApprox(q = breakR, mean = mu, oneOverSd = 1/sigma, oneOverH = oneOverH, pNormVal = pNormVal) - pnormApprox(q = breakL, mean = mu, oneOverSd = 1/sigma, oneOverH = oneOverH, pNormVal = pNormVal))
              countsPred[countsPred < 0] <- 0

  # precalculate sqrt of predicted counts
              sqrtCountsPred <- sqrt(countsPred)

  # selection of bins that can be approximated by normal distribution
              selectionCounts <- (countsPred >= 20)

  # selection of peak bins (without first and last bin)
              maxCountsPred   <- max(countsPred[2:length(countsPred)-1])
              maxIndexPred     <- which(abs(countsPred-maxCountsPred)<1e-20)[1]

              selectionPeak95 <- (countsPred >= 0.95*maxCountsPred)
              selectionPeak95[c(1, length(selectionPeak95))] <- FALSE

              selectionPeak80   <- (countsPred >= 0.80*maxCountsPred)
              selectionPeak80[c(1, length(selectionPeak80))] <- FALSE

              selectionPeak20L <- (countsPred <= 0.20*maxCountsPred & (1:length(countsPred))<maxIndexPred)
              selectionPeak20L[1] <- TRUE

              selectionPeak20R <- (countsPred <= 0.20*maxCountsPred & (1:length(countsPred))>maxIndexPred)
              selectionPeak20R[length(selectionPeak20R)] <- TRUE

  # lower 0.5% should not be negative, thus check if lowerBound is NA
              RRLowerBound <- invBoxCox(mu-qNormFactor*sigma, lambda = lambda)

              if (!is.na(sum(selectionCounts)) & sum(selectionCounts) > (HistData$NBins*HistData$overlapFactor)/16 &  !is.na(RRLowerBound)) {

  # perform pre-calculations for ratio count
                  ratioCounts20PreL <- sum(countsData[selectionPeak20L])/sum(countsPred[selectionPeak20L])
                  ratioCounts20PreR <- sum(countsData[selectionPeak20R])/sum(countsPred[selectionPeak20R])

  # get sum of data points in specified region next to the peak for the actual data and the prediction
                  peakArea  <- getSumForPArea(pLimitMin = pLimitMin, pLimitMax = pLimitMax, countsPred = countsPred, HistData = HistData, lambda = lambda, mu = mu, sigma = sigma, pCorr = pCorr)
                  sumDataPeak <- peakArea$sumDataPeak
                  sumPredPeak <- peakArea$sumPredPeak

                  ratio <- sumDataPeak/sumPredPeak
  #increase ratio until the majority exceed the observed data points
                  while (all(ratio <= 1.0) & sum(ratio*sumPredPeak-qNormFactor*sqrt(ratio*sumPredPeak) <= sumDataPeak)>5) {
                      rInd <- ratio*sumPredPeak-qNormFactor*sqrt(ratio*sumPredPeak) <= sumDataPeak
                      ratio[rInd] <- ratio[rInd] + 0.001
                  }
                  PMax <- min(max(0.401, min(ratio)), 1.000)

                  ratio <- sumDataPeak/sumPredPeak
  #decrease ratio until all but one fall below the observed data points
                  while(all(ratio >= 0.4) & sum(ratio*sumPredPeak+qNormFactor*sqrt(ratio*sumPredPeak) >= sumDataPeak)>1) {
                      rInd <- ratio*sumPredPeak+qNormFactor*sqrt(ratio*sumPredPeak) >= sumDataPeak
                      ratio[rInd] <- ratio[rInd] - 0.001
                  }
                  PMin <- min(max(0.400, min(ratio)), 0.999)

                  PVec  <- c(seq(from=PMin, to=PMax, by=min(max(0.002, (PMax-PMin)/8), 0.01)), PMax)

                  continuePLoop <- TRUE
                  PIndex           <- 1

  # test different values for P and calculate costs...
                  while (continuePLoop) {

                      P <- PVec[PIndex]

                      ratioCounts20L <- max(0.01, min(1, ratioCounts20PreL/P, na.rm=TRUE))^2
                      ratioCounts20R <- max(0.01, min(1, ratioCounts20PreR/P, na.rm=TRUE))^2

  # factors applied for relaxation of cost function
                      rfVec <- c(5, 3, 1)/1000
                      relaxFactors <- P*maxCountsPred / (rfVec*P*maxCountsPred + sqrt(P*maxCountsPred))^2
                      relaxFactors <- c(relaxFactors[relaxFactors>0.001 & relaxFactors<1], 1)

                      for(rf in relaxFactors)
                      {
  # calculate actual predicted counts and corresponding confidence region...
                          countsPredP <- countsPred*P*rf
                          confWidth     <- qNormFactor * sqrtCountsPred * sqrt(P*rf)
                          countsPredPLower <- countsPredP - confWidth
                          countsPredPUpper <- countsPredP + P*confWidth

  # select bins within confidence region
                          selectionPredBand <- (rf*countsData <= countsPredPUpper & rf*countsData >= countsPredPLower)
                          selectionEval     <- (selectionCounts & selectionPredBand)

  # count number of bins that are selected and not selected for evaluation
                          sumSelectionEval  <- sum(selectionEval)

  # count number of bins that are at the peak and selected for evaluation
                          sumSelectionPeakEval  <- sum(selectionEval & selectionPeak95)

  # decide if loop shall be continued...
                          continuePLoop <- (sum(rf*countsData < countsPredPLower) <= maxCountsBelow & PIndex < length(PVec))
  # if expected number of bins are selected...
                          if (sumSelectionEval > (HistData$NBins*HistData$overlapFactor)/32 & (sumSelectionPeakEval>=(0.2*(sum(selectionPeak95))) | sum(selectionPeak95)==1))
                          {
                              ratioCounts80  <- max(0.01, min(1, sum(countsData[selectionEval & selectionPeak80])/sum(countsData[!selectionEval & selectionPeak80]), na.rm=TRUE))^2

  # extract data of selection
                              cData <- rf*countsData[selectionEval]
                              cPred <- countsPredP[selectionEval]

                              costPreSum <- sum(log(sqrt(P*rf/(2*pi))*sqrtCountsPred[selectionEval]))

  # calculate costs (negative log likelihood with regularization term)

  # identical to the approximation with dnorm(), but faster
                              cost <- -(costPreSum + sumSelectionEval*(log(ratioCounts80)+log(ratioCounts20L)+log(ratioCounts20R)) + sum(-0.5*(cData-cPred)^2/cPred))/sqrt(sumSelectionEval)

  ### derivation of faster cost function calculation

  # original costs with poisson distribution based on statistical assumption of histogram data, with additional regularization
  # pois         <- dpois(x = round(cData*rf), lambda = cPred*rf)*cPred*rf*ratioCounts80*ratioCounts20L*ratioCounts20R

  # very good approximation of calculation of likelihood (with regularization term)
  # pois <- dnorm(x=cData*rf, mean=cPred*rf, sd=sqrt(cPred*rf))*cPred*rf*ratioCounts80*ratioCounts20L*ratioCounts20R

  # rearrange equation
  # cost <- -sum(log(dnorm(x=cData, mean=cPred, sd=sqrt(cPred))*cPred*ratioCounts80*ratioCounts20L*ratioCounts20R)/sqrt(sumSelectionEval))
  # cost <- -(  sumSelectionEval*(log(ratioCounts80)+ratioCounts20L+ratioCounts20R) + sum(log(sqrt((P*rf)/(2*pi))*sqrtCountsPred[selectionEval])) +sum((-0.5*(cData - cPred)^2)/cPred ))/sqrt(sumSelectionEval)
  # cost <- -(  sumSelectionEval*(log(ratioCounts80)+ratioCounts20L+ratioCounts20R) + costPreSum) +sum((-0.5*(cData - cPred)^2)/cPred ))/sqrt(sumSelectionEval)

  # update best cost
                              if (!is.na(cost) & cost < bestParam[5]) {
                                  bestParam <- c(lambda, mu, sigma, P, cost, rf)
                              }
                          }
                      }
                      PIndex <- PIndex + 1
                  }
              }
          }
      }

      bestParam
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  findPeaksAndValleys <- function(Dens) {

      diffDens <- diff(Dens$y)
      peakCrit <- valleyCrit <- 0
      peaksV   <- valleysV  <- NULL
      prev <- act <- diffDens[1]
      prevSgn <- ifelse(act < 0, "-", "+")

      lastP = FALSE
      lastV = FALSE

      for (i in 2:length(diffDens)) {
  # get current value and sign
          act <- diffDens[i]
          sgn <- ifelse(act < 0, "-", "+")

  # determine local minima and maxima depending on the change of signs (- to + is minima; + to - is maxima)
          if (prevSgn == '-' & sgn == '+') {
              valleyCrit <- 1
              prevSgn   <- sgn

          } else if(sgn == '+' & valleyCrit > 0) {
              valleyCrit <- valleyCrit + 1

          } else if (prevSgn == '+' & sgn == '-') {
              peakCrit <-  1
              prevSgn  <- sgn

          } else if (sgn == '-' & peakCrit > 0) {
              peakCrit <- peakCrit +1

          } else {
              prevSgn        <- sgn
              prev        <- act
              peakCrit     <- 0
              valleyCrit     <- 0
          }
  #after a minima or maxima, the same sign has to appear at least 5 times, otherwise it is considered noise
          if (peakCrit == 5) {
              if (!lastP) {
                  maxPeak <- which.max(Dens$y[(i-5):i])
                  peaksV  <- c(peaksV, (i-5+maxPeak-1))
                  lastP   <- TRUE
                  lastV   <- FALSE
              }
              peakCrit     <- 0
              valleyCrit     <- 0
              prevSgn     <- sgn
              prev         <- act

          }else if (peakCrit >= 1 & i == length(diffDens) & is.null(peaksV)){
              if (!lastP) {
                  maxPeak <- which.max(Dens$y[(i-5):i])
                  peaksV  <- c(peaksV, (i-5+maxPeak-1))
                  lastP   <- TRUE
                  lastV   <- FALSE
              }
              peakCrit     <- 0
              valleyCrit     <- 0
              prevSgn     <- sgn
              prev         <- act
          }

          if (valleyCrit == 5) {
              if (!lastV) {
                  minValley <- which.min(Dens$y[(i-5):i])
                  valleysV  <- c(valleysV, (i-5+minValley-1))
                  lastP     <- FALSE
                  lastV     <- TRUE
              }

              valleyCrit <- 0
              peakCrit   <- 0
              prevSgn    <- sgn
              prev        <- act
          }
      }

  # if distribution starts with a valley, delete it
      if (length(valleysV)!= 0 && !is.null(peaksV) && valleysV[1] < peaksV[1])
          valleysV <- valleysV[-1]

  # if there are no valleys found due to smoothed density function, add a valley at the end
      if(length(valleysV)==0)
          valleysV <- length(Dens$x)

      return (list("peaks" = peaksV, "valleys"= valleysV))
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  optimizeGrid <- function(currentBestParam, paramUnique, iter, sigmLimit = TRUE) {

      optInd <- which.min(abs(paramUnique - currentBestParam))

      if (!is.na(currentBestParam) && optInd == 1) {
  # param is at left border, adapted grid exceeds original search space
          diff <- paramUnique[2]-paramUnique[1]
  # in first iteration, adapted grid is larger than in following iterations
          if (iter == 2) {
  # account for lower limit of sigma, cannot be < 0
              if (sigmLimit)
                  tmp <- c(max(paramUnique[optInd]-2*diff, 0), max(paramUnique[optInd]-diff, 0), paramUnique)
              else
                  tmp <- c(paramUnique[optInd]-2*diff, paramUnique[optInd] - diff, paramUnique)
              optInd = optInd+2
          } else {
  # account for lower limit of sigma, cannot be < 0
              if (sigmLimit)
                  tmp <- c(max(paramUnique[optInd]-diff, 0), paramUnique)
              else
                  tmp <- c(paramUnique[optInd]-diff, paramUnique)
              optInd = optInd + 1
          }

          OptL     <- max(optInd - ifelse(iter == 2, 2, 1), 1)
          OptR      <- min(optInd + ifelse(iter == 2, 2, 1), length(tmp))
          paramVec <- c(currentBestParam, 0.5*(tmp[OptL] + tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))

      } else if (!is.na(currentBestParam) && optInd == length(paramUnique)) {
  #param is at right border /upper margin of vector, adapted grid exceeds original search space
          diff <- paramUnique[optInd] - paramUnique[optInd-1]
  # in first iteration, adapted grid is larger than in following iterations
          if (iter == 2) {
              tmp <- c(paramUnique, paramUnique[optInd]+diff, paramUnique[optInd]+2*diff)
          } else {
              tmp <- c(paramUnique, paramUnique[optInd]+diff)
          }

          OptL      <- max(optInd - ifelse(iter == 2, 2, 1), 1)
          OptR      <- min(optInd + ifelse(iter == 2, 2, 1), length(tmp))
          paramVec <- c(currentBestParam, 0.5*(tmp[OptL]+ tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))

      } else if (!is.na(currentBestParam) && optInd == 2 & iter == 2) {
  #param is at position 2 but iter is 2, adapted grid exceeds original search space
          diff <- paramUnique[2]- paramUnique[1]
  # account for lower limit of sigma, cannot be < 0
          if (sigmLimit)
              tmp <- c(max(paramUnique[optInd]-2*diff, 0), paramUnique)
          else
              tmp <- c(paramUnique[optInd]-2*diff, paramUnique)

          optInd = optInd+1

          OptL      <- max(optInd-2, 1)
          OptR      <- min(optInd+2, length(tmp))
          paramVec <- c(currentBestParam, 0.5*(tmp[OptL] + tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))

      } else if (!is.na(currentBestParam) && optInd == length(paramUnique)-1 & iter == 2) {
  #param is at position length -1 (penultimate position) and iter is 2, adapted grid exceeds original search space
          diff <- paramUnique[optInd]-paramUnique[optInd-1]
          tmp  <- c(paramUnique, paramUnique[optInd]+2*diff)
          OptL <- max(optInd-2, 1)
          OptR <- min(optInd+2, length(tmp))
          paramVec <- c(currentBestParam, 0.5*(tmp[OptL]+tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))

      } else {
  #param is in the middle of the vector/grid
          OptL <- max(optInd - ifelse(iter == 2, 2, 1), 1)
          OptR <- min(optInd + ifelse(iter == 2, 2, 1), length(paramUnique))
          paramVec <- c(currentBestParam, 0.5*(paramUnique[OptL] + paramUnique[OptL+1]), 0.5*(paramUnique[OptR] + paramUnique[OptR-1]))
      }

      return(paramVec)
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  getSumForPArea <- function(pLimitMin, pLimitMax, countsPred, HistData, lambda, mu, sigma, pCorr) {

      N <- length(countsPred)

      countsFullHist <- HistData$fullHist$counts
      breaks           <- HistData$fullHist$breaks
      nLimits         <- length(pLimitMin)

      sumDataPeak <- rep(1, times = 2*nLimits+1)

      borderLvec <- rep(0, times = 2*nLimits+1)
      borderRvec <- rep(0, times = 2*nLimits+1)

  # determine peak
      indexPeak <- which.max(countsPred)
      peak      <- countsPred[indexPeak]

      pLimitMin <- pLimitMin*peak
      pLimitMax <- pLimitMax*peak

  # selection of peak bins
      selectionPeak <- (countsPred >= 0.95*peak)

  #get leftmost and rightmost value of peak bins
      borderL <- min(HistData$breakL[selectionPeak])
      borderR <- max(HistData$breakR[selectionPeak])
      borderLvec[1] <- borderL
      borderRvec[1] <- borderR

  #get observed counts at the peak
      sumDataPeak[1] <- sum(countsFullHist[breaks >= borderL & breaks < borderR])

  # split data into left and right part
      countsPredLeft     <- countsPred[1:indexPeak]
      breakLLeft         <- HistData$breakL[1:indexPeak]
      breakRLeft         <- HistData$breakR[1:indexPeak]

      countsPredRight <- countsPred[indexPeak:N]
      breakLRight     <- HistData$breakL[indexPeak:N]
      breakRRight     <- HistData$breakR[indexPeak:N]

      for (i in 1:length(pLimitMin)) {

  #get index of area left from the peak
          selectionLeft <- (countsPredLeft >= pLimitMin[i] & countsPredLeft < pLimitMax[i])

          if (sum(selectionLeft) > 0) {
  #get leftmost and rightmost value of selected region
              borderL <- min(breakLLeft[selectionLeft])
              borderR <- max(breakRLeft[selectionLeft])
              borderLvec[2*i] <- borderL
              borderRvec[2*i] <- borderR

  #get observed counts in the region next to the peak
              sumDataPeak[2*i] <- sum(countsFullHist[breaks >= borderL & breaks < borderR])
          }

  #get index of area right from the peak
          selectionRight <- (countsPredRight >= pLimitMin[i] & countsPredRight < pLimitMax[i])

          if (sum(selectionRight) > 0) {
  #get leftmost and rightmost value of selected region
              borderL <- min(breakLRight[selectionRight])
              borderR <- max(breakRRight[selectionRight])
              borderLvec[2*i+1] <- borderL
              borderRvec[2*i+1] <- borderR

  #get observed counts in the region next to the peak
              sumDataPeak[2*i+1] <- sum(countsFullHist[breaks >= borderL & breaks < borderR])
          }
      }

  #get predicted counts at peak and in the regions next to the peak
      sumPredPeak <- pCorr*HistData$NData*(pnorm(q=BoxCox(borderRvec, lambda), mean = mu, sd = sigma) - pnorm(q = BoxCox(borderLvec, lambda), mean = mu, sd = sigma))
      sumPredPeak[abs(sumPredPeak) < 1e-20] <- 1

      return(list(sumDataPeak = sumDataPeak, sumPredPeak = sumPredPeak))
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  findMainPeak <- function(x, ab, mStart, withHeight = FALSE, prevPeak =NULL) {

      qDensFullTr <- ashDensity(x = x, ab = ab, nbin = 512, m = mStart)

  # find peak with largest area
      listPV   <- findPeaksAndValleys(qDensFullTr)
      peaks      <- listPV$peaks
      valleys  <- listPV$valleys
      area      <- rep(0, times = length(peaks))
      height95 <- rep(1, times = length(peaks))

  #get area under peak (between two valleys)
      for (counter in 1:length(peaks)) {
          if (counter ==1) {
              vLow <- 1
          } else {
              vLow <- valleys[counter-1] #-1 as the function returns the value after the peak/valley
          }
          p <- peaks[counter]

          if (counter > length(valleys)) {
              vHigh <- length(qDensFullTr$y)
          } else {
              vHigh <- valleys[counter]
          }

  # normalize the two measures used for peak detection
          area[counter]       <- sum(qDensFullTr$y[vLow:vHigh])

          if(withHeight)
              height95[counter] <- quantile(qDensFullTr$y[vLow:vHigh], probs=0.95)
      }

      area      <- area/max(area)
      height95 <- height95/max(height95)
      locPrevP <- 0

      if(withHeight & length(height95)>1)
      {
          maxIndex <- which.max(height95)

          ratioHeight <- height95[maxIndex]/max(height95[-maxIndex])
          wHeight95   <- min(1, max(0, (ratioHeight-2)*5 ))

          maxIndexA <- which.max(area)
          if(max(area[-maxIndexA]) >= 0.95 & wHeight95 == 0)
              wHeight95 <- min(1, max(0, (max(area[-maxIndexA])-0.95)*20)) #0.970 0.975 0.980 0.985 0.990 0.995 1.000 -->  0.4 0.5 0.6 0.7 0.8 0.9 1.0

          height95 <- height95*wHeight95
      }

      if(!is.null(prevPeak) & length(peaks) > 1){
          locPrevP <- abs(qDensFullTr$x[peaks]-prevPeak)
          locPrevP <- locPrevP/max(locPrevP)
      }

  #get index of peak in area with biggest surface
      peakInd <- min(length(qDensFullTr$x)-1, max(2, peaks[which.max(area+height95-locPrevP)]))
      modEst     <- sum(qDensFullTr$x[peakInd+c(-1:1)]*qDensFullTr$y[peakInd+c(-1:1)])/sum(qDensFullTr$y[peakInd+c(-1:1)])

  # select region around peak for stronger smoothing
      abPeak      <- NULL
      abPeak[1] <- max(ab[1], qDensFullTr$x[which(qDensFullTr$x < modEst & qDensFullTr$y <= 0.35*qDensFullTr$y[peakInd])])
      abPeak[2] <- min(ab[2], qDensFullTr$x[which(qDensFullTr$x > modEst & qDensFullTr$y <= 0.35*qDensFullTr$y[peakInd])])

  # apply more significant/stronger smoothing in selected region
      densMM <- ceiling(mStart/1.5)
      qDensTr <- ashDensity(x = x, ab = abPeak, nbin = 256, m = densMM)

  # find peak with largest area
      listPV   <- findPeaksAndValleys(qDensTr)
      peaks      <- listPV$peaks
      valleys  <- listPV$valleys
      area      <- rep(0, times=length(peaks))
      height95 <- rep(1, times=length(peaks))

  #get area under peak (between two valleys)
      for (counter in 1:length(peaks)) {
          if (counter ==1) {
              vLow <- 1
          } else {
              vLow <- valleys[counter-1] #-1 as the function returns the value after the peak/valley
          }
          p <- peaks[counter]

          if (counter > length(valleys)) {
              vHigh <- length(qDensTr$y)
          } else {
              vHigh <- valleys[counter]
          }

          if(length(vLow) > 0 & length(vHigh) > 0 ){

              area[counter] <- sum(qDensTr$y[vLow:vHigh])

              if(withHeight)
                  height95[counter] <- quantile(qDensTr$y[vLow:vHigh], probs=0.95)
          }
      }

  # normalize the two measures used for peak detection
      area      <- area/max(area)
      height95 <- height95/max(height95)
      location <- 0

      if(withHeight & length(height95)>1)
      {
          maxIndex <- which.max(height95)


          ratioHeight <- height95[maxIndex]/max(height95[-maxIndex])

          wHeight95 <- min(1, max(0, (ratioHeight-2)*5 ))

          maxIndexA <- which.max(area)
          if(max(area[-maxIndexA]) >= 0.95 & wHeight95 == 0)
              wHeight95 <- min(1, max(0, (max(area[-maxIndexA])-0.95)*20)) #0.970 0.975 0.980 0.985 0.990 0.995 1.000 -->  0.4 0.5 0.6 0.7 0.8 0.9 1.0

          height95 <- height95*wHeight95
      }

      if(length(peaks) > 1){
          location <- abs(qDensTr$x[peaks]-modEst)
          location <- location/max(location)
      }

  #get index of peak in area with biggest surface
      peakInd <- min(length(qDensTr$x)-1, max(2, peaks[which.max(area+height95-location)]))
      modEst     <- sum(qDensTr$x[peakInd+c(-1:1)]*qDensTr$y[peakInd+c(-1:1)])/sum(qDensTr$y[peakInd+c(-1:1)])

      return(list(peakInd=peakInd, modEst=modEst, Dens=qDensTr))
  }

  # Tatjana Ammer \email{tatjana.ammer@@roche.com}
  estimateAB <- function(x) {

      QQ <- as.numeric(quantile(x, probs=c(0, 0.01, 0.95, 0.99, 0.999), na.rm=TRUE))

  # default ab used as a starting point for optimization
      abOpt    <- QQ[c(2, 4)]
      rangeAB  <- diff(abOpt)
      abOpt[1] <- max(abOpt[1]-0.1*rangeAB, QQ[1], 1e-20)

      abOpt[2] <- abHistMax <- min(abOpt[2]+0.1*rangeAB, QQ[5])

  # default ab for histogram
      abHist <- abSearchReg <- abOpt

  # ab for peak detection
      abPeakDet    <- QQ[c(2, 4)]
      abPeakDet[1] <- max(abPeakDet[1], 1e-20)
      abPeakDet[1] <- abPeakDet[1] - 0.01*diff(abPeakDet)

      stepSize <- 0.5*diff(abOpt)
      peakM      <- 24

      iter <- 1
      continueABLoop <- TRUE

      while(continueABLoop & iter <= 30 & abOpt[2]>abOpt[1] & abOpt[2]>abPeakDet[1] )    {

          abPeakTemp         <- c(max(abPeakDet[1], 1e-20), min(abOpt[2], abPeakDet[2]))
          abPeakTemp[1]     <- abPeakTemp[1] - 0.01*diff(abPeakTemp)

  # find main peak of distribution
          mainPeak  <- findMainPeak(x=x, ab=c(abPeakTemp[1], min(abOpt[2], abPeakDet[2])), mStart=peakM, withHeight=TRUE)
          Dens       <- ashDensity(x=x, ab=abOpt, nbin=512, m=24)
          indexPeak <- which.min(abs(Dens$x-mainPeak$modEst))

  # left side of the peak
          yLeft <- rev(Dens$y[1:indexPeak])

  # pre-process density that it becomes monotonic function
          if(length(yLeft) > 1)
          {
              for(i in 2:length(yLeft))
              {
                  if(yLeft[i] > yLeft[i-1])
                      yLeft[i] <- yLeft[i-1]
              }
          }

  # right side of the peak
          yRight <- Dens$y[indexPeak:length(Dens$x)]

  # pre-process density that it becomes monotonic function
          if(length(yRight) > 1)
          {
              for(i in 2:length(yRight))
              {
                  if(yRight[i] > yRight[i-1])
                      yRight[i] <- yRight[i-1]
              }
          }

  # transform y into weights and calculate weighted proportion of low density bins
          y <- c(yLeft, yRight)
          y <- y/max(y)
          y <- 2-y*5
          y <- pmin(pmax(y, 0), 1)

          propLowDens <- sum(y)/length(y)

          if((iter==1 | propLowDens>0.85) & propLowDens<0.86)    {
              continueABLoop <- FALSE
              abHist            <- abOpt
              abHist[2]        <- min(abHist[2] + 0.5*diff(abHist), abHistMax)

          } else if(propLowDens>=0.86) {
              abOpt[2]    <- abOpt[2] - stepSize

          } else if(propLowDens<=0.85) {
              abOpt[2]    <- abOpt[2] + stepSize
          }

          stepSize <- stepSize/2
          iter      <- iter + 1
      }

      abOpt <- abSearchReg <- abHist

  # find main peak of distribution
      mainPeak  <- findMainPeak(x=x, ab=abOpt, mStart=peakM, withHeight=TRUE, prevPeak = mainPeak$modEst)
      Dens       <- ashDensity(x=x, ab=abOpt, nbin=512, m=24)
      indexPeak <- which.min(abs(Dens$x-mainPeak$modEst))

      yMin <- yMax <- Dens$y[indexPeak]

      indexL <- max(which((1:length(Dens$x))<indexPeak & Dens$y<0.5*yMax), 2)
      indexR <- min(which((1:length(Dens$x))>indexPeak & Dens$y<0.5*yMax), length(Dens$y)-1)

      area50 <- sum(Dens$y[indexL:indexR])

  # left side of the peak
      xLeft <- rev(Dens$x[1:(indexL-1)])
      yLeft <- rev(Dens$y[1:(indexL-1)])

      if(length(yLeft) > 1) {

          continueABLoop <- TRUE

          for(i in 1:length(yLeft)) {

              if(yLeft[i]<0.5*yMax & yLeft[i]<yMin)
                  yMin <- yLeft[i]

              if(yLeft[i] > yMin+0.5*yMax & continueABLoop) {

                  abOpt[1] <- xLeft[i]
                  abSearchReg[1] <- xLeft[i]
                  continueABLoop <- FALSE
              }

              if(sum(yLeft[1:i])<1.3*area50 & continueABLoop)
                  abOpt[1] <- xLeft[i]

              if(sum(yLeft[1:i])<1.1*area50 & continueABLoop)
                  abSearchReg[1] <- xLeft[i]
          }
      }

      yMin <- yMax

  # right side of the peak
      xRight <- Dens$x[(indexR+1):length(Dens$x)]
      yRight <- Dens$y[(indexR+1):length(Dens$x)]

      if(length(yRight) > 1)    {

          continueABLoop <- TRUE

          for(i in 1:length(yRight)) {

              if(yRight[i]<0.5*yMax & yRight[i]<yMin)
                  yMin <- yRight[i]

              if(yRight[i] > yMin+0.5*yMax & continueABLoop) {
                  abOpt[2] <- xRight[i]
                  abSearchReg[2] <- xRight[i]
                  continueABLoop <- FALSE
              }

              if(sum(yRight[1:i])<1.3*area50 & continueABLoop)
                  abOpt[2] <- xRight[i]


              if(sum(yRight[1:i])<1.1*area50 & continueABLoop)
                  abSearchReg[2] <- xRight[i]
          }
      }

      abHist <- abOpt

      return(list(abHist = abHist, abSearchReg = abSearchReg))
  }

  # The unmodified Box-Cox search path from findRI with NBootstrap = 0.
  # Build the point fit only, with the same two future seed calls as upstream.
  findBoxCoxPoint <- function(x, seed) {
    raw <- stage2_routine(x)
    n <- length(raw)
    if (n > 50000000)
      message(" Data has extremely large sample size (N > 50'000'000). This may lead to performance and memory issues.")
    if (n < 1000)
      message(" Data has small sample size (N < 1000). Evaluate results carefully.")
    stage2_job(seed, function() {
      data <- rep(raw, times = ceiling(250 / n)) - 0
      rounding <- findRoundingBase(data)
      lambda <- (1 / 12 * (0:12))^1.8170595
      starts <- defineSearchRegions(data, lambda, rounding)
      hist_data <- generateHistData(data, starts$abEst$abHist, rounding)
      initial <- c(NA, NA, NA, NA, 1e9, NA)
      best <- testParam(lambda, initial, data, hist_data, starts, 7, 0.01, 0.1)

      middle <- min(length(lambda) - 2, max(3, which.min(abs(best[1L] - lambda))))
      lambda <- unlist(lapply((middle - 2):(middle + 1), function(i) {
        seq(lambda[i], lambda[i + 1], length.out = 5)[2:4]
      }), use.names = FALSE)
      starts <- defineSearchRegions(data, lambda, rounding, starts$abEst)
      best <- testParam(lambda, best, data, hist_data, starts, 7, 0.01, 0.1)
      if (best[5L] >= initial[5L]) best <- initial

      structure(list(Data = raw, Method = "refineR", PkgVersion = package_version("2.0.0"),
                     Model = "BoxCox", NBootstrap = 0, Lambda = best[1L], Mu = best[2L],
                     Sigma = best[3L], P = best[4L], Shift = 0, Cost = best[5L],
                     abOr = hist_data$abOr, roundingBase = rounding, rf = best[6L]),
                class = "RWDRI")
    })
  }

  environment()
})
