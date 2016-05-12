# $Id: gstat.cv.q,v 1.9 2009-10-30 16:11:21 edzer Exp $

"gstat.cv" <-
  function (object, nfold = nrow(object$data[[1]]$data), remove.all = FALSE, verbose = interactive(), 
            all.residuals = FALSE, ...) {
    
    if (!inherits(object, "gstat")) 
      stop("first argument should be of class gstat")
    
    # ASR: Select first list of data/model, get model frame, and get model formula
    var1 <- object$data[[1]]
    data <- var1$data
    formula <- var1$formula
    
    if (all.residuals) {
      # ASR: Get number of variables and create data.frame to store the residuals of all variables
      nc <- length(object$data)
      ret <- data.frame(matrix(NA, nrow(data), nc))
    } else {
      
      # ASR: Create SpatialPointsDataFrame to store the residuals of the first variable
      cc <- coordinates(data)
      rownames(cc) <- NULL
      df <- data.frame(matrix(as.numeric(NA), nrow(data), 2))
      ret <- SpatialPointsDataFrame(cc, df)
    }
    
    # ASR: Define the number of folds
    if (missing(nfold)) 
      nfold <- 1:nrow(data)                             # ASR: loo-cv
    if (length(nfold) == nrow(data))
      fold <- nfold                                     # ASR: loo-cv
    else if (nfold < nrow(data)) 
      fold <- sample(nfold, nrow(data), replace = TRUE) # k-fold cv
    else fold <- 1:nrow(data)                           # ASR: loo-cv
    
    # ASR: Check if the residuals of all variables should be returned or 
    #      if the test points should be removed from all variables and there are more than one variable.
    #      If TRUE, copy all data to a list with as many itens as there are variables.
    if (all.residuals || (remove.all && length(object$data) > 1)) {
      all.data <- list()
      for (v in 1:length(object$data))
        all.data[[v]] <- object$data[[v]]$data
    }
    
    # ASR: Print status
    if (verbose)
      pb <- txtProgressBar(1, length(unique(fold)), style=3)
    
    # ASR: Start main loop
    # 
    # ASR: Check how many folds we have
    for (i in sort(unique(fold))) {
      if (verbose)
        setTxtProgressBar(pb, i)
      
      # ASR: Identify test observations
      sel <- which(fold == i)
      
      # ASR: Select calibration observations
      object$data[[1]]$data <- data[-sel, ]
      
      # ASR: Check if the test points should be removed from all variables and there are more than one 
      #      variable.
      #      If TRUE, reformulate the model frame in the data object.
      if (remove.all && length(object$data) > 1) {
        
        # ASR: Reformulate the model frame in the data object. Skip the first because it was reformulated 
        #      above (default).
        for (v in 2:length(object$data)) {
          varv <- object$data[[v]]
          varv$data <- all.data[[v]] # ASR: I THINK THIS LINE IS UNECESSARY BECAUSE IT SEEMS TO REPEAT
          #                                 THE LINE ABOVE!

          #atv = gstat.formula(varv$formula, varv$data)$locations
          #at1 = gstat.formula(formula, data[sel, ])$locations
          
          atv <- coordinates(varv$data)  # asr: get coordinates of all observations
          at1 <- coordinates(data[sel,]) # asr: get coordinates of selected test observations
          cc <- rbind(atv, at1)          # asr: stack coordinates
          rownames(cc) <- NULL # as there will be duplicates
          all <- SpatialPoints(cc)       # asr: create spatial points
          zd <- zerodist(all)            # asr: identify duplicates
          skip <- zd[, 1]                # ASR: I THINK 'skip' IS EXACTLY EQUAL TO 'sel'.
          object$data[[v]]$data <- varv$data[-skip, ] # asr: keep non-duplicated points for calibration
          
          # ASR: Perhaps the previous lines of code could be replaced with:
          # varv <- object$data[[v]]
          # object$data[[v]]$data <- varv$data[-sel, ]
          # However, the code seems to be optimized to deal with non-collocated cokriging.
        }
      }
      
      # asr: this is where the re-fit of the variogram comes in!!!
      
      # asr: make predictions at test locations using the object reformulated above
      x <- predict(object, newdata = data[sel, ], ...)
      
      # asr: get residuals of the test set
      if (all.residuals) {
        for (i in 1:length(object$data)) {
          var.i <- object$data[[i]] # asr: get first list of data/models
          data.i <- all.data[[i]]   # asr: get all observed values
          formula.i <- var.i$formula # asr: get formula
          observed <- gstat.formula(formula.i, data.i)$y[sel] # asr: get observed values at test locations
          pred.name <- paste(names(object$data)[i], "pred", sep = ".")
          residual <- as.numeric(observed - x[[pred.name]]) # asr: compute residuals
          ret[sel, i] <- residual # asr: store residuals in the specific object (data.frame)
        }
      } else {
        # asr: store prediction and prediction error variance of first variable in 
        #      SpatialPointsDataFrameObject
        ret[[1]][sel] <- x[[1]]
        ret[[2]][sel] <- x[[2]]
      }
    }
    
    if (verbose)
      cat("\n")
    
    # asr: prepare output when interested only in the first variable
    if (!all.residuals) {
      names(ret) <- names(x)[1:2]
      ret$observed <- gstat.formula(formula, data)$y
      pred.name <- paste(names(object$data)[1], "pred", sep = ".")
      ret$residual <- ret$observed - ret[[pred.name]]
      var.name <- paste(names(object$data)[1], "var", sep = ".")
      ret$zscore <- ret$residual/sqrt(ret[[var.name]])
      ret$fold <- fold
    } else
      names(ret) <- names(object$data)
    
    if (!is.null(object$locations))
      ret <- as.data.frame(ret)
    
    ret
  }
