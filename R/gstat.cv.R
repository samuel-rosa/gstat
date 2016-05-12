# $Id: gstat.cv.q,v 1.9 2009-10-30 16:11:21 edzer Exp $

"gstat.cv" <-
  function (object, nfold = nrow(object$data[[1]]$data), remove.all = FALSE, verbose = interactive(), 
            all.residuals = FALSE, 
            # New arguments to re-fit the variogram model
            refit = FALSE, # should the variogram model be re-fitted?
            ...) {
    
    if (!inherits(object, "gstat")) 
      stop("first argument should be of class gstat")
    
    # ASR: Select first list of data/model, get model frame, and get model formula
    var1 <- object$data[[1]]
    data <- var1$data # keep the original data
    formula <- var1$formula
    
    # Create object to store cross-validation results
    if (all.residuals) {
      nc <- length(object$data)
      ret <- data.frame(matrix(NA, nrow(data), nc))
    } else {
      cc <- coordinates(data)
      rownames(cc) <- NULL
      df <- data.frame(matrix(as.numeric(NA), nrow(data), 2))
      ret <- SpatialPointsDataFrame(cc, df)
    }
    
    # Define cross-validation folds
    if (missing(nfold)) 
      nfold <- 1:nrow(data)
    if (length(nfold) == nrow(data))
      fold <- nfold
    else if (nfold < nrow(data)) 
      fold <- sample(nfold, nrow(data), replace = TRUE)
    else fold <- 1:nrow(data)
    
    # ASR: Check if the residuals of all variables should be returned or 
    #      if the test points should be removed from all variables and there are more than one variable.
    #      If TRUE, copy all data to a list with as many itens as there are variables.
    if (all.residuals || (remove.all && length(object$data) > 1)) {
      all.data <- list()
      for (v in 1:length(object$data))
        all.data[[v]] <- object$data[[v]]$data
    }
    
    if (verbose)
      pb <- txtProgressBar(1, length(unique(fold)), style=3)
    
    # Cross-validation loop
    for (i in sort(unique(fold))) {
      if (verbose)
        setTxtProgressBar(pb, i)
      sel <- which(fold == i)               # test observations
      object$data[[1]]$data <- data[-sel, ] # calibration observations
      
      # ASR: Check if the test points should be removed from all variables and there are more than one 
      #      variable.
      #      If TRUE, reformulate the model frame in the data object.
      if (remove.all && length(object$data) > 1) {
        
        # Reformulate model frames (skip the first because it was reformulated above).
        for (v in 2:length(object$data)) {
          varv <- object$data[[v]]
          varv$data <- all.data[[v]]
          #atv = gstat.formula(varv$formula, varv$data)$locations
          #at1 = gstat.formula(formula, data[sel, ])$locations
          atv <- coordinates(varv$data)
          at1 <- coordinates(data[sel,])
          cc <- rbind(atv, at1)
          rownames(cc) <- NULL # as there will be duplicates
          all <- SpatialPoints(cc)
          zd <- zerodist(all)
          skip <- zd[, 1]
          object$data[[v]]$data <- varv$data[-skip, ]
        }
      }
      
      # Re-fit variogram models when needed.
      # NOTE: This is working only for multivariate variogram models.
      if (refit) {
        
        # Start getting the fitted variogram models that will be used as 
        # initial values when re-fitting the variogram models.
        data_names <- names(object$data)
        fitted_models <- object$model[data_names]
        
        # Create new gstat object with initial values for variogram models
        for (v in 1:length(object$data)) {
          if (v == 1) {
            g <- gstat(
              g = NULL, id = data_names[v], data = object$data[[v]]$data, formula = object$data[[v]]$formula)
          } else {
            g <- gstat(
              g = g, id = data_names[v], data = object$data[[v]]$data, formula = object$data[[v]]$formula)
          }
        }
        g$model <- object$model
        
        # Compute empirical variograms
        # How do we get arguments to be passed to 'variogram'? I think we have to create new arguments or use
        # a 'control' function, something like 'refit.control'.
        v <- variogram(g, boundaries = attr(vario, "boundaries"))
        
        # Re-fit the variogram models:
        # - 'fit.lmc' is used to fit multivariate variogram models
        # - 'fit.variogram' is used to fit univariate variogram models
        if (length(object$data) > 1) {
          object_refit <- fit.lmc(v = v, g = g, correct.diagonal = 1.01)
        } else {
          # object_refit <- fit.variogram(object = v, ...)
        }
        
      }
      
      # ASR: Here we make predictions at the test locations ('sel') using the
      #      (1) reformulated and re-fited or (2) reformulated model object.
      if (refit) {
        x <- predict(object_refit, newdata = data[sel, ], ...)
      } else {
        x <- predict(object, newdata = data[sel, ], ...)
      }
      
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
