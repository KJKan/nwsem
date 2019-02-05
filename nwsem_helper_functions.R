###################################################################
#                                                                 #
#      Kan, K. J., van der Maas, L. J. & Levine, S. Z. (2019).    #
#           Extending psychometric network analysis:              #
#      Empirical evidence against g in favor of mutualism?        #
#                   Intelligence, 73, 52-62.                      #
#          https://doi.org/10.1016/j.intell.2018.12.004           #
#                                                                 #
# This script contains helper functions                           #
#                                                                 #
# Author: Kees-Jan Kan                                            #
#                                                                 #
###################################################################

# This function reads the factor loading pattern of a measurement model.
readpat<-function(ne){
  matrix( scan( what = 'character' ), , ne, byrow = TRUE )
}


# This function labels for elements of (e.g. OpenMx) matrices.
label <- function(lab, nrow, ncol = NA, add = NA, type = NA, triangular = NA ){
  if( is.na( ncol ) ) {
    ncol <- nrow
    if( is.na( type ) ) type <- 'Symm'
    if ( is.na( triangular ) ) { triangular = FALSE }
  } else if (is.na( triangular ) ) { triangular = FALSE
  type = 'Full' }
  labels <- matrix( paste (lab,
                           rep( 1:nrow, each = ncol ),
                           rep( 1:ncol, nrow ),
                           sep = '_' ),
                    nrow, ncol, byrow = TRUE )

  if (!is.na(add)) labels <- matrix( paste( labels, add, sep = '_'), nrow, ncol, byrow = TRUE )

  if ( type == 'Symm' ) labels[ upper.tri( labels ) ] <- t( labels )[ upper.tri( t( labels ) ) ]
  if ( type == 'Diag' ) labels <- diag( labels )
  if ( triangular) labels <- labels[ lower.tri( labels, TRUE ) ]
  return( labels )
}

# This function pastes a string to another string, e.g. the labels of an OpenMx matrix.
add2labels <- function( labels, add, sep = '_' ) paste( labels, add, sep = sep )

#' This function pastes a string to the labels of an OpenMx matrix within an MxModel.
relabel <- function( Model, Matrix, add, sep = '_' ) {
  labels <- attributes( Model )[[ 'matrices' ]] [[ Matrix ]]$labels
  omxSetParameters(Model,
                   labels = labels,
                   newlabels = paste( labels, add, sep = sep ) )
}

# This function relabels all OpenMx matrices within an MxModel.
relabelAll <- function( Model, add, sep = '_' ) {
  oldlabels <- names( omxGetParameters( Model ) )
  newlabels <- paste( oldlabels, add, sep = sep )
  omxSetParameters( Model, oldlabels, newlabels = newlabels )
}

#' This function extracts the fit statistiscs from the summary of a fitted OpenMx model.
extrFitStats<-function( mxSummary, nd = 2 ){
  # name of the model
  Name                <- mxSummary[ 'modelName' ][[ 1 ]]

  # loglikelihood results
  Min2LL              <- mxSummary[ 'Minus2LogLikelihood' ] [[ 1 ]]
  df                  <- mxSummary[ 'degreesOfFreedom' ] [[ 1 ]]
  nep                 <- mxSummary[ 'estimatedParameters' ] [[ 1 ]]
  dMin2LL             <- round( mxSummary[ 'Chi' ] [[ 1 ]], nd )
  ddf                 <- mxSummary[ 'ChiDoF' ] [[ 1 ]]
  p                   <- round( mxSummary[ 'p' ] [[ 1 ]], nd )
  lltest              <- data.frame( Name, Min2LL, df, nep, dMin2LL, ddf, p )

  # Information criteria
  infCrit             <- t( as.data.frame( unlist ( mxSummary[ 'informationCriteria' ] )[ c( 3, 4, 6 ) ] ) )
  colnames( infCrit ) <- c( 'AIC.nep', 'BIC.nep', 'SABIC' )

  # Alternative Fit Measures
  CFI                 <- round( mxSummary[ 'CFI' ] [[ 1 ]], nd )
  NNFI                <- round( mxSummary[ 'TLI' ] [[ 1 ]], nd )
  rmsea               <- round( mxSummary[ 'RMSEA' ] [[ 1 ]], nd )
  rmsea.CI95L         <- round( mxSummary[ 'RMSEACI' ] [[ 1 ]] [ 1 ], nd+1 )
  rmsea.CI95U         <- round( mxSummary[ 'RMSEACI' ] [[ 1 ]] [ 2 ], nd+1 )
  rmsea.p             <- round( mxSummary[ 'RMSEAClose' ] [[ 1 ]], nd+1 )
  altfit              <- data.frame( CFI, NNFI, rmsea, rmsea.CI95L, rmsea.CI95U, rmsea.p )

  # Bind together and return
  all                 <- data.frame( lltest, infCrit, altfit )
  rownames(all)       <- NULL
  return(all)
}



