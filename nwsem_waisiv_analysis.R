###################################################################
#                                                                 #
#      Kan, K. J., van der Maas, L. J. & Levine, S. Z. (2019).    #
#           Extending psychometric network analysis:              #
#      Empirical evidence against g in favor of mutualism?        #
#                   Intelligence, 73, 52-62.                      #
#          https://doi.org/10.1016/j.intell.2018.12.004           #
#                                                                 #
# This script (re-)analyzes the WAIS IV (2008) correlation matrix #
#                                                                 #
# In OpenMx                                                       #
#  - an explicit saturated model is fitted                        #
#  - a measurement model is fitted                                #
#  - a higher order g factor model is fitted                      #
#  - a network extracted through qgraph is fitted                 #
#                                                                 #
# Author: Kees-Jan Kan                                            #
#                                                                 #
###################################################################


# load required R packages
require( OpenMx )
require( qgraph )

# clear work space, not run
rm(list = ls())

# load some helper functions
source("https://raw.githubusercontent.com/KJKan/nwsem/master/nwsem_helper_functions.R")

# The variable names
varnames <- c( 'SI', # Similarities
               'VC', # Vocabulary
               'IN', # Information
               'CO', # Comprehension
               'BD', # Block Design
               'MR', # Matrix Reasoning
               'VP', # Visual Puzzles
               'FW', # Figure Weights
               'PCm',# Picture Completion
               'DS', # Digit Span
               'AR', # Arithmetic
               'LN', # Letter-Number Sequencing
               'SS', # Symbol Search
               'CD', # Coding
               'CA') # Cancellation

# The number of variables
ny <- length ( varnames )

# Reading in the correlation matrix
WAISCorMat <- matrix( scan(), 
                          ny, ny, 
                          byrow = TRUE, 
                          dimnames = list( varnames, varnames) )
1.00 0.74 0.64 0.71 0.49 0.51 0.44 0.53 0.44 0.48 0.54 0.45 0.35 0.41 0.23
0.74 1.00 0.73 0.74 0.45 0.51 0.42 0.53 0.39 0.50 0.57 0.48 0.34 0.41 0.24
0.64 0.73 1.00 0.66 0.44 0.49 0.43 0.51 0.41 0.43 0.57 0.43 0.34 0.34 0.22
0.71 0.74 0.66 1.00 0.44 0.49 0.43 0.53 0.40 0.48 0.55 0.47 0.32 0.39 0.21
0.49 0.45 0.44 0.44 1.00 0.54 0.64 0.56 0.49 0.45 0.50 0.42 0.41 0.40 0.34
0.51 0.51 0.49 0.49 0.54 1.00 0.53 0.57 0.42 0.47 0.52 0.45 0.39 0.45 0.26
0.44 0.42 0.43 0.43 0.64 0.53 1.00 0.58 0.48 0.40 0.48 0.41 0.38 0.37 0.32
0.53 0.53 0.51 0.53 0.56 0.57 0.58 1.00 0.41 0.50 0.61 0.48 0.34 0.36 0.29
0.44 0.39 0.41 0.40 0.49 0.42 0.48 0.41 1.00 0.39 0.37 0.37 0.41 0.38 0.33
0.48 0.50 0.43 0.48 0.45 0.47 0.40 0.50 0.39 1.00 0.60 0.69 0.40 0.45 0.34
0.54 0.57 0.57 0.55 0.50 0.52 0.48 0.61 0.37 0.60 1.00 0.46 0.37 0.43 0.31
0.45 0.48 0.43 0.47 0.42 0.45 0.41 0.48 0.37 0.69 0.46 1.00 0.37 0.38 0.30
0.35 0.34 0.34 0.32 0.41 0.39 0.38 0.34 0.41 0.40 0.37 0.37 1.00 0.65 0.46
0.41 0.41 0.34 0.39 0.40 0.45 0.37 0.36 0.38 0.45 0.43 0.38 0.65 1.00 0.42
0.23 0.24 0.22 0.21 0.34 0.26 0.32 0.29 0.33 0.34 0.31 0.30 0.46 0.42 1.00

# Defining the mean vector
# Note: Not neccesary to model the means
#       If not modeled, subtract ny dfs of freedom
WAISMeansVec <- rep( 0, ny )
names( WAISMeansVec ) <- varnames

# The sample size
n <- 1800

# The WAIS factor pattern matrix, see function readpat ( HelperFunctions )
WAIS.pat <- readpat( 4 )
1 0 0 0
* 0 0 0
* 0 0 0
* 0 0 0
0 1 0 0
0 * 0 0
0 * 0 0
0 * 0 0
0 * 0 0
0 0 1 0
0 0 * 0
0 0 * 0
0 0 0 1
0 0 0 *
0 0 0 *

# number of first order factors
ne <- ncol( WAIS.pat )


# ------------------------------ OpenMx Explicit Saturated Model

# The data as OpenMx object
Data   <- mxData( WAISCorMat, 
                  type = 'cov', 
                  means = WAISMeansVec,
                  numObs = n )

# Mean matrix
ExpMean <- mxMatrix( name = 'ExpMean',
                     type = 'Full',
                     nrow = 1,
                     ncol = ny,
                     free = TRUE,
                     values = rep( 0, ny ) ,
                     labels = label( 'm', ny ) )

# Covariance matrix
ExpCovSat <- mxMatrix( name = 'ExpCov',
                       type = 'Symm',
                       nrow = ny,
                       free = TRUE,
                       values = WAISCorMat,
                       labels = label( 's', ny, tri = TRUE ) )

# ny by ny identity matrix
I <- mxMatrix( name = 'I',
               type = 'Iden',
               nrow = ny )

# Correlation matrix
ExpCor <- mxAlgebra( name = 'ExpCor',
                     expression = solve( sqrt( I * ExpCov ) ) %*% ExpCov %*% t( solve( sqrt( I * ExpCov ) ) ) )

# Objective
Obj <- mxExpectationNormal( covariance = 'ExpCov',
                            means = 'ExpMean', # Remove this line when not modeling the mean
                            dimnames = varnames )

# The complete model
SatModel <- mxModel( name = 'Saturated Model',
                     Data, 
                     ExpMean, # Remove this line when not modeling the mean
                     ExpCovSat,
                     I,
                     ExpCor,
                     Obj,
                     mxFitFunctionML() )

# Fit the model
SatFit <- mxRun( SatModel )

# Result 
SatRes <- summary( SatFit )  


# ------------------------ Build Measurement Model

# Mean matrix
ExpMean <- mxMatrix( name = 'ExpMean',
                     type = 'Full',
                     nrow = 1,
                     ncol = ny,
                     free = FALSE,
                     values = rep( 0, ny ) ,
                     labels = label( 'm', ny ) )

# MxMatrix containing the factor loadings 
Lambda <- mxMatrix( name = 'Lambda',
                    type = 'Full',
                    nrow = ny,
                    ncol = ne,
                    free = WAIS.pat == '*',
                    values = ( WAIS.pat != '0' )*1,
                    labels = label( 'lambda', ny, ne ) )

# MxMatrix containing the latent variances and covariances
Psi <- mxMatrix( name = 'Psi',
                 type = 'Symm',
                 nrow = ne,
                 free = TRUE,
                 values = diag( ne ),
                 labels = label ( 'psi', ne ) )

# MxMatrix containing the residual variances of the observed variables
Theta <- mxMatrix( name = 'Theta',
                   type = 'Diag',
                   nrow = ny,
                   free = TRUE,
                   values = 1,
                   labels = label ( 'theta', ny ) )

# The factor model implied variance-covariance matrix
ExpCov <- mxAlgebra( name = 'ExpCov',
                     expression = Lambda %*% Psi %*% t( Lambda ) + Theta )

# The complete model
MeasModel <-mxModel( name = 'Measurement Model',
                     Data,
                     ExpMean, # Remove this line when mean is not considered (see Saturated Model)
                     Lambda,
                     Psi,
                     Theta,
                     ExpCov,
                     I,
                     ExpCor,
                     Obj,
                     mxFitFunctionML() )

# Fit the model 
MeasFit <- mxRun( MeasModel )

# Result
MeasRes <- summary( MeasFit )


# ------------------------ Build Second-Order g-Factor Model

# MxMatrix containing factor loadings (second order)
Gamma <- mxMatrix( name = 'Gamma',
                   type = 'Full',
                   nrow = ne,
                   ncol = 1,
                   free = TRUE,
                   values = 1,
                   labels = label( 'gamma', ne, 1 ) )

# MxMatrix containing the variance of the general factor
Phi <- mxMatrix( name = 'Phi',
                 type = 'Symm',
                 nrow = 1,
                 free = FALSE,
                 values = 1,
                 labels = 'phi' )

# MxMatrix containing the residual variances of the first order factors
Psi <- mxMatrix( name = 'Psi',
                 type = 'Diag',
                 nrow = ne,
                 free = TRUE,
                 values = diag(ne),
                 labels = label ( 'psi', ne ) )

# The factor model implied variance-covariance matrix
ExpCovg <- mxAlgebra( name = 'ExpCov',
                      expression = Lambda %*% ( Gamma %*% Phi %*% t( Gamma ) + Psi ) %*% t( Lambda ) + Theta )

# The complete model
gModel <-mxModel( name = 'g Model',
                  Data,
                  ExpMean, # Remove this line when mean is not considered (see Saturated Model)
                  Lambda,
                  Gamma,
                  Phi,
                  Psi,
                  Theta,
                  ExpCovg,
                  I,
                  ExpCor,
                  Obj,
                  mxFitFunctionML() )


# Fit the model 
gFit <- mxRun( gModel )

# Result
gRes <- summary( gFit )



# ------------- Build Network Model

pmat <- ggmModSelect ( WAISCorMat, n )$graph

# Matrix containing the scaling parameters
Delta <- mxMatrix( name = 'Delta',
                   type = 'Diag',
                   nrow = ny,
                   ncol = ny,
                   free = TRUE,
                   values = 1,
                   labels = label ('delta', ny ) )

# Matrix containing the partial relations (except the diagonal contains zeroes)
Omega <- mxMatrix( name = 'Omega',
                   type = 'Symm',
                   nrow = ny,
                   ncol = ny,
                   free = pmat!=0,
                   values = 0,
                   labels = label( 'omega', ny ) )

# Expected partial correlation Matrix
ExpPcor <- mxAlgebra( name = 'ExpPcor',
                      expression = Omega + I )

# Expected Covariance matrix
ExpCovNW <- mxAlgebra( name = 'ExpCov',
                       expression = Delta %*% solve( I - Omega ) %*% t( Delta ) )

# Model implied correlation matrix
ExpCor <- mxAlgebra( name = 'ExpCor',
                     expression = solve( sqrt( I * ExpCov ) ) %*% ExpCov %*% t( solve( sqrt( I * ExpCov ) ) ) )


# The complete model
NWModel <- mxModel(name = 'Network Model',
                   Data,
                   I,
                   Delta,
                   Omega,
                   ExpMean, # Remove this line when mean is not considered (see Saturated Model)
                   ExpCovNW,
                   ExpPcor,
                   ExpCor,
                   Obj,
                   mxFitFunctionML() )


# Fit the model 
NWFit <- mxRun( NWModel )

# Result
NWRes <- summary( NWFit )



# ------------- Collect results

# Fits statistics
fits    <- rbind( extrFitStats( SatRes  ),
                  extrFitStats( MeasRes ),
                  extrFitStats( gRes    ),
                  extrFitStats( NWRes   ) )

# Print
print( fits )                                
