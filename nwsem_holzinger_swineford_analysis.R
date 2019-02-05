###################################################################
#                                                                 #
#      Kan, K. J., van der Maas, L. J. & Levine, S. Z. (2019).    #
#           Extending psychometric network analysis:              #
#      Empirical evidence against g in favor of mutualism?        #
#                   Intelligence, 73, 52-62.                      #
#          https://doi.org/10.1016/j.intell.2018.12.004           #
#                                                                 #
# This script (re-)analyzes The Holzinger-Swineford (1939) data   #
#                                                                 #
# The sample is split in half (randomly), after which             #
# In OpenMx                                                       #
#  - an explicit saturated model is fitted                        #
#  - a measurement model is fitted                                #
#  - a higher order g factor model is fitted                      #
#  - a network extracted through qgraph is fitted                 #
#  in subsample 1                                                 #
#                                                                 #
#  Next the models are refitted in sample 2                       #
#                                                                 #
# Author: Kees-Jan Kan                                            #
#                                                                 #
###################################################################


# load required R packages
require( OpenMx )
require( qgraph )
require( lavaan )

# clear work space, not run
rm(list = ls())

# load some helper functions
source("https://raw.githubusercontent.com/KJKan/nwsem/master/nwsem_helper_functions.R")

# load data (from R package lavaan)
data( HolzingerSwineford1939 )

# sample size
n <- nrow( HolzingerSwineford1939 )

# split the sample 
# seed to replicate the results in our paper was 17
# N.B: Using other seeds may 
#      - sometimes result in better fits for the g model 
#      - sometimes result in better fits for the network model
#      Conclusion: the two are close (better use large samples)
set.seed( 17 ) 
sample1 <-  sample( 1:n, n/2 )
sample2 <- -sample1

# make selection of the data 
HZ.dat1 <- HolzingerSwineford1939[ sample1, -(1:6) ] # colums 1 to 6 contain irrelevant data for the analysis 
HZ.dat2 <- HolzingerSwineford1939[ sample2, -(1:6) ] # and are removed

# The HolzingerSwineford factor pattern matrix, see function readpat ( HelperFunctions )
HZ.pat <- readpat(3)
1 0 0
* 0 0
* 0 0
0 1 0
0 * 0
0 * 0
0 0 1
0 0 *
0 0 *

# number of first order factors
ne <- ncol( HZ.pat ) 

# The variable names
varnames <- colnames( HZ.dat1 )

# The number of variables
ny <- length( varnames )


# ------------------------------ OpenMx Explicit Saturated Model

# The data as OpenMx object
Data1   <- mxData( HZ.dat1, 
                   type = 'raw')

Data2   <- mxData( HZ.dat2, 
                   type = 'raw')

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
                       values = cov( HZ.dat1 ),
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
                            means = 'ExpMean',
                            dimnames = varnames )

# The complete model
SatModel <- mxModel( name = 'Saturated Model',
                     ExpMean, 
                     ExpCovSat,
                     I,
                     ExpCor,
                     Obj,
                     mxFitFunctionML() )

# Fit the model on the two data sets
SatFit1 <- mxRun( mxModel( name = 'Saturated Model Subsample 1', SatModel, Data1 ) )
SatFit2 <- mxRun( mxModel( name = 'Saturated Model Subsample 2', SatModel, Data2 ) )
                 
# Reference models
RefModel1 <- mxRefModels( SatFit1, run = TRUE )
RefModel2 <- mxRefModels( SatFit2, run = TRUE )

# Result 
SatRes1 <- summary( SatFit1, refModels = RefModel1 )  
SatRes2 <- summary( SatFit2, refModels = RefModel2 )


# ------------------------ Build Measurement Model

# MxMatrix containing the factor loadings 
Lambda <- mxMatrix( name = 'Lambda',
                    type = 'Full',
                    nrow = ny,
                    ncol = ne,
                    free = HZ.pat == '*',
                    values = ( HZ.pat != '0' )*1,
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
                     ExpMean,
                     Lambda,
                     Psi,
                     Theta,
                     ExpCov,
                     I,
                     ExpCor,
                     Obj,
                     mxFitFunctionML() )

# Fit the model on the two data sets
MeasFit1 <- mxRun( mxModel( name = 'Measurement Model Subsample 1', MeasModel, Data1 ) )
MeasFit2 <- mxRun( mxModel( name = 'Measurement Model Subsample 2', MeasModel, Data2 ) )

# Result 
MeasRes1 <- summary( MeasFit1, refModels = RefModel1 )  
MeasRes2 <- summary( MeasFit2, refModels = RefModel2 )


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
                  ExpMean, 
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

# Fit the model on the two data sets
gFit1 <- mxRun( mxModel( name = 'g Model Subsample 1', gModel, Data1 ) )
gFit2 <- mxRun( mxModel( name = 'g Model Subsample 2', gModel, Data2 ) )

# Result 
gRes1 <- summary( gFit1, refModels = RefModel1 )  
gRes2 <- summary( gFit2, refModels = RefModel2 )



# ------------- Build Network Model

# Extract Network in subsample 1, e.g. 
pmat <- ggmModSelect( cor( HZ.dat1 ), nrow( HZ.dat1 ) )$graph

# Apparently I forgot to use ggmModSelect and still used the old method ( EBICglasso )
# So to replicate the results in the paper, use
pmat <- EBICglasso( cor( HZ.dat1 ), nrow( HZ.dat1 ) )

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
                   I,
                   Delta,
                   Omega,
                   ExpMean, 
                   ExpCovNW,
                   ExpPcor,
                   ExpCor,
                   Obj,
                   mxFitFunctionML() )

# Fit the model on the two data sets
NWFit1 <- mxRun( mxModel( name = 'Network Subsample 1', NWModel, Data1 ) )
NWFit2 <- mxRun( mxModel( name = 'Network Subsample 2', NWModel, Data2 ) )

# Result 
NWRes1 <- summary( NWFit1, refModels = RefModel1 )  
NWRes2 <- summary( NWFit2, refModels = RefModel2 )



# ------------- Collect results

# Fits statistics
fits1 <- rbind( extrFitStats( SatRes1  ),
                extrFitStats( MeasRes1 ),
                extrFitStats( gRes1    ),
                extrFitStats( NWRes1   ) )
fits2 <- rbind( extrFitStats( SatRes2  ),
                extrFitStats( MeasRes2 ),
                extrFitStats( gRes2    ),
                extrFitStats( NWRes2   ) )

# Print
print( fits1 )                                
print( fits2 )    

