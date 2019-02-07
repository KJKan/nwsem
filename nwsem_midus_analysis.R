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
require(nwsem)

# clear work space, not run
rm(list = ls())

# load some helper functions
source("https://raw.githubusercontent.com/KJKan/nwsem/master/nwsem_helper_functions.R")

# load MIDUS data
con <- url( "https://raw.githubusercontent.com/KJKan/nwsem/master/DataYounger.Rdata" )
load( con )
close( con)
con <- url( "https://raw.githubusercontent.com/KJKan/nwsem/master/DataMiddle.Rdata" )
load( con )
close( con)
con <- url( "https://raw.githubusercontent.com/KJKan/nwsem/master/DataOlder.Rdata" )
load( con )
close( con)

# The MIDUS factor pattern matrix, see function readpat ( HelperFunctions )
MIDUS.pat <- readpat(2)
1 0
* 0
0 1
0 *
0 *
0 *
0 *
  
# number of first order factors
ne <- ncol( MIDUS.pat ) 

# The variable names
varnames <- colnames( DataYounger )

# The number of variables
ny <- length( varnames )


### Unigroup modeling

# ------------------------------ OpenMx Explicit Saturated Model

# The data as OpenMx object
DataY   <- mxData( DataYounger, 
                   type = 'raw')

DataM   <- mxData( DataMiddle, 
                   type = 'raw')

DataO   <- mxData( DataOlder, 
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
                       values = cov( DataYounger  ),
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
SatFitY <- mxRun( mxModel( name = 'Saturated Model Younger', SatModel, DataY ) )
SatFitM <- mxRun( mxModel( name = 'Saturated Model Middle',  SatModel, DataM ) )
SatFitO <- mxRun( mxModel( name = 'Saturated Model Older',   SatModel, DataO ) )

# Reference models
RefModelY <- mxRefModels( SatFitY, run = TRUE )
RefModelM <- mxRefModels( SatFitM, run = TRUE )
RefModelO <- mxRefModels( SatFitO, run = TRUE )

# Result 
SatResY <- summary( SatFitY, refModels = RefModelY )  
SatResM <- summary( SatFitM, refModels = RefModelM )  
SatResO <- summary( SatFitO, refModels = RefModelO )



# ------------------------ Build Measurement Model

# MxMatrix containing the factor loadings 
Lambda <- mxMatrix( name = 'Lambda',
                    type = 'Full',
                    nrow = ny,
                    ncol = ne,
                    free = MIDUS.pat == '*',
                    values = ( MIDUS.pat != '0' )*1,
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

# Fit the model on the three data sets
MeasFitY <- mxTryHard( mxModel( name = 'Measurement Model Younger', MeasModel, DataY ) )
MeasFitM <- mxTryHard( mxModel( name = 'Measurement Model Middle',  MeasModel, DataM ) )
MeasFitO <- mxTryHard( mxModel( name = 'Measurement Model Older',   MeasModel, DataO ) )

# Result 
MeasResY <- summary( MeasFitY, refModels = RefModelY )  
MeasResM <- summary( MeasFitM, refModels = RefModelM )
MeasResO <- summary( MeasFitO, refModels = RefModelO )



# ------------- Build Network Model

# Extract Network in subsample 1, e.g. 
pmat <- ggmModSelect( cor( DataYounger ), nrow( DataYounger) )$graph

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
NWFitY <- mxTryHard( mxModel( name = 'Network Younger', NWModel, DataY ) )
NWFitM <- mxTryHard( mxModel( name = 'Network Middel',  NWModel, DataM ) )
NWFitO <- mxTryHard( mxModel( name = 'Network Older',   NWModel, DataO ) )

# Result 
NWResY <- summary( NWFitY, refModels = RefModelY )  
NWResM <- summary( NWFitM, refModels = RefModelM )
NWResO <- summary( NWFitO, refModels = RefModelO )



# ------------- Collect results

# Fits statistics
fitsY <- rbind( extrFitStats( SatResY  ),
                extrFitStats( MeasResY ),
                extrFitStats( NWResY   ) )
fitsM <- rbind( extrFitStats( SatResM  ),
                extrFitStats( MeasResM ),
                extrFitStats( NWResM   ) )
fitsO <- rbind( extrFitStats( SatResO  ),
                extrFitStats( MeasResO ),
                extrFitStats( NWResO   ) )

# Print
print( fitsY )                                
print( fitsM ) 
print( fitsO )


### Multigroup modeling


# ------------- Build Network Model


# Copy from the unigroup model and give unique labels
NWFitY  <- relabelAll( NWFitY , add = 1 )
NWFitY  <- mxModel( name = 'Younger', NWFitY )

NWFitM <- relabelAll( NWFitM, add = 2 )
NWFitM <- mxModel( name = 'Middle',   NWFitM )

NWFitO <- relabelAll( NWFitO, add = 3 )
NWFitO <- mxModel( name = 'Older',    NWFitO )

# Multigroup fit function
mgFitFun <- mxFitFunctionMultigroup( c( "Younger", "Middle", "Older" ) )

# The multigroup network model
NWMGModel <- mxModel( name = 'Multigroup Network', 
                      mgFitFun, 
                      NWFitY, 
                      NWFitM, 
                      NWFitO )

# Fit the model
NWMGFit <- mxRun( NWMGModel )



# ------------- Build the (Omega) Invariant Network Model

# Copy the multigroup model
NWMGModel_Omega_Inv <- mxModel( NWMGFit, name = 'Omega invariant' )

# Add constraints by giving parameters in Omega the same label
NWMGModel_Omega_Inv <- omxSetParameters( NWMGModel_Omega_Inv, 
                                         labels = NWFitY$Omega$labels, 
                                         newlabels = label( 'omega', 7 ) )
NWMGModel_Omega_Inv <- omxSetParameters( NWMGModel_Omega_Inv, 
                                         labels = NWFitM$Omega$labels, 
                                         newlabels = label( 'omega', 7 ) )
NWMGModel_Omega_Inv <- omxSetParameters( NWMGModel_Omega_Inv, 
                                         labels = NWFitO$Omega$labels, 
                                         newlabels = label( 'omega', 7 ) )

# Prevent starting value problems
NWMGModel_Omega_Inv  <- omxAssignFirstParameters( NWMGModel_Omega_Inv  )

# Fit the multigroup model
NWMGFit_Omega_Inv   <- mxRun( NWMGModel_Omega_Inv  )



# ------------- Build a Saturated Model

# Copy the multigroup network model
# Let all parameters of the network be free
SatFitY <- NWFitY
SatFitY <- omxSetParameters( SatFitY, 
                             labels = SatFitY$Omega$labels, 
                             newlabels = label( 'omega', 7, add = 1 ), 
                             free = !diag( ny ) )

SatFitM <- NWFitM
SatFitM <- omxSetParameters( SatFitM, 
                             labels = SatFitM$Omega$labels, 
                             newlabels = label( 'omega', 7, add = 1 ), 
                             free = !diag( ny ) )

SatFitO <- NWFitO
SatFitO <- omxSetParameters( SatFitO, 
                             labels = SatFitO$Omega$labels, 
                             newlabels = label( 'omega', 7, add = 1 ), 
                             free = !diag( ny ) )

# Define multigroup model
SatMGModel <- mxModel( name = 'multigroup saturated', 
                       mgFitFun, 
                       SatFitY, 
                       SatFitM, 
                       SatFitO )

# Prevent starting value problems
SatMGModel  <- omxAssignFirstParameters( SatMGModel  )

# Fit the model
SatMGFit <- mxRun( SatMGModel )

# Define it the reference model
refModel <- mxRefModels(SatMGFit, run = TRUE )

# print fit statistics
rbind( extrFitStats( summary (SatMGFit,          refModels = refModel ) ),
       extrFitStats( summary (NWMGFit,           refModels = refModel ) ),
       extrFitStats( summary (NWMGFit_Omega_Inv, refModels = refModel ) ))

# Test of the invariance holds
mxCompare( NWMGFit, NWMGFit_Omega_Inv ) # no



### plot results

# Omega estimates in the three groups
OmegaY <- mxEval( Omega, NWFitY )
OmegaM <- mxEval( Omega, NWFitM )
OmegaO <- mxEval( Omega, NWFitO )

# Add labels
colnames( OmegaY ) <- colnames( DataYounger )
colnames( OmegaY ) <- colnames( DataMiddle  )
colnames( OmegaY ) <- colnames( DataOlder   )

# Use the multigroup invariant network for the layout
q <- mxEval( Younger.Omega, NWMGFit_Omega_Inv )

# Layout
lo <- qgraph( q, layout= 'spring', 
              DoNotPlot = TRUE )$layout

# group the variables (according to the two latent dimensions)
groups <- list( c(1:2), 
                c(3:7) )

# plot the networks using the above layout
layout( matrix(1:3,1,3))

qgraph( OmegaY,
        layout=lo, 
        groups = groups, 
        vsize = 14, 
        color = c('yellow','lightblue') )
title( 'Network in the group younger\n(age 35-54, n = 1876)')

qgraph( OmegaY,
        layout=lo, 
        groups = groups, 
        vsize = 14, 
        color = c('yellow','lightblue'))
title( 'Network in the group middle\n(age 55-64, n = 1000)')

qgraph( OmegaY,
        layout=lo, 
        groups = groups, 
        vsize = 14, 
        color = c('yellow','lightblue'))
title ( 'Network in the group older\n(age 65 and older, n = 903)')

