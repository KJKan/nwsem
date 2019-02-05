###################################################################
#                                                                 #
#      Kan, K. J., van der Maas, L. J. & Levine, S. Z. (2019).    #
#           Extending psychometric network analysis:              #
#      Empirical evidence against g in favor of mutualism?        #
#                   Intelligence, 73, 52-62.                      #
#          https://doi.org/10.1016/j.intell.2018.12.004           #
#                                                                 #
# This script (re-)analyzes Spearman (1904)'s correlation matrix  #
#                                                                 #
# In OpenMx                                                       #
#  - A one-factor model is fitted, after which                    #
#        the model implied correlation matrix (C) is calculated   #
#  - A saturated symmetrical mutualism model (M) is fitted to C   #
#        yielding perfect fit                                     #
#                                                                 #
# The two models have indistinguishable observable implications   #
# The first is more parsimonious in term of number of parameters  #
# the latter with respect to the number of variables in the model #
#                                                                 #
# Author: Kees-Jan Kan                                            #
#                                                                 #
###################################################################

# load required R packages
require( OpenMx )
require( qgraph )

# clear work space, not run
rm(list = ls())

# The variable names
varnames <- c( 'Classics', 
               'French',
               'English',
               'Math',
               'Pitch',
               'Music' )

# The number of variables
ny <- length ( varnames )

# Reading in the correlation matrix
SpearmanCorMat <- matrix( scan(), 
                          ny, ny, 
                          byrow = TRUE, 
                          dimnames = list( varnames, varnames) )
1 .83 .78 .70 .66 .63
.83 1 .67 .67 .65 .57
.78 .67 1 .64 .54 .51
.70 .67 .64 1 .45 .51
.66 .65 .54 .45 1 .40
.63 .57 .51 .51 .40 1

# The sample size
n <- 33

# ------------------------------ OpenMx Factor analysis

# The data as OpenMx object
Data   <- mxData( SpearmanCorMat, 
                  type = 'cov', 
                  numObs = n )

# Matrix containing the factor loadings ( estimated freely )
Lambda <- mxMatrix( name = 'L',  
                    type = 'Full',  
                    nrow = ny, 
                    ncol = 1, 
                    free = TRUE, 
                    values = .5 )

# Matrix containing the factor variance ( fixed, to 0 )
Psi    <- mxMatrix( name = 'Ps', 
                    type = 'Symm', 
                    nrow = 1, 
                    free = FALSE, 
                    values = 1 )

# Matrix containing the residual variances ( estimated freely )
Theta  <- mxMatrix( name = 'Th', 
                    type = 'Diag', 
                    nrow = 6, 
                    free = TRUE, 
                    values = .1 )

# Identity matrix
I      <- mxMatrix( name = 'I',  
                    type = 'Iden', 
                    nrow = ny )

# Expected covariance matrix
ExpCov   <- mxAlgebra( name = 'S',   
                       expression = L %*% Ps %*% t( L ) + Th )

# Standardization

# Diagonal matrix containing the standard devations
SD       <- mxAlgebra( name = 'SD',  expression = sqrt( I * S ) )

# Inverse of that  matrix
invSD    <- mxAlgebra( name = 'iSD', expression = solve( SD ) )

# Standardized factor loadings
st_Lamba <- mxAlgebra( name ='st_L', iSD %*% L )

# Standardized residual variances
st_Theta <- mxAlgebra( name ='st_Th', iSD %*% Th %*% t( iSD) )

# Expected correlation matrix
ExpCor <- mxAlgebra( name = 'ExpCor',
                     expression = iSD %*% S %*% t( iSD) )

# The complete model
Model <- mxModel( name = '1 Factor Model',
                  Data,
                  Lambda,
                  Psi,
                  Theta,
                  ExpCov,
                  I,
                  SD,
                  invSD,
                  st_Lamba,
                  st_Theta,
                  ExpCor,
                  mxExpectationNormal( covariance = 'S',
                                       dimnames = varnames ),
                  mxFitFunctionML() ) 

# Fit the factor model
Fit <- mxTryHard( Model )

# Standardized factor loadings
L  <- mxEval( st_L,  Fit ) 

# Standardized residual variances
Th <- mxEval( st_Th, Fit ) 

# Model implied correlation matrix
C  <- mxEval( ExpCor, Fit ) 
dimnames( C ) <- list( varnames, varnames )

# Print standardized factor loadings
print( round( L, 2 ) ) # As reported in Fig 1a, p. 53




# ------------------------------ OpenMx Mutualism Matrix

# The data as OpenMx object
# Make sure the correlation matrix is regarded as symmetric  
ExpData <- mxData( ( C + t( C ) )/2 , type = 'cov', numObs = n  )

# Mutualism matrix M, symmetrical, free ( except diagonal, which contains zeroes )
M       <- mxMatrix( name = 'M', 
                     type = 'Symm', 
                     nrow = ny, 
                     free = !diag( ny ), 
                     values = 0 )

# Residual variances 
Psi_mut <- mxMatrix( name = 'Ps', 
                     type = 'Diag', 
                     nrow = ny, 
                     free = TRUE, 
                     values = .1, 
                     ubound =  1 )

# Expected covariance matrix
ExpCov_mut   <- mxAlgebra( name = 'S',   expression = solve( I - M ) %*% Ps %*% t( solve( I - M ) ) )

# Standardized mutualism weights 
st_M         <- mxAlgebra( name = 'st_M', expression = iSD %*% M %*% SD ) 

# Standardized residual variances
st_Psi       <- mxAlgebra( name ='st_Ps', iSD %*% Ps %*% t( iSD ) )

# The complete model
MutModel <- mxModel( name = 'Mutualism Model',
                     ExpData,
                     M,
                     Psi_mut,
                     ExpCov_mut,
                     I,
                     SD,
                     invSD,
                     st_M,
                     st_Psi,
                     ExpCor,
                     mxExpectationNormal( covariance = 'S',
                                          dimnames = varnames ),
                     mxFitFunctionML() )

# Fit the mutualism model
MutFit <- mxTryHardWideSearch( MutModel )

# Standardized mutualism matrix
MutMat <- mxEval( st_M, MutFit )
dimnames( MutMat ) <- dimnames( C )

# Standardized residual variances
resVar <- mxEval( st_Ps, MutFit )

# Model implied correlation matrix
Cmut  <- mxEval( ExpCor, MutFit ) 
dimnames( Cmut ) <- list( varnames, varnames )

# Print standardized matrices
print( round( MutMat, 2 ) )          # As reported in Fig 1b, p. 53
print( diag( round( resVar, 2 ) ) ) 


# ------------------------------ Results

# Model fitting results
print( summary( Fit ) )    # -2loglikelihood: 61.195; chi-square ( df=9 ) = 2.82,  p = 0.97
print( summary( MutFit ) ) # -2loglikelihood: 61.195 

# print correlation matrices
print( C )    # 1-factor model implied 
print( Cmut ) # mutualism model implied

# Conclusion: Ignoring rounding errors, the two matrices are the same
#             In terms of the number of parameters, the 1-factor model is more parsimoneous
#             In terms of the number of variables, the mutualism model is more parsimoneous
