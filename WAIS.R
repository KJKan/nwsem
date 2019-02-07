###################################################################
#                                                                 #
#      Kan, K. J., van der Maas, L. J. & Levine, S. Z. (2019).    #
#           Extending psychometric network analysis:              #
#      Empirical evidence against g in favor of mutualism?        #
#                   Intelligence, 73, 52-62.                      #
#          https://doi.org/10.1016/j.intell.2018.12.004           #
#                                                                 #
# This script contains WAIS-IV details                            #
#                                                                 #
# Author: Kees-Jan Kan                                            #
#                                                                 #
###################################################################

# Sample size
n            <- 1800

# Raw data? (no; correlation matrix)
raw          <- FALSE

# The information is stored in the file WAISData.R
dat <- read.table( "WAISData.R", header = TRUE ) 

# Number of subtests/variables 
ny <- ncol( dat )

# Pattern of first order factor loadings 
# ( 0 = fixed at 0; 1 = fixed at 1; rest = free )
lambda_g <- lambdaPat_g <- matrix( c(
                                   1, 0, 0,  0,
                                   2, 0, 0,  0,
                                   3, 0, 0,  0,
                                   4, 0, 0,  0,
                                   0, 1, 0,  0,
                                   0, 5, 0,  0,
                                   0, 6, 0,  0,
                                   0, 7, 0,  0,
                                   0, 8, 0,  0,
                                   0, 0, 1,  0,
                                   0, 0, 9,  0,
                                   0, 0, 10, 0,
                                   0, 0, 0,  1,
                                   0, 0, 0,  11,
                                   0, 0, 0,  12
                                   ), nrow = ny, byrow = TRUE )

# Pattern of second order factor loadings
gamma_g <- gammaPat_g  <- matrix( c(
                                  1,
                                  13,
                                  14,
                                  15
                                  ),ncol = 1, byrow = TRUE )


                                  
                                  