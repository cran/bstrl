## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bstrl)

## ----showdata-----------------------------------------------------------------
length(geco_small)
head(geco_small[[1]])

## ----preparecomparisons-------------------------------------------------------
# Names of the columns on which to perform linkage
fieldnames <- c("given.name", "surname", "age", "occup", "extra1", "extra2", "extra3", "extra4", "extra5", "extra6")

# How to compare each of the fields
types <- c("lv", "lv", # First name and last name use normalized edit distance
           "bi", "bi", "bi", "bi", "bi", "bi", "bi", "bi") # All others binary equal/unequal
breaks <- c(0, 0.25, 0.5) # Break continuous difference measures into 4 levels using these split points


## ----streaming----------------------------------------------------------------
res.twofile <- bipartiteRL(geco_small[[1]], geco_small[[2]],
                           flds = fieldnames, types = types, breaks = breaks,
                           nIter = 600, burn = 100,
                           seed = 0)

## ----pprbupdate---------------------------------------------------------------
res.pprb3 <- PPRBupdate(res.twofile, geco_small[[3]], # Comparison details are stored with previous result
                       nIter = 600, burn = 100,
                       seed = 0,
                       refresh = 0.05)
res.pprb4 <- PPRBupdate(res.pprb3, geco_small[[4]], # Comparison details are stored with previous result
                       nIter = 600, burn = 100,
                       seed = 0,
                       refresh = 0.05)

## ----filtersamples------------------------------------------------------------
filtered <- thinsamples(res.twofile, 50) # Don't need 500, 50 for demonstration

## ----smcmcupdate--------------------------------------------------------------
res.smcmc3 <- SMCMCupdate(filtered, geco_small[[3]],
                         nIter.jumping=2, nIter.transition = 8,
                         proposals.jumping="component", proposals.transition="component", #Either can be LB, but increase corresponding number of iterations
                         cores = 2) # Parallel execution
res.smcmc4 <- SMCMCupdate(res.smcmc3, geco_small[[4]],
                         nIter.jumping=2, nIter.transition = 8,
                         proposals.jumping="component", proposals.transition="component", #Either can be LB, but increase corresponding number of iterations
                         cores = 2) # Parallel execution

## ----resultsobject------------------------------------------------------------
names(res.pprb4)

## ----objectsizes--------------------------------------------------------------
# All have 500 columns and different numbers of rows.
dim(res.pprb4$Z)
dim(res.pprb4$m)
dim(res.pprb4$u)

## ----examiningZ---------------------------------------------------------------
# The first post-burn MCMC sample of Z
Zexample <- res.pprb4$Z[,1]
Zexample

## ----examplelink--------------------------------------------------------------
Zexample[8]

## ----examplecluster-----------------------------------------------------------
Zexample[27]
Zexample[8]

## ----linkprocessing-----------------------------------------------------------
# Create a list of length 500, where each element is one streaming link object
# for each posterior sample.
samples <- extractlinks(res.pprb4)

# Are record 9 in file 1 and record 7 in file 4 linked in the first posterior sample?
islinked(samples[[1]], file1=1, record1=9, file2=4, record2=7)

# In what proportion of posterior samples are record 9 in file 1 and record 7 in file 4 linked?
mean(sapply(samples, islinked, file1=1, record1=9, file2=4, record2=7))

# In what proportion of posterior samples are record 8 in file 1 and record 1 in file 2 linked?
mean(sapply(samples, islinked, file1=1, record1=8, file2=2, record2=1))

## ----multifile----------------------------------------------------------------
# Link three files from scratch, returning 500 posterior samples
res.threefile <- multifileRL(geco_small[1:3],
                             flds = fieldnames, types = types, breaks = breaks,
                             nIter = 600, burn=100, # Number of iterations to run
                             proposals = "comp", # Change to "LB" for faster iterations, slower convergence
                             seed = 0,
                             refresh = 0.05) # Print progress every 5% of run.

