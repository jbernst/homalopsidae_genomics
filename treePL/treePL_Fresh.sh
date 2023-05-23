# To run treePL, you just run the following:
treePL Homalopsidae_Fresh_NoGEO.config

# Then run random subsample and replicate cross-validation (RSRCV) by adding 'randomcv' at the end of the configuration file.
# After you do this, you will get a cv.out file, and pick the best smoothing parameter (the one with the lowest error, i.e., the lowest value in the parentheses). 
# Finally, rerun the analysis above with prime commented out, randomcv commented out, and the smooth parameter adjusted to the best option in cv.out.
