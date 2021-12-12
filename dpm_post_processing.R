###READ .DPM SAMPLE FILE###
#---------------------------

#Converts the a .dpm file produced by sampling a surface following particle tracking in Ansys fluent into an R dataframe

read.dpm <- function(filepath){
  df <- readLines(filepath)
  df <- strsplit(gsub('[()]', '', df)[-1L], "\\s+")
  names <- df[[1]]
  df <- data.frame(t(sapply(df[2:length(df)], c)))
  df <- setNames(df, names)
  df <- data.frame(apply(df[, 2:13], 2, as.numeric))
  return(df)
}

#Groups data into n bins, operates on dataframe produced from above. X-variable (ordered) is split into bins, into which y variables are placed.

bin_data <- function(df, x.var, y.var, bins){
  df <- df[order(df[,x.var]),]
  df$b <- cut(df[,x.var], bins)
  u <- split(df, df$b)
  newy <-  unlist(lapply(u, function(X){median(X[,y.var])}))
  newx <-  unlist(lapply(u, function(X){median(X[,x.var])}))
  newxy <- data.frame(cbind(newx, newy))
  return(newxy)
}