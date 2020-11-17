
#' Cluster Resampling: resampling long format longitudinal or otherwise clustered data
#'
#' @param data a data frame containing the variables in the model.
#' @param cluster.name the name (as a character) of the column containing the cluster identifiers.
#' @param size a non-negative integer giving the number of items to choose, i.e. the number of clusters to resample.
#'
#' @return returns a new dataframe with resampled clusters
#' @export
#'
#' @examples
#' table(cfd.example.data$id)
#' # every ID appears 5 times
#' cfd.example.sample <- cluster.resample(cfd.example.data, cluster.name='id',
#' size=length(unique(cfd.example.data[,'id'])))
#' table(cfd.example.sample$id)
#' # some ID's now don't appear, and some appear more times (multiples of 5)
#' # the important part is that if a person (id) is resampled, all their rows of data are taken
#' # i.e. the function resamples clusters, rather than rows.
#' # this has produced 1 resample, so generally this function would be used
#' # inside a loop where it is used multiple times
#' # we use this function inside our decomposition functions when cluster.sample is
#' # set to TRUE in those functions.
#'
cluster.resample <- function(data, cluster.name, size) {
  # select a bunch of IDs
  IDs <- unique(data[,cluster.name])
  y <- sample(IDs,size,replace=TRUE)
  z <- table(table(y))

  # from there, select a group once
  selectID <- sample(IDs,size=z[1],replace=FALSE)
  newdata <- data[which(data[,cluster.name] %in% selectID),]

  if(length(z) > 1) {

    for(i in 2:length(z)) {

      # select a new group of IDs that was not yet selected
      IDs2 <- setdiff(IDs,selectID)

      # from there, randomly select a group of people of the right size
      selectID2 <- sample(IDs2,size=z[i],replace=FALSE)
      selectID <- c(selectID,selectID2) # so we don't re-select the newly
      # selected people either

      for(j in 1:i) {

        # copy the new dataset i number of times
        newdata <- rbind(newdata,data[which(data[,cluster.name] %in% selectID2),])

      }
    }
  }
  return(newdata)
}
