
#' Cluster Resampling: resampling long format longitudinal or otherwise clustered data
#'
#' @param data a data frame containing the variables in the model.
#' @param cluster.name the name (as a character) of the column containing the cluster identifiers.
#' @param size a non-negative integer giving the number of items to choose, i.e. the number of clusters to resample. If not specified, takes the original data cluster size.
#' @param newID if set to TRUE, gives each replicant a new ID, rather than having the old ID appear multiple times. This might be needed if you work with e.g. cluster fixed effects
#'
#' @return returns a new dataframe with resampled clusters
#' @export
#'
#' @examples
#' table(cfd.example.data$id)
#' # every ID appears 5 times
#' cfd.example.sample <- cluster.resample(cfd.example.data, cluster.name='id')
#' table(cfd.example.sample$id)
#' # some ID's now don't appear, and some appear more times (multiples of 5)
#' # the important part is that if a person (id) is resampled, all their rows of data are taken
#' # i.e. the function resamples clusters, rather than rows.
#' # this has produced 1 resample, so generally this function would be used
#' # inside a loop where it is used multiple times
#' # we use this function inside our decomposition functions when cluster.sample is
#' # set to TRUE in those functions.
#'
cluster.resample <- function(data,cluster.name,size=NA,newID=FALSE) {

  # 1. sample cluster ids with replacement
  IDs <- unique(data[cluster.name])[,1]
  # if size has not been specified, use original cluster size
  if(is.na(size)) {size <- length(IDs)}
  y <- table(sample(IDs,size,replace=T))

  # 2. determine how many times a cluster id appears, save cluster id and how often it appeared
  id.rep <- data.frame(t(y))[,-1]
  names(id.rep) <- c(cluster.name,'freq')

  # 3. get all the observations for each cluster id in id.rep - # then right join id.rep onto data
  data.select <- data[which(data[cluster.name][,1] %in% id.rep[,1]),]
  data.merge <- merge(data.select,id.rep,all.x=TRUE)

  # remove unnecessary data from memory
  rm(data.select)

  # 4. use lapply to copy the rows
  data.sample <- as.data.frame(lapply(data.merge, rep, data.merge$freq))

  if(newID==TRUE) {
    # 5. label the replicates
    data.sample$replicate <- unlist(lapply(data.merge$freq,function(x) if(x > 0) 1:x))

    # 6. each clusterid & replication combination is the new unique id
    data.sample$newID <- paste(data.sample[cluster.name][,1],data.sample$replicate,sep='')

    # 7. order by newid
    data.sample[order(data.sample$newID),]

    # remove unneccesary data from memory
    data.sample$freq <- NULL
    data.sample$replicate <- NULL
  }
  rm(data.merge)

  # output the result
  return(data.sample)
}
