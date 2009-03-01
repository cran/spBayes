hexGrid <- function(nodes.per.layer, b.box, ...){
   
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(nodes.per.layer)){stop("error: nodes.per.layer must be specified")}
  if(nodes.per.layer < 1){stop("error: nodes.per.layer must be >=1")}

  if(missing(b.box)){stop("error: bbox must be specified")}  
  if(!is.vector(b.box) || length(b.box) != 4) {stop("error: bbox must be a vector of length 4")}

  storage.mode(nodes.per.layer) <- "integer"
  storage.mode(b.box) <- "double"

  ##check the bounding box
  if(b.box[1] >= b.box[3]){stop("error: b.box x.min >= x.max")}
  if(b.box[2] >= b.box[4]){stop("error: b.box y.min >= y.max")}
  
  out <- .Call("hexGrid", nodes.per.layer, b.box);
  out$hex.centroids <- matrix(out$hex.centroids, ncol=2, byrow=TRUE)

  ##get hex edges
  n <- nrow(out$hex.centroids)
  hex.polygons <- vector("list", n)

  for(i in 1:n){
    dx <- out$hx.hy[1]/2
    dy <- dx/sqrt(3)
    hex.polygons[[i]] <- cbind(out$hex.centroids[i,1]+c(dx, dx,     0, -dx, -dx,    0),
                               out$hex.centroids[i,2]+c(dy,-dy, -2*dy, -dy,  dy, 2*dy))
  }

  out$hex.polygons <- hex.polygons
  out
}
