# moving<-read.table('/Users/danielfurth/Downloads/cpd-master 2/tests/fixtures/fish_distorted.csv', sep=',')
# fixed<-read.table('/Users/danielfurth/Downloads/cpd-master 2/tests/fixtures/fish.csv', sep=',')
# 
# moving<-cbind(as.integer( (5+moving[,1])*1000) , as.integer((5+moving[,2])*1000) ) 
# fixed<-cbind(as.integer((5+fixed[,1])*1000 ), as.integer((5+fixed[,2])*1000) ) 
# 
# plot(rbind(fixed, moving))
# points(moving, pch=3)
# 
# 
# regi<-cpd(fixed, moving)

#' Coherent Point Drift 
#'
#' Perform registration between two sets of points using coherent point drift equation.
#' @param fixed fixed set of points. 
#' @param moving moving set of points to be matched with the fixed.
#' @param plot boolean, if the result should be displayed.
#' @param beta numeric, .
#' @param lambda numeric, .
#' @param gamma numeric, .
#' @param sigma numeric, .
#' @param max.iter integer, maximum number of iterations to perform.
#' @keywords registration
#' @export
#' @examples
#' regi<-cpd(fish, fish.dist)
cpd<-function(fixed, moving, plot = TRUE, beta = 3.0, lambda = 3.0, gamma = 0.7, sigma = 0.0, max.iter = 150){
  
  output<-.Call("CoherentPointDriftRegistration", fixed[,1], fixed[,2], moving[,1], moving[,2], beta, lambda, gamma, sigma^2, max.iter)
  output$points<-matrix(output$points, nrow = length(output$points)/2)
  return(output)
}


# for(i in 1:83){
#   regi<-cpd(fixed, moving, max.iter=i, gamma=0.2)
#   plot(fixed, col='red')
#   points(regi$points)
#   Sys.sleep(0.5)
# }