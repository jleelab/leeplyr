#' Update leeplyr package
#'
#' This function updates the leeplyr packae.
#' @keywords update
#' @export
#' @examples
#' update.leeplyr
update.leeplyr<-function(){
  detach('package:leeplyr', unload=TRUE)
  remove.packages("leeplyr")
  devtools::install_github("jleelab/leeplyr")
}