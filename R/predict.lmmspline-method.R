# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the predict function from the Stats package the predict function of the nlme package
# and functions from the lmeSplines, parallel and snow packages
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Moleculesral Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Moleculesral Public License for more details.
#
# You should have received a copy of the GNU Moleculesral Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#' Predicts fitted values of an \code{lmmspline} Object
#' 
#' Predicts the fitted values of an \code{lmmspline} object for time points of interest.
#' 
#' @import parallel
#' @import nlme
#' @import lmeSplines
#' @param object an object inheriting from class \code{lmmspline}.
#' @param timePredict an optional \code{numeric} vector. Vector of time points to predict fitted values. If \code{missing} uses design points. 
#' @param numCores alternative \code{numeric} value indicating the number of CPU cores to be used for parallelization. By default estimated automatically.
#' @param ... Additional arguments which are passed to \code{predict}.
#' @return \code{matrix} containing predicted values for the requested time points from argument \code{timePredict}. 
#' @examples
#' \dontrun{
#' data(kidneySimTimeGroup)
#' G1 <- which(kidneySimTimeGroup$group=="G1")
#' testLMMSpline<- lmmSpline(data=kidneySimTimeGroup$data[G1,],
#'                  time=kidneySimTimeGroup$time[G1],
#'                  sampleID=kidneySimTimeGroup$sampleID[G1])
#' mat.predict <- predict(testLMMSpline, timePredict=c(seq(1,4, by=0.5)))}

#' @export
predict.lmmspline<- function(object, timePredict, numCores, ...){

  if(missing(timePredict)){
    return(object@pred.spline)
  }else{
    
  models <- object@models
  t <- na.omit(object@data$time)
  if(min(timePredict)<min(t) | max(timePredict)>max(t))
    stop(cat('Can only predict values within the time range',range(t)[1],'to',range(t)[2]))
  
  pred.spline <- rep(NA,length(timePredict))
  pred.df <- data.frame(all=rep(1,length(timePredict)), time=timePredict)
  pred.df$Zt = approx.Z(object@data$Zt, object@data$time, timePredict)
  
  if(missing(numCores)){
    num.Cores <- detectCores()
  }else{
    num.Cores <- detectCores()
    if(num.Cores<numCores){
      warning(paste('The number of cores is bigger than the number of detected cores. Using the number of detected cores',num.Cores,'instead.'))
    }else{
      num.Cores <- numCores
    }
    
  }
  lme <- nlme::lme
 cl <- makeCluster(num.Cores,"SOCK")
  clusterExport(cl, list('models','pred.spline','pred.df','predict'),envir=environment())

  new.data <- parLapply(cl,1:length(models),fun = function(i){
  # library(nlme)
    cl <- class(models[[i]])
    pred.spline <- switch(cl,
                          lm=predict.lm(models[[i]], newdata=pred.df, level=1, na.action=na.exclude),
                          lme=predict(models[[i]], newdata=pred.df, level=1, na.action=na.exclude)          
    )
    return(pred.spline)
  
  })
  
  stopCluster(cl)
  pred.spl <- matrix(unlist(new.data),nrow=length(models),ncol=length(timePredict),byrow=T)
  return(pred.spl)}
}