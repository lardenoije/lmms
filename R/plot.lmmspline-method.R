# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed the graphics package.
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

#' Plot of \code{lmmspline} object
#' 
#' Plots the raw data, the mean and the fitted or derivative information of the \code{lmmspline} object.
#' 
#' @import graphics
#' @param x An object of class \code{lmmspline}.
#' @param y \code{character} or \code{numeric} value. Determining which feature should be plotted can be either the index or the name of the feature. 
#' @param smooth an optional \code{logical} value. Default \code{FALSE}, if \code{TRUE} smooth representation of the fitted values. 
#' @param \ldots Additional arguments which are passed to \code{plot}.
#' @return xyplot showing raw data, mean profile and fitted profile. 
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' # running for samples in group 1
#' G1 <- which(kidneySimTimeGroup$group=="G1")
#' testLmmspline <- lmmSpline(data=kidneySimTimeGroup$data[G1,],
#'                  time=kidneySimTimeGroup$time[G1],
#'                  sampleID=kidneySimTimeGroup$sampleID[G1])
#' plot(testLmmspline, y=2)
#' plot(testLmmspline, y=2, smooth=TRUE)}

#' @export
plot.lmmspline <- function(x, y, smooth, ...){
  if(length(y)>1)
    stop('Can just plot a single feature.')
  name <- ""
  if(class(y)=='numeric'){
    model <- x@models[[y]]
    name <- rownames(x@predSpline)[y]
  }
  if(class(y)=='character'){
    nam <- rownames(x@predSpline)
    if(sum(nam%in%y)>0){
      name <- y
      y <- which(nam%in%y)
      model <- x@models[[y]]
   
    }else{
      stop(paste('Could not find feature',y,'in rownames(x@pred.spline).'))
    }
  }
  if(x@derivative){
   plotLmmsDeriv(model, smooth=smooth,data=x@predSpline[y,],name,...) 
  }else{
     plotLmms(model,smooth=smooth,name,...)
  }

}


plotLmms <- function(model,smooth,name,...){
#library(graphics)

  if(missing(smooth))
    smooth <- F
  cl <- class(model)
  p <- NULL

  if(cl=="lm"){
    yl <- range(na.omit(model$model$Expr))
    yl[2] <- yl[2]+0.2
    
    plot(model$model$Expr~model$model$time,xlab='Time',ylab='Intensity',main=name,pch=16,col="blue",ylim=yl,...)
    ext <- tapply(model$model$Exp,model$model$time,function(x)mean(x,na.rm=T))
    ut <- unique(model$model$time)
     lines(ext~ut, type='l',col='grey',lty=2)
    if(smooth){
      s <- spline(x = model$model$time, y = fitted(model), n = 500, method = "natural")

    }else{
      s <- data.frame(y=fitted(model),x=model$model$time)
    }
  }
  
  if(cl=='lme'){
    yl <- range(na.omit(model$data$Expr))
    yl[2] <- yl[2]+0.2

    plot(model$data$Expr~model$data$time,xlab='Time',ylab='Intensity',pch=16,col='blue',main=name,ylim=yl,...)
    ext <- tapply(model$data$Exp,model$data$time,function(x)mean(x,na.rm=T))
    ut <- unique(model$data$time)
   
    f <- fitted(model,level=1)
    if(smooth){
      s <- spline(x = model$data$time, y = f, n = 500, method = "natural")
     
    }else{
      s <- data.frame(y=na.omit(f),x=model$data$time[!is.na(f)])
    }
  }
  lines(ext~ut, col='grey',lty=2)
  lines(s$y~s$x,col="black",lwd=2)
  legend('topleft',legend=c("Raw","Fitted","Mean"),pch=c(16,NA,NA),lty=c(NA,1,2),col=c('blue','black','grey'),cex=0.8,bty='n',ncol=3)

  
}

plotLmmsDeriv <- function(model,smooth,data,name,...){
  #library(graphics)
  
  if(missing(smooth))
    smooth <- F
  cl <- class(model)
  p <- NULL
  
  if(cl=="lm"){
    yl <- range(na.omit(unlist(c(model$model$Expr,data))))
    yl[2] <- yl[2]+0.2
    plot(model$model$Expr~model$model$time,xlab='Time',ylab='Intensity',ylim=yl,pch=16,col="blue",main=name,...)
    ext <- tapply(model$model$Exp,model$model$time,function(x)mean(x,na.rm=T))
    ut <- unique(model$model$time)
    if(smooth){
      s <- spline(x = as.numeric(colnames(data)), y = data, n = 500, method = "natural")
    }else{
      s <- data.frame(y=as.numeric(as.character(data)),x=as.numeric(colnames(data)))
    }
  }
  
  if(cl=='lme'){
    yl <- range(na.omit(unlist(c(model$data$Expr,data))))
    yl[2] <- yl[2]+0.2
    plot(model$data$Expr~model$data$time,xlab='Time',ylab='Intensity',ylim=yl,pch=16,col="blue",main=name,...)
    ext <- tapply(model$data$Exp,model$data$time,function(x)mean(x,na.rm=T))
    ut <- unique(model$data$time)

    f <- fitted(model,level=1)
    if(smooth){
      s <- spline(x = as.numeric(colnames(data)), y = data, n = 500, method = "natural")
    }else{
      s <- data.frame(y=as.numeric(as.character(data)),x=as.numeric(colnames(data)))
    }
   
  }
  
  lines(ext~ut, type='l',col='grey',lty=2)
  lines(s$y~s$x,type='l',col="black",lwd=2)
  legend('topleft',legend=c("Raw","Deriv","Mean"),pch=c(16,NA,NA),lty=c(NA,1,2),col=c('blue','black','grey'),cex=0.8,bty='n',ncol=3)
  
}
