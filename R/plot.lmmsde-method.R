# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the graphics and stats package.
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

#' Plot of \code{lmmsde} objects
#' 
#' Plot of the raw data the mean and the fitted \code{lmmsde} profile.
#' 
#' @import graphics
#' @param x An object of class \code{lmmsde}.
#' @param y \code{numeric} or \code{character} value. Either the row index or the row name determining which feature should be plotted. 
#' @param type a \code{character} indicating what model to plot. Default  \code{'all'}, options: \code{'time'}, \code{'group'},\code{'time*group'}.
#' @param smooth an optional \code{logical} value.By default set to \code{FALSE}. If \code{TRUE} smooth representation of the fitted values. 
#' @param \ldots Additional arguments which are passed to \code{plot}.
#' @return plot showing raw data, mean profile and fitted profile. 
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' lmmsDEtestl1 <-lmmsDE(data=kidneySimTimeGroup$data,time=kidneySimTimeGroup$time,
#'                 sampleID=kidneySimTimeGroup$sampleID,
#'                 group=kidneySimTimeGroup$group,
#'                 experiment="longitudinal1",basis="p-spline") 
#' plot(lmmsDEtestl1,y=2,type="all")
#' plot(lmmsDEtestl1,y=2,type="time")
#' plot(lmmsDEtestl1,y=2,type="group")
#' plot(lmmsDEtestl1,y=2,type="time*group",smooth=TRUE)}  

#' @export
plot.lmmsde <- function(x, y, type, smooth, ...){
 # library(graphics)
  if(missing(type)){
    type <- c()
    if(!is.na(x@modelGroup))
      type <- c(type,'group')
    if(!is.na(x@modelTime))
      type <- c(type,'time')
    if(!is.na(x@modelTimeGroup))
      type <- c(type,'time*group')
  }
   
  if(type=="all"){
    type <- c('time','group','time*group') 
  }
  
  if(class(y)=='numeric')
    name <- x@DE$Molecule[y]
    
  if(class(y)=='character'){
    nam <- x@DE$Molecule
    if(sum(nam%in%y)>0){
      name <- y
      y <-which(nam%in%y)
    }else{
      stop(paste('Could not find feature',y,'in rownames(x@pred.spline).'))
    }
  }
  if(length(type)==3)
    par(mfrow=c(2,2))
  
  if(sum(type%in%'time')>0){
    name2 <- paste(name,'time')
    plotLmms(x@modelTime[[y]],smooth=smooth,name2,...)
  }
  if(sum(type%in%"group")>0){
   name2 <- paste(name,'group')
   plotLmmsdeFunc(x@modelGroup,index=y,smooth=smooth,name2,...)
  }
  if(sum(type%in%"time*group")>0){
    name2 <- paste(name,'time*group')
    plotLmmsdeFunc(x@modelTimeGroup,index=y,smooth=smooth,name2,...)
  }
  par(mfrow=c(1,1))
    
  }


plotLmmsdeFunc <- function(object,index,smooth,name,...){
  if(missing(smooth))
    smooth <- F
  model <- object[[index]]

  if(is.null(model))
    stop("Requested model not available")
  cl <- class(model)
  p <- NULL
  yl <- range(na.omit(model$data$Expr))
  yl[2] <- yl[2]+1

    group <- model$data$Group
    g1 <-which(group==unique(group)[1])
    g2 <- which(group==unique(group)[2])
    plot(model$data$Expr~model$data$time,xlab='Time',ylab='Intensity',ylim=yl,main=name,col=ifelse(group==group[g1[1]],"blue","red"),pch=ifelse(group==group[g1[1]],16,17),...)
    ext1 <- tapply(model$data$Exp[g1],model$data$time[g1],function(x)mean(x,na.rm=T))
    ext2 <- tapply(model$data$Exp[g2],model$data$time[g2],function(x)mean(x,na.rm=T))
    ut1 <- unique(model$data$time[g1])
    ut2 <- unique(model$data$time[g2])
    lines(ext1~ut1, col='grey',lty=2)
    lines(ext2~ut2, type='l',col='rosybrown',lty=2)
  

    f <- fitted(model,level=1)
    f1 <- f[g1]
    f2<- f[g2]
  
    if(smooth){
      s1 <- spline(x = model$data$time[g1], y = f1, n = 500, method = "natural")
      s2 <- spline(x = model$data$time[g2], y = f2, n = 500, method = "natural")
      lines(s1$y~s1$x,type='l',col="black",lwd=2)
      lines(s2$y~s2$x,type='l',col="brown",lwd=2)
    }else{
      lines(na.omit(f1)~model$data$time[intersect(which(!is.na(f)),g1)],col="black",lwd=2)
      lines(na.omit(f2)~model$data$time[intersect(which(!is.na(f)),g2)],col="brown",lwd=2)
    }
 
    legend('topleft',legend=paste(group[c(g1[1],g2[1])],rep(c("Raw","Fitted","Mean"),each=2)),col=c('blue','red','black','brown','lightgrey','rosybrown'),lty=c(NA,NA,1,1,2,2),pch=c(16,17,NA,NA,NA,NA),ncol=3,cex=0.8,bty='n')
}
