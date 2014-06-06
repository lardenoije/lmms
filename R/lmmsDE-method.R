# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the lm function from the Stats package the lme function of the nlme package
# and functions from the lmeSpline, reshape, parallel, snow and gdata packages
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

#' Differential expression analysis using linear mixed effect model splines.
#' 
#' Function to fit a linear mixed effect model splines to perform differential expression analysis. The \code{\link{lmmsDE}} function fits LMM models with either a \code{cubic}, \code{p-spline} or \code{cubic p-spline} basis and compares the models to the null models. The type of basis to use is specified with the \code{basis} argument.
#' 
#' @import nlme
#' @import lmeSplines
#' @import gdata
#' @import reshape
#' @import parallel
#' @import methods
#' @usage lmmsDE(data, time, sampleID, group, type,
#' experiment, basis, knots, numCores)
#' @param data \code{data.frame} or \code{matrix} containing the samples as rows and features as columns
#' @param time \code{numeric} vector containing the sample time point information.
#' @param sampleID \code{character}, \code{numeric} or \code{factor} vector containing information about the unique identity of each sample
#' @param group \code{character}, \code{numeric} or \code{factor} vector containing information about the group (or class) of each sample
#' @param type \code{character} indicating what type of analysis is to be performed. Options are \code{"time"} to identify differential expression over time, \code{"group"} to identify profiles with different baseline levels (intercepts), and \code{"time*group"} an interaction between these two . Use \code{"all"} to calculate all three types.
#' @param experiment \code{character} describing the experiment performed. \code{"timecourse"} for replicated experiments with less variation in individual expression values (e.g. model organism, cell culture), \code{"longitudinal1"} for different intercepts and \code{"longitudinal2"} for different intercepts and slopes. 
#' @param basis \code{character} string. What type of basis to use, matching one of \code{"cubic"} smoothing spline as defined by Verbyla \emph{et al.} 1999, \code{"p-spline"} Durb\'an \emph{et al.} 2005 or a \code{"cubic p-spline"}.
#' @param knots can take an integer value corresponding to the number of knots for the chosen basis or by default calculated as  in Ruppert 2002. Not in use for the 'cubic' smoothing spline basis.
#' @param numCores alternative \code{numeric} value indicating the number of CPU cores to be used for parallelization. Default value is automatically estimated.
#' @details
#' lmmsDE extends the LMMS modelling framework to permit tests between groups, across time, and for interactions between the two. 
#' Suppose we have \eqn{R} different groups of individuals, with \eqn{h_i} denoting the group for each individual \eqn{i}. 
#' Further we define \eqn{h_{ir}} to be the indicator for the \eqn{r^{th}} group, that is, \eqn{h_{ir}=1} if \eqn{h_i=r} and 0 otherwise. 
#' The mean curve for each group \eqn{f_{h_i}} in the full LMMSDE \code{experiment}="timecourse" is given by:
#' \deqn{f_{h_i}(t_{ij})= \beta_0+ \beta_1 t_{ij} + \sum\limits_{k=1}^{K}u_k(t_{ij}-\kappa_k)_+ }
#' \deqn{+\sum\limits_{r=2}^{R}h_{ir}(\alpha_{0r} + \alpha_{1r}t_{ij}) + \sum\limits_{r=2}^{R} h_{ir}{ \sum\limits_{k=1}^{K}v_{rk}(t_{ij}-\kappa_k)_+}.}
#' Here \eqn{\mathbf{\alpha_0}={\alpha_{0r}}} are the differences between the intercepts between each group and the first group; \eqn{\mathbf{\alpha_1}={\alpha_{1r}}} are the differences in slope between each group and the first group; and \eqn{v_{rk}} are the differences in spline coefficients between each group and the first group. 
#' We can use this model to test the effects of time, group and interactions between the two. 
#' For a single group (\code{type}="time"), all \eqn{h_{ir}=0}, and time effects will be detected if the null hypothesis \eqn{\beta_1=0} is rejected. 
#' To detect differences between groups (\code{type}="group"), we set \eqn{\mathbf{\alpha_1}=0}, and test the null hypotheses \eqn{\mathbf{\alpha_0=0}}. 
#' Finally, including all parameters allows us to test for time * group interactions (\code{type}="time*group"), by permitting different slopes in different groups and different intercepts. 
#' For \code{experiment}="longitudinal1" we include subject-specific random effects and assume them to be parallel to the mean curve. Finally, for \code{experiment}="longitudinal2" we assume subject-specific random effects to be straight lines and assuming independence between the random intercept and slope, so the covariance matrix for the random effects \eqn{\Sigma} is diagonal.
#' In each case we compare the model fit of the expanded model with the respective null model using the function \code{\link{anova}}.
#' @return lmmsDE returns an object of class \code{lmmsde} containing the following components:
#' \item{DE}{\code{data.frame} returning p-values and adjusted p-values using Benjamini-Hochberg correction for multiple testing of the differential expression testing over time, group or their interaction.}
#' \item{modelTime}{a \code{list} of class \code{\link{lme}}, containing the models for every feature modelling the time effect.} 
#' \item{modelGroup}{a \code{list} of class \code{\link{lme}}, containing the models for every feature modelling group effect. }
#' \item{modelTimeGroup}{a \code{list} of class \code{\link{lme}}, containing the models for every feature modelling time and group interaction effect. }
#' \item{type}{an object of class \code{character}, describing the test performed either time, group, time*group or all. }
#' \item{experiment}{an object of class \code{character} describing the model used to perform differential expression analysis.}
#' @references  Durbàn, M., Harezlak, J., Wand, M. P., & Carroll, R. J. (2005). \emph{Simple fitting of subject-specific curves for longitudinal data.} Statistics in medicine, 24(8), 1153-67.
#' @references  Ruppert D. (2002). \emph{Selecting the number of knots for penalized splines.} Journal of Computational and Graphical Statistics 11, 735-757
#' @references  Verbyla, A. P. Cullis, B. R., & Kenward, M. G. (1999). \emph{The analysis of designed experiments and longitudinal data by using smoothing splines.} Appl.Statist.(1999), 18(3), 269-311.
#' @seealso \code{\link{summary.lmmsde}}, \code{\link{plot.lmmsde}}
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' lmmsDEtest <-lmmsDE(data=kidneySimTimeGroup$data,time=kidneySimTimeGroup$time,
#'               sampleID=kidneySimTimeGroup$sampleID,group=kidneySimTimeGroup$group)
#' summary(lmmsDEtest)}
#' @docType methods
#' @rdname lmmsDE-methods
#' @export
setGeneric('lmmsDE',function(data,time,sampleID,group,type,experiment,basis,knots,numCores){standardGeneric('lmmsDE')})
setClassUnion("missingOrnumeric", c("missing", "numeric"))
setClassUnion("missingOrcharacter", c("missing", "character"))
setClassUnion("factorOrcharacterOrnumeric", c( "factor","character","numeric"))
setClassUnion("matrixOrframe",c('matrix','data.frame'))
#' @rdname lmmsDE-methods
#' @aliases lmmsDE lmmsDE,matrixOrframe,numeric,factorOrcharacterOrnumeric,
#' factorOrcharacterOrnumeric,missingOrcharacter,missingOrcharacter,missingOrcharacter,
#' missingOrnumeric,missingOrnumeric-method
#' @exportMethod lmmsDE

setMethod('lmmsDE',c(data="matrixOrframe",time="numeric",sampleID="factorOrcharacterOrnumeric",group="factorOrcharacterOrnumeric",type="missingOrcharacter",experiment="missingOrcharacter",basis="missingOrcharacter",knots="missingOrnumeric",numCores="missingOrnumeric"), function(data,time,sampleID,group,type,experiment,basis,knots,numCores){
  lmmsDEPara(data=data,time=time,sampleID=sampleID,group=group,type=type,experiment=experiment,basis=basis,numCores=numCores)
})

lmmsDEPara <- function(data, sampleID, time, group, type,experiment,basis,knots,numCores){

  model.time <- list()
  model.time.group <- list()
  model.group <- list()

  if(missing(basis))
    basis<-'cubic'
  if(missing(experiment))
    experiment<-'timecourse'
  if(missing(type))
    type<-'all'
  
  if(type=="time*group")
    type <- "grouptime"
  
  experiment.collection <-  c("timecourse","longitudinal1","longitudinal2")
  if(!experiment%in% experiment.collection)
    stop(paste("Chosen experiment is not available. Choose ", paste(experiment.collection,collapse=', ')))
  
  
  basis.collection <-  c("cubic","p-spline","cubic p-spline")
  if(!basis%in% basis.collection)
    stop(paste("Chosen basis is not available. Choose ",paste(basis.collection ,collapse=", ")))
  
  
  type.collection <-  c("all","time","group","grouptime")
  if(!type%in% type.collection)
    stop(paste("Chosen type is not available. Choose ",paste(type.collection,collapse=", ")))
  
  if(diff(range(c(length(sampleID),length(time),length(group),nrow(data))))>0)
    stop("Size of the input vectors sampleID, time, group and ncol(data) are not equal")
  if(missing(knots)& (basis=="p-spline"|basis=='cubic p-spline'))
    warning("The number of knots is automatically estimated")

  data <- other.reshape.group(Rep=sampleID,Time=time,Data=data,Group=group)
  
  if(length(unique(group))==1){
    warning("Only single group detected! Performing only time effect test.")
    type <- 'time'
  }
  
  
  data$Group <- as.factor(data$Group)
  data$time = as.numeric(as.character(data$Time))
  data$Expr = as.numeric(as.character(data$Expr))
  data$all= rep(1,nrow(data))
  

  #### CUBIC SPLINE BASIS ####
  if(basis=="cubic"){
    data$Zt <- smspline(~ time, data=data)
    knots <- sort(unique(data$time))[-c(1,length(unique(data$time)))]
  }
  #### PENALIZED SPLINE BASIS#####
  if(basis%in%c("p-spline","cubic p-spline")){
    if(missing(knots)){
      K <- max(5,min(floor(length(unique(data$time))/4),40))
      knots <- quantile(na.omit(unique(data$time)),seq(0,1,length=K+2))[-c(1,K+2)]
    }
    
    PZ <- outer(data$time,knots,"-")
    if(basis=="cubic p-spline"){
      PZ <- PZ^3
    }
    PZ <- PZ *(PZ>0)
    data$Zt <- PZ 
    
  }
  
  nMolecules <- NULL
  nMolecules <- length(unique(data$Molecule))
  
  levels.genotype <- unique(data$Group)
  pvals <- c()
  k<-0
    
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

  clusterExport(cl, list('lme','data','try','anova','type','pvals'),envir=environment())
  models <-list()
  
  new.data <- parLapply(cl,1:nMolecules,fun = function(i){
    
 
    options(contrasts=c('contr.treatment','contr.poly'))
    tmp.data <- data[which(data$Molecule==unique(data$Molecule)[i]),]
    p1 <- NA
    p2 <- NA
    p3 <- NA
    
    ###########group effect################# 
    
    if(type=="group"|type=="all"){
      fit0 <- NULL
      fit2 <- NULL
      

      if(experiment=="timecourse"){
        fit0 <- try(lme(Expr ~1, data=tmp.data,
                        random=list(all=pdIdent(~Zt - 1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        fit2<- try(lme(Expr ~as.factor(Group), data=tmp.data,
                       random=list(Group=pdIdent(~Zt - 1)),
                       na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
      }
      
      if(experiment=='longitudinal1'){      
             fit0 <- try(lme(Expr ~1, data=tmp.data,
                                   random=list(all=pdIdent(~Zt - 1),Rep=pdIdent(~1)),
                                  na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
      
             fit2<- try(lme(Expr ~as.factor(Group), data=tmp.data,
                            random=list(Group=pdIdent(~Zt - 1), Rep=pdIdent(~1)),
                            na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        }
             
             if(experiment=="longitudinal2"){
               fit0 <- try(lme(Expr ~1, data=tmp.data,
                               random=list(all=pdIdent(~Zt - 1),Rep=pdDiag(~time)),
                               na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
               
               fit2<- try(lme(Expr ~as.factor(Group), data=tmp.data,
                              random=list(Group=pdIdent(~Zt - 1), Rep=pdSymm(~time)),
                              na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))

             }

      
      
      model.group <- fit2
      if(class(fit0) != 'try-error' & class(fit2) != 'try-error'){
        p1 <- anova(fit0,fit2)$`p-value`[2][1]
        
      }
    }
    
    ###########time effect#################
    if(type=="time"|type=="all"){
      fit0 <- NULL
      fit3 <- NULL
      if(experiment=="timecourse"){
        fit0<- try(lme(Expr ~1, data=tmp.data,
                       random=list(all=pdIdent(~ 1)),
                       na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        fit3<- try(lme(Expr ~time, data=tmp.data,
                       random=list(all=pdIdent(~Zt -1)),
                       na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
      }
      if(experiment=="longitudinal1"){
        fit0<- try(lme(Expr ~1, data=tmp.data,
                       random=list(Rep=pdIdent(~1)),
                       na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        fit3<- try(lme(Expr ~time, data=tmp.data,
                     random=list(all=pdIdent(~Zt - 1),Rep=pdIdent(~1)),
                     na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
      }
      if(experiment=="longitudinal2"){
        fit0<- try(lme(Expr ~1, data=tmp.data,
                       random=list(all=pdIdent(~1),Rep=pdDiag(~time)),
                       na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit3<- try(lme(Expr ~time, data=tmp.data,
                       random=list(all=pdIdent(~Zt - 1),Rep=pdDiag(~time)),
                       na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
      }
      
      model.time <- fit3
    
      if(class(fit0) != 'try-error'& class(fit3) != 'try-error'){
        p2 <- anova(fit0,fit3)$`p-value`[2][1]
         
      }
    }
    
    ###########group*time interaction#################
    if(type=="grouptime"|type=="all"){
      fit0 <- NULL
      if(experiment=="timecourse"){
        fit0 <- try(lme(Expr ~ as.factor(Group)+time,data=tmp.data,
                        random=list(all=pdIdent(~Zt - 1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit5 <- try(lme(Expr ~ as.factor(Group)*time,data=tmp.data,
                        random=list(Group=pdIdent(~Zt - 1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))

      }  
      if(experiment=="longitudinal1") {
        fit0 <- try(lme(Expr ~ as.factor(Group)+time, data=tmp.data,
                        random=list(all=pdIdent(~Zt - 1), Rep=pdIdent(~1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit5 <- try(lme(Expr ~ as.factor(Group)*time,data=tmp.data,
                        random=list(Group=pdIdent(~Zt - 1), Rep=pdIdent(~1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
      
      }
      
      if(experiment=="longitudinal2") {
        fit0 <- try(lme(Expr ~ as.factor(Group)+time, data=tmp.data,
                        random=list(all=pdIdent(~Zt - 1), Rep=pdDiag(~time)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit5 <- try(lme(Expr ~ as.factor(Group)*time,data=tmp.data,
                        random=list(Group=pdIdent(~Zt - 1), Rep=pdDiag(~time)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
      }

      model.time.group <- fit5
      if(class(fit0) != 'try-error'& class(fit5) != 'try-error'){
        p3 <- anova(fit0,fit5)$`p-value`[2][1]
        
      }
    }
    
    pvals<-c(p1,p2,p3)
    return(list(pvals=pvals,model.time=model.time,model.group=model.group,model.time.group=model.time.group))
})
 

                  
                        new.df <- matrix(sapply(new.data,'[[','pvals'),nrow=nMolecules,ncol=3,byrow=T)
                        model.time <- sapply(new.data,'[','model.time')
                        model.group <- sapply(new.data,'[','model.group')
                        model.time.group <- sapply(new.data,'[','model.time.group')

  df <- data.frame(Molecule=as.character(unique(data$Molecule)),Time=new.df[,2],adj.Time=signif(p.adjust(new.df[,2],method="BH"),2), Group=new.df[,1],adj.Group=signif(p.adjust(new.df[,1],method="BH"),2),Group_Time=new.df[,3],adj.Group_Time=signif(p.adjust(new.df[,3],method="BH"),2))

  l <- new('lmmsde',DE=df, modelTime=model.time, modelGroup=model.group, modelTimeGroup=model.time.group,knots=knots,basis=basis,type=type,experiment=experiment)
  return(l)

}

other.reshape.group <- function(Rep, Time, Data,Group){
  lme.data<-NULL
  if(sum(table(Rep,Time)>1)!=0)
    stop('Make sure you have one time point per individual')
  lme.data <- data.frame(Time=Time,Rep=Rep,Group=Group,as.matrix(Data))
  
  lme.data$Time = factor(drop.levels(lme.data$Time))
  lme.data$Rep = factor(drop.levels(lme.data$Rep))
  lme.data$Group = factor(drop.levels(lme.data$Group))
  
  melt.lme.data <-NULL
  melt.lme.data <- melt(lme.data)
  cast.lme.data  <- NULL
  cast.lme.data <- cast(melt.lme.data, variable+ Group+Rep ~ Time)
  melt.lme.data2 <- NULL
  melt.lme.data2 <-  melt(data.frame(cast.lme.data))
  
  names(melt.lme.data2) <- c("Molecule",  "Group","Rep", "Time", "Expr")
  melt.lme.data2$Time <- factor(gsub("^X", "", as.character(melt.lme.data2$Time)))
  return(as.data.frame(melt.lme.data2))
}