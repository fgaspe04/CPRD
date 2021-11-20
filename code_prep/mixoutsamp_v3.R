
#ruth.keogh@lshtm.ac.uk
#mixoutsamp is a function which gives out-of sampling predictions from a mixed model fitted using lme. 
#'@param model: this is a mixed model fitted using lme
#'@param newdata: this is a new data set for which it is desired to make predictions. The columns in the data should have the same names as those used to fit the model. The variables should also be of the same type as in the data used to fit the mixed model (numeric, factor etc)

mixoutsamp=function(model,newdata){
  
  n=dim(newdata)[1]
  
  #----error message
  if(class(model)!="lme"){stop("Error: model is not of type lme")}
  #----end of error message
  
  #----error message
  if(is.null(newdata)==T){stop("Error: no data was provided in newdata")}
  #----end of error message
  
  #name of the id (group) variable
  id.name=names(model$groups)
  
  #name of the response variable
  y.name=as.character(attr(model$terms,"variables")[[2]])
  
  #name of variable which distinguished between different response variables (for multivariate situation)
  if(is.null(model$modelStruct$varStruct)==0){
    ytype.name=as.character(formula(model$modelStruct$varStruct)[[2]][[3]])
  }else{ytype.name=NULL}
  
  re.names = attr(model$modelStruct$reStr[[1]],"Dimnames")[[1]]
  
  #type of within-person correlation structure
  if(is.null(attributes(model$modelStruct$corStruct)$class[1])==1){
    corr.struct.type="indep"
  } else {
    corr.struct.type=attributes(model$modelStruct$corStruct)$class[1]
  }
  
  #----error message
  
  #names of all variables used in any part of the model
  ## Note this is a new way of calculating this. It hasn't been fully tested yet. Let me know if it breaks anything. David
  data.names = unique(c(all.vars(model$terms),all.vars(attr(model$modelStruct$reStruct$patid,"formula")),all.vars(attr(model$modelStruct$varStruct,"formula"))))
  
  
  #   if(is.null(ytype.name)==1){
  #     data.names=c(id.name,
  #     as.character(attributes(model$terms)$variables)[-1],#response variable and variables with fixed effects
  #     attributes(getVarCov(model))$dimnames[[1]])#variables with random effects
  #   }else{
  #       data.names=c(id.name,
  #                    as.character(attributes(model$terms)$variables)[-1],#response variable and variables with fixed effects
  #                    attributes(getVarCov(model))$dimnames[[1]]#variables with random effects
  #                    )}
  # 
  #   data.names=data.names[data.names!="(Intercept)"]
  #   data.names=unique(data.names)
  #   data.names=data.names[which(sapply(1:length(data.names),function(v)grepl(":",data.names[v]))==T)]
  #   # MS: code above removes all variable names when there are no interactions # DS: Maybe but the code bellow creates a bug for Example 5.4
  #   #data.names<-data.names[!sapply(1:length(data.names),function(v)grepl(":",data.names[v]))]
  #   #this omits interactions between variables already listed (this is necessary to avoid an error occurring in the error message where I check whether newdata contains all the variables in the model)
  
  #check that newdata contains all variable names in data.names
  if(sum(names(newdata)%in%data.names)!=length(data.names)){stop("Error: newdata does not contain the variables used in the model")}
  #----end of error message
  
  #----error message
  for(v in 1:length(data.names)){
    if(sum(is.na(newdata[,v]))==n){
      stop(paste0("Error: the variable ''", data.names[v], "'' has no non-missing values"))
    }
  }
  #----end of error message
  
  #-------------------------------------------------------------
  #model parameters
  
  coef.fixed=matrix(model$coef$fixed,nrow=length(model$coef$fixed),ncol=1) #vector of fixed coefficients
  G=as.matrix(getVarCov(model))  #random effect variance-covariance matrix
  #residual variance: this is a single value if the variance is homoscedastic but a vector of values in cases of heteroscedasticity (e.g. when we have a multivariable outcome)
  
  if(is.null(attributes(model$modelStruct$varStruct)$weights)==1){
    resid.var=(model$sigma)^2  
  }else if (is.null(attributes(model$modelStruct$varStruct)$weights)==0){
    #resid.var= matrix((model$sigma/unique(attributes(model$modelStruct$varStruct)$weights))^2,nrow=1,byrow=TRUE)
    #resid.var= matrix((model$sigma*coef(varFunc(model$modelStruct$varStruct),unconstrained = F,allCoef = T)[
    #                       as.numeric(attributes(model$modelStruct$varStruct)$groupNames)])^2,nrow=1,byrow=TRUE)
    resid.var= (model$sigma*coef(varFunc(model$modelStruct$varStruct),unconstrained = F,allCoef = T))^2
  }
  
  ## filter the dataset removing rows where there is missing data in needed variables other than the response
  observed = !logical(dim(newdata)[1]) # logical defaults to false so !
  for(v in 1:length(data.names)){
    if (data.names[v]==y.name)
      next
    observed = observed & (!is.na(newdata[,data.names[v]]))
  }
  
  #----error message
  if (sum(observed)==0)
  {
    stop(paste0("Error: there are no data rows with all needed variables observed"))
  }
  #----end of error message
  
  # remove rows that don't have all variables 
  # TODO: What do we want the behaviour to be in this case? Do we want to error out if not complete cases?
  newdata = newdata[observed,]
  
  #-------------------------------------------------------------
  #the observed data (rows of newdata where the response y is observed)
  
  ids <- unique(newdata[,id.name])  #vector of unique id numbers (aka the grouping variable)
  is.y.obs=!is.na(newdata[,y.name]) #indicator of whether y is observed (0: no, 1: yes)
  
  x.mat=model.matrix(model,data=newdata[is.y.obs,])  #fixed effects design matrix
  z.mat=model.matrix(formula(model$modelStruct$reStr)[[1]],data=newdata[is.y.obs,]) #random effects design matrix
  y = newdata[,y.name][is.y.obs] #the outcome
  xb = x.mat%*%coef.fixed   #fitted values for the fixed part of the model
  if (is.null(model$modelStruct$varStruct)==0){
    ytype=newdata[,ytype.name][is.y.obs] #response variable type (for multivariate situation)
  }
  
  obs.ids=newdata[,id.name][is.y.obs] #id numbers
  
  #vector of times used in the correlation structure 
  if(corr.struct.type!="indep"){
    vec.corr.times=model.matrix(as.formula(paste("~",formula(model$modelStruct$corStruct)[[2]][[2]])),data=newdata[is.y.obs,])[,-1]
  }
  
  #-------------------------------------------------------------
  #within-person correlation parameters
  
  num.raneff=dim(z.mat)[2]   #number of random effects
  #parameters associated with the within-person correlation
  if(corr.struct.type=="corExp"){
    range=coef(model$modelStruct$corStruct, unconstrained = F)
  } else if(corr.struct.type=="corLin"){
    range=coef(model$modelStruct$corStruct, unconstrained = F)
  }else if(corr.struct.type=="corGaus"){
    range=coef(model$modelStruct$corStruct, unconstrained = F)
  }else if(corr.struct.type=="corRatio"){
    range=coef(model$modelStruct$corStruct, unconstrained = F)
  }else if(corr.struct.type=="corSpher"){
    range=coef(model$modelStruct$corStruct, unconstrained = F)
  }else if(corr.struct.type=="corCompSymm"){
    Rho=coef(model$modelStruct$corStruct, unconstrained = F)
  } else if(corr.struct.type=="corAR1"|corr.struct.type=="corCAR1"){
    Phi=coef(model$modelStruct$corStruct, unconstrained = F)
  } 
  
  
  
  #-------------------------------------------------------------
  #the data for which predictions are required - this is for all individuals in newdata
  
  xstar.mat=model.matrix(model, model.frame(~ ., newdata, na.action=na.pass)) #fixed effects design matrix
  zstar.mat=model.matrix(formula(model$modelStruct$reStr)[[1]],data=newdata) #random effects design matrix
  xbstar<-xstar.mat%*%coef.fixed   #fitted values for the fixed part of the model
  
  pred.ids=newdata[,id.name]  #id numbers 
  ids=unique(pred.ids) #unique id numbers
  ## MS: Why do we need to order the IDs? ## DS: Not sure. Don't think we need to here apart from making it easier to debug.
  # ids=ids[order(ids)]
  #-------------------------------------------------------------
  #loop over individuals to give random effect part of the fitted values and individual random effects
  
  reffects=matrix(nrow=length(pred.ids),ncol=1)
  reffects.individual=matrix(nrow=length(unique(obs.ids)),ncol=num.raneff)
  
  vec.zeros=as.vector(rep(0,num.raneff))
  
  # order.pred.ids=order(pred.ids)
  #sorted.pred.ids=pred.ids[order.pred.ids]
  
  # order.obs.ids=order(obs.ids)
  
  # pred.index = 1
  # obs.index = 1
  
  split.obs.ids = split(1:length(obs.ids),newdata[is.y.obs,id.name])
  split.pred.ids = split(1:length(pred.ids),newdata[,id.name])
  
  pb = txtProgressBar(0,length(split.obs.ids), style=3)
  for(i in 1:length(split.obs.ids)){
    setTxtProgressBar(pb,i)
    
    
    # ## The following loop enables us to go through the data only once to find the observations for each person.
    # ## On a sorted list of ids it iterates until the id changes and records how many there are for that id.
    # ## It then indexes into the unsorted data. This is far faster than testing whether each id==the current id each loop.
    # count=0
    # while ((pred.index+count<=length(order.pred.ids)) & ids[i]==pred.ids[order.pred.ids[pred.index+count]]) {
    #   count = count + 1
    # }
    # count = count - 1
    # 
    # pids.list = order.pred.ids[pred.index:(pred.index+count)]
    # #pids = invwhich(pids.list,length(pred.ids))
    # pred.index = pred.index+count+1
    # 
    # count=0
    # while ((obs.index+count<=length(order.obs.ids)) & ids[i]==obs.ids[order.obs.ids[obs.index+count]]) {
    #   count = count + 1
    # }
    # count = count - 1
    # 
    # if (count==-1)
    #   next
    # oids.list = order.obs.ids[obs.index:(obs.index+count)]
    # #oids = invwhich(oids.list,length(obs.ids))
    # obs.index = obs.index+count+1
    # 
    # stopifnot(oids.list==split.obs.ids[[i]])
    # stopifnot(pids.list==split.pred.ids[[i]])
    
    oids.list=split.obs.ids[[i]]
    pids.list=split.pred.ids[[names(split.obs.ids)[i]]]
    
    #if (names(split.obs.ids)[i]==24805519)
    #  browser()
    ## It is really easy to mess this up so enable these to test each version.
    #stopifnot(oids==obs.ids%in%ids[i]) 
    #stopifnot(pids==pred.ids%in%ids[i])
    
    if(corr.struct.type!="indep"){
      t.rep=matrix(rep(vec.corr.times[oids.list],length(oids.list)),nrow=length(oids.list),ncol=length(oids.list))
      distances=abs(t.rep-t(t.rep))  
      distances.2=t.rep-t(t.rep)
    }
    if(corr.struct.type=="indep"){
      corr.mat=diag(1,length(oids.list)) 
    }else if(corr.struct.type=="corExp"){
      corr.mat=exp(-distances/range)
    }else if(corr.struct.type=="corLin"){
      corr.mat=ifelse(distances<range,1-distances/range,0)
    }else if(corr.struct.type=="corGaus"){
      corr.mat=exp(-(distances/range)^2)
    }else if(corr.struct.type=="corRatio"){
      corr.mat=1/(1+(distances/range)^2)
    }else if(corr.struct.type=="corSpher"){
      corr.mat=ifelse(distances<range,1-1.5*(distances/range)+0.5*(distances/range)^3,0)
    } else if(corr.struct.type=="corCompSymm"){
      corr.mat=matrix(Rho,nrow=nrow(t.rep),ncol=ncol(t.rep))
    } else if(corr.struct.type=="corAR1"|corr.struct.type=="corCAR1"){
      corr.mat=Phi^distances
    } else if(corr.struct.type=="corARMA" & p==1 & q==0){
      corr.mat=Phi^distances
    } 
    
    if(is.null(model$modelStruct$varStruct)==1){
      # Is this next line overcomplicated? Is dim(corr.mat)[1] always 1 in this case?
      resid.var.mat=matrix(resid.var,nrow=dim(corr.mat)[1],ncol=dim(corr.mat)[1])
    } else if(is.null(model$modelStruct$varStruct)==0){
      resid.var.mat=as.matrix(diag(resid.var[as.character(ytype[oids.list])],nrow=dim(corr.mat)[1],ncol=dim(corr.mat)[1]))
    }
    z.mat.oids=z.mat[oids.list,,drop=FALSE]
    t.z.mat.oids=t(z.mat.oids)
    
    reffects.individual[i,] = (G%*%t.z.mat.oids)%*%
      solve((z.mat.oids%*%G%*%t.z.mat.oids+resid.var.mat*corr.mat))%*%
      (y[oids.list]-xb[oids.list])
    reffects[pids.list]=zstar.mat[pids.list,]%*%t(reffects.individual[i,,drop=FALSE])
  }
  
  #-------------------------------------------------------------
  #calculate fitted values
  
  fitted=xbstar+reffects
  
  #--------------
  #things to output
  
  preddata=cbind(newdata,xbstar,reffects,fitted)
  names(preddata)=c(names(newdata),"fixed","random","fitted")
  random=data.frame(cbind(unique(obs.ids),reffects.individual))
  #names(random)=c(id.name,  paste0("reff",1:num.raneff))
  names(random)=c(id.name,re.names)
  
  list(preddata=preddata,random=random)
  
}


## Not actually used any more and could be removed at some point.
invwhich<-function(indices, outlength, useNames = TRUE)
{
  rv<-logical(outlength)
  if(length(indices) > 0)
  {
    rv[indices]<-TRUE
    if(useNames) names(rv)[indices]<-names(indices)
  }
  return(rv)
}
