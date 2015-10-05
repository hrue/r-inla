library(shiny)
library(INLA)

logit <- function(x){return(log(x/(1-x)))}
invlogit <- function(x){return(1/(1+exp(-x)))}

summary.plot <- function(original.data,model,prev.flag,ppv.flag){
  data = analysis.data(original.data)
  I = dim(data)[1]
  # intervals
  f = pf(0.95, 2, I-2)
  c = sqrt(f)
  mu = model$summary.fixed[,1]
  mean.logitse = mu[1]
  mean.logitsp = mu[2]
  t = seq(0,2*pi,by = 2*pi/100)
  # confidence
  sd = model$summary.fixed[,2]
  sd.Elogitse = sd[1]
  sd.Elogitsp = sd[2]
  rho.Elogit = model$misc$lincomb.derived.correlation.matrix[1,2]
  mu.A = mean.logitse + sd.Elogitse*c*cos(t)
  mu.B = mean.logitsp + sd.Elogitsp*c*cos(t + acos(rho.Elogit))
  
  confidence.se = invlogit(mu.A)
  confidence.sp = invlogit(mu.B)
  
  # predict
  prec = model$summary.hyperpar[,1]
  var.logitse = 1/prec[1]
  var.logitsp = 1/prec[2]
  if (prev.flag==TRUE){
    rho.logit = prec[4]
  }else{
    rho.logit = prec[3]
  }
  
  varM.mulogit = model$misc$lincomb.derived.covariance.matrix
  sd.predict.logitse = sqrt(var.logitse + varM.mulogit[1,1])
  sd.predict.logitsp = sqrt(var.logitsp + varM.mulogit[2,2])
  cov.predict = rho.logit*sqrt(var.logitse*var.logitsp) + varM.mulogit[1,2]
  rho.predict = cov.predict/(sd.predict.logitse*sd.predict.logitsp)
  mu.A = mean.logitse + sd.predict.logitse*c*cos(t)
  mu.B = mean.logitsp + sd.predict.logitsp*c*cos(t + acos(rho.predict))
  
  predict.se = invlogit(mu.A)
  predict.sp = invlogit(mu.B)
  
  # plot
  if (ppv.flag==TRUE){
  tt = seq(0,1,by=0.0001)
  logit.tt = logit(tt)
  logit.ww = mu[1] + rho.logit*sd.predict.logitse/sd.predict.logitsp*(logit.tt-mu[2])
  ww = invlogit(logit.ww)
  symbols(data$PPV,data$NPV,circles=data$N,inches=0.35,fg=rgb(0,0,100,50,maxColorValue=255),bg=rgb(200,250,30,80,maxColorValue=255))
  lines(confidence.sp,confidence.se,col="red", lwd=3,lty=3)
  lines(predict.sp,predict.se,col="green",lty=2, lwd=3)
  lines(tt,ww,col="blue",lty=1, lwd=3)
  }else{
    tt = seq(0,1,by=0.0001)
    logit.tt = logit(tt)
    logit.ww = mu[1] + rho.logit*sd.predict.logitse/sd.predict.logitsp*(logit.tt-mu[2])
    ww = invlogit(logit.ww)
    symbols(1-data$Sp,data$Se,circles=data$N,inches=0.35,fg=rgb(0,0,100,50,maxColorValue=255),bg=rgb(200,250,30,80,maxColorValue=255))
    lines(1-confidence.sp,confidence.se,col="red", lwd=3,lty=3)
    lines(1-predict.sp,predict.se,col="green",lty=2, lwd=3)
    lines(1-tt,ww,col="blue",lty=1, lwd=3)
  }
}

analysis.data <- function(data){
  I = dim(data)[1]
  TP = data$TP
  TN = data$TN
  FP = data$FP
  FN = data$FN
  n1 = TP+FN
  n0 = FP+TN
  Se = TP/n1
  Sp = TN/n0
  N = n0+n1
  m1 = TP+FP
  m0 = TN+FN
  PPV = TP/m1
  NPV = TN/m0
  
  analydata = data.frame(id=(1:I),N=N,TP=TP,TN=TN,FP=FP,FN=FN,
                          n1=n1,n0=n0,m1=m1,m0=m0,
                          Se=Se,Sp=Sp,PPV=PPV,NPV=NPV)
  
}

 inlabiv <- function(original.data,ppv.flag){
   data = analysis.data(original.data)
   I = dim(data)[1]

   if (ppv.flag==FALSE){
     Ntrials = matrix(rbind(data$n1,data$n0),2*I,1)
   } else{
     Ntrials = matrix(rbind(data$m1,data$m0),2*I,1)
   }
   Y = matrix(rbind(data$TP,data$TN),2*I,1)
   x1 = matrix(rbind(rep(1,I),rep(0,I)),2*I,1)
   x2 = matrix(rbind(rep(0,I),rep(1,I)),2*I,1)
   
   # include a linear combination to get (tp + tn) or (ppv + npv)
   lc1 = inla.make.lincomb(x1=1)
   names(lc1) = "lc1"
   lc2 = inla.make.lincomb(x2=1)
   names(lc2) = "lc2"
   lc = c(lc1, lc2)
   
   run.data = data.frame(Ntrials=c(Ntrials), Y=c(Y), id=c(1:(2*I)), x1=c(x1), x2=c(x2))
   
   formula <- Y~f(id, model="2diid", param=c(0.25, 0.025, 0.25, 0.025, 0, 0.2), n=2*I)  + x1 + x2 -1
   
   model = inla(formula, family="binomial", data=run.data, Ntrials=Ntrials, verbose=T,
                lincomb = lc, control.inla=list(strategy="gaussian", lincomb.derived.correlation.matrix = TRUE))
   
   return(model)
 }

inlativ <- function(original.data,ppv.flag){
  data = analysis.data(original.data)
  I = dim(data)[1]
  if (ppv.flag==FALSE){
    Ntrials = c(data$n1,data$n0,data$N)
    Y = c(data$TP,data$TN,data$n1)
  } else{
    Ntrials = c(data$m1,data$m0,data$N)
    Y = c(data$TP,data$TN,data$m1)
  }
  x1 = c(rep(1,I),rep(NA,2*I))
  x2 = c(rep(NA,I),rep(1,I),rep(NA,I))
  x3 = c(rep(NA,2*I),rep(1,I))
  
  lc1 = inla.make.lincomb(x1=1)
  names(lc1) = "lc1"
  lc2 = inla.make.lincomb(x2=1)
  names(lc2) = "lc2"
  lc3 = inla.make.lincomb(x3=1)
  names(lc3) = "lc3"
  lc = c(lc1, lc2, lc3)
  
  run.data = data.frame(Ntrials=c(Ntrials), Y=c(Y), id=c(1:(3*I)), x1=c(x1), x2=c(x2), x3=c(x3))
  formula <- Y~f(id, model="iid3d", n=3*I)  + x1 + x2 + x3 -1
  model = inla(formula, family="binomial", data=run.data, Ntrials=Ntrials, verbose=T,
               lincomb = lc, control.inla=list(strategy="gaussian", lincomb.derived.correlation.matrix = TRUE))
}


################# Server code #######################
shinyServer(function(input, output) {
  
  
  # Return the requested dataset
  datasetInput <- reactive({
    in.file = input$iFile
    if (!is.null(in.file)){
      read.table(in.file)
    } else{
      switch(input$dataset,
             "Telomerase" = read.table("data/Telomerase.txt", header=TRUE),
             "ScheidlerLag" = read.table("data/scheidler_LAG.txt", header=TRUE),
             "NMP22" = read.table("data/NMP22.txt", header=TRUE))
  }
  })
  
  
  
  # output data review
  output$view <- renderTable({
    data <- datasetInput()
    head(data,dim(data)[1])
  })
  
  runmodel <- reactive({
    data <- datasetInput()
    if (input$prev==TRUE) {
      model <- inlativ(data,input$ppvnpv)
    } else {
      model <- inlabiv(data,input$ppvnpv)
    }
    })
  # output model summary
  output$summary <- renderPrint({
    summary(runmodel())
  })
  # output summary plot
  output$sumplot <- renderPlot({
    model = runmodel()
    data = datasetInput()
    summary.plot(data,model,input$prev,input$ppvnpv)
    #summary.plot(data,model,FALSE,FALSE)
  })
  # output plot
  output$plot <- renderPlot({
    model <- runmodel()
    if (input$prev){
      par(mfrow=c(3,3))
      plot(inla.smarginal(model$marginals.hyperpar[[1]]),type="l",main="Precision for logit(se)",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.hyperpar[[2]]),type="l",main="Precision for logit(sp)",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.hyperpar[[3]]),type="l",main="Precision for logit(pi)",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.hyperpar[[4]]),type="l",main="cor(logit(se), logit(sp))",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.hyperpar[[5]]),type="l",main="cor(logit(se), logit(pi))",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.hyperpar[[6]]),type="l",main="cor(logit(sp), logit(pi))",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.fixed[[1]]),type="l",main="mean for logit(se)",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.fixed[[2]]),type="l",main="mean for logit(sp)",xlab="",ylab="") 
      plot(inla.smarginal(model$marginals.fixed[[3]]),type="l",main="mean for logit(pi)",xlab="",ylab="") 
    } else{
      par(mfrow=c(2,3))
      plot(inla.smarginal(model$marginals.hyperpar[[1]]),type="l",main="Precision for logit(se)",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.hyperpar[[2]]),type="l",main="Precision for logit(sp)",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.hyperpar[[3]]),type="l",main="cor(logit(se), logit(sp))",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.fixed[[1]]),type="l",main="mean for logit(se)",xlab="",ylab="")
      plot(inla.smarginal(model$marginals.fixed[[2]]),type="l",main="mean for logit(sp)",xlab="",ylab="") 
    }
    })
  
})
