require(data.table)
require(randomForestSRC)
require(survival)
require(rms)
require(riskRegression)
require(pec)

library(survAUC)
library(riskRegression)


data(pbc,package="survival")
setDT(pbc)
pbc <- na.omit(pbc[,.(time,status,edema,age,bili,protime,albumin)])
pbc[,time:=time/365.25]
pbc[,log.bili:=log(bili)]
pbc[,log.protime:=log(protime)]
pbc[,log.albumin:=log(albumin)]
data(follic) 
setDT(follic)
data(cost,package="pec")
setDT(cost)
cost[,time:=time/365.25]
data(GBSG2,package="pec")
setDT(GBSG2)
GBSG2[,time:=time/365.25]
##
## Table 1
## 
ccc <- function(formula,data,horizon,tau){
  #
  # 0. fit Cox model and evaluate risk prediction at prediction time horizon
  #
  fit <- coxph(formula,data=data,x=TRUE,y=TRUE)
  S <- fit$y
  prediction <- predictRisk(fit,newdata=data,times=horizon)
  #
  # 1. IPCW estimate of AUC(t)
  # 
  AUC.IPCW <- Score(list(cox=prediction),
                    formula=formula,
                    data=data,
                    times=horizon,
                    metrics="AUC",
                    split.method="none",
                    conf.int=FALSE,
                    cens.model="km")$AUC$score[model=="cox",AUC]
  # 
  # 2. IPCW estimate of C(t,tau,infty)
  C.IPCW = pec::cindex(list(cbind(1-prediction,1-prediction)),
                        formula,
                        data=data,
                        eval.times=c(horizon,tau))$AppC[[1]]
  # 
  # 3. Compute Harrell's C-index with outcome stopped at time horizon 
  #
  harrellC.t = rcorr.cens(x=1-prediction,S=stopTime(S,horizon))["C Index"]
  # 
  # 4. Compute C-index from coxph
  #
  C.coxph = summary(fit)$concordance["C"]
  # 
  # 5. Collect results output
  #
  output = data.table(matrix(round(100*c(AUC.IPCW,C.IPCW,harrellC.t,C.coxph),1),nrow=1))
  output = cbind("Data set"=as.character(substitute(data)),output,round(tau,1))
  setnames(output,c("Data set",
                    "\\widehat{AUC}_{\\mathrm{IPCW}}(t)",
                    "\\widehat{C}_{\\mathrm{IPCW}}(t)",
                    "\\widehat{C}_{\\mathrm{IPCW}}(\\tau)",
                    "\\widehat{C}(t)",
                    "\\widehat{C}(\\tau)",
                    "\\tau (yrs)"))
  output
}
table1 <- rbindlist(
  list(ccc(formula=Surv(time,status!=0)~edema+age+log.bili+log.protime+log.albumin,
           data=pbc,
           horizon=5,
           tau=max(pbc$time)),
       ccc(Surv(time,status!=0)~age+hgb+clinstg+ch,
           data=follic,
           horizon=5,
           tau=max(follic$time)),
       ccc(Surv(time, cens!=0) ~ horTh+age+menostat+tsize+tgrade+pnodes+progrec+estrec,
           data=GBSG2,
           horizon=5,
           tau=max(GBSG2$time)),
       ccc(Surv(time, status!=0) ~ age+sex+hypTen+ihd+prevStroke+cholest+atrialFib+strokeScore,
           data=cost,
           horizon=5,
           tau=max(cost$time))))
table1













TR = ovarian[1:16,]
TE = ovarian[17:26,]

summary(TR$futime)
table(TR$fustat)

train.fit = coxph(Surv(futime, fustat) ~ age,
                   x=TRUE, y=TRUE, method="breslow", data=TR)


train.fit$y

#fit <- coxph(formula,data=data,x=TRUE,y=TRUE)
S <- with(TE, Surv(futime, fustat))
prediction = predictRisk(train.fit,newdata=TE,times=700)
#
# 1. IPCW estimate of AUC(t)
# 
AUC.IPCW <- Score(list(cox=prediction),
                  formula=Surv(futime, fustat) ~ age,
                  data=TE,
                  times=c(400,500,700,800),
                  metrics="AUC",
                  split.method="none",
                  conf.int=FALSE,
                  cens.model="km")$AUC$score[model=="cox",AUC]
# 
# 2. IPCW estimate of C(t,tau,infty)
C.IPCW <- pec::cindex(list(cbind(1-prediction,1-prediction)),
                      formula,
                      data=TE,
                      eval.times=c(400,500,700,800))$AppC[[1]]
# 
# 3. Compute Harrell's C-index with outcome stopped at time horizon 
#
harrellC.t <- rcorr.cens(x=1-prediction,S=stopTime(S,700))["C Index"]
# 
# 4. Compute C-index from coxph
#
C.coxph <- summary(fit)$concordance["C"]

#cindex PEC
cindex_pec = pec::cindex( train.fit,
            formula = Surv(futime, fustat) ~ age,
            data=TE,
            eval.times=700)

#cindex Harrell
p1 = predictSurvProb( train.fit, newdata = TE, times = 700 )#c(400,500,700,855))    ###Determine concordance
harrelC1 = rcorrcens( formula = Surv(futime, fustat) ~ age, data = TE )

lpnew    = predict(train.fit, newdata=TE)
Surv.rsp = Surv(TR$futime, TR$fustat)
Surv.rsp.new = Surv(TE$futime, TE$fustat)

Cstat = UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat
