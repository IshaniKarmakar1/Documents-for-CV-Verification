rm(list=ls())
setwd("C:/Users/ISHANI/Desktop/STSA/sem6/Dissertation/Comparison of ICC for 1st and 2nd waves")
library(lme4)

#for phase 1
data.p1 = read.csv("phase 1 data.csv")
View(data.p1)

clus.p1 = data.p1$Cluster
y.p1 = data.p1$Affected
age.p1 = data.p1$Age
sex.p1 = data.p1$Sex
comr.p1 = data.p1$Comorbidity
expo.p1 = data.p1$Exposure

model.p1 = glmer(y.p1 ~ age.p1 + sex.p1 + comr.p1 + expo.p1 + (1|clus.p1), family=binomial(link="logit"), nAGQ=25)
summary(model.p1)

rho.hat.p1=2.612/(2.612 + (pi**2 / 3))

fam.p1 = aggregate(rep(1,times = length(clus.p1)), by = list(clus.p1), sum)     #family sizes
k.p1 = dim(fam.p1)[1]
p.p1 = fam.p1 $ x
a.p1 = 1 - (1/p.p1)
p.p1.bar = mean(p.p1)
c.sq.p1 = 1 - (2 * ((1 - rho.hat.p1)**2) * sum(a.p1)/k.p1) + ((1 - rho.hat.p1)**2) * (sum(a.p1)/k.p1 + ((sum(a.p1)/k.p1) ** 2) / (p.p1.bar - 1))

var.p1 = 2 * (1 - rho.hat.p1)**2 * ( (1/(p.p1.bar - 1)) + c.sq.p1 - (2 * (1 - rho.hat.p1) * sum(a.p1) / k.p1 / (p.p1.bar - 1)) )

#for phase 2
data.p2 = read.csv("phase 2 data.csv")
View(data.p2)

clus.p2 = data.p2$Cluster
y.p2 = data.p2$Affected
age.p2 = data.p2$Age
sex.p2 = data.p2$Sex
comr.p2 = data.p2$Comorbidity
expo.p2 = data.p2$Exposure

model.p2 = glmer(y.p2 ~ age.p2 + sex.p2 + comr.p2 + expo.p2 + (1|clus.p2), family=binomial(link="logit"), nAGQ=25)
summary(model.p2)

rho.hat.p2=3.445/(3.445 + (pi**2 / 3))

fam.p2 = aggregate(rep(1,times = length(clus.p2)), by = list(clus.p2), sum)     #family sizes
k.p2 = dim(fam.p2)[1]
p.p2 = fam.p2 $ x
a.p2 = 1 - (1/p.p2)
p.p2.bar = mean(p.p2)
c.sq.p2 = 1 - (2 * ((1 - rho.hat.p2)**2) * sum(a.p2)/k.p2) + ((1 - rho.hat.p2)**2) * (sum(a.p2)/k.p2 + ((sum(a.p2)/k.p2) ** 2) / (p.p2.bar - 1))

var.p2 = 2 * (1 - rho.hat.p2)**2 * ( (1/(p.p2.bar - 1)) + c.sq.p2 - (2 * (1 - rho.hat.p2) * sum(a.p2) / k.p2 / (p.p2.bar - 1)) )

#testing Ho: rho1 = rho2 against H1: rho1 != rho2
Z = (rho.hat.p1 - rho.hat.p2) / sqrt((var.p1/k.p1) + (var.p2/k.p2))

aplha = c(0.01,0.05,0.10)
tau.aplha = qnorm(aplha/2, lower.tail = FALSE)

for(i in 1:length(tau.aplha))
{
  print(100*aplha[i])
  if(abs(Z) > tau.aplha[i])
    print("Reject Ho at the above level of significance") else
      print("Not enough evidence to reject Ho at the above level of significance")
}

