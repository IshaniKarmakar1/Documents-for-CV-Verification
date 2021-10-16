rm(list=ls())
setwd("C:/Users/ISHANI/Desktop/STSA/sem6/Dissertation/Comparison of ICC for 1st and 2nd waves")

#WAVE 1

data.p1 = read.csv("phase 1 data.csv")
View(data.p1)

clus.p1 = unique(data.p1$Cluster)

fam.size.p1 = aggregate(rep(1,times = length(data.p1$Cluster)), by = list(data.p1$Cluster), sum)     #family sizes

library(Compositional)
Q.func = function(n)
{
  first.row = rep((1/sqrt(n)), times = n)
  A = rbind(first.row, helm(n))
  return(A)
}

data.new.x = list()

for(i in 1:(dim(fam.size.p1)[1]))
  data.new.x[[i]] = data.p1$Affected[data.p1$Cluster == clus.p1[i]]

data.new.u = list()
for(i in 1:(dim(fam.size.p1)[1]))
  data.new.u[[i]] = Q.func(length(data.new.x[[i]])) %*% data.new.x[[i]]

px = fam.size.p1$x
ax = 1 - (1/px)
kx = dim(fam.size.p1)[1]

mu.hat.x = 0
for(i in 1:kx)
  mu.hat.x = mu.hat.x + data.new.u[[i]][1]
mu.hat.x = mu.hat.x/kx

gamma.hat.sq.x = 0
for(i in 1:kx)
  for(r in 2:px[i])
    gamma.hat.sq.x = gamma.hat.sq.x + (data.new.u[[i]][r])**2
gamma.hat.sq.x = gamma.hat.sq.x / sum(px-1)

sigma.hat.sq.x = 0
for(i in 1:kx)
  sigma.hat.sq.x = sigma.hat.sq.x + (data.new.u[[i]][1] - mu.hat.x)**2
sigma.hat.sq.x = sigma.hat.sq.x/(kx-1) + gamma.hat.sq.x*sum(ax)/kx

rho.hat.x = 1 - (gamma.hat.sq.x/sigma.hat.sq.x)


#WAVE 2

data.p2 = read.csv("phase 2 data.csv")
View(data.p2)

clus.p2 = unique(data.p2$Cluster)

fam.size.p2 = aggregate(rep(1,times = length(data.p2$Cluster)), by = list(data.p2$Cluster), sum)     #family sizes

data.new.y = list()

for(i in 1:(dim(fam.size.p2)[1]))
  data.new.y[[i]] = data.p2$Affected[data.p2$Cluster == clus.p2[i]]

data.new.v = list()
for(i in 1:(dim(fam.size.p2)[1]))
  data.new.v[[i]] = Q.func(length(data.new.y[[i]])) %*% data.new.y[[i]]

py = fam.size.p2$x
ay = 1 - (1/py)
ky = dim(fam.size.p2)[1]

mu.hat.y = 0
for(i in 1:ky)
  mu.hat.y = mu.hat.y + data.new.v[[i]][1]
mu.hat.y = mu.hat.y/ky

gamma.hat.sq.y = 0
for(i in 1:ky)
  for(r in 2:py[i])
    gamma.hat.sq.y = gamma.hat.sq.y + (data.new.v[[i]][r])**2
gamma.hat.sq.y = gamma.hat.sq.y / sum(py-1)

sigma.hat.sq.y = 0
for(i in 1:ky)
  sigma.hat.sq.y = sigma.hat.sq.y + (data.new.v[[i]][1] - mu.hat.y)**2
sigma.hat.sq.y = sigma.hat.sq.y/(ky-1) + gamma.hat.sq.y*sum(ay)/ky

rho.hat.y = 1 - (gamma.hat.sq.y/sigma.hat.sq.y)


#UNDER Ho

mu.hat = (kx*mu.hat.x + ky*mu.hat.y) / (kx + ky)

gamma.hat.sq = ((gamma.hat.sq.x * sum(px - 1)) + (gamma.hat.sq.y * sum(py -1))) / (sum(px - 1) + sum(py - 1))

sg.p1 = 0
for(i in 1:kx)
  sg.p1 = sg.p1 + (data.new.u[[i]][1] - mu.hat)**2

sg.p2 = 0
for(i in 1:ky)
  sg.p2 = sg.p2 + (data.new.v[[i]][1] - mu.hat)**2

sigma.hat.sq = (sg.p1 + sg.p2)/(kx + ky - 1) + (gamma.hat.sq * (sum(ax) + sum(ay)) / (kx + ky))

rho.hat = 1 - (gamma.hat.sq / sigma.hat.sq)


#TEST STATISTIC AND DECISION

part1 = sum(log((1/px) * (1 + (px - 1)*rho.hat)))
part2 = sum((px - 1) * log(1 - rho.hat))

part3 = sum(log((1/py) * (1 + (py - 1)*rho.hat)))
part4 = sum((py - 1) * log(1 - rho.hat))

part5 = 0
for(i in kx)
  part5 = part5 + (px[i] * (data.new.u[[i]][1] - mu.hat.x)**2 / (1 + (px[i] - 1)*rho.hat))

part6 = 0
for(i in kx)
  for(r in 2:px[i])
    part6 = part6 + (data.new.u[[i]][r])**2
part6 = part6 / (1 - rho.hat)

part7 = 0
for(i in ky)
  part7 = part7 + (py[i] * (data.new.v[[i]][1] - mu.hat.y)**2 / (1 + (py[i] - 1)*rho.hat))

part8 = 0
for(i in ky)
  for(r in 2:py[i])
    part8 = part8 + (data.new.v[[i]][r])**2
part8 = part8 / (1 - rho.hat)

part9 = sum(log((1/px) * (1 + (px - 1)*rho.hat.x)))
part10 = sum((px - 1) * log(1 - rho.hat.x))

part11 = sum(log((1/py) * (1 + (py - 1)*rho.hat.y)))
part12 = sum((py - 1) * log(1 - rho.hat.y))

part13 = 0
for(i in kx)
  part13 = part13 + (px[i] * (data.new.u[[i]][1] - mu.hat.x)**2 / (1 + (px[i] - 1)*rho.hat.x))

part14 = 0
for(i in kx)
  for(r in 2:px[i])
    part14 = part14 + (data.new.u[[i]][r])**2
part14 = part14 / (1 - rho.hat.x)

part15 = 0
for(i in ky)
  part15 = part15 + (py[i] * (data.new.v[[i]][1] - mu.hat.y)**2 / (1 + (py[i] - 1)*rho.hat.y))

part16 = 0
for(i in ky)
  for(r in 2:py[i])
    part16 = part16 + (data.new.v[[i]][r])**2
part16 = part16 / (1 - rho.hat.y)

test.stat = (part1 + part2 + part3 + part4) + (part5 + part6 + part7 + part8)/sigma.hat.sq - 
  (part9 + part10 + part11 + part12) - (part13 + part14 + part15 + part16)/sigma.hat.sq

#Ho is rejected

alpha = c(0.05, 0.01, 0.10)
w1 = qchisq(alpha/2, 1, lower.tail = TRUE)
w2 = qchisq(alpha/2, 1, lower.tail = FALSE)
