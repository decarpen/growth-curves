library(ggplot2)
library(reshape2)

dat = read.csv("120807_4ZEV-snf1-GC.csv")

strain = scan("strain.csv", sep=",")
dose = scan("best_dose.csv", sep=",")
media = scan("media.csv", sep=",", what="char")
replicate=scan("replicate.csv",sep=",")

time = seq(15,15*4*48+15,by=15)/60

dat = dat[2:73]
dat = as.data.frame(t(dat))
dat$strain = strain
dat$dose = dose
dat$media = media
dat$rep = replicate

rm(strain)
rm(dose)
rm(media)
rm(replicate)

dat.m = melt(dat, id.vars=c("dose", "media", "strain", "rep"))
dat.m$strain = as.factor(dat.m$strain)
time = rep(time,each=72)
dat.m$time = time
dat.m$strainrep = paste(dat.m$strain, dat.m$rep)

p = ggplot(dat.m, aes(x=time, y=value,group=c(strainrep)))+geom_line(aes(colour=strain))+facet_grid(media~dose)+theme_bw()
print(p)

dat.m.sub = dat.m[-which(dat.m$strain==12424),]
p = ggplot(dat.m.sub, aes(x=time, y=value,group=c(strainrep)))+geom_line(aes(colour=strain))+facet_grid(media~dose)+theme_bw()
print(p)


