# EMPIRICAL DATA

# poilog computes Poisson log normal distribution statistics
# the source functions fit distributions, provide random numbers, and get some stats

library(poilog)

source('R/cegsML.R')
source('R/sadrad.R')
source('R/rcegs.R')
source('R/fisher.R')
source('R/shannon.R')
source('R/simpson.R')
source('R/negbin.R')

# used by cegsML

trim<-function(n)	{
  ns <- sort(n, decreasing=TRUE)
  ns[1:2] <- sort(unique(n))[2]
	ns
}


# the data are drawn from Dryad
# DOI: 10.5061/dryad.brv15dvdc

if(!file.exists('data/Ecological_Register_data.txt.gz')) stop('data/Ecological_Register_data.txt.gz does not exist. You need to download it from https://doi.org/10.5061/dryad.brv15dvdc')

data <- read.delim(gzfile('data/Ecological_Register_data.txt.gz')) #, exdir = "data")


# prepare the data
data[is.na(data[,'count']),'count'] <- 0
data[is.na(data[,'count.2']),'count.2'] <- 0
data$co <- data$count + data$count.2


# compute basic descriptive variables
# shannon and simpson transform the values into diversity equivalents by default
# see Hill (1973) for definitions

CalcBasicStats <- function(n, removeDom=FALSE, vegan=FALSE) {
  if(removeDom & length(unique(n))>1) n <- n[n!=max(n)]
  
  if(vegan) {
    res <- c(N = sum(n), S = length(n), alpha = vegan::fisher.alpha(n), 
             H = vegan::diversity(n, index="shannon"), 
             D = 1/(1-vegan::simpson.unb(n)), 
             UniqueCounts = length(unique(n)))
  } else {
    res <- c(N = sum(n), S = length(n), alpha = fisher(n), H = shannon(n), 
             D = simpson(n), UniqueCounts = length(unique(n)))
  }
  # compute Pielou's J & the log-transformed Hill ratio of H and D
  res <- c(res, J=log(res['H'])/log(res['S']), HD = log(res['H']/res['D']))
  names(res) <- gsub("\\.H", "", names(res))
  
  res
}

CalcBasicStats(data$co[data$sample.no==data$sample.no[1]], removeDom=TRUE)
CalcBasicStats(data$co[data$sample.no==data$sample.no[1]], removeDom=TRUE, vegan=TRUE)


ExctractMetadata <- function(df) {
  meta <- c(eco=df$ecozone[1], lat=df$latitude[1], 
            long=df$longitude[1], 
            form=df$life.form[1]
  )
  meta
}

# Note: the original code skipped samples where there were <3 unique counts. 
# Here we include them but have UniqueCounts so we can remove later if needed
BasicStats.by <- tapply(data$co, list(data$sample.no), CalcBasicStats, simplify=TRUE)
BasicStats.m <- do.call(rbind, BasicStats.by)
BasicStatsNoDom.by <- tapply(data$co, list(data$sample.no), CalcBasicStats, simplify=TRUE, removeDom=TRUE)
BasicStatsNoDom.m <- do.call(rbind, BasicStatsNoDom.by)
Meta.by <- by(data, list(data$sample.no), ExctractMetadata, simplify=TRUE)
MetaData.m <- do.call(rbind, Meta.by)

BasicStats <- cbind(as.data.frame(BasicStats.m), as.data.frame(MetaData.m))
BasicStatsNoDom <- cbind(as.data.frame(BasicStatsNoDom.m), as.data.frame(MetaData.m))

write.csv(BasicStats, file='data/basic_statistics.txt')
write.csv(BasicStatsNoDom, file='data/basic_statisticsNoDom.txt')

# correlations of S and J and of S and the ratio, which are reported in the text


par(mfrow=c(1,2))
plot(log(BasicStats$S),BasicStats$J, xlab="log(S)", ylab="J")
plot(log(BasicStats$S),BasicStats$HD, xlab="log(S)", ylab="HD")

cor.test(log(BasicStats$S),BasicStats$J)
cor.test(log(BasicStats$S),BasicStats$HD)

cor.test(BasicStats$S,BasicStats$J,method='s')$estimate
cor.test(BasicStats$S,BasicStats$J,method='s')$p.value
cor.test(BasicStats$S,BasicStats$HD,method='s')$estimate
cor.test(BasicStats$S,BasicStats$HD,method='s')$p.value


# DISTRIBUTION PARAMETERS

# use poilog and provided functions to compute distribution parameters


# Calculate Parametric Richness Estimates

CalcParSR <- function(n, trim=TRUE, removeDom=FALSE) {
  if(removeDom & length(unique(n))>1) n <- n[n!=max(n)]
  if(trim) n <- trim(n)

# These are in case anything crashes
  ln <- list(p=as.numeric(NA), par=as.numeric(c(NA,NA)), logLval=as.numeric(NA))
  eg <- list(richness=as.numeric(NA), shape=as.numeric(NA))
  nb <- list(richness=as.numeric(NA), size=as.numeric(NA), 
             prob=as.numeric(NA), loglik=as.numeric(NA))

  try(ln <- poilog::poilogMLE(n), silent = TRUE)
  try(nb <- FitNB(n), silent = TRUE)
  try(eg <- cegsML(n, fitted = FALSE), silent = TRUE)
  cegs_loglik <- ifelse(is.na(eg$scale), NA, sum(dcegs(x=n, l=eg$scale,g=eg$shape, log=TRUE)))
  
  res <- c(
    pln_richness = length(n)/ln$p, pln_sigma = ln$par[2], pln_loglik = ln$logLval,
    cegs_richness = eg$richness, cegs_gamma = eg$shape, cegs_loglik = cegs_loglik, 
    nb_richness = nb$richness, nb_size=nb$size, nb_loglik = nb$loglik,
    nb_prob=nb$prob
  )
  res
}

n <- data$co[data$sample.no==1254]
CalcParSR(n, trim=TRUE)
CalcParSR(n, trim=FALSE)

system.time(
ParSR.by <- tapply(data$co, list(data$sample.no), CalcParSR, 
                   simplify=TRUE, trim=FALSE)
)
ParSR <- data.frame(do.call(rbind, ParSR.by))
write.csv(ParSRtrim,file='data/richness_estimates.txt')


thing <- unlist(lapply(ParSR.by, length))
which (thing!=10)

ParSR.by[[1200]]

ParSRNoDom.by <- tapply(data$co, list(data$sample.no), CalcParSR, 
                   simplify=TRUE, trim=FALSE, removeDom=TRUE)
ParSRNoDom <- data.frame(do.call(rbind, ParSRNoDom.by))
write.csv(ParSRtrim,file='data/richness_estimatesNoDom.txt')

ParSRtrim.by <- tapply(data$co, list(data$sample.no), CalcParSR, 
                       simplify=TRUE, trim=TRUE)
ParSRtrim <- as.data.frame(do.call(rbind, ParSRtrim.by))
write.csv(ParSRtrim,file='data/richness_estimatestrimmed.txt')


pairs(cbind(log(ParSR$pln_richness), log(ParSRtrim$pln_richness), log(ParSRNoDom$pln_richness)))
pairs(cbind(log(ParSR$cegs_richness), log(ParSRtrim$cegs_richness), log(ParSRNoDom$cegs_richness)))
pairs(cbind(log(ParSR$nb_richness), log(ParSRtrim$nb_richness), log(ParSRNoDom$nb_richness)))


pairs(cbind(ParSR$pln_loglik, ParSR$nb_loglik, ParSR$cegs_loglik))

sum(ParSR$pln_loglik<ParSR$cegs_loglik, na.rm = TRUE)
mean(ParSR$pln_loglik<ParSR$nb_loglik, na.rm = TRUE)

range(ParSR$pln_loglik-ParSR$cegs_loglik, na.rm = TRUE)
hist(ParSR$nb_loglik-ParSR$cegs_loglik)



pairs(ParSR)
plot(log(ParSR$pln_richness), log(ParSR$cegs_richness))
plot(log(ParSR$pln_richness), log(ParSR$nb_richness))
plot(log(ParSR$cegs_richness), log(ParSR$nb_richness))



# FIGURE 1: STATISTICS WITHOUT DOMINANTS

x <- read.csv('data/richness_estimates.txt')

pln_richness <- x[,1]
cegs_richness <- x[,2]

x <- read.delim('basic_statistics.txt')

alpha <- x[,3]
H <- x[,4]
D <- x[,5]
J <- x[,6]

x <- read.delim('subdominant_diversity_estimates.txt')

sub_pln_richness <- x[,1]
sub_cegs_richness <- x[,2]
sub_alpha <- x[,3]
sub_H <- x[,4]
sub_D <- x[,5]
sub_J <- x[,6]

# CEGS sample size reported at the beginning of the Results

length(sort(sub_cegs_richness))


pdf(file='Evenness_Fig_1.pdf',height=11.25,width=7.5)

par(mfrow=c(3,2))
par(mai=c(0.8,0.8,0.3,0.3))
par(cex=1)

plot(pln_richness,sub_pln_richness,pch=21,cex=0.85,col='white',bg='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Full sample PLN',ylab='PLN without dominant',xlim=c(1,1e4),ylim=c(1,1e4),log='xy',axes=F)
axis(side=1,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e4,text='A',cex=1.3,line=2.65,las=2)

plot(cegs_richness,sub_cegs_richness,pch=21,cex=0.85,col='white',bg='dodgerblue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Full sample CEGS',ylab='CEGS without dominant',xlim=c(1,1e4),ylim=c(1,1e4),log='xy',axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e4,text='B',cex=1.3,line=2.65,las=2)

plot(alpha,sub_alpha,pch=21,cex=0.85,col='white',bg='blue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste("Full sample Fisher's ",alpha)),ylab=expression(paste("Fisher's ",alpha," w/o dominant")),xlim=c(0.1,1e3),ylim=c(0.1,1e3),log='xy',axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e3,text='C',cex=1.3,line=2.65,las=2)

plot(H,sub_H,pch=21,cex=0.85,col='white',bg='orange',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('Full sample exp(',italic(H),')')),ylab=expression(paste('exp(',italic(H),') without dominant')),xlim=c(0.1,1000),ylim=c(0.1,1000),log='xy',axes=F)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1000,text='D',cex=1.3,line=2.75,las=2)

plot(D,sub_D,pch=21,cex=0.85,col='white',bg='violet',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('Full sample 1/',italic(D))),ylab=expression(paste('1/',italic(D),' without dominant')),xlim=c(0.1,1000),ylim=c(0.1,1000),log='xy',axes=F)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1000,text='E',cex=1.3,line=2.65,las=2)

plot(J,sub_J,pch=21,cex=0.85,col='white',bg='green3',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('Full sample ',italic(J))),ylab=expression(paste(italic(J),' without dominant')),xlim=c(0,1),ylim=c(0,1))
abline(0,1,col='gold')
mtext(side=2,at=1,text='F',cex=1.3,line=2.65,las=2)

dev.off()


# statistics reported in the figure caption

HD <- log(H / D)
sub_HD <- log(sub_H / sub_D)

cor.test(pln_richness,sub_pln_richness,method='s')$estimate
cor.test(cegs_richness,sub_cegs_richness,method='s')$estimate
cor.test(alpha,sub_alpha,method='s')$estimate
cor.test(H,sub_H,method='s')$estimate
cor.test(D,sub_D,method='s')$estimate
cor.test(J,sub_J,method='s')$estimate
cor.test(HD,sub_HD,method='s')$estimate


# BOOTSTRAP ANALYSIS

boot_pln_sigma <- array()
boot_pln_richness <- array()
boot_cegs_gamma <- array()
boot_cegs_richness <- array()
boot_alpha <- array()
boot_H <- array()
boot_D <- array()
boot_J <- array()

d <- 0

for (i in sample(nos))	{
	n <- co[sample == i]
	if (length(table(n)) < 3)
		next
	d <- d + 1
	if (d %% 100 == 0)
		cat('\r',d)
	h <- sample(1:length(n),ceiling(length(n) / 2))
	n <- sample(n,replace=T)
	if (length(table(n)) < 3)
		next
	p <- try(poilogMLE(n))
	if (length(p) > 1)	{
		boot_pln_sigma[i] <- p$par[2]
		boot_pln_richness[i] <- length(n) / p$p
	}
	eg <- try(cegsML(trim(n)))
	if (length(eg) > 1)	{
		boot_cegs_gamma[i] <- eg$shape
		boot_cegs_richness[i] <- eg$richness
	}
	boot_alpha[i] <- fisher(n)
	boot_H[i] <- shannon(n)
	boot_D[i] <- simpson(n)
	boot_J[i] <- log(boot_H[i]) / log(length(n))
}

x <- cbind(boot_pln_sigma,boot_pln_richness,boot_cegs_gamma,boot_cegs_richness,boot_alpha,boot_H,boot_D,boot_J)

colnames(x) <- c('bootstrapped PLN sigma','bootstrapped PLN richness','bootstrapped CEGS gamma','bootstrapped CEGS richness','bootstrapped alpha','bootstrapped H','bootstrapped D','bootstrapped J')

write.table(x,file='bootstrapped_estimates.txt',sep='\t',quote=F)


# FIGURE 2: BOOTSTRAP ANALYSIS

x <- read.delim('richness_estimates.txt')

pln_richness <- x[,1]
cegs_richness <- x[,2]

x <- read.delim('basic_statistics.txt')

alpha <- x[,3]
H <- x[,4]
D <- x[,5]
J <- x[,6]

x <- read.delim('bootstrapped_estimates.txt')

while (nrow(x) < length(H))
	x <- rbind(x,NA)

boot_pln_richness <- x[,2]
boot_cegs_richness <- x[,4]
boot_alpha <- x[,5]
boot_H <- x[,6]
boot_D <- x[,7]
boot_J <- x[,8]


pdf(file='Evenness_Fig_2.pdf',height=11.25,width=7.5)

par(mfrow=c(3,2))
par(mai=c(0.8,0.8,0.3,0.3))
par(cex=1)

plot(pln_richness,boot_pln_richness,pch=21,cex=0.85,col='white',bg='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Original PLN',ylab='Bootstrapped PLN',xlim=c(1,1e4),ylim=c(1,1e4),log='xy',axes=F)
axis(side=1,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e4,text='A',cex=1.3,line=2.65,las=2)

plot(cegs_richness,boot_cegs_richness,pch=21,cex=0.85,col='white',bg='dodgerblue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Original CEGS',ylab='Bootstrapped CEGS',xlim=c(1,1e4),ylim=c(1,1e4),log='xy',axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e4,text='B',cex=1.3,line=2.65,las=2)

plot(alpha,boot_alpha,pch=21,cex=0.85,col='white',bg='blue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste("Original Fisher's ",alpha)),ylab=expression(paste("Bootstrapped Fisher's ",alpha)),xlim=c(0.1,1e3),ylim=c(0.1,1e3),log='xy',axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e3,text='C',cex=1.3,line=2.7,las=2)

plot(H,boot_H,pch=21,cex=0.85,col='white',bg='orange',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('Original exp(',italic(H),')')),ylab=expression(paste('Bootstrapped exp(',italic(H),')')),xlim=c(0.1,1000),ylim=c(0.1,1000),log='xy',axes=F)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1000,text='D',cex=1.3,line=2.75,las=2)

plot(D,boot_D,pch=21,cex=0.85,col='white',bg='violet',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('Original 1/',italic(D))),ylab=expression(paste('Bootstrapped 1/',italic(D))),xlim=c(0.1,1000),ylim=c(0.1,1000),log='xy',axes=F)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1000,text='E',cex=1.3,line=2.8,las=2)

plot(J,boot_J,pch=21,cex=0.85,col='white',bg='green3',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('Original ',italic(J))),ylab=expression(paste('Bootstrapped ',italic(J))),xlim=c(0,1),ylim=c(0,1))
abline(0,1,col='gold')
mtext(side=2,at=1,text='F',cex=1.3,line=2.8,las=2)

dev.off()


# statistics reported in the figure caption

HD <- log(H / D)
boot_HD <- log(boot_H / boot_D)

cor.test(boot_pln_richness,pln_richness,method='s')$estimate
cor.test(boot_cegs_richness,cegs_richness,method='s')$estimate
cor.test(boot_alpha,alpha,method='s')$estimate
cor.test(boot_H,H,method='s')$estimate
cor.test(boot_D,D,method='s')$estimate
cor.test(boot_J,J,method='s')$estimate
cor.test(boot_HD,HD,method='s')$estimate

sum(boot_pln_richness > pln_richness,na.rm=T) / sum(! is.na(boot_pln_richness))
sum(boot_cegs_richness > cegs_richness,na.rm=T) / sum(! is.na(boot_cegs_richness))
sum(boot_alpha > alpha,na.rm=T) / sum(! is.na(boot_alpha))
sum(boot_H > H,na.rm=T) / sum(! is.na(boot_H))
sum(boot_D > D,na.rm=T) / sum(! is.na(boot_D))
sum(boot_J > J,na.rm=T) / sum(! is.na(boot_J))
sum(boot_HD > HD,na.rm=T) / sum(! is.na(boot_J))


# MATCHED SAMPLE ANALYSIS

# matches are based on sums of reciprocals of counts

ls <- array()
nn <- array()

for (i in which(! is.na(S)))
	ls[i] <- sum(1 / co[sample == i])

for (i in which(! is.na(S)))	{
	d <- abs(ls[i] - ls)
	d[i] <- NA
	nn[i] <- which(d == min(d,na.rm=T))[1]
}

# save the data

x <- cbind(1:length(nn),nn)
colnames <- c('sample number','matched.sample')

write.table(x,file='matched_samples.txt',sep='\t',quote=F)


# FIGURE 3: MATCHED SAMPLES

x <- read.delim('richness_estimates.txt')

pln_richness <- x[,1]
cegs_richness <- x[,2]

x <- read.delim('basic_statistics.txt')

alpha <- x[,3]
H <- x[,4]
D <- x[,5]
J <- x[,6]

nn <- read.delim('matched_samples.txt')[,2]
nn[nrow(x)] <- NA


pdf(file='Evenness_Fig_3.pdf',height=11.25,width=7.5)

par(mfrow=c(3,2))
par(mai=c(0.8,0.8,0.3,0.3))
par(cex=1)

plot(pln_richness,pln_richness[nn],pch=21,cex=0.85,col='white',bg='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='PLN richness estimate',ylab='Matched PLN estimate',xlim=c(1,1e4),ylim=c(1,1e4),log='xy',axes=F)
axis(side=1,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e4,text='A',cex=1.3,line=2.65,las=2)

plot(cegs_richness,cegs_richness[nn],pch=21,cex=0.85,col='white',bg='dodgerblue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='CEGS richness estimate',ylab='Matched CEGS estimate',xlim=c(1,1e4),ylim=c(1,1e4),log='xy',axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e4,text='B',cex=1.3,line=2.65,las=2)

plot(alpha,alpha[nn],pch=21,cex=0.85,col='white',bg='blue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste("Fisher's ",alpha)),ylab=expression(paste("Matched Fisher's ",alpha)),xlim=c(0.1,1e3),ylim=c(0.1,1e3),log='xy',axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1e3,text='C',cex=1.3,line=2.65,las=2)

plot(H,H[nn],pch=21,cex=0.85,col='white',bg='orange',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('exp(',italic(H),')')),ylab=expression(paste('Matched exp(',italic(H),')')),xlim=c(0.1,1000),ylim=c(0.1,1000),log='xy',axes=F)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1000,text='D',cex=1.3,line=2.75,las=2)

plot(D,D[nn],pch=21,cex=0.85,col='white',bg='violet',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste('1/',italic(D))),ylab=expression(paste('Matched 1/',italic(D))),xlim=c(0.1,1000),ylim=c(0.1,1000),log='xy',axes=F)
axis(side=1,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
abline(0,1,col='gold')
mtext(side=2,at=1000,text='E',cex=1.3,line=2.75,las=2)

plot(J,J[nn],pch=21,cex=0.85,col='white',bg='green3',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste("Pielou's ",italic(J))),ylab=expression(paste('Matched ',italic(J))),xlim=c(0,1),ylim=c(0,1))
abline(0,1,col='gold')
mtext(side=2,at=1,text='F',cex=1.3,line=2.75,las=2)

dev.off()


# statistics reported in the figure caption

HD <- log(H / D)

cor.test(pln_richness,pln_richness[nn],method='s')$estimate
cor.test(cegs_richness,cegs_richness[nn],method='s')$estimate
cor.test(alpha,alpha[nn],method='s')$estimate
cor.test(H,H[nn],method='s')$estimate
cor.test(D,D[nn],method='s')$estimate
cor.test(J,J[nn],method='s')$estimate
cor.test(HD,HD[nn],method='s')$estimate


# FIGURE 4: DISTRIBUTION PARAMETERS

x <- read.delim('distribution_parameters.txt')

pln_sigma <- x[,1]
cegs_gamma <- x[,2]

x <- read.delim('basic_statistics.txt')

H <- x[,4]
D <- x[,5]

HD <- log(H / D)


pdf(file='Evenness_Fig_4.pdf',height=7.5,width=7.5)

par(mfrow=c(2,2))
par(mai=c(0.8,0.8,0.3,0.3))
par(cex=1)

plot(J,pln_sigma,pch=21,cex=0.85,col='white',bg='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste("Pielou's ",italic(J))),ylab=expression(paste('PLN parameter ',sigma)),xlim=c(0,1),ylim=c(0,10),axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,mgp=c(2.5,0.75,0))
mtext(side=2,at=10,text='A',cex=1.3,line=2.75,las=2)

plot(HD,pln_sigma,pch=21,cex=0.85,col='white',bg='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Hill ratio',ylab=expression(paste('PLN parameter ',sigma)),xlim=c(0,1.2),ylim=c(0,10),axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=seq(0,10,2),mgp=c(2.5,0.75,0))
mtext(side=2,at=10,text='B',cex=1.3,line=2.75,las=2)

plot(J,cegs_gamma,pch=21,cex=0.85,col='white',bg='dodgerblue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab=expression(paste("Pielou's ",italic(J))),ylab=expression(paste('CEGS parameter ',gamma)),xlim=c(0,1),ylim=c(0,10),axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,mgp=c(2.5,0.75,0))
mtext(side=2,at=10,text='C',cex=1.3,line=2.75,las=2)

plot(HD,cegs_gamma,pch=21,cex=0.85,col='white',bg='dodgerblue',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Hill ratio',ylab=expression(paste('CEGS parameter ',gamma)),xlim=c(0,1.2),ylim=c(0,10),axes=F)
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=seq(0,10,2),mgp=c(2.5,0.75,0))
mtext(side=2,at=10,text='D',cex=1.3,line=2.75,las=2)

dev.off()


# statistics reported in the figure caption

cor.test(J,pln_sigma,method='s')$estimate
cor.test(HD,pln_sigma,method='s')$estimate
cor.test(J,cegs_gamma,method='s')$estimate
cor.test(HD,cegs_gamma,method='s')$estimate


# SIMULATION ANALYSIS

# richness is fixed and counts derive from two distributions
# statistics are reported in the text

# s_ means observed richness and h_ means observed Shannon's H
# the distributions are:
# ln = Poisson-sampled log normal
# eg = CEGS
# the distribution parameters are explained in the text

s_ln <- array()
s_eg <- array()
h_ln <- array()
h_eg <- array()
d_ln <- array()
d_eg <- array()

# generate 10000 count distributions and compute S and H

for (i in 1:10000)	{
	n <- rpois(100,exp(rnorm(100,mean=0,sd=2)))
	s_ln[i] <- sum(n > 0)
	h_ln[i] <- shannon(n[n > 0])
	d_ln[i] <- simpson(n)
	n <- rcegs(100)
	s_eg[i] <- sum(n > 0)
	h_eg[i] <- shannon(n[n > 0])
	d_eg[i] <- simpson(n)
}

# save the data

x <- cbind(s_ln,h_ln,d_ln,s_eg,h_eg,d_eg)
colnames(x) <- c('S and PLN','H and PLN','D and PLN','S and CEGS','H and CEGS','D and CEGS')

write.table(x,'S_H_and_D_simulation.txt',sep='\t',quote=F)

x <- read.delim('S_H_and_D_simulation.txt')

s_ln <- x[,1]
h_ln <- x[,2]
d_ln <- x[,3]
s_eg <- x[,4]
h_eg <- x[,5]
d_eg <- x[,6]

# summary statistics reported in the text

exp(mean(log(s_ln)))
exp(mean(log(s_eg)))

sd(log(s_ln))
sd(log(s_eg))

mean(h_ln)
mean(h_eg)

sd(log(h_ln))
sd(log(h_eg))

sd(log(d_ln))
sd(log(d_eg))


# FIGURE 5: SIMULATION RESULTS

x <- read.delim('S_H_and_D_simulation.txt')

s_ln <- x[,1]
h_ln <- x[,2]
d_ln <- x[,3]
s_eg <- x[,4]
h_eg <- x[,5]
d_eg <- x[,6]


pdf(file='Evenness_Fig_5.pdf',height=8,width=4.5)

par(mfrow=c(2,1))
par(mai=c(0.8,0.8,0.3,0.3))
par(cex=1)

hist(log(s_ln),xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),border=NA,col='grey80',main=NA,breaks=10,xlim=log(c(1,100)),ylim=c(0,2500),cex.lab=1.2,xlab='Diversity metric',ylab='Number of trials',axes=F,xpd=NA)
axis(side=1,at=log(c(1,3,10,30,100)),labels=c('1','3','10','30','100'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0,500,1000,1500,2000,2500,3000),labels=c('0','500','1000','1500','2000','2500','3000'),mgp=c(2.5,0.75,0))
hist(log(h_ln),add=T,border=NA,col=hsv(h=0.1,alpha=0.5),main=NA,breaks=40,xlim=log(c(1,100)),axes=F)
hist(log(d_ln),add=T,border=NA,col=hsv(h=0.8,alpha=0.5),main=NA,breaks=40,xlim=log(c(1,100)),axes=F)
text(0.7,2300,labels='PLN data',cex=0.9)
text(0.45,1700,labels='richness',cex=0.9,pos=4)
text(0.45,1520,labels=expression(italic(H)),cex=0.9,pos=4)
text(0.45,1340,labels=expression(italic(D)),cex=0.9,pos=4)
rect(0.3,1680,0.45,1730,col='grey80',border=NA)
rect(0.3,1500,0.45,1550,col=hsv(h=0.1,alpha=0.5),border=NA)
rect(0.3,1320,0.45,1370,col=hsv(h=0.8,alpha=0.5),border=NA)
mtext(side=2,at=2500,text='A',cex=1.35,line=2.6,las=2)

hist(log(s_eg),xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),border=NA,col='grey80',main=NA,breaks=10,xlim=log(c(1,100)),ylim=c(0,2500),cex.lab=1.2,xlab='Diversity metric',ylab='Number of trials',axes=F,xpd=NA)
axis(side=1,at=log(c(1,3,10,30,100)),labels=c('1','3','10','30','100'),mgp=c(2.5,0.75,0))
axis(side=2,at=c(0,500,1000,1500,2000,2500),labels=c('0','500','1000','1500','2000','2500'),mgp=c(2.5,0.75,0))
hist(log(h_eg),add=T,border=NA,col=hsv(h=0.1,alpha=0.5),main=NA,breaks=40,xlim=log(c(1,100)),axes=F)
hist(log(d_eg),add=T,border=NA,col=hsv(h=0.8,alpha=0.5),main=NA,breaks=40,xlim=log(c(1,100)),axes=F)
text(0.7,2300,labels='CEGS data',cex=0.9)
mtext(side=2,at=2500,text='B',cex=1.35,line=2.6,las=2)

dev.off()


# FIGURE 6: LATITUDINAL GRADIENTS

x <- read.delim('richness_estimates.txt')

pln_richness <- x[,1]
cegs_richness <- x[,2]

x <- read.delim('basic_statistics.txt')

H <- x[,4]
D <- x[,5]
J <- x[,6]

x <- read.delim('metadata.txt')

lat <- x[,2]
long <- x[,3]
form <- x[,4]

t <- which(form == 'trees' & long < -30 & ! is.na(pln_richness) & ! is.na(cegs_richness))
b <- which(form == 'bats' & long < -30 & ! is.na(pln_richness) & ! is.na(cegs_richness))

ot <- order(lat[t])
ob <- order(lat[b])

# statistics reported in the caption
# there are no values for J because it is not in units analogous to species

sd(log(pln_richness[t]))
sd(log(pln_richness[b]))

sd(log(cegs_richness[t]))
sd(log(cegs_richness[b]))

sd(log(alpha[t]))
sd(log(alpha[b]))

sd(log(H[t]))
sd(log(H[b]))

sd(log(D[t]))
sd(log(D[b]))


pdf(file='Evenness_Fig_6.pdf',height=11.25,width=7.5)

par(mfrow=c(3,2))
par(mai=c(0.8,0.8,0.3,0.3))
par(cex=1)

plot(lat[t],pln_richness[t],pch=19,cex=0.4,col='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Latitude',ylab='PLN richness estimate',xlim=c(-60,60),ylim=c(1,1e4),log='y',axes=F)
points(lat[b],pln_richness[b],pch=19,cex=0.4,col='dodgerblue')
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
mtext(side=2,at=1e4,text='A',cex=1.3,line=2.65,las=2)

l <- exp(predict(loess(log(pln_richness[t]) ~ lat[t])))
lines(lat[t[ot]],l[ot],col='red')

l <- exp(predict(loess(log(pln_richness[b]) ~ lat[b])))
lines(lat[b[ob]],l[ob],col='dodgerblue')

points(-45,5000,col='red',pch=19,cex=0.9)
points(-45,3000,col='dodgerblue',pch=19,cex=0.9)
text(-45,5000,labels='trees',cex=0.9,pos=4)
text(-45,3000,labels='bats',cex=0.9,pos=4)


plot(lat[t],cegs_richness[t],pch=19,cex=0.4,col='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Latitude',ylab='CEGS richness estimate',xlim=c(-60,60),ylim=c(1,1e4),log='y',axes=F)
points(lat[b],cegs_richness[b],pch=19,cex=0.4,col='dodgerblue')
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(1,10,100,1000,1e4),mgp=c(2.5,0.75,0))
mtext(side=2,at=1e4,text='B',cex=1.3,line=2.65,las=2)

l <- exp(predict(loess(log(cegs_richness[t]) ~ lat[t])))
lines(lat[t[ot]],l[ot],col='red')

l <- exp(predict(loess(log(cegs_richness[b]) ~ lat[b])))
lines(lat[b[ob]],l[ob],col='dodgerblue')


plot(lat[t],alpha[t],pch=19,cex=0.4,col='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Latitude',ylab=expression(paste("Fisher's ",alpha)),xlim=c(-60,60),ylim=c(0.1,1e3),log='y',axes=F)
points(lat[b],alpha[b],pch=19,cex=0.4,col='dodgerblue')
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),mgp=c(2.5,0.75,0))
mtext(side=2,at=1e3,text='C',cex=1.3,line=2.65,las=2)

l <- exp(predict(loess(log(alpha[t]) ~ lat[t])))
lines(lat[t[ot]],l[ot],col='red')

l <- exp(predict(loess(log(alpha[b]) ~ lat[b])))
lines(lat[b[ob]],l[ob],col='dodgerblue')


plot(lat[t],H[t],pch=19,cex=0.4,col='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Latitude',ylab=expression(paste('exp(',italic(H),')')),xlim=c(-60,60),ylim=c(0.1,1e3),log='y',axes=F)
points(lat[b],H[b],pch=19,cex=0.4,col='dodgerblue')
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
mtext(side=2,at=1e3,text='D',cex=1.3,line=2.65,las=2)

l <- exp(predict(loess(log(H[t]) ~ lat[t])))
lines(lat[t[ot]],l[ot],col='red')

l <- exp(predict(loess(log(H[b]) ~ lat[b])))
o <- order(lat[b])
lines(lat[b[ob]],l[ob],col='dodgerblue')


plot(lat[t],D[t],pch=19,cex=0.4,col='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Latitude',ylab=expression(paste('1/',italic(D))),xlim=c(-60,60),ylim=c(0.1,1e3),log='y',axes=F)
points(lat[b],D[b],pch=19,cex=0.4,col='dodgerblue')
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,at=c(0.1,1,10,100,1000),labels=c('0.1','1','10','100','1000'),mgp=c(2.5,0.75,0))
mtext(side=2,at=1e3,text='E',cex=1.3,line=2.75,las=2)

l <- exp(predict(loess(log(D[t]) ~ lat[t])))
lines(lat[t[ot]],l[ot],col='red')

l <- exp(predict(loess(log(D[b]) ~ lat[b])))
o <- order(lat[b])
lines(lat[b[ob]],l[ob],col='dodgerblue')


plot(lat[t],J[t],pch=19,cex=0.4,col='red',bty='l',xaxs='i',yaxs='i',mgp=c(2.5,0.75,0),cex.lab=1.2,xlab='Latitude',ylab=expression(paste("Pielou's ",italic(J))),xlim=c(-60,60),ylim=c(0,1),axes=F)
points(lat[b],J[b],pch=19,cex=0.4,col='dodgerblue')
axis(side=1,mgp=c(2.5,0.75,0))
axis(side=2,mgp=c(2.5,0.75,0))
mtext(side=2,at=1,text='F',cex=1.3,line=2.65,las=2)

l <- predict(loess(J[t] ~ lat[t]))
lines(lat[t[ot]],l[ot],col='red')

l <- predict(loess(J[b] ~ lat[b]))
o <- order(lat[b])
lines(lat[b[ob]],l[ob],col='dodgerblue')

dev.off()

