library(fishtree)
library(rfishbase)
library(ape)
library(phytools)
library(geiger)
#library(ggplot2)
#library(rgbif)
#library(maps)
library(tidyverse)
#library(nlme)
#library(devtools)#
#library(parallel)
#library(doParallel)
library(brms)
#library(dplyr)
library(tidybayes)
library(MCMCglmm)
library(data.table)
#library(reshape2)
#library(plyr)
#library(foreach)
#library(raster)#
#library(ncdf4)#
#library(viridis)
#library(rgdal)#
#library(rnaturalearthdata)
#library(rgeos)# - geos-config
#library(modelr)
#library(knitr)

#fishtree species
rob.tree <- fishtree_phylogeny(type = "chronogram")
ape::is.ultrametric(rob.tree, tol = 0.00001)
  #retrieve tip labels (species) and remove underscore
rob.sp <- gsub("_"," ",rob.tree$tip.label)
rob.tree$tip.label <- gsub("_"," ",rob.tree$tip.label)

#species
fb.dat <- readRDS("fishbase.data.RDS")
gbif <- readRDS("GBIF.data.table.RDS")

# limit species
tax <- fishtree_taxonomy()
fams <- tax%>%filter(rank=="family")%>%dplyr::select(name)%>%pull
fish.sp.l <- list()
for(i in fams){
  tax.i <- fishtree_taxonomy(rank=i)
  fish.sp.l[[i]] <- tax.i[[1]]$species}
mean(sapply(fish.sp.l,function(x) length(x)),na.rm=T)
set.seed(123)
sample.sp <- function(x,prop=0.3) {
  x <- unlist(x)
  if(length(x)<=5) y <- x
  if(length(x)>5){
    n <- length(x)
    n.prop <- round(n*prop)
    y <- sample(x,size = n.prop,replace = FALSE)}}
sp.1 <- lapply(fish.sp.l, sample.sp)
sp.1 <- unlist(sp.1)
fam.sp <- fb.dat%>%
  filter(Species%in%sp.1)
rob.fam <- drop.tip(rob.tree,tip=rob.tree$tip.label[!rob.tree$tip.label%in%fam.sp$Species])
rob.fam <- phytools::force.ultrametric(rob.fam)
fam.dat <- data.frame(log.l=fam.sp$log.L)
rownames(fam.dat) <- fam.sp$Species 


#Covariance Matrix
phy=rob.fam
#phy <- root(phy,outgroup = last(rob.fam$tip.label),resolve.root = TRUE)
phy <- phytools::force.ultrametric(phy,method="nnls",message=FALSE)
#using nnls method for ultramecizing causes unrooting, using extend method keeps root
#phy <- phytools::reroot(phy,node.number=1618) - other option for rooting
phy$root.edge<-0 
#keep nnls ultrametric and rooting
is.rooted(phy) 
is.ultrametric(phy)


dat <- fb.dat%>%
  filter(Species%in%phy$tip.label)%>%
  mutate(phylo=Species)%>%
  as.data.frame()
#fix 0-length edges
phy$edge.length[phy$edge.length<0.000001] <- 0.00001
inv.phylo<-inverseA(phy,nodes="all",scale=TRUE,tol = 1e-7)
As <- solve(inv.phylo$Ainv)
rownames(As) <- rownames(inv.phylo$Ainv)


#Latsum
lat.sum <- gbif%>%
  group_by(species)%>%
  dplyr::summarise(abslat=mean(abs(decimalLatitude)), lat_sd=sd(decimalLatitude))
colnames(lat.sum) <- c('Species', 'lat', 'lat_sd')
lat.sum <- dat%>%
  left_join(lat.sum)%>%
  mutate(abslat = abs(lat))%>%
  filter(!is.na(abslat))
latsum <- lat.sum%>%
  mutate(abslat = abs(lat))%>%
  filter(!is.na(abslat))

#Latmodel
g.lat.parallel <- brm(
  log.L ~ abslat*mar.fresh2 + (1|gr(Species, cov = As)),   
  data = latsum, 
  family = student(link = "identity"),
  iter = 1000, 
  chains = 4,
  cores = 4,#comment in on andromeda run
  prior = c(
    prior(normal(0, 10), "b"),
    prior(student_t(3, 0, 20), "sd"),
    prior(student_t(3, 0, 20), "sigma")
  ),
  data2=list(As=As),
  control = list(adapt_delta = 0.99, max_treedepth = (15)),
#  cores = parallel::detectCores(),#comment out on andromeda run
  file = paste0("test_",gsub(" ","",as.character(Sys.time()))),
  backend = "cmdstanr",
  threads=threading(16)
)
saveRDS(g.lat.parallel,"final0.RDS")


