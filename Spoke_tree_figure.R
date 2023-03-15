library(ape)
library(geiger)
library(nlme)
library(tidyverse)
library(fishtree)
library(ggtree)
library(wesanderson)

#fishtree species
rob.tree <- fishtree_phylogeny(type = "chronogram")
ape::is.ultrametric(rob.tree, tol = 0.00001)
#retrieve tip labels (species) and remove underscore
rob.sp <- gsub("_"," ",rob.tree$tip.label)
rob.tree$tip.label <- gsub("_"," ",rob.tree$tip.label)


d<-  readRDS("data0301.RDS")
fb <- readRDS("fishbase.data.RDS")
gbif <- readRDS("GBIF.data.table.RDS")
#Latsum
lat.sum <- gbif%>%
  group_by(species)%>%
  dplyr::summarise(abslat=mean(abs(decimalLatitude)), lat_sd=sd(decimalLatitude))
colnames(lat.sum) <- c('Species', 'lat', 'lat_sd')
lat.sum <- fb%>%
  left_join(lat.sum)%>%
  mutate(abslat = abs(lat))%>%
  filter(!is.na(abslat))
latsum <- lat.sum%>%
  mutate(abslat = abs(lat))%>%
  filter(!is.na(abslat))

rob.tree2 <- keep.tip(rob.tree,latsum$Species)
rob.tree2$edge.length <- rob.tree2$edge.length*5000
t_rob <- system.time(
  pglsOU_rob<- gls(log.L~abslat*mar.fresh2,data = latsum,correlation = corMartins(1,phy = rob.tree2,form=~Species),method = "ML")
)


### tree figures



bl <- latsum$log.L
lat <- latsum$abslat
names(bl) <- names(lat) <- latsum$Species
lat.pic <- pic(lat,rob.tree2)
bl.pic <- pic(bl,rob.tree2)
plot(lat.pic,bl.pic)

tip_df <- tibble(id=latsum$Species,marfresh=latsum$mar.fresh2,bl=latsum$log.L,lat=latsum$abslat,y2=p$data %>% filter(isTip==T) %>% pull(y))


p<- ggtree(rob.tree2)
tip_df2 <- tibble(id=names(lat),bl=bl.pic,lat=loglat.pic) %>% na.omit

facet_plot(p,data=tip_df,geom=geom_point,panel="bl",mapping=aes(,x=bl,col=lat),size=0.3,stat="identity")+scale_colour_gradientn(colours = viridis(50, option = "D"))

gheatmap(p,data=tip_df)

p<- ggtree(rob.tree2,layout="fan",open.angle = 180,alpha=.3,size=0.1)

tip_df <- tibble(label=latsum$Species,marfresh=latsum$mar.fresh2,bl=latsum$log.L,lat=latsum$abslat) %>% left_join(p$data %>% filter(isTip==T))

pal <- wes_palette(21, name = "Zissou1", type = "continuous")

p2 <- p+geom_spoke(aes(x=x,y=y,angle=0,radius=bl*200000,col=lat),linewidth=0.01,data=tip_df,show.legend=TRUE)+scale_colour_gradientn(colours = pal)+theme(legend.position=c(0.475,0.5),legend.direction = "horizontal",legend.title.align=0.5,legend.title = element_blank())+geom_text(aes(Inf,Inf,label="Abs. Lat."),inherit.aes = F,size=4,vjust=-2.6,hjust=3.7)

ggsave("BL_Lat_phylo.pdf",p2)

