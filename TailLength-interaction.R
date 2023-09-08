library(MCMCglmm)
library(phangorn)

x<-read.csv("RodentData.csv") #read in data

#if you would like to run a model on only a suborder, take a subset at this point
#e.g., x<-x[which(x$Suborder=="Scuriomorpha"),]

t100<-read.tree("RodentTrees.tre") #read in trees -- these are 100 of the trees from Upham et al 2019, trimmed to only the rodents

tree<-t100[[1]] #select one tree for trimming purposes

if(sum(is.na(x$UphamTreeName.full)>0)){x<-x[-which(is.na(x$UphamTreeName.full)),]}
if(sum(is.na(x$A)>0)){x<-x[-which(is.na(x$A)),]}
if(sum(is.na(x$F)>0)){x<-x[-which(is.na(x$F)),]}
if(sum(is.na(x$Tail_length)>0)){x<-x[-which(is.na(x$Tail_length)),]}
if(sum(is.na(x$Litter_size)>0)){x<-x[-which(is.na(x$Litter_size)),]}
if(sum(is.na(x$Mean_Annual_Temp)>0)){x<-x[-which(is.na(x$Mean_Annual_Temp)),]}


t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full)) #trim out everything from the tree that's not in the dataset

#t100 and x should now match


######
#prepare the data
######

x$zLength<-scale(x$Tail_length)
x$zLS<-scale(x$Litter_size)
x$zTemp<-scale(x$Mean_Annual_Temp)


######
#set up a dummy run
######

i=1 #this is arbitrary

tree<-t100[[i]]  
tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)

animalA<-inverseA(tree)$Ainv 

prior.PN<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

mod<-MCMCglmm(zLength~zLS*zTemp+zTemp+as.factor(A)*zTemp+as.factor(F)*zTemp,
              random=~UphamTreeName.full, 
              ginverse=list(UphamTreeName.full=animalA), 
              prior = prior.PN, 
              verbose=TRUE, 
              family="gaussian", 
              data = x,
              nitt=11000,
              thin=10,
              burnin=1000,
              pl=TRUE, 
              pr=TRUE, 
              slice=TRUE) 

Final.disp<-mod #set up a structure that we'll populate with the full model
Final.disp$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] 
Final.disp$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
Final.disp$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 

nsamp.l<-nrow(mod$VCV)
start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"]))

save(Final.disp,file="tail-length-inter.Rdata")

######
#run the real model
######

for(i in 1:100){ #loop through 100 trees
  tree<-t100[[i]] #select the ith tree 
  tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) #force ultrametric
  
  animalA<-inverseA(tree)$Ainv 
  
  mod<-MCMCglmm(zLength~zLS*zTemp+zTemp+as.factor(A)*zTemp+as.factor(F)*zTemp,
                random=~UphamTreeName.full, 
                ginverse=list(UphamTreeName.full=animalA), 
                prior = prior.PN, 
                verbose=FALSE, 
                family="gaussian", 
                start= start1.l,
                data = x,
                nitt=30000,
                thin=2000,
                burnin=10000,
                pl=TRUE,
                pr=TRUE,
                slice=TRUE)
  
  print(i)
  
  Final.disp$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,]  #put this tree's run into the overall structure
  Final.disp$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
  Final.disp$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 
  
  nsamp.l<-nrow(mod$VCV)
  start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"]))
  
  save(Final.disp,file="tail-length-inter.Rdata") 
  
}

save(Final.disp,file="tail-length-inter.Rdata")