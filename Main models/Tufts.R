library(MCMCglmm)
library(phangorn)

x<-read.csv("RodentData.csv") #read in data

#if you would like to run a model on only a suborder, take a subset at this point
#e.g., x<-x[which(x$Suborder=="Scuriomorpha"),]

t100<-read.tree("RodentTrees.tre") #read in trees -- these are 100 of the trees from Upham et al 2019, trimmed to only the rodents

#need to get the dataset and phylogeny to match
tree<-t100[[1]] #this is arbitrary

if(sum(is.na(x$UphamTreeName.full)>0)){x<-x[-which(is.na(x$UphamTreeName.full)),]}
if(sum(is.na(x$Tufts)>0)){x<-x[-which(is.na(x$Tufts)),]}
if(sum(is.na(x$A)>0)){x<-x[-which(is.na(x$A)),]}
if(sum(is.na(x$T)>0)){x<-x[-which(is.na(x$T)),]}
if(sum(is.na(x$OwlThreat)>0)){x<-x[-which(is.na(x$OwlThreat)),]}
if(sum(is.na(x$RaptorThreatNonOwl)>0)){x<-x[-which(is.na(x$RaptorThreatNonOwl)),]}
if(sum(is.na(x$Openness.Score)>0)){x<-x[-which(is.na(x$Openness.Score)),]}
if(sum(is.na(x$Mass)>0)){x<-x[-which(is.na(x$Mass)),]}

t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full)) #trim out everything from the tree that's not in the dataset

#t100 and x should now match


######
#prepare the data
######

x$zMass<-scale(x$Mass)
x$zOwl<-scale(x$OwlThreat)
x$zHawk<-scale(x$RaptorThreatNonOwl)
x$zOpen<-scale(x$Openness.Score)



######
#set up a dummy run
######

i=1 #this is arbitrary

tree<-t100[[i]]  
tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) 

animalA<-inverseA(tree)$Ainv 


gelmanprior<-list(B=list(mu=rep(0,7), #set up priors for this model
                         V=gelman.prior(~as.factor(T)+as.factor(A)+zMass+zOwl+zHawk+zOpen, 
                                        data = x,  scale=1+pi^2/3)), 
                  R=list(V=1,fix=1),G=list(G1=list(V=1E-10,nu=-1)))

mod<-MCMCglmm(as.factor(Tufts)~as.factor(T)+as.factor(A)+zMass+zOwl+zHawk+zOpen,
              random=~UphamTreeName.full, 
              ginverse=list(UphamTreeName.full=animalA), 
              prior = gelmanprior, 
              verbose=TRUE, 
              family="categorical", 
              data = x,
              nitt=11000,
              thin=10,
              burnin=1000,
              pl=TRUE, 
              pr=TRUE) 







Final.mod<-mod #set up a structure that we'll populate with the full model
Final.mod$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] 
Final.mod$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
Final.mod$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 

nsamp.l<-nrow(mod$VCV)
start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"]))

save(Final.mod,file="tufts.Rdata")

######
#run the full model
######

for(i in 1:100){ #loop through 100 trees
  tree<-t100[[i]] #select the ith tree 
  tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) #force ultrametric
  
  animalA<-inverseA(tree)$Ainv 
  
  mod<-MCMCglmm(as.factor(Tufts)~as.factor(T)+as.factor(A)+zMass+zOwl+zHawk+zOpen,
                random=~UphamTreeName.full, 
                ginverse=list(UphamTreeName.full=animalA), 
                prior = gelmanprior, 
                verbose=FALSE, 
                family="categorical", 
                start= start1.l,
                data = x,
                nitt=30000,
                thin=2000,
                burnin=10000,
                pl=TRUE,
                pr=TRUE)
  
  print(i) 
  
  Final.mod$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,]  #put this tree's run into the overall structure
  Final.mod$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
  Final.mod$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 
  
  nsamp.l<-nrow(mod$VCV)
  start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"]))
  
  save(Final.mod,file="tufts.Rdata") 
  
}

save(Final.mod,file="tufts.Rdata")