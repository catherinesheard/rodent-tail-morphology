library(MCMCglmm)
library(phangorn)

x<-read.csv("RodentData.csv") #read in data

#if you would like to run a model on only a suborder, take a subset at this point
#e.g., x<-x[which(x$Suborder=="Scuriomorpha"),]

t100<-read.tree("RodentTrees.tre") #read in trees -- these are 100 of the trees from Upham et al 2019, trimmed to only the rodents

#need to get the dataset and phylogeny to match
tree<-t100[[1]] #this is arbitrary

if(sum(is.na(x$UphamTreeName.full)>0)){x<-x[-which(is.na(x$UphamTreeName.full)),]}
if(sum(is.na(x$Fatty)>0)){x<-x[-which(is.na(x$Fatty)),]}
if(sum(is.na(x$Mean_Annual_Temp )>0)){x<-x[-which(is.na(x$Mean_Annual_Temp)),]}
if(sum(is.na(x$Hib)>0)){x<-x[-which(is.na(x$Hib)),]}
if(sum(is.na(x$PrecipVar)>0)){x<-x[-which(is.na(x$PrecipVar)),]}


t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full)) #trim out everything from the tree that's not in the dataset

#t100 and x should now match


######
#prepare the data
######

x$zTemp<-scale(x$Mean_Annual_Temp)
x$zPV<-scale(x$PrecipVar)



######
#set up a dummy run
######

i=1 #this is arbitrary

tree<-t100[[i]]  
tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) 

animalA<-inverseA(tree)$Ainv 


gelmanprior<-list(B=list(mu=rep(0,4), #set up priors for this model
                         V=gelman.prior(~as.factor(Hib)+zTemp+zPV, 
                                        data = x,  scale=1+pi^2/3)), 
                  R=list(V=1,fix=1),G=list(G1=list(V=1E-10,nu=-1)))

prior.PN<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

mod<-MCMCglmm(as.factor(Fatty)~as.factor(Hib)+zTemp+zPV,
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
              pr=TRUE, 
              slice=TRUE) 







Final.disp<-mod #set up a structure that we'll populate with the full model
Final.disp$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] 
Final.disp$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
Final.disp$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 

nsamp.l<-nrow(mod$VCV)
start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"]))

save(Final.disp,file="FattyDeposits.Rdata")

######
#run the full model
######

for(i in 1:100){ #loop through 100 trees
  tree<-t100[[i]] #select the ith tree 
  tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) #force ultrametric
  
  animalA<-inverseA(tree)$Ainv 
  
  mod<-MCMCglmm(as.factor(Fatty)~as.factor(Hib)+zTemp+zPV,
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
                pr=TRUE,
                slice=TRUE)
  
  print(i) 
  
  Final.disp$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,]  #put this tree's run into the overall structure
  Final.disp$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
  Final.disp$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 
  
  nsamp.l<-nrow(mod$VCV)
  start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"]))
  
  save(Final.disp,file="FattyDeposits.Rdata")
  
}

save(Final.disp,file="FattyDeposits.Rdata")