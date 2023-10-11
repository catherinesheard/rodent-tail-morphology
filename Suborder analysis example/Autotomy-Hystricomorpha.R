#Note: works with R version 4.1.3 (and others). 
#Does NOT work in R version 4.3.1.
#If you get an error after the line "animalA<-inverseA(tree)$Ainv", this is due to a problem with the command "nnls.tree"
#The way to fix this is to load the library "phytools" and then replace the line "tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)" with "tree <- force.ultrametric(tree, method=c("extend"))"


#Hystricomorpha

library(MCMCglmm)
library(phangorn)

x<-read.csv("RodentData.csv") #read in data

#if you would like to run a model on only a suborder, take a subset at this point
x<-x[which(x$Suborder=="Hystricomorpha"),]

t100<-read.tree("RodentTrees.tre") #read in trees -- these are 100 of the trees from Upham et al 2019, trimmed to only the rodents

#need to get the dataset and phylogeny to match
tree<-t100[[1]] #this is arbitrary

#remove rows with NAs in the variables that we're about to analyse
if(sum(is.na(x$UphamTreeName.full)>0)){x<-x[-which(is.na(x$UphamTreeName.full)),]}
if(sum(is.na(x$Autotomy)>0)){x<-x[-which(is.na(x$Autotomy)),]}
if(sum(is.na(x$A)>0)){x<-x[-which(is.na(x$A)),]}
if(sum(is.na(x$T)>0)){x<-x[-which(is.na(x$T)),]}
if(sum(is.na(x$Di)>0)){x<-x[-which(is.na(x$Di)),]}
if(sum(is.na(x$Tail_length)>0)){x<-x[-which(is.na(x$Tail_length)),]}
if(sum(is.na(x$Openness.Score)>0)){x<-x[-which(is.na(x$Openness.Score)),]}
if(sum(is.na(x$Mass)>0)){x<-x[-which(is.na(x$Mass)),]}

t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full)) #trim out everything from the tree that's not in the dataset

#t100 and x should now match


######
#prepare the data
######

#scale the continuous variables, for ease of coefficient interpretation
x$zLength<-scale(x$Tail_length)
x$zOpen<-scale(x$Openness.Score)
x$zMass<-scale(log(x$Mass))


######
#set up a dummy run
######

#this is a short run that will determine a start point for the phylogenetic variance and will create a structure of the correct dimensions for us to save our results into

i=1 #this is arbitrary

tree<-t100[[i]]  
tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) #the tree is very slightly non-ultrametric -- we're just forcing it to be ultrametric
#PLEASE NOTE -- this line breaks in newer versions of R -- the tree won't root
#the function from phytools " tree <- force.ultrametric(tree, method=c("extend")) " should work in that case

animalA<-inverseA(tree)$Ainv #this is line is breaking for you, PLEASE read the above note, and either use a previous version of R (4.1.3 works fine) or replace the command using 'nnls.tree' with the phytools command 'force.ultrametric' (make sure you load in phytools as well!)


#set priors for the model
#Gelman priors (Gelman, A. et al. (2008) The Annals of Appled Statistics 2 4 1360-1383) provide a standardised covariance matrix for the fixed effects
#the residual variance, R, is typically fixed in a logistic/categorical regression, conventionally to 1
#the phylogenetic variance, G, is set with an extremely uninformative inverse Wishart prior
#mu needs to match the number of coefficients estimated by your fixed effects formula -- that's 1 for the intercept, 1 for every continuous variable, and k-1 for every categorical variable of size k

gelmanprior<-list(B=list(mu=rep(0,7), #set up priors for this model
                         V=gelman.prior(~as.factor(T)+as.factor(A)+as.factor(Di)+zLength+zOpen+zMass, #you'll have to adjust this equation to match your linear model, and mu above may have to change
                                        data = x,  scale=1+pi^2/3)), 
                  R=list(V=1,fix=1),G=list(G1=list(V=1E-10,nu=-1)))

#this model considers the presence-absence of tail autonomy as it relates to terrestriality (T), arboreality (A), diurnality (Di), tail length (zLength), habitat openness (zOpen), and body mass (zMass)

mod<-MCMCglmm(as.factor(Autotomy)~as.factor(T)+as.factor(A)+as.factor(Di)+zLength+zOpen+zMass,
              random=~UphamTreeName.full, 
              ginverse=list(UphamTreeName.full=animalA), 
              prior = gelmanprior, 
              verbose=TRUE, 
              family="categorical", 
              data = x,
              nitt=11000, #number of iterations -- for a dummy run, this and the thin/burn-in rate don't matter very much, as long as the total sample number is equal to the number that you're aiming for in the final model
              thin=10, #sampling ("thinning") rate
              burnin=1000, #burn-in rate (number of iterations initially discarded)
              pl=TRUE, 
              pr=TRUE) 


Final.mod<-mod #set up a structure that we'll populate with the full model
Final.mod$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] #random effects (phylogeny)
Final.mod$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] #fixed effects (and node values)
Final.mod$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] #latent variables

nsamp.l<-nrow(mod$VCV)
start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"])) #determine the start point for residual and phylogenetic variance the main model

save(Final.mod,file="autonomy-H.Rdata")

######
#run the full model with start points from dummy run
######

for(i in 1:100){ #loop through 100 trees
  tree<-t100[[i]] #select the ith tree 
  tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) #force ultrametric
  #again, if you're using a newer version of R, your code might break here. see above: either use an older version of R (4.1.3 works fine) or replace this line with tree <- force.ultrametric(tree, method=c("extend"))
  
  animalA<-inverseA(tree)$Ainv 
  
  mod<-MCMCglmm(as.factor(Autotomy)~as.factor(T)+as.factor(A)+as.factor(Di)+zLength+zOpen+zMass,
                random=~UphamTreeName.full, 
                ginverse=list(UphamTreeName.full=animalA), 
                prior = gelmanprior, 
                verbose=FALSE, 
                family="categorical", 
                start= start1.l,
                data = x,
                nitt=30000, #number of iterations per tree -- for this, we were aiming for a posterior sample of 10/tree, with a thin rate large enough to avoid autocorrelation and a burn-in rate large enough to ensure convergence
                thin=2000,
                burnin=10000,
                pl=TRUE,
                pr=TRUE)
  
  print(i) #print which tree you're on (so that you know if the loop's running)
  
  Final.mod$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,]  #puts this tree's run into the overall structure
  Final.mod$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
  Final.mod$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 
  
  nsamp.l<-nrow(mod$VCV)
  start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"UphamTreeName.full"]))
  
  save(Final.mod,file="autonomy-H.Rdata") #this will slow down your model; you can turn it off
  
}

save(Final.mod,file="autonomy-H.Rdata")
#to interact with this file once it's done, load the workspace, make sure the MCMCglmm library is loaded, and then use the "summary" command on the model ("Final.mod")
#you can also interact directly with the posterior distributions, i.e. "Final.mod$Sol"