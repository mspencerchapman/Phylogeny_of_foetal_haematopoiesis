my_working_directory="~/R Work/Fetal HSPCs/Phylogeny_of_foetal_haematopoiesis/"
setwd(my_working_directory)

#Polytomies simulation
#This is a bayesian framework to estimate the upper limit of mutation rate at the 8-cell stage
#(where the mutation rate abruptly drops) using the observed polytomies at this stage
library(ape)
#The data vector is the observed polytomy orders from the 8-cell stage divisions (for which we have resolution from the blood tree)
data=c(2,5,5,10,12)
#Define the plausible mutation rates to test
mutation_rates=seq(0.5,1.0, by=0.05)
#Set up a named "probability" vector to fill
prob=rep(NA,length(mutation_rates))
names(prob)=mutation_rates

#Now test each mutation rate in turn
for(k in 1:length(mutation_rates)) {
  mutation_rate=mutation_rates[k]
  print(mutation_rate)
  poly=vector()
  #Simulate 1000 trees which go through 12 rounds of co-ordinated division and see the observed baseline polytomy
  for(i in 1:1000) {
    #Create the tree
    tree=ape::stree(2**12,type="balanced")
    #Assign the mutations using a poisson distribution with lambda of the mutation rate
    tree$edge.length=rpois(nrow(tree$edge),lambda = mutation_rate)
    #Convert to polytomous tre
    tree.multi=di2multi(tree)
    #Store the observed order of the baseline polytomy/dichotomy
    poly[i]=sum(tree.multi$edge[,1]==(length(tree$tip.label)+1))
  }
  hist(poly,breaks=40)
  
  #Set up empty vector to fill
  pass=vector()
  #Re-sample 5 values from this vector 1000 times
  #For each sample, sort the values from small to large (as with the data)
  #If all values are at least as big as their true data counterparts, the sample is a "PASS" and TRUE is stored
  #Note: the real data may be underestimates of the "true" underlying polytomy, therefore this test only
  #checks if the values are bigger than the data to get an upper limit on the mutation rate at this point
  for(j in 1:1000) {
    sort_samp=sort(sample(poly,size = 5,replace = TRUE))
    pass[j]=all(sort_samp>=data)
  }
  #Store the proportino of draws that are "PASS" draws for a given mutation rate
  prob[k]<-sum(pass)/length(pass)
}

plot(prob~names(prob))

#Display with ggplot
library(ggplot2)
df=data.frame(mutation_rate=as.numeric(names(prob)),p_value=prob)
df %>%
  ggplot(aes(x=mutation_rate,y=p_value)) +
  geom_point()+
  geom_smooth(se = FALSE) +
  scale_y_continuous(breaks=seq(0,0.6,by=0.05))
