#Workflow
#Loading dataset
absolute_dataset <- read.csv("C:/Users/sneha/OneDrive/Desktop/new_scripts/coda4_absolute_2.csv",header=TRUE) #absolute abundance
head(absolute_dataset) #but this package does not read the dataset if the OTU column is there
#Working on zeroes
#zCompositions package
library(zCompositions)
zero_fixed_dataset.GBM <- cmultRepl(absolute_dataset,method="GBM") #Bayesian-Multiplicative replacement of count zeroes
#if in the above command we choose option method= "p-counts", then pseudocounts works only columns holding 0, the columns not having a single 0, remain unchanged. But if we dont select this option then it works on entire dataset. So we are not selecting it.
#Clr-transformation
#Package-compositions
library('compositions')
clr_data_compositions <- compositions::clr(zero_fixed_dataset.GBM) #applying clr transformation
#If the sum of all the values are the clr transformed matrix becomes zero, then it works)
data2 <- clr_data_compositions
head(clr_data_compositions) #result is correct clr trnaformation since the values in the matrix sum upto zero
write.csv(data2, "C:/Users/sneha/OneDrive/Desktop/new_scripts/data2.csv")
#Trying clr transformation with hotelling package
library(Hotelling)
data("bottle.df") #test data set for hotelling packge- includes decimals, but no zeroes
clr_test_data <- clr(bottle.df) #trying clr transformation from hotelling package on its test dataset
#the clr transformed matrix of hotelling test dataset, the matrix summed upto 0.000000000000341061. Is it correct?
write.csv(clr_test_data,"C:/Users/sneha/OneDrive/Desktop/new_scripts/clr_test_data.csv")
#Trying hotelling package on our dataset (the data on which zCompositions was applied)
clr_data_hotelling <- Hotelling::clr(zero_fixed_dataset.GBM)
#clr transformation
#package- Tjazi
#Trying clr_lite function on our dataset of absolute count directly since
library(Tjazi)
head(absolute_dataset)
clr_tjazi_data <- Tjazi::clr_lite(absolute_dataset, samples_are = "cols")
write.csv(clr_data_Tjazi,"C:/Users/sneha/OneDrive/Desktop/new_scripts/clr_tjazi_data.csv") #gives output as 0.03, is that correct?

#Aitchison distance- dissimilarity index

library(vegan)
#test data
data(varespec) #data includes decimals as well as zeroes and only positive values
aitchison_testdata <- vegdist(varespec)
head(aitchison_testdata,10)
head(aitchison_testdata)
aitchison_clr_data <- vegdist(clr_data_compositions, method = "aitchison") #result might be meaningless, as clr transformed data have negative values
print(aitchison_clr_data)
print(clr_data_compositions)
write.csv(aitchison_clr_data,"C:/Users/sneha/OneDrive/Desktop/new_scripts/aitchison_clr_data.csv") 

library(microbiome)
data(dietswap)
x <- dietswap
x_clr <- microbiome::transform(x,'clr')
