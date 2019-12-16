#---------------------------------------------------------------------------------------
# R script for projection of new Procrustes landmark data into training space.
# You can then evaluate the network on the projected data. 
#---------------------------------------------------------------------------------------
# Import packages.
library(geomorph)
library(shapes)
#---------------------------------------------------------------------------------------
# Set working directory.
setwd("/path/to/landmarks/")

# We coerced the automated landmark data (.tag) into the .dta format using
# standard Bash commands. Please see an example of the .dta file in Data.
# Let's assume you have a new .dta file called Test_Landmarks.dta. 
Test_LMs <- readland.nts("Test_Landmarks.dta")
Test_LMs_GPA <- gpagen(Test_LMs)
Test_LMs <- Test_LMs_GPA$coords
Test_LMs_CS <- Test_LMs_GPA$Csize

# Original test set index used for subsetting training/testing datasets.
# We used this in our study. You probably may or may not need to do this.
Test_Index <- read.csv("Test_Set_Index.csv", header=T)

# Load in original training data. Be careful -- what non-linear registration 
# algorithm did you use for training? We're going to use the LUS (SyN/Unoptimized/Single Atlas)
# training data here. 
LUS_Train <- read.csv("Large_Single_and_Manual_Train.csv",header=T)
LUS_Train <- LUS_Train[,2:ncol(LUS_Train)]
LUS_Train <- LUS_Train[1:218,]

# Remove test specimens from sample to isolate training specimens.
LUS_Train <- LUS_Train[-c(Test_Index[,1]),]

# Compute mean shape of training data.
Train_Mean <- matrix(colMeans(LUS_Train), ncol = 3, byrow = T)
Train_Mean <- arrayspecs(Train_Mean,68,3)

# Create empty 3d array for new configurations?
# We have 68 landmarks in 3 dimensions.
Empty <- array(numeric(),c(68,3,0)) 
Empty <- replicate(dim(Test_LMs)[3], Empty)

# Build loop to OPA specimens to training mean. 
for (i in seq(1,dim(Test_LMs)[3])) {
  SUP <- procOPA(Train_Mean[,,1],Test_LMs[,,i], scale=FALSE)
  SUP <- SUP$Bhat
  Empty[[i]] <- SUP
}

# Let's say our new dataset, Test_LMs, has 100 configurations (nrow=100), 
# each with pk dimensions, e.g, (ncol=68*3=204). Obviously this may differ for you.
Test_LMs_OPA <- matrix(nrow=100,ncol=204)
COUNTER=0
for (i in Empty) {
  COUNTER=COUNTER+1
  i <- arrayspecs(i,68,3)
  i <- two.d.array(i)
  Test_LMs_OPA[COUNTER,] <- i
}

# We have specimens that have been procOPA aligned to our automated training mean. 
# Next, we want to PCA the training data and rotate the new data into this space.
# Note: The default for prcomp is centering without scale. 
LUS_PCA <- prcomp(LUS_Train)
Test_PCA <- prcomp(Test_LMs_OPA)

# Generate the projected LM data. Note that the maximum number of PCs below will differ.
Test_LMs_New <- t(t(Test_PCA$x[,1:171] %*% t(LUS_PCA$rotation)) + LUS_PCA$center)

# Write csv,
write.csv(Longshanks_new, "Test_Landmarks_Projected.csv")






