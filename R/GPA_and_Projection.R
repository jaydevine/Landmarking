#-----------------------------------------------------------------------------------------------------------------------
# GPA configurations onto manual training mean, then project into tangent space prior to training. 
#-----------------------------------------------------------------------------------------------------------------------
setwd("/path/to/wd/")
#-----------------------------------------------------------------------------------------------------------------------

# Load libraries.
library(Morpho)
library(geomorph)

# Read in the datasets. 
MAN_Raw_Data_Train_Specs <- read.csv("MAN_Raw_Data_Train.csv", header=T)
MAN_Raw_Data_Train <- MAN_Raw_Data_Train_Specs[,3:ncol(MAN_Raw_Data_Train_Specs)]
MAN_Raw_Data_Test_Specs <- read.csv("MAN_Raw_Data_Test.csv", header=T)
MAN_Raw_Data_Test <- MAN_Raw_Data_Test_Specs[,3:ncol(MAN_Raw_Data_Test_Specs)]
Large_Raw_Data_Train_Specs <- read.csv("Large_Raw_Data_Train.csv", header=T)
Large_Raw_Data_Train <- Large_Raw_Data_Train_Specs[,3:ncol(Large_Raw_Data_Train_Specs)]
Large_Raw_Data_Test_Specs <- read.csv("Large_Raw_Data_Test.csv", header=T)
Large_Raw_Data_Test <- Large_Raw_Data_Test_Specs[,3:ncol(Large_Raw_Data_Test_Specs)]

# Create 3D arrays. 
MAN_Train_3D_Arr <- arrayspecs(MAN_Raw_Data_Train,68,3)
MAN_Test_3D_Arr <- arrayspecs(MAN_Raw_Data_Test,68,3)
SyN_Train_3D_Arr <- arrayspecs(Large_Raw_Data_Train,68,3)
SyN_Test_3D_Arr <- arrayspecs(Large_Raw_Data_Test,68,3)

# GPA manual training data.
MAN_Only_procSym <- procSym(MAN_Train_3D_Arr, reflect = TRUE, CSinit = TRUE, orp = TRUE,
                            tol = 1e-05, pairedLM = NULL, sizeshape = FALSE,
                            use.lm = NULL, center.part = FALSE, weights = NULL,
                            centerweight = FALSE, pcAlign = TRUE, distfun = c("angle", "riemann"),
                            SMvector = NULL, outlines = NULL, deselect = FALSE, recursive = TRUE,
                            iterations = 0, initproc = FALSE, bending = FALSE)

# GPA data to manual training mean, including orthogonal projection.
MAN_Train_Align <- align2procSym(MAN_Only_procSym,MAN_Train_3D_Arr, orp = TRUE)
SyN_Train_Align <- align2procSym(MAN_Only_procSym,SyN_Train_3D_Arr, orp = TRUE)
MAN_Test_Align <- align2procSym(MAN_Only_procSym,MAN_Test_3D_Arr , orp = TRUE)
SyN_Test_Align <- align2procSym(MAN_Only_procSym,SyN_Test_3D_Arr, orp = TRUE)

# Create 2D arrays.
MAN_Train_Align_2d <- two.d.array(MAN_Train_Align)
SyN_Train_Align_2d <- two.d.array(SyN_Train_Align)
MAN_Test_Align_2d <- two.d.array(MAN_Test_Align)
SyN_Test_Align_2d <- two.d.array(SyN_Test_Align)

# Write .csv files for import into Julia. 
write.csv(MAN_Train_Align_2d, "MAN_Train_Orp_Revised.csv")
write.csv(SyN_Train_Align_2d, "SyN_Train_Orp_Revised.csv")
write.csv(MAN_Test_Align_2d, "MAN_Test_Orp_Revised.csv")
write.csv(SyN_Test_Align_2d, "SyN_Test_Orp_Revised.csv")

#---
