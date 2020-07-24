#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Python script for pairwise spatial normalization using ANIMAL's non-linear algorithm. It will generate a series of Bash scripts that will pairwise register your images to an atlas.
# This should happen on a compute cluster using the MINC toolkit module. Note that the parameters used throughout this script have been adapted for isotropic 35 micron uCT volumes.
# You'll need to adapt the blurring and registration values to your data.

# Citation: Percival, C.J., Devine, J., Darwin, B.C., Liu, W., van Eede, M., Henkelman, R.M. and Hallgrimsson, B., 2019. The effect of automated landmark identification on morphometric analyses. J Anat (2019). https://doi.org/10.1111/joa.12973
# Citation: Devine, J., Aponte, J.D., Katz, D.C. et al. A Registration and Deep Learning Approach to Automated Landmark Detection for Geometric Morphometrics. Evol Biol (2020). https://doi.org/10.1007/s11692-020-09508-8
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Import packages.
import os
import csv
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Some important notes:
# 1) Your local and remote directories should match the compute cluster paths below. Replace <PROJECT> with your project name;
# 2) Your local specimen list must be called test_list.txt (all test specimens).
# 3) Your landmark initialized test images MUST be called $spec.mnc, where $spec is the exact name of the specimen annotated in spec_list.txt;
# 4) The atlas files you intend to register/landmark your test images with must be called NL_4_average.mnc and NL_4_average_landmarks.tag. They must be
# sftp'd into your remote /home/$USER/<PROJECT>/nl/MNC directory on the cluster.
# 5) The mask and initialized source images must be sftp'd into your remote ~/<PROJECT>/Source/MNC directory on the cluster.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define local working directory.
os.chdir("/path/to/project/")
# Define list of all specimens. This should should be a single column of your specimen names (without the .mnc suffix). Make sure there aren't any hidden characters
# in the list (e.g., a space beside a name).
All_Specimens = "/path/to/project/test_list.txt"

# Create remote directory structure that matches your local structure. E.g.:
# mkdir -p <PROJECT\>{Scripts,Source/{aim,Blurred,MNC,Orig,Corr,Tag,Tiff,XFM},lsq6/{Blurred,MNC,XFM},lsq12/{Blurred,MNC,XFM},nl/{Ana_Test,Blurred,INIT,MNC,XFM}}

# Create and define remote directories (i.e., your compute cluster paths). We use the notation below because it is commonly seen in MINC.
Scripts_path = "~/<PROJECT>/Scripts/"
Source_XFM_path = "~/<PROJECT>/Source/XFM/"
Source_MNC_path = "~/<PROJECT>/Source/MNC/"
Source_Tag_path = "~/<PROJECT>/Source/Tag/"
lsq6_Blurred_path = "~/<PROJECT>/lsq6/Blurred/"
lsq6_XFM_path = "~/<PROJECT>/lsq6/XFM/"
lsq6_MNC_path = "~/<PROJECT>/lsq6/MNC/"
lsq12_Blurred_path = "~/<PROJECT>/lsq12/Blurred/"
lsq12_XFM_path = "~/<PROJECT>/lsq12/XFM/"
lsq12_MNC_path = "~/<PROJECT>/lsq12/MNC/"
nl_Init_path = "~/<PROJECT>/nl/INIT/"
nl_Blurred_path = "~/<PROJECT>/nl/Blurred/"
nl_XFM_path = "~/<PROJECT>/nl/XFM/"
nl_MNC_path = "~/<PROJECT>/nl/MNC/"

# Define average files.
LM_Avg = "~/<PROJECT>/Source/MNC/LM_average.mnc"
LM_Avg_Mask = "~/<PROJECT>/Source/MNC/LM_average_mask.mnc"
lsq6_Avg = "~/<PROJECT>/lsq6/<PROJECT>_lsq6_average.mnc"
lsq12_Avg = "~/<PROJECT>/lsq12/<PROJECT>_lsq12_average.mnc"
nl_4_Avg = "~/<PROJECT>/nl/MNC/NL_4_average.mnc"
nl_4_Avg_LM = "~/<PROJECT>/nl/MNC/NL_4_average_landmarks.tag"

# Define blur files without "_blur" suffix.
LM_Avg_Mask_352 = "~/<PROJECT>/Source/MNC/LM_average_mask_352"
LM_Avg_Mask_176 = "~/<PROJECT>/Source/MNC/LM_average_mask_176"
LM_Avg_Mask_098 = "~/<PROJECT>/Source/MNC/LM_average_mask_098"
LM_Avg_Mask_078 = "~/<PROJECT>/Source/MNC/LM_average_mask_078"
LM_Avg_Mask_064 = "~/<PROJECT>/Source/MNC/LM_average_mask_064"
LM_Avg_Mask_050 = "~/<PROJECT>/Source/MNC/LM_average_mask_050"
lsq12_Avg_030 = "~/<PROJECT>/nl/INIT/<PROJECT>_lsq12_average_030"
nl_4_Avg_352 = "~/<PROJECT>/lsq6/Blurred/NL_4_average_352"
nl_4_Avg_176 = "~/<PROJECT>/lsq6/Blurred/NL_4_average_176"
nl_4_Avg_098 = "~/<PROJECT>/lsq12/Blurred/NL_4_average_098"
nl_4_Avg_078 = "~/<PROJECT>/lsq6/Blurred/NL_4_average_078"
nl_4_Avg_064 = "~/<PROJECT>/lsq12/Blurred/NL_4_average_064"
nl_4_Avg_050 = "~/<PROJECT>/lsq12/Blurred/NL_4_average_050"

# Define blur files with "_blur" suffix.
LM_Avg_Mask_352_Blur = "~/<PROJECT>/Source/MNC/LM_average_mask_352_blur.mnc"
LM_Avg_Mask_176_Blur = "~/<PROJECT>/Source/MNC/LM_average_mask_176_blur.mnc"
LM_Avg_Mask_098_Blur = "~/<PROJECT>/Source/MNC/LM_average_mask_098_blur.mnc"
LM_Avg_Mask_078_Blur = "~/<PROJECT>/Source/MNC/LM_average_mask_078_blur.mnc"
LM_Avg_Mask_064_Blur = "~/<PROJECT>/Source/MNC/LM_average_mask_064_blur.mnc"
LM_Avg_Mask_050_Blur = "~/<PROJECT>/Source/MNC/LM_average_mask_050_blur.mnc"
lsq12_Avg_030 = "~/<PROJECT>/nl/INIT/<PROJECT>_lsq12_average_030_blur.MNC"
nl_4_Avg_352_Blur = "~/<PROJECT>/lsq6/Blurred/NL_4_average_352_blur.mnc"
nl_4_Avg_176_Blur = "~/<PROJECT>/lsq6/Blurred/NL_4_average_176_blur.mnc"
nl_4_Avg_098_Blur = "~/<PROJECT>/lsq12/Blurred/NL_4_average_098_blur.mnc"
nl_4_Avg_098_Dxyz = "~/<PROJECT>/lsq12/Blurred/NL_4_average_098_dxyz.mnc"
nl_4_Avg_078_Blur = "~/<PROJECT>/lsq6/Blurred/NL_4_average_078_blur.mnc"
nl_4_Avg_064_Blur = "~/<PROJECT>/lsq12/Blurred/NL_4_average_064_blur.mnc"
nl_4_Avg_050_Blur = "~/<PROJECT>/lsq12/Blurred/NL_4_average_050_blur.mnc"

# Define blur files with "dxyz" suffix.
LM_Avg_Mask_352_Dxyz = "~/<PROJECT>/Source/MNC/LM_average_mask_352_dxyz.mnc"
LM_Avg_Mask_176_Dxyz = "~/<PROJECT>/Source/MNC/LM_average_mask_176_dxyz.mnc"
LM_Avg_Mask_098_Dxyz = "~/<PROJECT>/Source/MNC/LM_average_mask_098_dxyz.mnc"
LM_Avg_Mask_078_Dxyz = "~/<PROJECT>/Source/MNC/LM_average_mask_078_dxyz.mnc"
LM_Avg_Mask_064_Dxyz = "~/<PROJECT>/Source/MNC/LM_average_mask_064_dxyz.mnc"
LM_Avg_Mask_050_Dxyz = "~/<PROJECT>/Source/MNC/LM_average_mask_050_dxyz.mnc"
lsq12_Avg_030_Dxyz = "~/<PROJECT>/nl/INIT/<PROJECT>_lsq12_average_030_dxyz.mnc"
nl_4_Avg_352_Dxyz = "~/<PROJECT>/lsq6/Blurred/NL_4_average_352_dxyz.mnc"
nl_4_Avg_176_Dxyz = "~/<PROJECT>/lsq6/Blurred/NL_4_average_176_dxyz.mnc"
nl_4_Avg_098_Dxyz = "~/<PROJECT>/lsq12/Blurred/NL_4_average_098_dxyz.mnc"
nl_4_Avg_078_Dxyz = "~/<PROJECT>/lsq6/Blurred/NL_4_average_078_dxyz.mnc"
nl_4_Avg_064_Dxyz = "~/<PROJECT>/lsq12/Blurred/NL_4_average_064_dxyz.mnc"
nl_4_Avg_050_Dxyz = "~/<PROJECT>/lsq12/Blurred/NL_4_average_050_dxyz.mnc"

# Open and read ('r') the all specimen list.
Specimen_List=open(All_Specimens,'r')
# Read specimen list as a single string.
Specimens=Specimen_List.read()
# Split the specimens by line to obtain a vector of specimen IDs.
Specimen_IDs=Specimens.split('\n')
# Remove the last blank entry caused by the final \n.
del Specimen_IDs[-1]
# Define the length of the specimen list.
Specimen_List_Length=len(Specimen_IDs)
# Divide the length of the list by the same number.
Specimen_Group=int(Specimen_List_Length/Specimen_List_Length)
# Calculate upper limit of specimen list to use in for loop sequences.
Specimen_Upper=(Specimen_List_Length)-1
# Close the specimen list.
Specimen_List.close()

# Code the strings for iterative blurring via mincblur. -clobber overwrites existing files; -no apodize turns off apodization, which is designed to reduce diffraction edge effects (e.g., detector noise); -gradient calculates the change in intensity of a single pixel in the source image to the target image; -fwhm stands for full-width at half-maximum. In other words, what we want to do is convolve a "smoothing kernel", or Gaussian function, over an input volume in order to average neighboring points. The full-width at half-maximum describes the width of the Gaussian function at half of its peak to the left and right;
MNC_Blur = "mincblur -clobber -no_apodize -gradient -fwhm "

# Code the hierarchical registration strings which call upon the minctracc command. -clobber overwrites existing files; -xcorr stands for cross-correlation, which is the similarity metric to be optimized (maximized); -lsq6 indicates that we want perform a rigid body transformation with six degrees of freedom (translation (z,y,x) and rotation (z,y,x)); -lsq12 indicates that we want to perform a 12-parameter affine transformation with twelve degrees of freedom (translation (z,y,x), rotation (z,y,x), scale (z,y,x), and shear (z,y,x)); -w_translations/rotations/scales/shear optimization weights along z,y,x; -step is the z,y,x resolution; -simplex is the optimizer; -tol is the value at which the optimization stops.
lsq6_Register_352_Blur = "minctracc -clobber -xcorr -lsq6 -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.352 0.352 0.352 -simplex 0.78 -use_simplex -tol 0.0001 "
lsq6_Register_176_Blur = "minctracc -clobber -xcorr -lsq6 -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.176 0.176 0.176 -simplex 0.54 -use_simplex -tol 0.0001 "
lsq6_Register_078_Blur = "minctracc -clobber -xcorr -lsq6 -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.078 0.078 0.078 -simplex 0.32 -use_simplex -tol 0.0001 "
lsq12_Register_098_Blur = "minctracc -clobber -xcorr -lsq12 -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.098 0.098 0.098 -simplex 0.980 -use_simplex -tol 0.0001 "
lsq12_Register_064_Blur = "minctracc -clobber -xcorr -lsq12 -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.064 0.064 0.064 -simplex 0.490 -use_simplex -tol 0.0001 "
lsq12_Register_050_Blur = "minctracc -clobber -xcorr -lsq12 -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.050 0.050 0.050 -simplex 0.333 -use_simplex -tol 0.0001 "
nl_1_Register_1400_Blur_Begin = "minctracc -clobber -xcorr -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 1.4 1.4 1.4 -simplex 3.5 -use_simplex -tol 0.0001 "
nl_1_Register_1400_Blur_End = "-iterations 40 -similarity 0.8 -weight 0.8 -stiffness 0.98 -nonlinear corrcoeff -sub_lattice 6 -lattice_diameter 4.2 4.2 4.2 -max_def_magnitude 1 -debug -xcorr -identity "
nl_2_Register_1000_Blur_Begin = "minctracc -clobber -xcorr -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 1.0 1.0 1.0 -simplex 3.5 -use_simplex -tol 0.0001 "
nl_2_Register_1000_Blur_End = "-iterations 20 -similarity 0.8 -weight 0.8 -stiffness 0.98 -nonlinear corrcoeff -sub_lattice 6 -lattice_diameter 3.0 3.0 3.0 -max_def_magnitude 1 -debug -xcorr -transform "
nl_3_Register_700_Blur_Begin = "minctracc -clobber -xcorr -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.7 0.7 0.7 -simplex 3.5 -use_simplex -tol 0.0001 "
nl_3_Register_700_Blur_End = "-iterations 15 -similarity 0.8 -weight 0.8 -stiffness 0.98 -nonlinear corrcoeff -sub_lattice 6 -lattice_diameter 2.1 2.1 2.1 -max_def_magnitude 1 -debug -xcorr -transform "
nl_4_Register_400_Blur_Begin = "minctracc -clobber -xcorr -w_translations 0.4 0.4 0.4 -w_rotations 0.0174533 0.0174533 0.0174533 -w_scales 0.02 0.02 0.02 -w_shear 0.02 0.02 0.02 -step 0.4 0.4 0.4 -simplex 3.5 -use_simplex -tol 0.0001 "
nl_4_Register_400_Blur_End = "-iterations 15 -similarity 0.8 -weight 0.8 -stiffness 0.98 -nonlinear corrcoeff -sub_lattice 6 -lattice_diameter 1.2 1.2 1.2 -max_def_magnitude 1 -debug -xcorr -transform "

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Begin writing .sh scripts that will be submitted to the cluster.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Open a file to write to; 'a' for append; lsq6_First is the first 6-parameter .sh script to submit, because it blurs the intended target (i.e., the landmark initialized average) as well as a mask of equivalent resolution to constrain our computation.
Atlas_Blur = open("Atlas_Blur.sh",'a')
# Write standard .sh header. This header is needed for SLURM job submission.
Atlas_Blur.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=15000M\n#SBATCH --time=07:00:00\n\nmodule load minc-toolkit/2016-11\n\ncd " + Scripts_path + "\n\necho \"The job started at $(date).\"\n\n")
# Write commands to blur the atlas file and LM_average_mask with isotropic blurring kernels. These blur values should be decided upon with respect to the original resolution of the image. Here, we're assuming 35 micron resolution.
Atlas_Blur.write(MNC_Blur + "1.4 " + nl_4_Avg + " " + nl_4_Avg_1400 + "\n")
Atlas_Blur.write(MNC_Blur + "1.0 " + nl_4_Avg + " " + nl_4_Avg_1000 + "\n")
Atlas_Blur.write(MNC_Blur + "0.7 " + nl_4_Avg + " " + nl_4_Avg_700 + "\n")
Atlas_Blur.write(MNC_Blur + "0.4 " + nl_4_Avg + " " + nl_4_Avg_400 + "\n")
Atlas_Blur.write(MNC_Blur + "0.352 " + nl_4_Avg + " " + nl_4_Avg_352 + "\n")
Atlas_Blur.write(MNC_Blur + "0.176 " + nl_4_Avg + " " + nl_4_Avg_176 + "\n")
Atlas_Blur.write(MNC_Blur + "0.098 " + nl_4_Avg + " " + nl_4_Avg_098 + "\n")
Atlas_Blur.write(MNC_Blur + "0.078 " + nl_4_Avg + " " + nl_4_Avg_078 + "\n")
Atlas_Blur.write(MNC_Blur + "0.064 " + nl_4_Avg + " " + nl_4_Avg_064 + "\n")
Atlas_Blur.write(MNC_Blur + "0.050 " + nl_4_Avg + " " + nl_4_Avg_050 + "\n")
Atlas_Blur.write(MNC_Blur + "1.4 " + LM_Avg_Mask + " " + LM_Avg_Mask_1400 + "\n")
Atlas_Blur.write(MNC_Blur + "1.0 " + LM_Avg_Mask + " " + LM_Avg_Mask_1000 + "\n")
Atlas_Blur.write(MNC_Blur + "0.7 " + LM_Avg_Mask + " " + LM_Avg_Mask_700 + "\n")
Atlas_Blur.write(MNC_Blur + "0.4 " + LM_Avg_Mask + " " + LM_Avg_Mask_400 + "\n")
Atlas_Blur.write(MNC_Blur + "0.352 " + LM_Avg_Mask + " " + LM_Avg_Mask_352 + "\n")
Atlas_Blur.write(MNC_Blur + "0.176 " + LM_Avg_Mask + " " + LM_Avg_Mask_176 + "\n")
Atlas_Blur.write(MNC_Blur + "0.098 " + LM_Avg_Mask + " " + LM_Avg_Mask_098 + "\n")
Atlas_Blur.write(MNC_Blur + "0.078 " + LM_Avg_Mask + " " + LM_Avg_Mask_078 + "\n")
Atlas_Blur.write(MNC_Blur + "0.064 " + LM_Avg_Mask + " " + LM_Avg_Mask_064 + "\n")
Atlas_Blur.write(MNC_Blur + "0.050 " + LM_Avg_Mask + " " + LM_Avg_Mask_050 + "\n")
Atlas_Blur.write("echo \"The job ended at $(date).\"")
# Close the file.
Atlas_Blur.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6-parameter (translation (z,y,x), rotation (z,y,x)) optimal rigid body registration stage.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Open a file to write to; 'a' for append; lsq6_Test is the only 6-parameter stage.
# Begin a for loop to blur all landmark initialized source files, register them to , and resample each image into the new translation and rotation invariant space. Note that we ideally want to resample the LM initialized images, rather than the original .mnc images. Hence, the LM .xfm is not included in the concatenation.
# Begin a counter.
lsq6_Counter=0
for SpecID in Specimen_IDs:
# Add 1 to the counter for every new file in the loop.
	lsq6_Counter += 1
	# Open a file to write to; 'a' for append.
	lsq6_Test = open("lsq6_Test_" + str(lsq6_Counter) + ".sh",'a')
	# Add header.
	lsq6_Test.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=15000M\n#SBATCH --time=07:00:00\n\nmodule load minc-toolkit/2016-11\n\ncd " + Scripts_path + "\n\necho \"The job started at $(date).\"\n\n")
	lsq6_Test.write("echo \"Begin the optimized 6-parameter registration for " + SpecID + ".\"\n\n")
	# Blur each image, as done with the average and average mask, with a 0.352, 0.176, and 0.078 isotropic Gaussian blurring kernel. These blur values should be decided upon with respect to the original resolution of the image.
	lsq6_Test.write(MNC_Blur + "0.352 " + Source_MNC_path + SpecID + ".mnc " + lsq6_Blurred_path + SpecID + "_352\n")
	lsq6_Test.write(MNC_Blur + "0.176 " + Source_MNC_path + SpecID + ".mnc " + lsq6_Blurred_path + SpecID + "_176\n")
	lsq6_Test.write(MNC_Blur + "0.078 " + Source_MNC_path + SpecID + ".mnc " + lsq6_Blurred_path + SpecID + "_078\n")
	# Call registration strings. We begin with the most blurred (e.g., 0.352) image. -model_mask specifies the mask we wish to use. Note that we use a mask with the same amount of blurring; -identity specifies an identity matrix that initializes the transformation matrix. Upon specifying a transformation matrix, we extract the relevant transformation parameters (here, rotation and translation) and optimize them to find the best transformation; -transformation specifies a file giving a previous source to target mapping, which is used as the new coordinate starting point for the optimization.
	lsq6_Test.write(lsq6_Register_352_Blur + lsq6_Blurred_path + SpecID + "_352_blur.mnc " + nl_4_Avg_352_Blur + " " + lsq6_XFM_path + SpecID + "_lsq6_0.xfm -model_mask " + LM_Avg_Mask_352_Blur + " -identity\n")
	lsq6_Test.write(lsq6_Register_176_Blur + lsq6_Blurred_path + SpecID + "_176_blur.mnc " + nl_4_Avg_176_Blur + " " + lsq6_XFM_path + SpecID + "_lsq6_1.xfm -model_mask " + LM_Avg_Mask_176_Blur + " -transformation " + lsq6_XFM_path + SpecID + "_lsq6_0.xfm\n")
	lsq6_Test.write(lsq6_Register_078_Blur + lsq6_Blurred_path + SpecID + "_078_blur.mnc " + nl_4_Avg_078_Blur + " " + lsq6_XFM_path + SpecID + "_lsq6_2.xfm -model_mask " + LM_Avg_Mask_078_Blur + " -transformation " + lsq6_XFM_path + SpecID + "_lsq6_1.xfm\n")
	# Resample the original image into the rotation and translation invariant space using the concatenated transformation. Be mindful of your "original" images and their naming convention.
	lsq6_Test.write("mincresample -like " + nl_4_Avg + " -clobber -transformation " + lsq6_XFM_path + SpecID + "_lsq6_2.xfm " + Source_MNC_path + SpecID + ".mnc " + lsq6_MNC_path + SpecID + "_lsq6.mnc\n\n")
	lsq6_Test.write("echo \"The job ended at $(date).\"")
	# Close the 6-parameter blur/register file.
	lsq6_Test.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 12-parameter (translation (z,y,x), rotation (z,y,x), scale (z,y,x), shear (z,y,x)) optimized affine registration stage.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Open a file to write to; 'a' for append; lsq12_Test is the only 12-parameter stage.
# Begin a for loop to blur all landmark initialized source files, register them to LM_Average, and resample each image into the new translation and rotation invariant space.
# Note that we ideally want to resample the LM initialized images, rather than the original .mnc images. Hence, the LM .xfm is not included in the concatenation.
# Begin a counter.
lsq12_Counter=0
for SpecID in Specimen_IDs:
# Add 1 to the counter for every new file in the loop.
	lsq12_Counter += 1
	# Open a file to write to; 'a' for append.
	lsq12_Test = open("lsq12_Test_" + str(lsq12_Counter) + ".sh",'a')
	# Add header.
	lsq12_Test.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=15000M\n#SBATCH --time=07:00:00\n\nmodule load minc-toolkit/2016-11\n\ncd " + Scripts_path + "\n\necho \"The job started at $(date).\"\n\n")
	lsq12_Test.write("echo \"Begin the optimized 12-parameter registration for " + SpecID + ".\"\n\n")
	# Blur each image, as done with the average and average mask, with a 0.098, 0.064, and 0.050 isotropic Gaussian blurring kernel. These blur values should be decided upon with respect to the original resolution of the image.
	lsq12_Test.write(MNC_Blur + "0.098 " + lsq6_MNC_path + SpecID + "_lsq6.mnc " + lsq12_Blurred_path + SpecID + "_098\n")
	lsq12_Test.write(MNC_Blur + "0.064 " + lsq6_MNC_path + SpecID + "_lsq6.mnc " + lsq12_Blurred_path + SpecID + "_064\n")
	lsq12_Test.write(MNC_Blur + "0.050 " + lsq6_MNC_path + SpecID + "_lsq6.mnc " + lsq12_Blurred_path + SpecID + "_050\n")
	# Call registration strings. We begin with the most blurred (e.g., 0.098) image. -model_mask specifies the mask we wish to use. Note that we use a mask with the same amount of blurring; -identity specifies an identity matrix that initializes the transformation matrix. Upon specifying a transformation matrix, we extract the relevant transformation parameters (here, rotation and translation) and optimize them to find the best transformation; -transformation specifies a file giving a previous source to target mapping, which is used as the new coordinate starting point for the optimization.
	lsq12_Test.write(lsq12_Register_098_Blur + lsq12_Blurred_path + SpecID + "_098_blur.mnc " + nl_4_Avg_098_Blur + " " + lsq12_XFM_path + SpecID + "_lsq12_0.xfm -model_mask " + LM_Avg_Mask_098_Blur + " -identity\n")
	lsq12_Test.write(lsq12_Register_064_Blur + lsq12_Blurred_path + SpecID + "_064_blur.mnc " + nl_4_Avg_064_Blur + " " + lsq12_XFM_path + SpecID + "_lsq12_1.xfm -model_mask " + LM_Avg_Mask_064_Blur + " -transformation " + lsq12_XFM_path + SpecID + "_lsq12_0.xfm\n")
	lsq12_Test.write(lsq12_Register_050_Blur + lsq12_Blurred_path + SpecID + "_050_blur.mnc " + nl_4_Avg_050_Blur + " " + lsq12_XFM_path + SpecID + "_lsq12_2.xfm -model_mask " + LM_Avg_Mask_050_Blur + " -transformation " + lsq12_XFM_path + SpecID + "_lsq12_1.xfm\n")
	# Concatenate .xfm files.
	lsq12_Test.write("xfmconcat -clobber " + lsq6_XFM_path + SpecID + "_lsq6_2.xfm " + lsq12_XFM_path + SpecID + "_lsq12_2.xfm " + lsq12_XFM_path + SpecID + "_origtolsq12.xfm\n")
	# Resample the original image into the rotation and translation invariant space using the concatenated transformation. Be mindful of your "original" images and their naming convention.
	lsq12_Test.write("mincresample -like " + nl_4_Avg + " -clobber -transformation " + lsq12_XFM_path + SpecID + "_origtolsq12.xfm " + Source_MNC_path + SpecID + ".mnc " + lsq12_MNC_path + SpecID + "_lsq12.mnc\n\n")
	lsq12_Test.write("echo \"The job ended at $(date).\"")
	# Close the 12-parameter blur/register file.
	lsq12_Test.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# First non-linear (most blurred) ANIMAL registration stage.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create .sh files for the first round of non-linear commands (nl_1).
for Element in list(range(Specimen_List_Length)):
	# Open a file to write to; 'a' for append.
	nl_Test_First = open("nl_Test_First_" + str(Element) + ".sh",'a')
	# Add header.
	nl_Test_First.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=15000M\n#SBATCH --time=11:00:00\n\nmodule load minc-toolkit/2016-11\n\ncd " + Scripts_path + "\n\necho \"The job started at $(date).\"\n\n")
	for SpecID in Specimen_IDs[(Specimen_Group*Element):(Specimen_Group*(Element+1))]:
		# Blur each image with a 1.4 isotropic Gaussian blurring kernel.
		nl_Test_First.write(MNC_Blur + "1.4 " + lsq12_MNC_path + SpecID + "_lsq12.mnc " + nl_Blurred_path + SpecID + "_1400\n")
		# Call the ANIMAL registration string. Details above.
		nl_Test_First.write(nl_1_Register_1400_Blur_Begin + nl_Blurred_path + SpecID + "_1400_blur.mnc " + nl_4_Avg_1400_Blur + " " + nl_XFM_path + SpecID + "_nl_1.xfm -model_mask " + LM_Avg_Mask_1400_Blur + " " + nl_1_Register_1400_Blur_End + "\n")
		# Concatenate .xfm files.
		nl_Test_First.write("xfmconcat -clobber " + lsq6_XFM_path + SpecID + "_lsq6_2.xfm " + lsq12_XFM_path + SpecID + "_lsq12_2.xfm " + nl_XFM_path + SpecID + "_nl_1.xfm " + nl_XFM_path + SpecID + "_origtonl_1.xfm\n")
		# Resample source image into the first on-linear space with concatenated .xfm.
		nl_Test_First.write("mincresample -like " + nl_4_Avg + " -clobber -transformation " + nl_XFM_path + SpecID + "_origtonl_1.xfm " + Source_MNC_path + SpecID + ".mnc "+ nl_MNC_path + SpecID + "_nl_1.mnc\n\n")
	nl_Test_First.write("echo \"The job ended at $(date).\"")
	# Close the file.
	nl_Test_First.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Second non-linear (second most blurred) ANIMAL registration stage.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create .sh files for the second round of non-linear commands (nl_2).
for Element in list(range(Specimen_List_Length)):
	# Open a file to write to; 'a' for append.
	nl_Test_Second = open("nl_Test_Second_" + str(Element) + ".sh",'a')
	# Add header.
	nl_Test_Second.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=15000M\n#SBATCH --time=11:00:00\n\nmodule load minc-toolkit/2016-11\n\ncd " + Scripts_path + "\n\necho \"The job started at $(date).\"\n\n")
	for SpecID in Specimen_IDs[(Specimen_Group*Element):(Specimen_Group*(Element+1))]:
		# Blur each image with a 1.0 isotropic Gaussian blurring kernel.
		nl_Test_Second.write(MNC_Blur + "1.0 " + lsq12_MNC_path + SpecID + "_lsq12.mnc " + nl_Blurred_path + SpecID + "_1000\n")
		# Call the ANIMAL registration string. Details above.
		nl_Test_Second.write(nl_2_Register_1000_Blur_Begin + nl_Blurred_path + SpecID + "_1000_blur.mnc " + nl_4_Avg_1000_Blur + " " + nl_XFM_path + SpecID + "_nl_2.xfm -model_mask " + LM_Avg_Mask_1000_Blur + " " + nl_2_Register_1000_Blur_End + nl_XFM_path + SpecID + "_nl_1.xfm\n")
		# Concatenate .xfm files.
		nl_Test_Second.write("xfmconcat -clobber " + lsq6_XFM_path + SpecID + "_lsq6_2.xfm " + lsq12_XFM_path + SpecID + "_lsq12_2.xfm " + nl_XFM_path + SpecID + "_nl_2.xfm " + nl_XFM_path + SpecID + "_origtonl_2.xfm\n")
		# Resample source image into the second non-linear space with concatenated .xfm.
		nl_Test_Second.write("mincresample -like " + nl_4_Avg + " -clobber -transformation " + nl_XFM_path + SpecID + "_origtonl_2.xfm " + Source_MNC_path + SpecID + ".mnc "+ nl_MNC_path + SpecID + "_nl_2.mnc\n\n")
	nl_Test_Second.write("echo \"The job ended at $(date).\"")
	# Close the file.
	nl_Test_Second.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Third non-linear (third most blurred) ANIMAL registration stage.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create .sh files for the third round of non-linear commands (nl_3).
for Element in list(range(Specimen_List_Length)):
	# Open a file to write to; 'a' for append.
	nl_Test_Third = open("nl_Test_Third_" + str(Element) + ".sh",'a')
	# Add header.
	nl_Test_Third.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=15000M\n#SBATCH --time=11:00:00\n\nmodule load minc-toolkit/2016-11\n\ncd " + Scripts_path + "\n\necho \"The job started at $(date).\"\n\n")
	for SpecID in Specimen_IDs[(Specimen_Group*Element):(Specimen_Group*(Element+1))]:
		# Blur each image with a 0.7 isotropic Gaussian blurring kernel.
		nl_Test_Third.write(MNC_Blur + "0.7 " + lsq12_MNC_path + SpecID + "_lsq12.mnc " + nl_Blurred_path + SpecID + "_700\n")
		# Call the ANIMAL registration string. Details above.
		nl_Test_Third.write(nl_3_Register_700_Blur_Begin + nl_Blurred_path + SpecID + "_700_blur.mnc " + nl_4_Avg_700_Blur + " " + nl_XFM_path + SpecID + "_nl_3.xfm -model_mask " + LM_Avg_Mask_700_Blur + " " + nl_3_Register_700_Blur_End + nl_XFM_path + SpecID + "_nl_2.xfm\n")
		# Concatenate .xfm files.
		nl_Test_Third.write("xfmconcat -clobber " + lsq6_XFM_path + SpecID + "_lsq6_2.xfm " + lsq12_XFM_path + SpecID + "_lsq12_2.xfm " + nl_XFM_path + SpecID + "_nl_3.xfm " + nl_XFM_path + SpecID + "_origtonl_3.xfm\n")
		# Resample source image into third non-linear space with concatenated .xfm.
		nl_Test_Third.write("mincresample -like " + nl_4_Avg + " -clobber -transformation " + nl_XFM_path + SpecID + "_origtonl_3.xfm " + Source_MNC_path + SpecID + ".mnc "+ nl_MNC_path + SpecID + "_nl_3.mnc\n\n")
	nl_Test_Third.write("echo \"The job ended at $(date).\"")
	# Close the script.
	nl_Test_Third.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fourth non-linear (fourth most blurred) ANIMAL registration stage.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create .sh files for the fourth round of non-linear commands (nl_4).
for Element in list(range(Specimen_List_Length)):
	# Open a file to write to; 'a' for append.
	nl_Test_Fourth = open("nl_Test_Fourth_" + str(Element) + ".sh",'a')
	# Add header.
	nl_Test_Fourth.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=15000M\n#SBATCH --time=11:00:00\n\nmodule load minc-toolkit/2016-11\n\ncd " + Scripts_path + "\n\necho \"The job started at $(date).\"\n\n")
	for SpecID in Specimen_IDs[(Specimen_Group*Element):(Specimen_Group*(Element+1))]:
		# Blur each image with a 0.4 isotropic Gaussian blurring kernel.
		nl_Test_Fourth.write(MNC_Blur + "0.4 " + lsq12_MNC_path + SpecID + "_lsq12.mnc " + nl_Blurred_path + SpecID + "_400\n")
		# Call the ANIMAL registration string. Details above.
		nl_Test_Fourth.write(nl_4_Register_400_Blur_Begin + nl_Blurred_path + SpecID + "_400_blur.mnc " + nl_4_Avg_400_Blur + " " + nl_XFM_path + SpecID + "_nl_4.xfm -model_mask " + LM_Avg_Mask_400_Blur + " " + nl_4_Register_400_Blur_End + nl_XFM_path + SpecID + "_nl_3.xfm\n")
		# Concatenate .xfm files.
		nl_Test_Fourth.write("xfmconcat -clobber " + lsq6_XFM_path + SpecID + "_lsq6_2.xfm " + lsq12_XFM_path + SpecID + "_lsq12_2.xfm " + nl_XFM_path + SpecID + "_nl_4.xfm " + nl_XFM_path + SpecID + "_origtonl_4.xfm\n")
		# Resample source image into fourth non-linear space with concatenated .xfm.
		nl_Test_Fourth.write("mincresample -like " + nl_4_Avg + " -clobber -transformation " + nl_XFM_path + SpecID + "_origtonl_4.xfm " + Source_MNC_path + SpecID + ".mnc "+ nl_MNC_path + SpecID + "_nl_4.mnc\n\n")
		nl_Test_Fourth.write("xfminvert -clobber " + nl_XFM_path + SpecID + "_origtonl_4.xfm " + nl_XFM_path + SpecID + "_origtonl_4_inverted.xfm\n")
		# To propagate the landmarks, we invert the entire transformation, then use "transformtags" to propagate the atlas landmarks along this path. The resulting landmarks are in "/home/$USER/<PROJECT>/Source/Tag/".
		nl_Test_Fourth.write("transformtags -vol1 -transformation " + nl_XFM_path + SpecID + "_origtonl_4_inverted.xfm " + nl_4_Avg_LM + " " + Source_Tag_path + SpecID + "_landmarks.tag\n")
	nl_Test_Fourth.write("echo \"The job ended at $(date).\"")
	# Close the script.
	nl_Test_Fourth.close()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create master job submission scripts for all stages. These scripts will automatically submit all .sh scripts to the queue and will be chained together. In other words, only the first script, Job_Submission_First.sh, needs to be submitted
# via $ sbatch Job_Submission_First.sh.
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Job_Submission_First (i.e., 6-parameter registrations).
Job_Submission_First = open("Job_Submission_First.sh",'a')
Job_Submission_First.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=2000M\n#SBATCH --time=05-00:00:00\n#SBATCH --job-name=Job_Submission_First.sh\n\necho \"The job started at $(date).\"\n\n")
Job_Submission_First.write("cd " + Scripts_path + "\n\n")
Job_Submission_First.write("sleep 10\n\n")
Job_Submission_First.write("Atlas_Blur=$(sbatch Atlas_Blur.sh)\n\n")
Job_Submission_First.write("until [[ $(squeue -t CD -u $USER --noheader -j ${Atlas_Blur##* } | wc -l) -eq 1 ]]; do\nsleep 5\ndone\n\n")
Job_Submission_First.write("# Generate a sequence of numbers.\nNUMBERS=$(seq 1 " + str(Specimen_List_Length) + ")\n\n")
Job_Submission_First.write("# For loop to automatically submit your jobs.\nfor NUM in $NUMBERS; do\nNAME=\"lsq6_Test_$NUM.sh\"\nLSQ6=\"sbatch $NAME\"\n$LSQ6\necho $LSQ6\n# Sleep script in 3 second intervals.\nsleep 3\ndone\n\n")
Job_Submission_First.write("# Sleep script for 5 seconds before while loop.\nsleep 5\n\n")
Job_Submission_First.write("# Create while loop to idle submission script.\nwhile [[ $(squeue -t R -u $USER --noheader | wc -l) -gt 1 ]]; do\nsleep 5\ndone\n\n")
Job_Submission_First.write("#----------------------------------------------------- Submit next job submission script and remove current submission script from the queue.\n\n")
Job_Submission_First.write("sbatch Job_Submission_Second.sh\n\n")
Job_Submission_First.write("echo \"The job ended at $(date).\"\n")
Job_Submission_First.write("scancel -u $USER --jobname=Job_Submission_First.sh")
# Job_Submission_Second (i.e., 12-parameter registrations).
Job_Submission_Second = open("Job_Submission_Second.sh",'a')
Job_Submission_Second.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=2000M\n#SBATCH --time=05-00:00:00\n#SBATCH --job-name=Job_Submission_Second.sh\n\necho \"The job started at $(date).\"\n\n")
Job_Submission_Second.write("cd " + Scripts_path + "\n\n")
Job_Submission_Second.write("sleep 10\n\n")
Job_Submission_Second.write("# Generate a sequence of numbers.\nNUMBERS=$(seq 1 " + str(Specimen_List_Length) + ")\n\n")
Job_Submission_Second.write("# For loop to automatically submit your jobs.\nfor NUM in $NUMBERS; do\nNAME=\"lsq12_Test_$NUM.sh\"\nLSQ12=\"sbatch $NAME\"\n$LSQ12\necho $LSQ12\n# Sleep script in 3 second intervals.\nsleep 3\ndone\n\n")
Job_Submission_Second.write("# Sleep script for 5 seconds before while loop.\nsleep 5\n\n")
Job_Submission_Second.write("# Create while loop to idle submission script.\nwhile [[ $(squeue -t R -u $USER --noheader | wc -l) -gt 1 ]]; do\nsleep 5\ndone\n\n")
Job_Submission_Second.write("#----------------------------------------------------- Submit next job submission script and remove current submission script from the queue.\n\n")
Job_Submission_Second.write("Job_Submission_Third=$(sbatch Job_Submission_Third.sh)\n\n")
Job_Submission_Second.write("echo \"The job ended at $(date).\"\n")
Job_Submission_Second.write("scancel -u $USER --jobname=Job_Submission_Second.sh")
# Job_Submission_Third (i.e., nl_First non-linear registrations).
Job_Submission_Third = open("Job_Submission_Third.sh",'a')
Job_Submission_Third.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=2000M\n#SBATCH --time=05-00:00:00\n#SBATCH --job-name=Job_Submission_Third.sh\n\necho \"The job started at $(date).\"\n\n")
Job_Submission_Third.write("cd " + Scripts_path + "\n\n")
Job_Submission_Third.write("sleep 10\n\n")
Job_Submission_Third.write("# Generate a sequence of numbers.\nNUMBERS=$(seq 0 " + str(Specimen_Upper) + ")\n\n")
Job_Submission_Third.write("# For loop to automatically submit your jobs.\nfor NUM in $NUMBERS; do\nNAME=\"nl_Test_First_$NUM.sh\"\nnl_First=\"sbatch $NAME\"\n$nl_First\necho $nl_First\n# Sleep script in 3 second intervals.\nsleep 3\ndone\n\n")
Job_Submission_Third.write("# Sleep script for 5 seconds before while loop.\nsleep 5\n\n")
Job_Submission_Third.write("# Create while loop to idle submission script.\nwhile [[ $(squeue -t R -u $USER --noheader | wc -l) -gt 1 ]]; do\nsleep 5\ndone\n\n")
Job_Submission_Third.write("#----------------------------------------------------- Submit next job submission script and remove current submission script from the queue.\n\n")
Job_Submission_Third.write("Job_Submission_Fourth=$(sbatch Job_Submission_Fourth.sh)\n\n")
Job_Submission_Third.write("echo \"The job ended at $(date).\"\n")
Job_Submission_Third.write("scancel -u $USER --jobname=Job_Submission_Third.sh")
# Job_Submission_Fourth (i.e., nl_Second non-linear registrations).
Job_Submission_Fourth = open("Job_Submission_Fourth.sh",'a')
Job_Submission_Fourth.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=2000M\n#SBATCH --time=05-00:00:00\n#SBATCH --job-name=Job_Submission_Fourth.sh\n\necho \"The job started at $(date).\"\n\n")
Job_Submission_Fourth.write("cd " + Scripts_path + "\n\n")
Job_Submission_Fourth.write("sleep 10\n\n")
Job_Submission_Fourth.write("# Generate a sequence of numbers.\nNUMBERS=$(seq 0 " + str(Specimen_Upper) + ")\n\n")
Job_Submission_Fourth.write("# For loop to automatically submit your jobs.\nfor NUM in $NUMBERS; do\nNAME=\"nl_Test_Second_$NUM.sh\"\nnl_Second=\"sbatch $NAME\"\n$nl_Second\necho $nl_Second\n# Sleep script in 3 second intervals.\nsleep 3\ndone\n\n")
Job_Submission_Fourth.write("# Sleep script for 5 seconds before while loop.\nsleep 5\n\n")
Job_Submission_Fourth.write("# Create while loop to idle submission script.\nwhile [[ $(squeue -t R -u $USER --noheader | wc -l) -gt 1 ]]; do\nsleep 5\ndone\n\n")
Job_Submission_Fourth.write("#----------------------------------------------------- Submit next job submission script and remove current submission script from the queue.\n\n")
Job_Submission_Fourth.write("Job_Submission_Fifth=$(sbatch Job_Submission_Fifth.sh)\n\n")
Job_Submission_Fourth.write("echo \"The job ended at $(date).\"\n")
Job_Submission_Fourth.write("scancel -u $USER --jobname=Job_Submission_Fourth.sh")
# Job_Submission_Fifth (i.e., nl_Third non-linear registrations).
Job_Submission_Fifth = open("Job_Submission_Fifth.sh",'a')
Job_Submission_Fifth.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=2000M\n#SBATCH --time=05-00:00:00\n#SBATCH --job-name=Job_Submission_Fifth.sh\n\necho \"The job started at $(date).\"\n\n")
Job_Submission_Fifth.write("cd " + Scripts_path + "\n\n")
Job_Submission_Fifth.write("sleep 10\n\n")
Job_Submission_Fifth.write("# Generate a sequence of numbers.\nNUMBERS=$(seq 0 " + str(Specimen_Upper) + ")\n\n")
Job_Submission_Fifth.write("# For loop to automatically submit your jobs.\nfor NUM in $NUMBERS; do\nNAME=\"nl_Test_Third_$NUM.sh\"\nnl_Third=\"sbatch $NAME\"\n$nl_Third\necho $nl_Third\n# Sleep script in 3 second intervals.\nsleep 3\ndone\n\n")
Job_Submission_Fifth.write("# Sleep script for 5 seconds before while loop.\nsleep 5\n\n")
Job_Submission_Fifth.write("# Create while loop to idle submission script.\nwhile [[ $(squeue -t R -u $USER --noheader | wc -l) -gt 1 ]]; do\nsleep 5\ndone\n\n")
Job_Submission_Fifth.write("#----------------------------------------------------- Submit next job submission script and remove current submission script from the queue.\n\n")
Job_Submission_Fifth.write("Job_Submission_Sixth=$(sbatch Job_Submission_Sixth.sh)\n\n")
Job_Submission_Fifth.write("echo \"The job ended at $(date).\"\n")
Job_Submission_Fifth.write("scancel -u $USER --jobname=Job_Submission_Fifth.sh")
# Job_Submission_Sixth (i.e., nl_Fourth non-linear registrations).
Job_Submission_Sixth = open("Job_Submission_Sixth.sh",'a')
Job_Submission_Sixth.write("#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --mem=2000M\n#SBATCH --time=05-00:00:00\n#SBATCH --job-name=Job_Submission_Sixth.sh\n\necho \"The job started at $(date).\"\n\n")
Job_Submission_Sixth.write("cd " + Scripts_path + "\n\n")
Job_Submission_Sixth.write("sleep 10\n\n")
Job_Submission_Sixth.write("# Generate a sequence of numbers.\nNUMBERS=$(seq 0 " + str(Specimen_Upper) + ")\n\n")
Job_Submission_Sixth.write("# For loop to automatically submit your jobs.\nfor NUM in $NUMBERS; do\nNAME=\"nl_Test_Fourth_$NUM.sh\"\nnl_Fourth=\"sbatch $NAME\"\n$nl_Fourth\necho $nl_Fourth\n# Sleep script in 3 second intervals.\nsleep 3\ndone\n\n")
Job_Submission_Sixth.write("# Sleep script for 5 seconds before while loop.\nsleep 5\n\n")
Job_Submission_Sixth.write("# Create while loop to idle submission script.\nwhile [[ $(squeue -t R -u $USER --noheader | wc -l) -gt 1 ]]; do\nsleep 5\ndone\n\n")
Job_Submission_Sixth.write("#----------------------------------------------------- Remove current submission script from the queue.\n\n")
Job_Submission_Sixth.write("echo \"The job ended at $(date).\"\n")
Job_Submission_Sixth.write("scancel -u $USER --jobname=Job_Submission_Sixth.sh")
