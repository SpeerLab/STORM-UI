Instructions for batch automated STORM analysis using Python3 sCMOS fitting code. 

# Colenso Speer 05_30_18

This folder contains scripts for running STORM analysis on data collected from Sperry using 640 x 640 parameter files.

Changes to data acquisition format will require:
A) modification of data path as appropriate
B) modification of queried image frame in bead, conv, ffc movies as needed
C) modification of image size parameters depending on acquisition
D) use of resliced calibration file for appropriate frame size

For data acquired using Dave and 640 x 640 frame size parameters, the following scripts should be run in order:

0) "0_fitting_parameters_evaluation":
	
	This script is the sCMOS portion of the following script (1). Running this will allow the user to quickly test fitting 
	parameters on STORM movies prior to running full analysis of the STORM files. To accomplish this, the user should
	change the appropriate XML files in "/Batch_analysis/XMLs/" to analyze only a small portion of the STORM movies for testing.
	Fitting parameters (sigma, threshold, etc.) can be adjusted as necessary. Be sure to change the XML file back to analyze
	the entire STORM movie once appropriate parameters have been determined. 
	
	# XMLs folder need to be copied under Z://xxxx// . 
	# The difference between 0_fitting_parameters_evaluation and 1_scmos_batch_analysis is whether beads are processed. 
	# After 0_fitting_parameters_evaluation, 2_master_python_analysis can be used to render the high-resolution .tiff images to test whether the parameters used are suitable. 

1) "1_scmos_batch_analysis":
	
	This script will run multifit scmos analysis on all STORM and bead movies. It will then generate warping transforms
	for bead movies and test the transforms generated from the second bead set on the first bead set. 

	The path to the expfolder should be changed for your experiment. All data files should be stored in a single folder
	entitled "Acquisition". The # of processes to run concurrently (max_processes) should be changed to suit machine 
	usage. The XML files for fitting and the calib.npy file for camera calibration should be copied to the expfolder.
	This file can be run by the .bat push file. Fitting should be tested on small portions of your data (see previous
	script) to confirm the best parameters prior to running batch analysis.


	This script will also convert all STORM .hdf5 fitting files into .TIF image files. 


3) 
4) "4_XY_alignment":

	This script builds folder structure for aligning images for each physical section and for saving out all processed
	and merged images. 


	This script will then run a subprocess calling the bead_processing.py code. This will load the FFC .dax files, 
	identify the appropriate frames of interest, stretch the histograms for
	each FFC image to normalize them across the entire set and save out the images. This script also saves the final
	FFC images for use in generating the lens distortion correction. The lens distortion
	correction using a 3x3 array of dense bead field FFC images is currently not being performed. 

	The frames and image size will need to be changed if acquisition parameters differ from 640 x 640 standard.
	If running the lens distortion correction, FIJI will throw a dialogue to adjust the lambda parameter (see lens distortion
	documentation for details on adjusting this)- press "cancel" to accept the lambda term and continue with the 
	distortion correction. (http://www.kaynig.de/downloads/DistortionCorrectionPlugin_Manual.pdf)

	The lens distortion correction Plugin is currently not reconfigured to avoid throwing dialogue. This makes it impossible
	to use when applying the distortion correction field to all subsequent conventional and STORM images. Going forward,
	this part of the error correction should either be fixed, or abandoned entirely. 
	
	# Fiji in accound 'Colenso' is used in this script. If FIji doesn't run, open C:/Users/Colenso/ in the file explorer 
	# manually to obtain the authority. 


	The code then finds the conventional images from .dax files, stretches the histograms, applies the FFC
	based on relative intensity values, and saves the images out. The code adjusts the contrast of the STORM images and saves
	them out to the "rawimages" folder for further processing. The code then runs Matlab to apply the polynomical warping
	transform to all conventional and STORM images. Lastly, the code runs an instance of FIJI to apply the lens distortion
	correction and save out merged and downsampled conventional and STORM images. 

	The path to the expfolder should be changed for your experiment. All data files should be stored in a single folder
	entitled "Acquisition". The # of processes to run concurrently (max_processes) should be changed to suit machine 
	usage. The frames queried and image size will need to be changed if acquisition parameters differ from 640 x 640 
	standard settings. 

	# The rendered channels might requires modified for different experiments. 
	
	