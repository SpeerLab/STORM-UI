# STORM_GUI  -- "STORMpro.py"


This is the Speer Lab STORM Analysis Application called "STORMpro" (beta). 

STORMpro is described in the following publication: Vatan et al., "Volumetric super-resolution imaging by serial ultrasectioning and STochastic Optical Reconstruction Microscopy (STORM) in neural tissue" Star Protocols 2021

Users should download and install FIJI (https://imagej.net/software/fiji/downloads).

Users should download and install XML files and macros folder containing Fiji Z alignment .py code for 3D serial section alignment. 

Macros from the STORM-UI repository should be added to the 'FIJI\macros' directory on the user's analysis machine. 

## How to Run? 

* Run `main.py` through command `python main.py`
* This will load "STORMpro", the GUI 


## Functionality

* This UI will implement a pipeline for processing STORM data from raw .dax files to 3D aligned .tiff image stacks. Important steps in the process are:
  1.) Test/visualize fitting parameters on small portion of STORM data prior to running full analysis 
  2.) 3D-DAOSTORM single molecule fitting on all STORM and bead movies 
  3.) XY alignment with drift/chromatic aberration correction for images in physical sections 
  4.) Z stack alignment by rigid and elastic registration of all serial sections

## Instructions 

### Visualizing DAX Files 

Load any chosen .dax STORM file into the visualizer. This enables interactive visualization of the raw STORM data using keyboard entries and mouse control to scroll through the movie. 

To load a DAX file, you can either drag and drop a file to the movie-viewer or simply load one using the 
`Load Movie` button at the top left. 

To move through the frames: 

* Press '.' to step forward one frame 
* Press ',' to step backward one frame 
* Press 'l' to step forward 200 frames 
* Press 'k' to step backward 200 frames 
* Press 'Home' to go to beginning of movie 
* Press 'End' to go to end of movie 

To change the contrast, adjust the vertical contrast bar to desired value. 

Additionally, you can load localizations which can overlay the movie frame.

To get the sub window to rescale as you resize the overall window, please maximize the sub window. 

### Evaluating Fitting Parameters 

First select the desired channels for analysis through the check boxes next to "Storm Channels". These selections reflect the specific wavelengths collected in the image acquisition. To adjust single-molecule fitting parameters, the user will modify the downloaded XML files and change settings within these files to specify SMLM parameters. Click on "Set Test Parameters" to open a window for multi-channel XML parameter adjustment. For each desired channel, adjust the input properties (e.g. frame range, background, etc.). Once properties are set, press "Update XMLs" to overwrite the XML file with the updated settings.

A description of 3D-DAOSTORM fitting parameters can be found here:
https://storm-analysis.readthedocs.io/en/latest/parameters.html

Set “Movie for fitting evaluation” to desired .dax file for analysis. This can be any movie that the user chooses. The goal is to fit a small portion of the movie to test the performance of the selected fit parameters for SMLM. Adjust start/end frames in the "Set test parameters" window to analyze a small sample portion of the selected STORM movie. Typically, 10-20 frames are selected in the middle or near the end of the .dax STORM movie. Select "Update XMLs" to update the desired test fitting range.

Select “Fit Parameters" to run test fits on the chosen .dax movie. 

Load the output .hdf5 single-molecule localization list by dragging/dropping into the visualizer or selecting “Load localization 1” to select the desired .hdf5 file. 

Note: to visualize localizations in correct frames, temporarily set the radius for matching peaks from frame to frame in .xml to 0. Or else all localizations will show up in the first frame. 

Inspect the fit quality by visualizing the fit beacons (green circles in Figure 9A, D) to determine if single-molecules are over-fit (fits with no single-molecule PSFs present) or under-fit (raw PSFs present with no corresponding fits).  

Based on the fit quality determined in step 12, repeat steps 9-12 to optimize the fitting parameters.  

### Final image corrections and 3D alignment

Select the appropriate “Alignment channel”, which is the imaging channel to be used for Z-section alignment. This should be the channel with the greatest overall biological structure in the dataset (greater overall signal intensity). Typically this is either a global neuropil stain (e.g. WGA or similar), a fluorescent-labeled neuron, or large synaptic structures (e.g. vesicle labels). 

Select the appropriate “Number of processes” to run concurrently during the STORM fitting analysis. This number is limited by the total threads available for computation on the analysis machine. 
You can choose how many processes to run by either typing in the field directly or using the slider. 

Select the start and end frames for each of the desired channels by pressing the "Set Frame Range" button and completing the form. This inputs the desired full frame range for final single molecule fitting of the STORM data. 

Note: The selected start frame for full STORM analysis should be set to a number where the majority of single molecules have begun to photoswitch.  

Once frame ranges have been set, the desired alignment channel is selected, and the number of processes has been specified, press "Multifit Analysis" to run the automated processing of all data into 3D stacks. 

When the first FIJI subprocess launches, input rigid registration parameters to continue the alignment. Information on rigid alignment in FIJI is found here: https://imagej.net/imaging/registration. 

A Python subprocess opens displaying the compressed final image field. Drag cursor within the displayed image to set a sample ROI for cropping the image stack. Release the cursor to identify the ROI. Close the window and press “Escape” on the keyboard to crop the dataset. 

When the second FIJI subprocess launches, input elastic registration parameters to continue the alignment. Elastic alignment is implemented in TrakEM2 based on the original publication (Saalfeld et al., 2012). Information on elastic alignment parameters is found here: https://imagej.net/plugins/elastic-alignment-and-montage. 
