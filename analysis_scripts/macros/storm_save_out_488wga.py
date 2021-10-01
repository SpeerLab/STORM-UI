import os,sys
import glob
import shutil
from ij import IJ  
from ij.process import ImageStatistics as IS
from ij import Prefs
import time


# get experiment folder from python
expfolder = getArgument()
#expfolder = "Z:\\Chenghang\\7.1.20.WT_P2_C\\"
# set paths
storm_folder = expfolder + "stormtiffs\\"
print storm_folder
s_analysisfolder = expfolder + "analysis\\"
ISanalysisfolder = s_analysisfolder + "individual_sections\\"
slicenum = len(glob.glob(ISanalysisfolder + "0*"))
print slicenum
print ISanalysisfolder

channels = ["750storm_","647storm_","561storm_","488storm_"]


for i in range(slicenum):
    for channel in channels:
        base = str(channel)
        if os.path.exists(ISanalysisfolder + "%04d" % i + "\\rawimages\\for_matlab\\" + base + "%03d" % i + ".tiff"):
            print("STORM images ready...beginning Matlab XY alignment")
            IJ.run("Quit")
        else:
            print("saving out 16-bit STORM image #")
            print(ISanalysisfolder + "%04d" % i + "\\rawimages\\for_matlab\\" + base + "%03d" % i + ".tif")
            imp = IJ.openImage((storm_folder + base + "%03d" % i + "_mlist.tiff"))
            IJ.run(imp, "16-bit", "")
            imp.show()
            IJ.run("Enhance Contrast", "saturated=0.3")
            IJ.run("Apply LUT")
            IJ.run(imp, "16-bit", "")
            IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "\\rawimages\\for_matlab\\" + base + "%03d" % i + ".tiff"))
            imp.close();
        IJ.run("Close All", "");

IJ.run("Quit");
