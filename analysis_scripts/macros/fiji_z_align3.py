import os, re, time
import glob
#from ij import IJ
from ini.trakem2 import ControlWindow, Project
from ini.trakem2.display import Display, Patch, LayerSet
from java.awt import Color
from ij import IJ  
import math
from ij.process import ImageStatistics as IS
from ij import Prefs

argv = getArgument()

arg_array = argv.split(" ")


exp_folder = arg_array[0]
print (exp_folder)
IJ.run("Memory & Threads...", "maximum=50000 parallel=24 run");
storm_merged_folder = exp_folder + "cropped/storm_merged/"
conv_merged_folder = exp_folder + "cropped/conv_merged/"
conv_align_folder = exp_folder + "cropped/for_align/"
conv_561_folder = exp_folder + "cropped/conv_561/"
conv_561_folder_ds = exp_folder + "cropped/for_align_ds/"
out_folder = exp_folder + "elastic_align/"
imlist = glob.glob(conv_align_folder + '/*')
num_sections = len(imlist)
print (num_sections)
if not os.path.exists(out_folder):
	os.mkdir (out_folder)
	os.mkdir (out_folder + "/for_align_ds/")
	os.mkdir (out_folder + "/for_align/")
	os.mkdir (out_folder + "/conv_merged_ds/")
	os.mkdir (out_folder + "/conv_merged/")
	os.mkdir (out_folder + "/conv_561_ds/")
	os.mkdir (out_folder + "/conv_561/")
	os.mkdir (out_folder + "/storm_merged_ds/")
	os.mkdir (out_folder + "/storm_merged/")
		
##determine image canvas size (query each image)

ControlWindow.setGUIEnabled(False);  
# 1. Create a TrakEM2 project
IJ.run("TrakEM2 XML...", "select=" + out_folder + "elastic_alignment.xml");
#  2. Create 10 layers (or as many as you need)
project = Project.getProjects().get(0)
#project = Project.newFSProject("blank", None, base_folder)
layerset = project.getRootLayerSet()

for i in range(0,num_sections):  

	layer = layerset.getLayers().get(i)
	rectangle = layerset.get2DBounds()
	backgroundColor = Color.black

	for j in range(4):
		patch_num = j
		scale = 1
		scalesmall = 0.1
		patches = layer.getAll(Patch)
		current_patch = patches[patch_num:(patch_num + 1)]
		name = current_patch[0].getImageFilePath() 
		basename = os.path.basename(name)

		#print current_patch[0].getImageFilePath() 

		ip = Patch.makeFlatImage(
         ImagePlus.COLOR_RGB,
         layer,
         rectangle,
         scale,
         current_patch,
         backgroundColor,
         True)  # use the min and max of each tile
		imp = ImagePlus("img", ip)
		imp2 = imp.flatten()
		ip2 = Patch.makeFlatImage(
         ImagePlus.COLOR_RGB,
         layer,
         rectangle,
         scalesmall,
         current_patch,
         backgroundColor,
         True)  # use the min and max of each tile
		imp3 = ImagePlus("img", ip2)
		imp4 = imp3.flatten()
		#imp2.show()
		if patch_num == 0:
			IJ.saveAs(imp2, "TIF", (exp_folder + "elastic_align/storm_merged/" + basename[:-4] + ".tif"))
			IJ.saveAs(imp4, "TIF", (exp_folder + "elastic_align/storm_merged_ds/" + basename[:-4] + ".tif"))
		if patch_num == 1:
			IJ.saveAs(imp2, "TIF", (exp_folder + "elastic_align/conv_merged/" + basename[:-4] + ".tif"))
			IJ.saveAs(imp4, "TIF", (exp_folder + "elastic_align/conv_merged_ds/" + basename[:-4] + ".tif"))
			IJ.run("Close All", ""); 
		
		if patch_num == 2:
			IJ.run(imp2, "8-bit", "")
			IJ.saveAs(imp2, "TIF", (exp_folder + "elastic_align/conv_561/" + basename[:-4] + ".tif"))
			IJ.run(imp4, "8-bit", "")
			IJ.saveAs(imp4, "TIF", (exp_folder + "elastic_align/conv_561_ds/" + basename[:-4] + ".tif"))
			IJ.run("Close All", ""); 

		if patch_num == 3:
			IJ.run(imp2, "8-bit", "")
			IJ.saveAs(imp2, "TIF", (exp_folder + "elastic_align/for_align/" + basename[:-4] + ".tif"))
			IJ.run(imp4, "8-bit", "")
			IJ.saveAs(imp4, "TIF", (exp_folder + "elastic_align/for_align_ds/" + basename[:-4] + ".tif"))
			IJ.run("Close All", ""); 
			
	print (basename)     

time.sleep(60)
IJ.run("Quit")	

