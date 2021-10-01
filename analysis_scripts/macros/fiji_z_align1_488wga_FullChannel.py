import os, re, time
import glob
from ini.trakem2 import ControlWindow, Project
from ini.trakem2.display import Display, Patch, LayerSet
from java.awt import Color
from ij import IJ  
import math
from ij.process import ImageStatistics as IS
from ij import Prefs
from ij import ImagePlus

argv = getArgument()

arg_array = argv.split(" ")


exp_folder = arg_array[0]
#exp_folder = "Z:/Colenso/10_02_18_SAC/analysis/"
print exp_folder
IJ.run("Memory & Threads...", "maximum=300000 parallel=24 run");

storm_merged_folder = exp_folder + "unaligned/storm_merged/"
conv_merged_folder = exp_folder + "unaligned/conv_merged/"
conv_align_folder = exp_folder + "unaligned/for_align/"
conv_488_folder = exp_folder + "unaligned/conv_488/"
extra_folder = exp_folder + "unaligned/storm_561_488/"
out_folder = exp_folder + "rigid_align/"
imlist = glob.glob(conv_align_folder + '/*')
num_sections = len(imlist)
print num_sections
if not os.path.exists(out_folder):
	os.mkdir (out_folder)
	os.mkdir (out_folder + "/for_align_ds/")
	os.mkdir (out_folder + "/for_align/")
	os.mkdir (out_folder + "/conv_merged_ds/")
	os.mkdir (out_folder + "/conv_merged/")
	os.mkdir (out_folder + "/conv_488_ds/")
	os.mkdir (out_folder + "/conv_488/")
	os.mkdir (out_folder + "/storm_merged_ds/")
	os.mkdir (out_folder + "/storm_merged/")
	os.mkdir (out_folder + "/storm_561_488/")
	os.mkdir (out_folder + "/storm_561_488_ds")
		
##determine image canvas size (query each image)
im_width = 0
im_height =0
for i in range(num_sections):
	imp = IJ.openImage((imlist[i]))
	#print imp.width
	if imp.width>im_width:
		im_width = imp.width
		print im_width

	if imp.height>im_height:
		im_height = imp.height
		print im_height

ControlWindow.setGUIEnabled(False);  
# 1. Create a TrakEM2 project
project = Project.newFSProject("blank", None, out_folder)
loader = project.getLoader()
loader.setMipMapsRegeneration(False) # disable mipmaps
layerset = project.getRootLayerSet()
 
#  2. Create 10 layers (or as many as you need)
for i in range(num_sections):
  layerset.getLayer(i, 1, True)
  layerset.setDimensions(im_width,im_height,LayerSet.NORTHWEST)
# ... and update the LayerTree:
project.getLayerTree().updateList(layerset)
# ... and the display slider
Display.updateLayerScroller(layerset)
filenames = sorted(os.listdir(storm_merged_folder))
filenames1 = sorted(os.listdir(conv_488_folder))
filenames2 = sorted(os.listdir(conv_merged_folder))
filenames3 = sorted(os.listdir(conv_align_folder))
filenames4 = sorted(os.listdir(extra_folder))
#3 Load in images
for i,layer in enumerate(layerset.getLayers()):
    #i = i +150 #use this to start in the middle of a dataset
    filepath = os.path.join(storm_merged_folder, filenames[i])
    patch = Patch.createPatch(project, filepath)
#    patch.updateMipMaps()
    layer.add(patch, False)
#
    filepath = os.path.join(conv_merged_folder, filenames2[i])
    patch = Patch.createPatch(project, filepath)
#    patch.updateMipMaps()
    layer.add(patch, False)

    filepath = os.path.join(conv_488_folder, filenames1[i])
    patch = Patch.createPatch(project, filepath)
#    patch.updateMipMaps()
    layer.add(patch, False)

    filepath = os.path.join(extra_folder, filenames4[i])
    patch = Patch.createPatch(project, filepath)
#    patch.updateMipMaps()
    layer.add(patch, False)

    filepath = os.path.join(conv_align_folder, filenames3[i])
    patch = Patch.createPatch(project, filepath)
#    patch.updateMipMaps()
    layer.add(patch, False)
 
    
Display.update(layer);  

from mpicbg.trakem2.align import Align, AlignTask

Align.alignLayersLinearly(layerset.getLayers(),1) #rigid alignment

# 5. Resize width and height of the world to fit the montages
layerset.setMinimumDimensions()

project.saveAs(project.getLoader().getStorageFolder() + "rigid.xml", True); # overwriting any existing xml file with that name  

for i in range(num_sections):
	layer = layerset.getLayers().get(i)
	rectangle = layerset.get2DBounds()
	backgroundColor = Color.black

	for j in range(5):
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
			IJ.saveAs(imp2, "PNG", (exp_folder + "rigid_align/storm_merged/" + basename[:-4] + ".png"))
			IJ.saveAs(imp4, "PNG", (exp_folder + "rigid_align/storm_merged_ds/" + basename[:-4] + ".png"))
		if patch_num == 1:
			IJ.saveAs(imp2, "PNG", (exp_folder + "rigid_align/conv_merged/" + basename[:-4] + ".png"))
			IJ.saveAs(imp4, "PNG", (exp_folder + "rigid_align/conv_merged_ds/" + basename[:-4] + ".png"))
			IJ.run("Close All", ""); 

		if patch_num == 2:
			IJ.run(imp2, "8-bit", "")
			IJ.saveAs(imp2, "PNG", (exp_folder + "rigid_align/conv_488/" + basename[:-4] + ".png"))
			IJ.run(imp4, "8-bit", "")
			IJ.saveAs(imp4, "PNG", (exp_folder + "rigid_align/conv_488_ds/" + basename[:-4] + ".png"))
			IJ.run("Close All", ""); 

		if patch_num == 4:
			IJ.run(imp2, "8-bit", "")
			IJ.saveAs(imp2, "PNG", (exp_folder + "rigid_align/for_align/" + basename[:-4] + ".png"))
			IJ.run(imp4, "8-bit", "")
			IJ.saveAs(imp4, "PNG", (exp_folder + "rigid_align/for_align_ds/" + basename[:-4] + ".png"))
			IJ.run("Close All", ""); 
		
		if patch_num == 3:
			IJ.saveAs(imp2, "PNG", (exp_folder + "rigid_align/storm_561_488/" + basename[:-4] + ".png"))
			IJ.saveAs(imp4, "PNG", (exp_folder + "rigid_align/storm_561_488_ds/" + basename[:-4] + ".png"))
			
	print basename     

time.sleep(60)
IJ.run("Quit")	

