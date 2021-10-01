import sys 
import glob
import os
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.interpolate import make_interp_spline,BSpline
#from scipy.misc import imresize
#from scipy.misc import imsave
from PIL import Image
from imageio import imwrite

def cal_hist(wgafile,num_images):
    print(wgafile)
    A = mpimg.imread(wgafile)
    hist,bins = np.histogram(A.ravel(),255,[1,255])
    return hist

def new_hist(varuse,wgafile,path4,path4a):
    print(wgafile)
    hgram4 = varuse / sum(varuse)
    A = mpimg.imread(wgafile)
    hist,bins = np.histogram(A.ravel(),256,[0,255])
    hist_cum = np.cumsum(hist)
    a = np.array(hist[0])
    b = hgram4*(sum(hist)-hist[0])
    hgram4a = np.concatenate((a,b),axis=None)
    
    constructed = np.array([])
    for i in range(256):
        constructed = np.concatenate((constructed,np.array(i*np.ones(int(np.round(hgram4a[i]))))),axis=None)
    InCaseShort = 255 * np.ones(300)
    constructed = np.concatenate((constructed,InCaseShort),axis=None)
    
    new = np.array(range(256))
    for i in range(256):
        if i == 0:
            first = 0
        else:
            first = hist_cum[i-1]
        second = hist_cum[i]
        result = np.average(constructed[first:second])
        new[i] = np.round(result)

    B = A.ravel()
    C = new[B]
    C = C.reshape(6400,6400)
    C = np.uint8(C)
    #D = imresize(C,(640,640))
    D = np.array(Image.fromarray(C).resize((640,640)))
    name = wgafile.split('\\')
    name = name[-1]
    #!!!!!!
    #Important: need to double check whether the file name should start with 0 or 1. 
    imwrite(path4 + name , C)
    imwrite(path4a + name , D)

    return 0

    
#def main(argment):
#if "__name__" == "__main__":
expfolder = 'Y:\\Chenghang\\06_Testing\\GUI_test\\analysis\\'
#expfolder = argment
path3 = expfolder + 'unaligned\\conv_488\\'
wgafiles = glob.glob(path3 + '*.tif')
num_images = len(wgafiles)

hy4c = np.zeros([num_images,255])
Orders = range(num_images)

i = 0
for x in wgafiles:
    hy4c[i] = cal_hist(x,num_images)
    i = i+1
#pool = mp.Pool(processes = 12)
#print('Starting histogram extraction')
#hy4c = [pool.apply(cal_hist,args = (x,num_images)) for x in wgafiles]
#print('Hitogram extraction done!')

x_hist = np.arange(1,256,1)
x_sec = np.arange(1,num_images+1,1)
cnt = 0
hy4cb = np.zeros([num_images,255])
for file in wgafiles:
    hy4cb[cnt,:] = hy4c[cnt]/sum(hy4c[cnt])
    cnt = cnt + 1

chan = hy4cb
Number = range(255)
varuse4 = np.zeros([num_images,255])
list_x = range(num_images)
#Some warnings here when using test num_images = 4.
for i in Number:
    poly = np.polyfit(list_x, chan[:,i],4)
    poly_y = np.poly1d(poly)(list_x)
    varuse4[:,i] = poly_y

path4 = expfolder + 'unaligned\\for_align\\'
path4a = expfolder + 'unaligned\\for_align_ds\\'
if not os.path.exists(path4):
    os.mkdir(path4)
    os.mkdir(path4a)

print('Starting rendering images')

i = 0
for wgafile in wgafiles:
    new_hist(varuse4[i,:],wgafile,path4,path4a)
    i = i+1

#temp = [pool.apply_async(new_hist,args = (x,y,path4,path4a)) for x,y in zip(varuse4,wgafiles)]
print('Done!, close the parallel processing pool')
#pool.close()
