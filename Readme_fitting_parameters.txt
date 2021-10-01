
Notes on XML parameters for STORM movies
06_01_2018_colenso

Parameters for testing:
Pixel size 153nm
Background 8
Foreground smoothing 0
Iterations 4
Threshold 6
PSF Sigma 1.2

Experiment: the parameters above were held constant and each one was systematically varied 
to examine impacts on single molecule fitting.

Background: tested sigma of 8 and sigma of 10 on 647/750 movies and this parameter did not 
make a difference in the number of localizations identified. 

Foreground smoothing: tested 1.2 and 0 on 647/750 movies (this essentially tested smoothing
on versus off). For 647, this lead to a small number of spurious localizations. For 750, 
the smoothing really helped to pick up dim localizations and led to better overall fitting. 

Iterations: tested 4 and 20 on 647/750 movies. This parameter did not have an impact on the
number of localizations identified. 

Threshold: tested 6 and 8 on 647/750 movies. As expected, increasing the threshold leads to
the detection of fewer molecules (missed localizations at 8). 

PSF Sigma: tested 1.0, 1.2, 1.4, and 1.6 on 647/750. For 647, 1.0 and 1.4 had slightly fewer 
correct localizations compared to 1.2. Additionally, the higher PSF values tended to shift 
some fits away from the center of mass (as determined by eye) and also led to some small 
number of spurios localizations. Similar results were seen for localizations in the 750 
channel. Overall, 1.2 seems reasonable. Perhaps 1.1 or 1.3 will yield better results once 
threshold and smoothing parameters are concurrently modified. 
----
Experiment: after testing each of the parameters individually, I began combining changes 
to multiple parameters to produce an optimal fitting result. I held the sigma constant at 
1.2 and engaged the smoothing function with the same sigma (1.2). I increased the iterations 
to 20 and held the background at 8. Then, the parameter that had the greatest impact was the 
threshold. I increased the threshold for each channel to remove spurious locations in the 
background. This approach worked very well and I was able to achieve parameters sets that 
led to excellent fitting results in 647/750. 
----
The caveat to this work is that I tested the fitting parameters at the end of each 647/750 
STORM movie (last 10 frames), where the background is typically lower and the individual 
molecules are more sparse. It is always a good idea to confirm the quality of your fitting 
at earlier stages of your STORM movies and make adjustments to fitting parameters accordingly. 


