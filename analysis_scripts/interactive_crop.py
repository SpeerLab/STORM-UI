import cv2
import numpy as np

from functools import partial 

""" Adapted from https://www.life2coding.com/crop-image-using-mouse-click-movement-python/""" 

cropping = False
x_start, y_start, x_end, y_end = 0, 0, 0, 0

def mouse_crop(ori_image, event, x, y, flags, param):
    # grab references to the global variables
    global x_start, y_start, x_end, y_end, cropping
    # if the left mouse button was DOWN, start RECORDING
    # (x, y) coordinates and indicate that cropping is being
    if event == cv2.EVENT_LBUTTONDOWN:
        x_start, y_start, x_end, y_end = x, y, x, y
        cropping = True
    # Mouse is Moving
    elif event == cv2.EVENT_MOUSEMOVE:
        if cropping == True:
            x_end, y_end = x, y
    # if the left mouse button was released
    elif event == cv2.EVENT_LBUTTONUP:
        # record the ending (x, y) coordinates
        x_end, y_end = x, y
        cropping = False # cropping is finished
        refPoint = [(x_start, y_start), (x_end, y_end)]
        if len(refPoint) == 2: #when two points were found
            roi = ori_image[refPoint[0][1]:refPoint[1][1], refPoint[0][0]:refPoint[1][0]]
            cv2.imshow("Cropped", roi)
            
def interactive_crop(image):

    ori_image = image.copy()
    scale = 2
    
    print(ori_image.shape) 
    
    window_width = int(ori_image.shape[1] * scale)
    window_height = int(ori_image.shape[0] * scale)
        
    cv2.namedWindow("Interactive Image Cropping", cv2.WINDOW_NORMAL)
    cv2.resizeWindow('Interactive Image Cropping', window_width, window_height)

    cv2.setMouseCallback("Interactive Image Cropping", partial(mouse_crop, ori_image))

    while True:
        i = image.copy()
        if not cropping:
            cv2.imshow("Interactive Image Cropping", image)
        elif cropping:
            cv2.rectangle(i, (x_start, y_start), (x_end, y_end), (255, 0, 0), 2)
            cv2.imshow("Interactive Image Cropping", i)
        key = cv2.waitKey(1)
        
        if key == 27: 
            break 
        
    # close all open windows
    cv2.destroyAllWindows()
    return [x_start, y_start, x_end, y_end]
    
if __name__ == '__main__':
    image = cv2.imread('test.jpg')
    out = interactive_crop(image)
    print(out)