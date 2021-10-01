import interactive_crop 
import cv2

image = cv2.imread('test.jpg')
out = interactive_crop.interactive_crop(image)
print(out)