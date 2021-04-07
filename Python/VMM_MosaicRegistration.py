# =============================================================================
# Created By  : Jackson Coole, Brady Hunt, David Brenes, and Patrick Chickey
# Created Date: March 2nd, 2021
# =============================================================================
"""The Module Has Been Built for the manuscript titled "A high frame rate video mosaicking
microendoscope to image large regions of intact tissue with subcellular resolution" 
by Hunt et al"""
# =============================================================================
# Import packages
import cv2
import os
import numpy as np

# Initalize global variables that will be passed from loop to loop
prevFrame = []
kpOld = []
desOld = []
xshiftLast = 0
yshiftLast = 0
detector = cv2.AKAZE_create(2,threshold= 0.00001)


def featureDetection(newFrame,maskFrame):
    # Declare global variables to enable overwriting for the next loop
    global prevFrame 
    global kpOld
    global desOld
    global xshiftLast
    global yshiftLast
    global detector

    # Convert new image into python format 
    newFrame = np.array(newFrame).astype(np.uint8)
    maskFrame = np.array(maskFrame).astype(np.uint8)

    # Detect and Computer AKAZE features
    kpNew, desNew = detector.detectAndCompute(newFrame, maskFrame)
    
    try:
        # Check if there is a prevFrame. If not, this is the first iteration of the process and the prevFrame will be the newFrame
        if prevFrame == []:
            prevFrame = newFrame
            kpOld = kpNew
            desOld = desNew

        # Match the keypoint and descriptor pairs
        ratio = .65

        bf = cv2.BFMatcher(cv2.NORM_L2,False)
        matches = bf.knnMatch(desNew, desOld, k=2)
        good = []

        #status1 = ', len matches = ' + str(len(matches))
        for mat in matches:
            # apply ratio test to find good matches
            try:
                if mat[0].distance < ratio * mat[1].distance:
                    good.append(mat[0])
            except:
                # pass match if only 1 closest match found, knnMatch searches for 2 closest
                pass

        # initialize coordinate lists
        list_kpNew = []
        list_kpOld = []
        for mat in good:
            # Get indices for the matching feature points
            img1_idx = mat.queryIdx
            img2_idx = mat.trainIdx
            # x, y - columns, rows
            # Get the coordinates by indexing the keypoint objects
            (x1, y1) = kpNew[img1_idx].pt
            (x2, y2) = kpOld[img2_idx].pt
            # Append to each list
            list_kpNew.append((x1, y1))
            list_kpOld.append((x2, y2))
            #status3 = ', Mat in good worked'
        try:
            # calculates shifts from the average change in coordinates of the matched feature points
            xshift = (sum(a[0] for a in list_kpNew) - sum(a[0] for a in list_kpOld))/len(list_kpNew)
            yshift = (sum(a[1] for a in list_kpNew) - sum(a[1] for a in list_kpOld))/len(list_kpNew)
            status = 'Passed'
        except:
            # if no good matches are found, prints error statement, new shift is same as last shift to keep trajectory
            xshift = xshiftLast
            yshift = yshiftLast
            #print('err: no matches found, mosaic is off', mat, good)
            status = 'Failed'
        #status = 'Calculated'
    except:
        status = 'Failed'
        xshift = xshiftLast
        yshift = yshiftLast

    # Now that the shift has been calculated, prepare the global variables for the next loop
    prevFrame = newFrame
    kpOld = kpNew
    desOld= desNew
    xshiftLast = xshift
    yshiftLast = yshift

    return([xshift,yshift])