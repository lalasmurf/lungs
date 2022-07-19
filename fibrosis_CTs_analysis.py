# This is a sample Python script.

import time

import matplotlib.pyplot as plt
import numpy as np
from pydicom import dcmread
from scipy.spatial import distance


def process_cts():

    start_time = time.time()

    isExperiment = False
    ctNumber = 1

    fibrosisLungPath = r"PATH/TO/DICOM/FILE.DCM";

    ds = dcmread(fibrosisLungPath)

    plt.imshow(ds.pixel_array, cmap=plt.cm.bone)
    plt.show()

    ct = ds.pixel_array
    threshold = 50 #gradient threshold
    n = 65 # sample area width/length

    adjacencyMatrix = []

    sampleCT = ct[295:360, 140:205]  # sampled area coordinates

    # calculate ggo/emphisema/consolidation matrix
    convertedImageGgo = []
    convertedImageGgo_e = []
    convertedImageGgo_g = []
    convertedImageGgo_c = []
    convertToHu(convertedImageGgo, n, ds, sampleCT)
    convertToHuEmphysema(convertedImageGgo_e, n, ds, sampleCT)
    convertToHuGgo(convertedImageGgo_g, n, ds, sampleCT)
    convertToHuConsolid(convertedImageGgo_c, n, ds, sampleCT)

    plt.imshow(convertedImageGgo, cmap=plt.cm.bone)
    plt.show()
    plt.clf()


    if isExperiment == False :
        # generate_histo(ctNumber, sampleCT, ds.RescaleIntercept, ds.RescaleSlope)

        # Se presupune ca folderul e deja creat de tine in prealabil
        imagePath = "C:\\PATH\\TO\\graph_adj\\" + str(ctNumber) + "\\separated_stripes\\" + "CT" + str(ctNumber)

        generate_separated_images(imagePath, convertedImageGgo, convertedImageGgo_e, convertedImageGgo_g, convertedImageGgo_c)

        adjFilePath = "C:\\PATH\\TO\\graph_adj\\" + str(ctNumber) + "\\separated_stripes\\edgeList_CT_"+ str(ctNumber)+"-t50_rd4"
        generate_all_adj_matrices_at_once(adjFilePath,threshold, n, convertedImageGgo, convertedImageGgo_e, convertedImageGgo_g, convertedImageGgo_c)

    print("--- %s seconds ---" % (time.time() - start_time))


def generate_histo(folderNo, sampleCT, intercept, slope):
    first_patient_pixels = get_pixels_hu(sampleCT, intercept, slope)
    plt.hist(first_patient_pixels.flatten(), color='c')
    plt.xlabel("Hounsfield Units (HU)")
    plt.ylabel("Frequency")
    plt.savefig(str(folderNo) + "PATH/TO/HISTO/FILE.png")
    plt.clf()

def generate_separated_images(imagePath, convertedImageGgo, convertedImageGgo_e, convertedImageGgo_g, convertedImageGgo_c):

    plt.imshow(convertedImageGgo, cmap=plt.cm.bone)
    plt.savefig(imagePath + "_snapshot_e_g_c" +".png")

    plt.imshow(convertedImageGgo_e, cmap=plt.cm.bone)
    plt.savefig(imagePath + "_snapshot_e" + ".png")

    plt.imshow(convertedImageGgo_g, cmap=plt.cm.bone)
    plt.savefig(imagePath + "_snapshot_g" + ".png")

    plt.imshow(convertedImageGgo_c, cmap=plt.cm.bone)
    plt.savefig(imagePath + "_snapshot_c" + ".png")

def generate_all_adj_matrices(adjFilePath, threshold,n, convertedImageGgo,convertedImageGgo_e, convertedImageGgo_g, convertedImageGgo_c ):
    radialDistance = 4;
    # ALL
    adjacencyMatrix1 = []
    adjFileName = adjFilePath + "_e_ggo_consolid_65x65_1.csv"
    adjacencyMatrix1 = calculateAdjacencyMatrix3(adjacencyMatrix1, radialDistance, convertedImageGgo, threshold,
                                                 adjFileName, n)
    adjacencyMatrix1 = []

    # EMPHYSEMA

    adjacencyMatrix2 = []
    adjFileName = adjFilePath + "_e_65x65_1.csv"
    adjacencyMatrix2 = calculateAdjacencyMatrix3(adjacencyMatrix2, radialDistance, convertedImageGgo_e, threshold,
                                                 adjFileName, n)
    adjacencyMatrix2 = []

    # GGO

    adjacencyMatrix3 = []
    adjFileName = adjFilePath + "_ggo_65x65_1.csv"
    adjacencyMatrix3 = calculateAdjacencyMatrix3(adjacencyMatrix3, radialDistance, convertedImageGgo_g, threshold,
                                                 adjFileName, n)
    adjacencyMatrix3 = []

    # CONSOLIDATION

    adjacencyMatrix4 = []
    adjFileName = adjFilePath + "_c_65x65_1.csv"
    adjacencyMatrix4 = calculateAdjacencyMatrix3(adjacencyMatrix4, radialDistance, convertedImageGgo_c, threshold,
                                                 adjFileName, n)
    adjacencyMatrix4 = []


def generate_all_adj_matrices_at_once(adjFilePath, threshold,n, convertedImageGgo,convertedImageGgo_e, convertedImageGgo_g, convertedImageGgo_c):
    radialDistance = 4

    # ALL
    adjacencyMatrix1 = []
    adjFileName1 = adjFilePath + "_e_ggo_consolid_65x65_1.csv"

    adjFileName2 = adjFilePath + "_e_65x65_1.csv"

    adjFileName3 = adjFilePath + "_ggo_65x65_1.csv"

    adjFileName4 = adjFilePath + "_c_65x65_1.csv"

    adjacencyMatrix1 = calculateAdjacencyMatrix4(radialDistance, convertedImageGgo, convertedImageGgo_e,
                                                 convertedImageGgo_g, convertedImageGgo_c, threshold,
                                                 adjFileName1, adjFileName2, adjFileName3, adjFileName4, n)


def calculateAdjacencyMatrix4(radialDistance, sampleCT1, sampleCT2, sampleCT3, sampleCT4, threshold, fileName1,
                              fileName2, fileName3, fileName4, n):
    f1 = open(fileName1, "w")
    f2 = open(fileName2, "w")
    f3 = open(fileName3, "w")
    f4 = open(fileName4, "w")

    adjacencyMatrix1 = []
    adjacencyMatrix2 = []
    adjacencyMatrix3 = []
    adjacencyMatrix4 = []

    adjList1 = []
    adjList2 = []
    adjList3 = []
    adjList4 = []

    for i in range(0, n * n - 1):
        # print(i)
        currentLine1 = []
        currentLine2 = []
        currentLine3 = []
        currentLine4 = []
        for j in range(0, n * n - 1):
            xI = i // n
            yI = i % n
            xJ = j // n
            yJ = j % n
            a = (xI, yI)
            b = (xJ, yJ)
            distIJ = distance.euclidean(a, b)

            delta1 = abs(sampleCT1[xI][yI] - sampleCT1[xJ][yJ])
            delta2 = abs(sampleCT2[xI][yI] - sampleCT2[xJ][yJ])
            delta3 = abs(sampleCT3[xI][yI] - sampleCT3[xJ][yJ])
            delta4 = abs(sampleCT4[xI][yI] - sampleCT4[xJ][yJ])
            # print(xI)
            # print(yI)

            # e_ggo_c
            if radialDistance >= distIJ > 0 and delta1 < threshold and sampleCT1[xI][yI] != 0 and sampleCT1[xJ][yJ] != 0:
                currentLine1.append(1)
                edge = str(i) + ' ' + str(j)
                altEdge1 = str(j) + ' ' + str(i)
                if edge not in adjList1:
                    adjList1.append(edge)
                    adjList1.append(altEdge1)
                    f1.write(str(i) + ' ' + str(j))
                    f1.write('\n')
            else:
                currentLine1.append(0)

            # emphysema
            if radialDistance >= distIJ > 0 and delta2 < threshold and sampleCT2[xI][yI] != 0 and sampleCT2[xJ][yJ] != 0:
                currentLine2.append(1)
                edge = str(i) + ' ' + str(j)
                altEdge2 = str(j) + ' ' + str(i)
                if edge not in adjList2:
                    adjList2.append(edge)
                    adjList2.append(altEdge2)
                    f2.write(str(i) + ' ' + str(j))
                    f2.write('\n')
            else:
                currentLine2.append(0)

            # ggo
            if radialDistance >= distIJ > 0 and delta3 < threshold and sampleCT3[xI][yI] != 0 and sampleCT3[xJ][
                yJ] != 0:
                currentLine3.append(1)
                edge = str(i) + ' ' + str(j)
                altEdge3 = str(j) + ' ' + str(i)
                if edge not in adjList3:
                    adjList3.append(edge)
                    adjList3.append(altEdge3)
                    f3.write(str(i) + ' ' + str(j))
                    f3.write('\n')
            else:
                currentLine3.append(0)

            # consolidation
            if radialDistance >= distIJ > 0 and delta4 < threshold and sampleCT4[xI][yI] != 0 and sampleCT4[xJ][
                yJ] != 0:
                currentLine4.append(1)
                edge = str(i) + ' ' + str(j)
                altEdge4 = str(j) + ' ' + str(i)
                if edge not in adjList4:
                    adjList4.append(edge)
                    adjList4.append(altEdge4)
                    f4.write(str(i) + ' ' + str(j))
                    f4.write('\n')
            else:
                currentLine4.append(0)

        adjacencyMatrix1.append(np.array(currentLine1))
        adjacencyMatrix2.append(np.array(currentLine2))
        adjacencyMatrix3.append(np.array(currentLine3))
        adjacencyMatrix4.append(np.array(currentLine4))
        # print(currentLine)
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    adjacencyMatrix = np.array(adjacencyMatrix1)
    return adjacencyMatrix

def get_pixels_hu(crop, intercept, slope):
    image = np.stack([crop])
    # Convert to int16 (from sometimes int16),
    # should be possible as values should always be low enough (<32k)
    image = image.astype(np.int16)

    # Set outside-of-scan pixels to 0
    # The intercept is usually -1024, so air is approximately 0
    image[image == -2000] = 0

    # Convert to Hounsfield units (HU)
    # intercept = scans[0].RescaleIntercept
    # slope = scans[0].RescaleSlope

    if slope != 1:
        image = slope * image.astype(np.float64)
        image = image.astype(np.int16)

    image += np.int16(intercept)

    return np.array(image, dtype=np.int16)

def convertToHu(convertedImageGgo, n, ds, sampleCT):

    rescaleSlope = ds.RescaleSlope
    rescaleIntercept = ds.RescaleIntercept

    for i in range(0, n):
        currentLine = []
        for j in range(0, n):
            gradient = sampleCT[i, j]
            nodeHu = rescaleSlope * gradient + rescaleIntercept
            if (-634 <= nodeHu <= -368 or -109 <= nodeHu <= 9 or -1024 < nodeHu < -984) and (i >10 or j < 47):
                currentLine.append(gradient)
            else:
                currentLine.append(0)
        convertedImageGgo.append(np.array(currentLine))


def convertToHuAll(convertedImageGgo, n, ds, sampleCT):

    rescaleSlope = ds.RescaleSlope
    rescaleIntercept = ds.RescaleIntercept

    for i in range(0, n):
        currentLine = []
        for j in range(0, n):
            gradient = sampleCT[i, j]
            nodeHu = rescaleSlope * gradient + rescaleIntercept
            if -634 <= nodeHu <= -368 or -109 <= nodeHu <= 9:
                currentLine.append(gradient)
            else:
                currentLine.append(0)
        convertedImageGgo.append(np.array(currentLine))

def convertToHuEmphysema(convertedImageGgo, n, ds, sampleCT):

    rescaleSlope = ds.RescaleSlope
    rescaleIntercept = ds.RescaleIntercept

    for i in range(0, n):
        currentLine = []
        for j in range(0, n):
            gradient = sampleCT[i, j]
            nodeHu = rescaleSlope * gradient + rescaleIntercept
            if (-1024 < nodeHu < -984) and (i >10 or j < 47):
                currentLine.append(gradient)
            else:
                currentLine.append(0)
        convertedImageGgo.append(np.array(currentLine))

def convertToHuGgo(convertedImageGgo, n, ds, sampleCT):

    rescaleSlope = ds.RescaleSlope
    rescaleIntercept = ds.RescaleIntercept

    for i in range(0, n):
        currentLine = []
        for j in range(0, n):
            gradient = sampleCT[i, j]
            nodeHu = rescaleSlope * gradient + rescaleIntercept
            if (-634 <= nodeHu <= -368) and (i >10 or j < 47):
                currentLine.append(gradient)
            else:
                currentLine.append(0)
        convertedImageGgo.append(np.array(currentLine))

def convertToHuConsolid(convertedImageGgo, n, ds, sampleCT):

    rescaleSlope = ds.RescaleSlope
    rescaleIntercept = ds.RescaleIntercept

    for i in range(0, n):
        currentLine = []
        for j in range(0, n):
            gradient = sampleCT[i, j]
            nodeHu = rescaleSlope * gradient + rescaleIntercept
            if (-109 <= nodeHu <= 9) and (i >10 or j < 47):
                currentLine.append(gradient)
            else:
                currentLine.append(0)
        convertedImageGgo.append(np.array(currentLine))

def calculateAdjacencyMatrix3(adjacencyMatrix, radialDistance, sampleCT, threshold, fileName, n):
    f = open(fileName,"w")
    adjList = []
    for i in range(0, n*n - 1):
        # print(i)
        currentLine = []
        for j in range(0, n * n - 1):
            xI = i // n
            yI = i % n
            xJ = j // n
            yJ = j % n
            a = (xI, yI)
            b = (xJ, yJ)
            distIJ = distance.euclidean(a, b)
            delta = abs(sampleCT[xI][yI] - sampleCT[xJ][yJ])
            # print(xI)
            # print(yI)
            if radialDistance >= distIJ > 0 and delta < threshold and sampleCT[xI][yI] != 0 and sampleCT[xJ][yJ] != 0:
                currentLine.append(1)
                edge = str(i) + ' ' + str(j)
                altEdge = str(j) + ' ' + str(i)
                if edge not in adjList:
                    adjList.append(edge)
                    adjList.append(altEdge)
                    f.write(str(i) + ' ' + str(j))
                    f.write('\n')
            else:
                currentLine.append(0)
        adjacencyMatrix.append(np.array(currentLine))
        # print(currentLine)
    f.close()
    adjacencyMatrix = np.array(adjacencyMatrix)
    return adjacencyMatrix


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    process_cts()
