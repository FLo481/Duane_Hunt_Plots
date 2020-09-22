import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def linear_fit(x, y0, a):

    return a*x+y0

def reader(dirName):
    imps = []
    angle = []
    wavelength = []
    intensity = []
    tau = 90*E-6
    int_time = 2
    d = 201.4E-12
    j = 0 #for skipping the first few lines
    
    for root, dirs , files in os.walk(dirName):
        for file in files:
            if file.endswith(".txt"):
                print("Reading in file " + os.path.join(root, file))
                with open(os.path.join(root, file)) as f:
                    reader = csv.reader(f, delimiter = '\t')
                    for row in reader:
                        if j > 2:
                            #print(row[0])
                            #print(row[1])
                            imps.append(row[1])
                            angle.append(row[0])
                        j += 1
                j = 0

    f.close()

    for i in range(0,len(angle)):
        intensity.append(imps[i]*int_time/(1-tau*imps[i]*int_time))
        wavelength.append(2*d*np.sin(angle[i]*np.pi/180))

    #for i in range(len(intensity)):
    #    print(intensity[i] + " ; " + wavelength[i])

    return intensity, wavelength

def main():
    dirName = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\DuaneHunt"
    reader(dirName)

if __name__ == "__main__" :
    main()