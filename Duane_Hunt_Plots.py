import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def linear_fit(x, y0, a):

    return a*x+y0

def Multi_Gauss_fit(x, *p):

    return (p[0] + p[1]*np.exp(-(x-p[2])**2/(2*p[3]**2)) + p[4] + p[5]*np.exp(-(x-p[6])**2/(2*p[3]**2)))

def reader(dirName):
    imps = [] #recorded counts per seconds
    angle = [] #measured angle of the goniometer
    energy = []
    intensity = []
    tau = 90*10**(-6) #dead time of the Geiger counter
    int_time = 2 #integration time of the detektor in seconds
    d = 201.4*10**(-12)
    h = 6.582119569*10**(-16)
    c = 299792458
    j = 0 #for skipping the first few lines
    
    for root, dirs , files in os.walk(dirName):
        for file in files:
            #if '.txt' in file:
                print("Reading in file " + os.path.join(root, file))
                with open(os.path.join(root, file)) as f:
                    reader = csv.reader(f, delimiter = '\t')
                    for row in reader:
                        if j > 2:
                            angle.append(row[0].replace(',','.'))
                            imps.append(row[1])
                        j += 1
                j = 0
    
    f.close()

    imps = list(map(float, imps ))
    angle = list(map(float, angle ))

    for i in range(0,len(angle)):
        intensity.append(imps[i]*int_time/(1-tau*imps[i]*int_time))
        energy.append(h*c/(2*d*np.sin(angle[i]*np.pi/180)))

    #for i in range(len(intensity)):
    #    print(intensity[i] + " ; " + wavelength[i])

    return intensity, energy

def plot_spectrum():

    Cu = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu"
    Test = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Test"
    intensity , energy = reader(Test)
    x_plt = np.empty(len(intensity), dtype = float)
    y_plt = np.empty(len(intensity), dtype = float)
    y_err = np.empty(len(intensity), dtype = float)
    x_min = 0.138
    x_min1 = 0.153
    x_max = 0.1575

    y_max = max(intensity)

    for i in range(0, len(intensity)):
        x_plt[i] = energy[i]/1000
        y_plt[i] = intensity[i]/y_max

    for i in range(0, len(y_plt)):
        if y_plt[i] == 0.0:
            y_err[i] == 0
        else:
            y_err[i] = y_plt[i]/(np.sqrt(y_plt[i]))

    initial_values = [0,0.01,x_min,0.04,0,0.01,x_min1]
    bounds = ([0,0.01,x_min,0.04,0,0,x_min1],[0.07,0.1,0.141,0.09,0.1,0.5,x_max])

    plt.errorbar(x_plt, y_plt, yerr = None, fmt = 'x', markersize = 3 )
    params, params_cov = scipy.optimize.curve_fit(Multi_Gauss_fit, x_plt, y_plt, p0 = initial_values, bounds = bounds, sigma = None, absolute_sigma = True)
    plt.plot(x_plt, Multi_Gauss_fit(x_plt, *params))

    print(params[2],"\n" ,params[6])

    plt.xlabel("Energy [keV]")
    plt.ylabel("Intensity [arbitrary units]")
    plt.xlim(1.2,1.8)
    plt.grid()
    plt.show()
    plt.clf()

    return 0

def main():
    plot_spectrum()

if __name__ == "__main__" :
    main()