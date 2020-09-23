import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.signal import find_peaks

def linear_fit(x, y0, a):

    return a*x+y0

def Multi_Gauss_fit(x, *p):
    #[y0_1,a_1,xc_1,sigma,y0_2,a_2,xc_2]
    return (p[0] + p[1]*np.exp(-(x-p[2])**2/(2*p[3]**2)) + p[4] + p[5]*np.exp(-(x-p[6])**2/(2*p[3]**2)))

def reader(dirName):
    imps = {}
    angle = {}
    lines = [] #for counting the lines of the read in files
    tau = 90*10**(-6) #dead time of the Geiger counter
    int_time = 2 #integration time of the detektor in seconds
    d = 201.4*10**(-12) #lattice constant for LiF
    h = 6.582119569*10**(-16) #h in eV*s
    c = 299792458 #c in m/s
    j = 0 #for skipping the first few lines
    file_numb = 0 #counts number of files in specified directory
    temp = []
    temp1 = []
    
    for root, dirs , files in os.walk(dirName):
        for file in files:
            if '.txt' not in file:
                print("Reading in file " + os.path.join(root, file))
                with open(os.path.join(root, file)) as f:
                    reader = csv.reader(f, delimiter = '\t')
                    for row in reader:
                        if j > 2:
                            #print(j-3, " ", row[0].replace(',', '.'))
                            imps[file_numb, j - 3] = float(row[1])
                            angle[file_numb, j - 3] = float(row[0].replace(',', '.'))
                            
                        j += 1
                lines.append(j-3)
                j = 0
                file_numb += 1 
    
    
    f.close()

    energy = np.empty([], dtype = float)
    intensity = np.empty([], dtype = float)
    y_err = np.empty([], dtype = float)
    y_max = np.empty((0,file_numb), dtype = float)

    for k in range(0, file_numb):
        for i in range(0, lines[k]):
            temp.append(float(imps[k, i]*int_time/(1-tau*imps[k, i]*int_time)))
            temp1.append(float(h*c/(2*d*np.sin(angle[k, i]*np.pi/180)*1000)))
        y_max = np.append(y_max, max(temp))
        #print(max(temp), "numpy values ", y_max[k])
        intensity = np.append(intensity, temp/y_max[k])
        energy = np.append(energy, temp1)
        temp.clear()
        temp1.clear()

    #for testing purposes

    #for k in range(0, file_numb):
    #    for i in range(0, lines[k]):
    #        temp.append(imps[k, i])
    #        temp1.append(angle[k, i])
    #    y_max = np.append(y_max, max(temp))
    #    intensity = np.append(intensity, temp/y_max[k])
    #    #print(max(temp))
    #    #print(temp)
    #    energy = np.append(energy, temp1)
    #    temp.clear()
    #    temp1.clear()

    del temp
    del temp1
    del imps
    del angle
   
    #for l in range(0, file_numb):
    #    for i in range(0, lines[l]):
    #        intensity[l, i] = intensity[l, i]/y_max[l]
    #        if intensity[l, i] == 0.0:
    #            y_err[l, i] == 0
    #        else:
    #            y_err[l, i] = np.sqrt(intensity[l, i])


    return energy, intensity, y_err, y_max, file_numb, lines

def fit_spectrum(dirName):

    energy, intensity, y_err, file_numb, lines = reader(dirName)

    perr = np.empty([file_numb, 7], dtype = float)
    params = [0]*file_numb
    params_cov = [0]*file_numb

    x_min = 1.253
    x_max = 1.3
    x_min1 = 1.39
    x_max1 = 1.43

    #[y0_1,a_1,xc_1,sigma,y0_2,a_2,xc_2]
    initial_values = [0,0.01,x_min,0.008,0,0.01,x_min1]
    bounds = ([0,0,x_min,0,0,0,x_min1],[0.06,1,x_max,0.01,0.1,0.5,x_max1])

    for i in range(0, file_numb):
        params[i], params_cov[i] = scipy.optimize.curve_fit(Multi_Gauss_fit, energy[i], intensity[i], p0 = initial_values, bounds = bounds, sigma = None, absolute_sigma = True)
        perr[i] = np.sqrt(np.diag(params_cov[i]))/np.sqrt(len(energy[i]))
    
    #plt.plot(x_data, Multi_Gauss_fit(x_data, *params))

    for i in range(0, file_numb):
        print("Center pos. 1 = ", params[i][2], "+/-", perr[i][2], "\n" , "Center pos. 2 = ", params[i][6], "+/-", perr[i][6])

    return params, params_cov, energy, intensity, y_err, file_numb

def plot_spectrum(x_plt, y_plt, y_err, n):

    plt.figure(n)
    plt.errorbar(x_plt, y_plt, yerr = y_err, fmt = 'x', markersize = 3 )
    #plt.title()

    plt.xlabel("Energy [keV]")
    plt.ylabel("Intensity [arbitrary units]")
    plt.grid()

    return 0
        
def plot_spectrum_w_fit(x_plt, y_plt, y_err, params, n):

    plt.figure(n)
    plt.errorbar(x_plt, y_plt, yerr = y_err, fmt = 'x', markersize = 3 )
    #plt.title()
    plt.plot(x_plt, Multi_Gauss_fit(x_plt, *params))

    plt.xlabel("Energy [keV]")
    plt.ylabel("Intensity [arbitrary units]")
    plt.grid()

    return 0

def Duane_Hunt():

    current_altered = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu\Current_altered"
    voltage_altered = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu\Voltage_altered"
    Test = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Test"

    #x_plt, y_plt, y_err, file_numb, line_nums = reader(Test)
    #plot_spectrum(x_plt, y_plt, None)
    #params, params_cov, x_plt, y_plt, y_err, file_numb = fit_spectrum(voltage_altered)
    
    #for n in range(0, file_numb):
    #    plot_spectrum_w_fit(x_plt[n], y_plt[n], None, params[n], n)

    x_plt, y_plt, y_err, y_max, file_numb, lines = reader(voltage_altered)

    sum = 0

    #print(x_plt[0:lines[0]+1])
    #print(len(y_plt[:lines[0]+1]))
    #print(len(x_plt[:lines[0]+1]))

    #plot spectra for different voltages or currents

    for i in range(0, file_numb):
        if i == 0:
            plot_spectrum(x_plt[1:lines[i]+1], y_plt[1:lines[i]+1], None, i)
        elif i < file_numb - 1:
            plot_spectrum(x_plt[sum:sum+lines[i+1]+1], y_plt[sum:sum+lines[i+1]+1], None, i)
        elif i == file_numb - 1:
            plot_spectrum(x_plt[sum:sum+lines[i]+1], y_plt[sum:sum+lines[i]+1], None, i)
        sum += lines[i] + 1

    peaks, _ = find_peaks(y_plt[1:lines[0]+1], height = (0.01, 1), threshold=0.05)
    for i in range(0,len(peaks)):
        print(x_plt[peaks[i]])

      

    plt.show()

    return 0

def main():
    Duane_Hunt()

if __name__ == "__main__" :
    main()