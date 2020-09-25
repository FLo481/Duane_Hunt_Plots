import os
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.signal import find_peaks

def linear_fit(x, y0, a):

    return a*x+y0

def Multi_Gauss_fit(x, *p):
    #[y0_1,a_1,xc_1,sigma,y0_2,a_2,xc_2]
    return (p[0] + p[1]*np.exp(-(x-p[2])**2/(2*p[3]**2)) + p[4] + p[5]*np.exp(-(x-p[6])**2/(2*p[3]**2)))

def three_halfs_fit(x, *p):

    return p[0]+p[1]*(x-p[2])**(p[3])

def reader(dirName):
    imps = {}
    angle = {}
    current_temp = []
    voltage_temp = []
    lines = [] #for counting the lines of the read in files
    tau = 90*10**(-6) #dead time of the Geiger counter
    d = 201.4*10**(-12) #lattice constant for LiF
    j = 0 #for skipping the first few lines
    file_numb = 0 #counts number of files in specified directory
    temp = []
    temp1 = []
    
    for root, dirs , files in os.walk(dirName):
        for file in files:
            if '.txt' not in file:
                print("Reading in file " + os.path.join(root, file))
                #print(str(file))
                voltage_temp.append(re.findall(r'\d+', str(file))[4])
                current_temp.append(re.findall(r'\d+', str(file))[5])
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

    voltage = np.empty(len(voltage_temp), dtype = float)
    current = np.empty(len(current_temp), dtype = float)

    voltage[:] = voltage_temp
    current[:] = current_temp

    energy = np.empty([], dtype = float)
    intensity = np.empty([], dtype = float)
    y_err = np.empty([], dtype = float)

    for k in range(0, file_numb):
        for i in range(0, lines[k]):
            temp.append(float(2*imps[k, i]/(1-tau*imps[k, i])))
            temp1.append(float(10**(10)*2.0*d*np.sin(angle[k, i]*np.pi/180.0)))
        intensity = np.append(intensity, temp)
        energy = np.append(energy, temp1)
        temp.clear()
        temp1.clear()
    
    #for testing purposes

    #for k in range(0, file_numb):
    #    for i in range(0, lines[k]):
    #        temp.append(imps[k, i])
    #        temp1.append(angle[k, i])
    #    intensity = np.append(intensity, temp)
    #    energy = np.append(energy, temp1)
    #    temp.clear()
    #    temp1.clear()

    del temp
    del temp1
    del imps
    del angle
    del voltage_temp
    del current_temp

    return energy, intensity, file_numb, lines, voltage, current

def fit_spectrum(dirName):

    energy, intensity, file_numb, lines, voltage, current = reader(dirName)

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

    return params, params_cov, energy, intensity, y_err, file_numb, voltage, current

def plot_spectrum(x_plt, y_plt, y_err, label, n):

    plt.figure(n)
    plt.errorbar(x_plt, y_plt, yerr = y_err, fmt = 'x', markersize = 3, label = label)
    #plt.title()

    if n == 1024:
        params, params_cov = scipy.optimize.curve_fit(three_halfs_fit, x_plt, y_plt, sigma = y_err, p0 = [-100,0.1,0,1], bounds = ([-500,0,-np.inf,0.5],[1000,100,np.inf,np.inf]), absolute_sigma = True, maxfev = 99999)
        perr = np.sqrt(np.diag(params_cov))/np.sqrt(x_plt.shape[0])
        plt.plot(x_plt, three_halfs_fit(x_plt, *params))
        print("U_K = ", params[2], "+/-", perr[2], "Exponent = ", params[3], "+/-", perr[3])
        plt.xlabel("Voltage [kV]", fontsize = 16)
        plt.ylabel("Intensity [arbitrary units]", fontsize = 16)
    elif n == 2048:
        params, params_cov = scipy.optimize.curve_fit(three_halfs_fit, x_plt, y_plt, sigma = y_err, p0 = [-100,0.1,0,1], bounds = ([-500,0,-np.inf,0.5],[1000,100,np.inf,np.inf]), absolute_sigma = True, maxfev = 99999)
        perr = np.sqrt(np.diag(params_cov))/np.sqrt(x_plt.shape[0])
        plt.plot(x_plt, three_halfs_fit(x_plt, *params))
        print("U_K = ", params[2], "+/-", perr[2], "Exponent = ", params[3], "+/-", perr[3])
        plt.xlabel("Current [mA]", fontsize = 16)
        plt.ylabel("Intensity [arbitrary units]", fontsize = 16)
    else:
        plt.xlabel("Energy [$A^°$]", fontsize = 16)
        plt.ylabel("Intensity [arbitrary units]", fontsize = 16)

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

def find_maxima(x_plt, y_plt, file_numb, lines):

    sum = 0

    lambda_min = []
    peaks = [0]*2*file_numb
    n = [0]*file_numb
    params1 = [0]*(file_numb-2)
    params_cov1 = [0]*(file_numb-2)
    ints_of_max = []

    #print(x_plt[0:lines[0]+1])
    #print(y_plt[0:lines[0]+1])
    #print(len(y_plt[:lines[0]+1]))
    #print(len(x_plt[:lines[0]+1]))

    #plot spectra for different voltages or currents and find K_alpha and K_beta maxima
    print(lines)

    for i in range(0, file_numb):
        #print(lines[i])
        if i == 0:
            #peak position calculation

            #plot_spectrum(x_plt[1:lines[i]+1], y_plt[1:lines[i]+1], np.sqrt(y_plt[1:lines[i]+1]), "File {}" .format(i), i)
            peaks[i], _ = find_peaks(y_plt[1:lines[i]+1], height = 100, threshold=20)
            ints_of_max.append(y_plt[peaks[i][0]+1])
            ints_of_max.append(y_plt[peaks[i][1]+1])
            #print("Peak pos. : ", x_plt[peaks[i]+1], "Intensity at peak pos. :", y_plt[peaks[i]+1])

            #print("File number ", i+1, " : " ,y_plt[1:lines[i]+1]) ### correct counting for the indices

            #lambda_min calculation

            for k in range(1, x_plt[1:lines[i]+1].shape[0]+1):
                if 1.65 > x_plt[k:lines[i]+1][0] > 1.595:
                    n[i] += 1
            params, params_cov = scipy.optimize.curve_fit(linear_fit, x_plt[lines[i]-n[i]+1:lines[i]+1] , y_plt[lines[i]-n[i]+1:lines[i]+1], sigma = np.sqrt( y_plt[lines[i]-n[i]+1:lines[i]+1]), absolute_sigma = True)
            #plt.plot(x_plt[lines[i]-n[i]+1:lines[i]+1], linear_fit(x_plt[lines[i]-n[i]+1:lines[i]+1], params[0], params[1]))
            lambda_min.append(-linear_fit(0, params[0], params[1])/params[1])
        elif i < file_numb - 1:
            #peak position calculation

            #plot_spectrum(x_plt[sum+1:sum+lines[i]+1], y_plt[sum+1:sum+lines[i]+1], np.sqrt(y_plt[sum+1:sum+lines[i]+1]), "File {}" .format(i), i)
            peaks[i], _ = find_peaks(y_plt[sum:sum+lines[i+1]+1], height = 120, threshold=80)
            ints_of_max.append(y_plt[sum+peaks[i][0]])
            ints_of_max.append(y_plt[sum+peaks[i][1]])
            #print("Peak pos. : ", x_plt[sum+peaks[i]], "Intensity at peak pos. :", y_plt[sum+peaks[i]])

            #print("File number ", i+1, " : " ,y_plt[sum+1:sum+lines[i]+1]) ### correct counting for the indices

            #lambda_min calculation

            for k in range(sum, sum + x_plt[sum:sum+lines[i+1]+1].shape[0]+1):
                if 1.65 > x_plt[k:sum+lines[i+1]+2][0] > 1.595:
                    n[i] += 1
            #print(sum+lines[i]-n[i]+2, " ", sum+lines[i]+2, n[i])
            #print(y_plt[sum+lines[i]-n[i]+2:sum+lines[i]+2])
            params1[i-1], params_cov1[i-1] = scipy.optimize.curve_fit(linear_fit, x_plt[sum+lines[i]-n[i]+2:sum+lines[i]+2] , y_plt[sum+lines[i]-n[i]+2:sum+lines[i]+2], sigma = np.sqrt( y_plt[sum+lines[i]-n[i]+2:sum+lines[i]+2]), absolute_sigma = True)
            #plt.plot(x_plt[sum+lines[i]-n[i]+2:sum+lines[i]+2], linear_fit(x_plt[sum+lines[i]-n[i]+2:sum+lines[i]+2], params1[i][0], params1[i][1]))
            #lambda_min.append(-linear_fit(0, params1[i][0], params1[i][1])/params1[i][1])
        elif i == file_numb - 1:
            #peak position calculation

            #plot_spectrum(x_plt[sum:sum+lines[i]+1], y_plt[sum:sum+lines[i]+1], np.sqrt( y_plt[sum:sum+lines[i]+1]), "File {}" .format(i), i)
            peaks[i], _ = find_peaks(y_plt[sum:sum+lines[i]+1], height = 5000, threshold = 1000)
            ints_of_max.append(y_plt[sum+peaks[i][0]])
            ints_of_max.append(y_plt[sum+peaks[i][1]])
            #print("Peak pos. : ", x_plt[sum+peaks[i]], "Intensity at peak pos. :", y_plt[sum+peaks[i]])

            #print("File number ", i+1, " : " ,y_plt[sum+1:sum+lines[i]+1]) ### correct counting for the indices

            #lambda_min calculation


        sum += lines[i]
        
    del n

    y = np.empty(len(ints_of_max), dtype = float)
    y[:] = ints_of_max
    r = len(ints_of_max)
    r = int(r/2)
    y_K_alpha = np.empty(r, dtype = float)
    y_K_beta = np.empty(r, dtype = float)
    y_K_alpha_err = np.empty(r, dtype = float)
    y_K_beta_err = np.empty(r, dtype = float)

    for i in range(0, r*2):
        if i%2 == 0.0:
            y_K_beta[int(1/2*i)] = y[i]
            y_K_beta_err[int(1/2*i)] = np.sqrt(y[i])
        else:
            y_K_alpha[int(1/2*(i-1))] = y[i]
            y_K_alpha_err[int(1/2*(i-1))] = np.sqrt(y[i])

    #normalize y axis to one here, otherwise all K_alpha and K_beta lines would have about the same intensity and thus the plot would make no sense

    #y_K_alpha = y_K_alpha/max(y_K_alpha)
    #y_K_beta = y_K_beta/max(y_K_beta)

    #for i in range(0, r):
    #    y_K_alpha_err[i] = y_K_alpha[i]*np.sqrt((np.sqrt(max(y_K_alpha))/max(y_K_alpha))**2+(y_K_alpha_err[i]/y_K_alpha[i])**2)
    #    y_K_beta_err[i] = y_K_beta[i]*np.sqrt((np.sqrt(max(y_K_beta))/max(y_K_beta))**2+(y_K_beta_err[i]/y_K_beta[i])**2)


    return y_K_alpha, y_K_alpha_err, y_K_beta, y_K_beta_err

def Duane_Hunt():

    current_altered = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu\Current_altered"
    voltage_altered = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu\Voltage_altered"
    Test = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Test"

    #x_plt, y_plt, y_err, file_numb, line_nums = reader(Test)
    #plot_spectrum(x_plt, y_plt, None)
    #params, params_cov, x_plt, y_plt, y_err, file_numb = fit_spectrum(voltage_altered)
    
    #for n in range(0, file_numb):
    #    plot_spectrum_w_fit(x_plt[n], y_plt[n], None, params[n], n)

    x_plt, y_plt, file_numb, lines, voltage, current = reader(voltage_altered)
    #x_plt1, y_plt1, file_numb1, lines1, voltage1, current1 = reader(current_altered)

    y_K_alpha, y_K_alpha_err, y_K_beta, y_K_beta_err =  find_maxima(x_plt, y_plt, file_numb, lines)
    #y_K_alpha1, y_K_alpha_err1, y_K_beta1, y_K_beta_err1 =  find_maxima(x_plt1, y_plt1, file_numb1, lines1)

    #plot_spectrum(voltage, y_K_alpha, y_K_alpha_err, "k alpha", 1024)
    #plot_spectrum(voltage, y_K_beta, y_K_beta_err, "k beta", 1024)
    #plot_spectrum(current1, y_K_alpha1, y_K_alpha_err1, "k alpha", 2048)
    #plot_spectrum(current1, y_K_beta1, y_K_beta_err1, "k beta", 2048)

    plt.legend()
    plt.show()

    return 0

def main():
    Duane_Hunt()

if __name__ == "__main__" :
    main()