import os
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.signal import find_peaks

def linear_fit(x, *p):

    return p[0] + p[1]*(x-p[2])

def lambda_min_fit(x, y0, a):

    return y0+a*x

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

def plot_spectrum(x_plt, y_plt, y_err, n):

    plt.figure(n)
    figure = 0
    data = plt.errorbar(x_plt, y_plt, yerr = y_err, fmt = 'x', markersize = 3)
    #plt.title()

    if n == 1024:
        params, params_cov = scipy.optimize.curve_fit(three_halfs_fit, x_plt, y_plt, sigma = y_err,  p0 = [-100,0.01,0.01,1], bounds = ([-500,0,0,0.5],[500,100,50,5]), absolute_sigma = True, maxfev = 99999)
        perr = np.sqrt(np.diag(params_cov))/np.sqrt(x_plt.shape[0])
        figure, = plt.plot(x_plt, three_halfs_fit(x_plt, *params))
        print("Voltage dependent data :")
        print("U_K = ", params[2], "+/-", perr[2], "Exponent = ", params[3], "+/-", perr[3])
        plt.xlabel("Voltage [kV]", fontsize = 16)
        plt.ylabel("Intensity [arbitrary units]", fontsize = 16)
    elif n == 2048:
        #params, params_cov = scipy.optimize.curve_fit(three_halfs_fit, x_plt, y_plt, sigma = y_err, p0 = [-100,0.1,0.1,1], bounds = ([-500,0,0.1,1],[4000,np.inf,50,50]), absolute_sigma = True, maxfev = 99999)
        params, params_cov = scipy.optimize.curve_fit(linear_fit, x_plt, y_plt, sigma = y_err, p0 = [-100,0.1,0.1], bounds = ([-100,0,0.1],[4000,np.inf,100]), absolute_sigma = True, maxfev = 99999)
        perr = np.sqrt(np.diag(params_cov))/np.sqrt(x_plt.shape[0])
        figure, = plt.plot(x_plt, linear_fit(x_plt, *params))
        print("Current dependent data :")
        #print("U_K = ", params[2], "+/-", perr[2], "Exponent = ", params[3], "+/-", perr[3])
        print("U_K = ", params[2], "+/-", perr[2], "Slope = ", params[1], "+/-", perr[1])
        plt.xlabel("Current [mA]", fontsize = 16)
        plt.ylabel("Intensity [arbitrary units]", fontsize = 16)
    else:
        plt.xlabel("Energy [$A^°$]", fontsize = 16)
        plt.ylabel("Intensity [arbitrary units]", fontsize = 16)

    plt.grid()

    return data, figure

def find_maxima(x_plt, y_plt, file_numb, lines):

    sum = 0
    numb = 0
    n = 0

    lambda_min = []
    peaks = [0]*2*file_numb
    params = [0]*(file_numb)
    params_cov = [0]*(file_numb)
    ints_of_max = []

    #plot spectra for different voltages or currents and find K_alpha and K_beta maxima

    for i in range(0, file_numb):
        #print(lines[i])
        if i == 0:
            #peak position calculation

            #plot_spectrum(x_plt[1:lines[i]+1], y_plt[1:lines[i]+1], np.sqrt(y_plt[1:lines[i]+1]), i)
            peaks[i], _ = find_peaks(y_plt[1:lines[i]+1], height = 100, threshold=20)
            ints_of_max.append(y_plt[peaks[i][0]+1])
            ints_of_max.append(y_plt[peaks[i][1]+1])
            #print("Peak pos. : ", x_plt[peaks[i]+1], "Intensity at peak pos. :", y_plt[peaks[i]+1])

            #print("File number ", i+1, " : " , y_plt[1:lines[i]+1]) ### correct counting for the indices

            #lambda_min calculation

            for k in range(sum + 1 , sum + lines[i] + 1):
                if 1.2 > x_plt[k] > 1.1:
                    if y_plt[k]-y_plt[k-1] >= 4:
                        lambda_min.append(x_plt[k])
                        break
          
        elif i < file_numb - 1:
            #peak position calculation

            #plot_spectrum(x_plt[sum+1:sum+lines[i]+1], y_plt[sum+1:sum+lines[i]+1], np.sqrt(y_plt[sum+1:sum+lines[i]+1]), i)
            peaks[i], _ = find_peaks(y_plt[sum+1:sum+lines[i]+1], height = 120, threshold=80)
            ints_of_max.append(y_plt[sum+peaks[i][0]])
            ints_of_max.append(y_plt[sum+peaks[i][1]])
            #print("Peak pos. : ", x_plt[sum+peaks[i]+1], "Intensity at peak pos. :", y_plt[sum+peaks[i]+1])

            #print("File number ", i+1, " : " , y_plt[sum+1:sum+lines[i]+1]) ### correct counting for the indices

            #lambda_min calculation

            for k in range(sum + 1 , sum + lines[i] + 1):
                if 1.05 > x_plt[k] > 0.35:
                    if y_plt[k]-y_plt[k-1] >= 10:
                        lambda_min.append(x_plt[k])
                        break

        elif i == file_numb - 1:
            #peak position calculation

            #plot_spectrum(x_plt[sum+1:sum+lines[i]+1], y_plt[sum+1:sum+lines[i]+1], np.sqrt(y_plt[sum+1:sum+lines[i]+1]), i)
            peaks[i], _ = find_peaks(y_plt[sum+1:sum+lines[i]+1], height = 1000, threshold = 1500)
            ints_of_max.append(y_plt[sum+peaks[i][0]])
            ints_of_max.append(y_plt[sum+peaks[i][1]])
            #print("Peak pos. : ", x_plt[sum+peaks[i]+1], "Intensity at peak pos. :", y_plt[sum+peaks[i]+1])

            #print("File number ", i+1, " : " , y_plt[sum+1:sum+lines[i]+1]) ### correct counting for the indices

            #lambda_min calculation

            for k in range(sum + 1, sum + lines[i] + 1):
                if 0.438 > x_plt[k] > 0.35:
                    n += 1
                    #print(y_plt[k])
                if x_plt[k] < 0.35:
                    numb += 1
            params[i], params_cov[i] = scipy.optimize.curve_fit(lambda_min_fit, x_plt[sum+numb:sum+numb+n] , y_plt[sum+numb:sum+numb+n], sigma = np.sqrt( y_plt[sum+numb:sum+numb+n]), absolute_sigma = True)
            #plt.plot(x_plt[sum+numb:sum+numb+n], lambda_min_fit(x_plt[sum+numb:sum+numb+n], params[i][0], params[i][1]))
            lambda_min.append(-lambda_min_fit(0, params[i][0], params[i][1])/params[i][1])
            n = 0
            numb = 0

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

    return y_K_alpha, y_K_alpha_err, y_K_beta, y_K_beta_err, lambda_min

def current_dep_spectrum():

    current_altered = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu\Current_altered"

    x_plt, y_plt, file_numb, lines, voltage, current = reader(current_altered)

    y_K_alpha, y_K_alpha_err, y_K_beta, y_K_beta_err =  find_maxima(x_plt, y_plt, file_numb, lines)

    data, figure = plot_spectrum(current, y_K_alpha, y_K_alpha_err, 2048)
    data1, figure1 = plot_spectrum(current, y_K_beta, y_K_beta_err, 2048)

    plt.legend(handles=[data, data1, figure, figure1],labels=[r"$K_{\alpha}$ lines",r"$K_{\beta}$ lines",r"Linear fit",r"Linear fit"], prop={'size': 12})

    return 0

def voltage_dep_spectrum():

    voltage_altered = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu\Voltage_altered"

    x_plt, y_plt, file_numb, lines, voltage, current = reader(voltage_altered)

    y_K_alpha, y_K_alpha_err, y_K_beta, y_K_beta_err, _ =  find_maxima(x_plt, y_plt, file_numb, lines)

    #data, figure = plot_spectrum(voltage, y_K_alpha, y_K_alpha_err, 1024)
    #data1, figure1 = plot_spectrum(voltage, y_K_beta, y_K_beta_err, 1024)

    #plt.legend(handles=[data, data1, figure, figure1],labels=[r"$K_{\alpha}$ lines",r"$K_{\beta}$ lines",r"$a+b(U_A-c)^{3/2}$ fit",r"$a+b(U_A-c)^{3/2}$ fit"], prop={'size': 12})

    return 0

def lambda_min_plot():

    voltage_altered = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Cu\Voltage_altered"

    x_plt, y_plt, file_numb, lines, voltage, current = reader(voltage_altered)

    y_K_alpha, y_K_alpha_err, y_K_beta, y_K_beta_err, lambda_min =  find_maxima(x_plt, y_plt, file_numb, lines)

    #plt.errorbar(1/voltage,  lambda_min, yerr = [0.5]*len(lambda_min), fmt = 'x', markersize = 3)

    c0 = 299792458
    qe = 1.602176634*10**(-19)

    y_err = []

    for i in range(0,len(lambda_min)):
        y_err.append(0.5*10**(-10))
        lambda_min[i] = lambda_min[i]*10**(-10)
        voltage[i] = voltage[i]*1000

    plt.errorbar(1/(voltage), lambda_min, yerr = y_err, fmt = 'x', markersize = 3)
    params, params_cov = scipy.optimize.curve_fit(lambda_min_fit, 1/voltage , lambda_min, sigma = y_err, absolute_sigma = True)
    plt.plot(1/voltage, lambda_min_fit( 1/voltage, params[0], params[1]))
    perr = np.sqrt(np.diag(params_cov))/np.sqrt(len(voltage))

    print("h = ", params[1]*qe/(c0), "+/-", perr[1]*qe/(c0))

    plt.xlabel(r"$U_A$ [$(kV)^{-1}$]", fontsize = 16)
    plt.ylabel(r"Wavelength [$A^°$]", fontsize = 16)

    return 0

def Duane_Hunt():

    Test = r"C:\Users\Flo\Desktop\F Praktikum\X Ray\Data\2020-09-22\data\Test"

    #x_plt, y_plt, y_err, file_numb, line_nums = reader(Test)

    #voltage_dep_spectrum()
    #current_dep_spectrum()
    lambda_min_plot()

    plt.show()

    return 0

def main():
    Duane_Hunt()

if __name__ == "__main__" :
    main()