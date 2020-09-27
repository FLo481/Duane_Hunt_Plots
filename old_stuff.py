#def Multi_Gauss_fit(x, *p):
#    #[y0_1,a_1,xc_1,sigma,y0_2,a_2,xc_2]
#    return (p[0] + p[1]*np.exp(-(x-p[2])**2/(2*p[3]**2)) + p[4] + p[5]*np.exp(-(x-p[6])**2/(2*p[3]**2)))

#def plot_spectrum_w_fit(x_plt, y_plt, y_err, params, n):

#    plt.figure(n)
#    plt.errorbar(x_plt, y_plt, yerr = y_err, fmt = 'x', markersize = 3 )
#    #plt.title()
#    plt.plot(x_plt, Multi_Gauss_fit(x_plt, *params))

  
#    plt.xlabel("Energy [keV]")
#    plt.ylabel("Intensity [arbitrary units]")

#    plt.grid()

#    return 0

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
