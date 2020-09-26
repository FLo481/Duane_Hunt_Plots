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
