"Mytools2 contains functions which are used to produce plots in the report. It also creates functions that are called in mytools."

"Contents"
"1. Functions called in mytools or that are run for additional analysis"
"2. Functions that produce plots"




import matplotlib.pyplot as plt
import numpy as np


"1.1 Function to define curve-fitting function"
#define fitting function
def quartic_transit(t, a, transit_loc, depth, baseline):
    return (a*(t-transit_loc)**4 + 0.5*depth) - abs((a*(t-transit_loc)**4 + 0.5*depth)) + baseline

"1.2 Function to create density curves"
#Compares to nasa data
def den_curve_cal(density):
    #density in [kg/m3]
    
    mass_val=[]
    min_rad =0.27 * 6371e3 #[m] found from nasa planet table
    max_rad = 60 *6371e3 #[m]
  
    rad= np.linspace(min_rad, max_rad, 250)
    mass= 4/3 * np.pi * (rad)**3 *density
    
    mass_val = mass/5.972e24 #back in earth masses #[M_E]
    rad_val = rad/6371e3 #back in earth rad #[M_R]
    
    return(mass_val, rad_val)





"2.1 Plot full lightcurve"
#full lightcurve
def full_lc(time_sort,flux_sort):
    
    plt.figure(figsize=(12, 6))
    plt.plot(time_sort, flux_sort, ls='None', marker='.', c='black')
    plt.xlabel('Time [days]')
    plt.ylabel('Normalised Flux')
    plt.xlim(100,1600)
    
    return


"2.2 Plot folded lightcurves"
#Plotting Exoplanet Transits using Subplots
def folded(phase_f, flux_f,ph_f, func_f, phase_e, flux_e, ph_e, func_e, phase_l, flux_l, ph_l, func_l,phase_s, flux_s, ph_s, func_s):
    fig, axs = plt.subplots(2,2, figsize= (15,15))
    
    parameters = {'axes.labelsize': 13,'axes.titlesize': 13}
    plt.rcParams.update(parameters)
    
    #planet F
    axs[0, 0].scatter(phase_f, flux_f, marker='.', c='black', label = 'original data')
    axs[0, 0].scatter(ph_f, func_f, ls='-', c='red', label= 'fitted data')
    axs[0, 0].set_xlim(0.885, 0.93)
    axs[0, 0].set_title('A) Planet-F')
    axs[0, 0].legend(loc="upper right")
    
    #planet E
    axs[0, 1].scatter(phase_e, flux_e, marker='.', c='black', label= 'original data')
    axs[0, 1].scatter(ph_e, func_e, ls='-', c='red', label = 'fitted data')
    axs[0, 1].set_xlim(0.695,0.72)
    axs[0, 1].set_title('B) Planet-E')
    axs[0, 1].legend(loc="upper right")
    
    #planet L
    axs[1, 1].scatter(phase_l, flux_l, marker='.', c='black', label='original data')
    axs[1, 1].scatter(ph_l, func_l, ls='-', c='red', label = 'fitted data')
    axs[1, 1].set_xlim(0.265,0.285)
    axs[1, 1].set_title('C) Planet-S')
    axs[1, 1].legend(loc="upper right")
    
    #planet S
    axs[1, 0].scatter(phase_s, flux_s, marker='.', c='black', label = 'original data')
    axs[1, 0].scatter(ph_s, func_s, ls='-', c='red', label='fitted data')
    axs[1, 0].set_xlim(0.197,0.220)
    axs[1, 0].set_title('D) Planet-S')
    axs[1, 0].legend(loc="upper right")
    
    for ax in axs.flat:
        ax.set(xlabel='Phase [days]', ylabel='Normalised Flux')
    
    return


"2.3 Plot Residuals"
def residuals(x_13, y_13, x_21, y_21, x_31, y_31, x_41, y_41):
    
    #plotting residuals
    fig, axs = plt.subplots(2,2, figsize= (15,15))
    
    parameters = {'axes.labelsize': 13,'axes.titlesize': 13}
    plt.rcParams.update(parameters)
    
    #planet F
    axs[0, 0].scatter(x_13,y_13, label='Planet-F', c='black', marker='.')
    axs[0, 0].axhline(y=0, color='r', linestyle='--')
    #axs[0, 0].legend(loc="upper right")
    axs[0, 0].set_title('A) Planet-F')
   
    #planet E
    axs[0, 1].scatter(x_21,y_21, label='Planet-E', c='black', marker='.')
    axs[0, 1].axhline(y=0, color='r', linestyle='--')
    #axs[0, 1].legend(loc="upper right")
    axs[0, 1].set_title('B) Planet-E')
    
    #planet L
    axs[1, 0].scatter(x_31,y_31, label='PLanet-L', c='black', marker='.')
    axs[1, 0].axhline(y=0, color='r', linestyle='--')
    #axs[1, 0].legend(loc="upper right")
    axs[1, 0].set_title('C) Planet-L')
    
    #planet S
    axs[1, 1].scatter(x_41,y_41, label='Planet-S', c='black', marker='.')
    axs[1, 1].axhline(y=0, color='r', linestyle='--')
    #axs[1, 1].legend(loc="upper right")
    axs[1, 1].set_title('D) Planet-S')
    
    for ax in axs.flat:
        ax.set(xlabel='Index', ylabel='Residual Period [days]')
    return


"2.4 Plot orbital length against radius "
#compares with nasa data

def len_rad(orbital_length,planet_radius, nasarad, nasalen):
    plt.figure(figsize=(7, 5))
    plt.scatter(nasalen,nasarad, marker='.',c='black', label='Exoplanet Population')
    plt.scatter(orbital_length[0], planet_radius[0], label='Planet-F', marker='.', c='brown')
    plt.scatter(orbital_length[1], planet_radius[1], label='Planet-E', marker='.', c='orangered')
    plt.scatter(orbital_length[2], planet_radius[2], label='Planet-L', marker='.', c='orange')
    plt.scatter(orbital_length[3], planet_radius[3], label='Planet-S', marker='.', c='darkkhaki')
    
    #plt.annotate('Hot Jupiter', xy=(10**-2, 10**1), xytext=(10**-2, 10**1))

    # Set the axis titles
    plt.ylabel(r'Exoplanet Radius [$R_{earth}$]',fontsize=14)
    plt.xlabel(r'Semi-Major Axis [Au]',fontsize=14)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
    plt.xlim(right =10**1)
    plt.ylim(top = 10**2)
    
    return



"2.5 Plot period against mass"
#Compares to nasa data
def mas_per(planet_mass,planet_period, nasamas, nasaper):
    
    
    plt.figure(figsize=(7, 5))
    plt.scatter(nasaper,nasamas, marker='.', c='black', label='Exoplanet Population') #add colour bar for mass
    plt.scatter(planet_period[0], planet_mass[0], marker='.', c='brown', label= 'Planet-F')
    plt.scatter(planet_period[1], planet_mass[1], marker='.', c='orangered', label= 'Planet-E')
    plt.scatter(planet_period[2], planet_mass[2], marker='.', c='orange', label= 'Planet-L')
    plt.scatter(planet_period[3], planet_mass[3], marker='.', c='blue', label= 'Planet-S')
    
    plt.xlabel(r'Orbital Period [days]',fontsize=14)
    plt.ylabel(r'Mass [$M_{earth}$]',fontsize=14)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
    plt.xlim(right =10**5)

    return


"2.6 Plot mass against radius with density curves plot"
def rad_mas2(planet_mass,planet_radius, earth_val, jup_val, water_val, nep_val):
    
    
    plt.figure(figsize=(5, 8))
    plt.scatter(planet_mass[0],planet_radius[0], marker='o', c='brown', label='Planet-F')
    plt.scatter(planet_mass[1],planet_radius[1], marker='o', c='orangered', label='Planet-E')
    plt.scatter(planet_mass[2],planet_radius[2], marker='o', c='orange', label='Planet-L')
    plt.scatter(planet_mass[3],planet_radius[3], marker='o', c='darkkhaki', label='Planet-S')
    
    plt.plot(earth_val[0], earth_val[1], c='darkgreen', label='Earth-Like')
    plt.plot(jup_val[0], jup_val[1], c='crimson', label='Jupiter-Like')
    plt.plot(water_val[0], water_val[1], c='deepskyblue', label='Water-like')
    plt.plot(nep_val[0], nep_val[1], c='blue', label ='Neptune-Like')
    
    plt.ylabel(r'Exoplanet Radius [$R_{earth}$]',fontsize=14)
    plt.xlabel(r'Mass [$M_{earth}$]',fontsize=14)
    plt.xlim(0,10)
    plt.ylim(0,7)
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
    
    return


"2.7 Plot habitable zone"
def hab_plot(uplim, lolim):
    
    theta = np.linspace(0, 2*np.pi, 100)

    rf= 0.113 #orbit radius of Planet-F
    f_1 = rf*np.cos(theta)
    f_2 = rf*np.sin(theta)

    re= 0.158 #orbit radius of Planet-E
    e_1 = re*np.cos(theta)
    e_2 = re*np.sin(theta)

    rl= 0.203 #orbit radius of Planet-L
    l_1= rl*np.cos(theta)
    l_2 = rl*np.sin(theta)

    rs= 0.240 #orbit radius of Planet-S
    s_1= rs*np.cos(theta)
    s_2 = rs*np.sin(theta)

    ru=uplim
    u_1= ru*np.cos(theta)
    u_2 = ru*np.sin(theta)

    rl= lolim
    lo_1= rl*np.cos(theta)
    lo_2 = rl*np.sin(theta)


    r_star= (1.66 * 696340)*6.68459e-9 #plotting star
    star=plt.Circle((0, 0), r_star, color='black')


    fig, ax = plt.subplots(1, figsize=(10,5))

    ax.plot(f_1, f_2, 'brown', label='Planet-F')
    ax.plot(e_1, e_2, 'orangered', label= 'Planet-E')
    ax.plot(l_1, l_2, 'orange', label= 'Planet-L')
    ax.plot(s_1, s_2, 'darkkhaki', label = 'Planet-S')
    ax.plot(lo_1, lo_2, 'g--')
    ax.plot(u_1, u_2, 'g--')

    ax.add_artist(star)
    plt.annotate('Habitable Zone', xy=(0.8, 0), xytext=(0.9, 0))
    plt.annotate('Exoplanet Orbits', xy=(0.2, 0.4), xytext=(0.03, 0.3))
    ax.set_aspect(1)

    plt.xlim(0,1.8)
    plt.ylim(-0.5,0.5)
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
    ax.set(xlabel='Distance from Star [Au]', ylabel='Distance from Star [Au]')
    plt.show()

    return
