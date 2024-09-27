"Mytools contains functions which are used for lightcurve analysis and exoplanet calculations"

"Contents"
"1. Lightcurve analysis functions"
"2. Functions to calculate charateristics of exoplanets"
"3. Functions to calculate error"
"4. Other"

#import libraries
import matplotlib.pyplot as plt
import numpy as np
import statistics

from scipy.signal import savgol_filter
from scipy.signal import medfilt
import scipy.signal
from scipy.signal import lombscargle
import scipy
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

import pandas as pd
from CompSample import density_scatter
import seaborn as sns
import pylab
import math
import mytools2



"1) Lightcurve Analysis Functions"

"1.1 Function to run normalising, either SG or median method "
#Option: used to determine which normalising method can be used option = 1 = SG method
def norm(Flux, Time, Error, option, Wd_length, kernel_s):
    
    #remove nans
    missing=np.isfinite(Flux) 
    flux_output =Flux[missing]
    time_output= Time[missing]
    error_output=Error[missing]
  
    
    #If option = 1 run SG else run = Median method
    if option ==1:
        
        interp_savgol = savgol_filter(flux_output, window_length=Wd_length, polyorder=3)
        #normalise
        flux_norm= flux_output/interp_savgol
        error_norm= error_output/interp_savgol
        
    else:
        testfil = medfilt(flux_output, kernel_size = kernel_s)
        #normalise
        flux_norm = flux_output/testfil 
        error_norm = error_output/testfil
    
    
    #using standard deviation to remove outliers
    stdlc=np.std(flux_norm)
    stdmask=(flux_norm < 1+3*stdlc)
    
    #Remove outliers
    time_final = time_output[stdmask]
    flux_final = flux_norm[stdmask]
    error_final = error_norm[stdmask]
    
    #SNR
    SNR=(flux_norm/error_norm)
    meanSNR = statistics.mean(SNR)
    
    return(time_final, flux_final, error_final, meanSNR)


"1.2 Function to run first lomb scargle"
#lombscargle
def lomb1(t_1st, t_dif, t_last, n_points, Time, Flux):
  
    t_2 = t_last - t_1st #time difference between data points 

    #lomb 1
    freqs = np.linspace((1/t_2),(1/t_dif), n_points)
    lomb = scipy.signal.lombscargle(Time, Flux, freqs, precenter=True)
    
     #plot both graphs separatly
    plot1 = plt.figure(1)
    plt.figure(figsize=(6, 4))
    plt.plot(freqs, lomb, c='black', label='Lomb Scargle 1')
    plt.xlabel('Frequency [1/days]')
    plt.ylabel('Power')
    
    return(lomb, freqs, t_2, t_dif)


"1.3 Function to run second lomb scargle" 
def lomb2(lomb, freqs):
    period = np.linspace(1,800, 10000) 
    lomb2 = scipy.signal.lombscargle(freqs, lomb, period, precenter=True)
    
    plot2 = plt.figure(2)
    plt.figure(figsize=(6, 4))
    plt.plot(period, lomb2, label='Lomb Scargle 2', c='black')
    #plt.title("Lomb Scargle 2")
    plt.xlabel('Period [days]')
    plt.ylabel('Power')
    
    return(period, lomb2)



"1.4 Function to fold lightcurve"
def fold_lc(time, flux, error, period):
    
    #Create a pandas dataframe from the inputted data
    data = pd.DataFrame({'time': time, 'flux': flux, 'error': error})
    
    #create the phase 
    data['phase'] = data.apply(lambda x: ((x.time/ period) - np.floor(x.time / period)), axis=1)
    
    phase_long = np.concatenate((data['phase'], data['phase'] + 1.0, data['phase'] + 2.0))
    flux_long = np.concatenate((flux, flux, flux))
    err_long = np.concatenate((error, error, error))
    
    return(phase_long, flux_long, err_long)
             
    
"1.5 Function to find peaks"
def peak_find(y, period):
    
    peaks = find_peaks(y, height=0.03e-12)
    time_peaks = period[peaks[0]] #output peak locations
    
    return(time_peaks)


"1.6 Function to iterate to find exact peak value via inspection"
#A function that runs through a range of exoplanet period values to detect an exoplanet transit
def find_period(max_value, min_value, increments):
    
    period_results = []
    difference= max_value - min_value
    period_test = min_value
    
    for i in range(increments):
        points= difference/increments
        period_test = period_test + points
        period_results.append(period_test)
    
    return(period_results)




"2) Functions to calculate the properties of the exoplanets"

"2.1 Planet Radius"
#formula: R_planet = sqrt(depth) * R_star
def exoplanet_rad(depth):
    star_rad=1.66 * 696340 #in km
    earth_rad=6371 #in km
    
    exo_rad= (np.sqrt(depth*-1)*star_rad)
    #convert to units of earth_mass
    exo_rad_Em = exo_rad/earth_rad
    
    return(exo_rad, exo_rad_Em)


"2.2 Semi Major Axis"
#formula: a = Cbrt( (G*Mtot)/4pi^2) * T^2/3
def orb_len(planet_period, planet_mass):
    
    star_mass=1.1 #[solar masses]
    G = 6.67e-11 #[m3 kg-1 s-2]
    
    #covert days to seconds
    period = planet_period*86400
    #covnert mass into kg
    mass= planet_mass *5.972e24
    #convert star mass into kg
    s_mass = star_mass*1.989e30
    
    pol = np.cbrt(((period**2)*G*(mass+s_mass))/(4*(np.pi)**2))#[m]
    
    #convert planet orbital length into AU
    length= pol*6.6846e-12
    
    return(length)


"2.3 Density"
#density = mass/volume
def density(mass, radius):
    
    #coverting values to SI (currently in earth mass and radius)
    mass_con = mass * 5.972e24 #[kg]
    radius_con = radius * 6371e3 #[m]
    
    volume = 4/3 * np.pi * (radius_con)**3
    
    density = mass_con/volume
    return(density)


"2.4 Escape Velocity"
#vesc = sqrt(2GM/r)
def esc(G, earth_mass, rad_f, rad_e, rad_l, rad_s):
    
    v_f = np.sqrt((2*G*0.390*earth_mass)/(rad_f*1000)) #[m/s]
    v_e = np.sqrt((2*G*4.1*earth_mass)/(rad_e*1000))
    v_l = np.sqrt((2*G*5.5*earth_mass)/(rad_l*1000))
    v_s= np.sqrt((2*G*9.6*earth_mass)/(rad_s*1000))
    return(v_f, v_e, v_l, v_s)


"2.5 Habitable Zone"
#formula: d = sqrt(L / (4pi(T_eff)^4) )
def hab(star_lum, SB):
    
    up_limit = (np.sqrt((star_lum*3.83e26)/(SB*4*np.pi*(373**4))))*6.6846e-12
    lo_limit = (np.sqrt((star_lum*3.83e26)/(SB*4*np.pi*(272**4))))*6.6846e-12
    print('uplimit, lowlimit [Au]:', up_limit, lo_limit)
    return(up_limit, lo_limit)


"2.6 Effective Temperature at Orbital Lengths"
#formula: T= (L / (4pi r^2 SB)) ^ 1/4
def teff(star_lum, SB, orbit_length):
    
    
    L = star_lum * 3.83e26 #converting to W
    OL = orbit_length * 1.496e+11 #converting AU to m
    
    
    teff = ( L / (4*np.pi*(OL**2)*SB) ) **(1/4)
    print('Effective Temp [K]', teff)
    return(teff)

"2.7 Finding RMS of common Greenhouse Gases"
#formaula: vrms = sqrt(3RT/m)
def rms(teff):
    
    Rg = 8.314 #J/mol K
    m_co2 = 0.04401 #kg/mol
    m_h20 = 0.018015 #kg/mol
    #formula: v = sqrt(3RT/m)
    
    rms_co2 = np.sqrt( (3*Rg*teff)/m_co2 )
    rms_h20 = np.sqrt( (3*Rg*teff)/m_h20 )
    
    print('RMS CO2 [m/s]', rms_co2, 'RMS H20 [m/s]', rms_h20)
    return(rms_co2, rms_h20)






"3) Functions to calculate error "

"3.1 Radius Error"
#Formula: R_planet = R_star * sqrt(depth)
#error in sqrt(depth) = 0.5 * (error_depth/ depth) *sqrt(depth)
#error in R_planet = sqrt((error_rad_star/rad_star)^2 + (error_sqrt_depth/sqrt_depth)^2) * exoplanet_radius
def rad_err(optimised_matrix, depth1, rad):
    
    depth = depth1 * -1
    
    error = np.sqrt(np.diag(optimised_matrix))
    error_depth = error[2]
    error_sqrt_depth = 0.5 * (error_depth/depth) * np.sqrt(depth)
    error_rad = np.sqrt( (0.07/1.66)**2 + ( (error_sqrt_depth)/np.sqrt(depth) )**2) * rad #if rad is in [earth rad] the error is
    
    return(error_rad)


"3.2 Semi Major Axis Error"
#formula: a^3 = (P^2 * G(M_star + M_Planet)) / 4pi^2
#error in P^2 => error_p^2/p^2 = (2 * (error_p/p)
#error in (Mtotal) => error_Mtotal = sqrt(error_m_star ^2 + error_m_planet ^2)

#error in a^3 => err_a^3/a^3 = sqrt( err_Mtotal/Mtotal ** 2 + err_p^2/p^2 **2 )
#error in a => err_a / a = 1/3(err_a^3/a^3)

def err_len(er_period, period, er_mass_planet_up, er_mass_planet_low, orb_len, mass_planet):
    
    #Define constants
    mass_star = 1.1 * 1.989e30 #[kg]
    er_mass_star = 0.1 * 1.989e30 #[kg]
    
    Mtot = (mass_star) + (mass_planet*5.972e24) #converting into kg
   
    err_P2 = 2 * (er_period/period) #days
    
    #Upper and Lower limits of error
    error_Mtot_up= np.sqrt((er_mass_star)**2 + (er_mass_planet_up*5.972e24)**2) #converted into kg (planet error)
    error_Mtot_low= np.sqrt((er_mass_star)**2 + (er_mass_planet_low*5.972e24)**2) #converted into kg (planet error)
    
    error_a3_up = np.sqrt( (error_Mtot_up/Mtot)**2 + (err_P2/(period**2))**2) * (orb_len**3)
    error_a_up = 1/3 * (error_a3_up/(orb_len**3))*orb_len
    
    error_a3_low = np.sqrt( (error_Mtot_low/Mtot)**2 + (err_P2/(period**2))**2) * (orb_len**3)
    error_a_low = 1/3 * (error_a3_low/(orb_len**3))*orb_len

    return(error_a_up, error_a_low)



"3.3 Density Error"
#formula: density = mass/((4/3)piR^3)
#err_R^3/R^3 = 3 (err_R/R)

#error in density => error_density/density = sqrt( err_mass/mass **2 + err_R^3/R^3 **2 )
def den_err(err_mass_up, err_mass_lo, err_rad, mass, radius, density):
    
    #error in R^3
    error_R3 = 3 *(err_rad/radius) * radius**3 #[earth radi]
    
    #error in density
    err_den_up = np.sqrt( (err_mass_up/mass)**2 + (error_R3/(radius**3))**2 ) * density
    err_den_lo = np.sqrt( (err_mass_lo/mass)**2 + (error_R3/(radius**3))**2 ) * density

    return(err_den_up, err_den_lo)


"3.4 Escape Velocity Error"
#formula: esc_v = sqrt( (2GM)/R )
#err_esc_v/ esc_v = sqrt( ({err_M^1/2}/{M^1/2})^2 + ({err_R^1/2}/ {R^1/2})^2 )
def err_esc(esc_v, M, R, err_M_up, err_M_lo, err_R):
    
    #error in M^1/2
    #upper
    err_M12_up = 0.5 * (err_M_up/M) * (M**0.5)
    #low
    err_M12_lo = 0.5 * (err_M_lo/M) * (M**0.5)
    
    #error in R^1/2
    err_R12 = 0.5 * (err_R/R) * (R**0.5)
    
    err_v_up = np.sqrt( (err_M12_up/M**0.5)**2 + (err_R12/R**0.5)**2 ) *esc_v
    err_v_lo = np.sqrt( (err_M12_lo/M**0.5)**2 + (err_R12/R**0.5)**2 ) *esc_v
    
    print('escape velocity error:', err_v_up, err_v_lo)
    return(err_v_up, err_v_lo)


"3.5 Period Error"

"3.51 Locating transit locations and applying curvefit"
#residuals method

#Functions to find when the exoplanet periods occur on the original lightcurve
def period_loc(time, flux,error, period,phase, depth, av_period ):
    #time: time from original lightcurve
    #flux flipped: for peak fitting
    #period: period found for the planet
    #phase: phase location of the transit
    #depth: depth of transit, used to create limits
    #av_period: average period to set a peak find limit
    

    full_time = time[-1] #finding the total time length to create a limit for, for loop 
    final_periods = [] #empty array to save results
    
    limit = math.ceil(full_time/period)
    start_point = phase * period #when we will observe the first transit
   
    
    #code that loops round to find peak locations and the upper and lower limit of the time [days], then finds the exact location of the peaks
    for i in range (int(limit)):
        
        loc = start_point +(period *i)
        uplim, lolim = (loc +2, loc -2)
        
        mask_period = (uplim>time) & (time>lolim)
        flux_1 = flux[mask_period]
        time_1 = time[mask_period]
        error_1 = error[mask_period]
        
       
        if len(flux_1) == 0: continue
        if len(flux_1) < 4: continue 
        
      
        para, err = curve_fit(mytools2.quartic_transit, time_1, flux_1, p0=[1e1, loc,depth,1])
        #func = mytools2.quartic_transit(time_1, para[0],para[1],para[2], para[3])
        
        final_periods.append(para[1])
   
        
    return(final_periods)


"3.52 Measuring the time difference between neighbouring transits"
#function to calculate the difference in the period locations to determine the time [days] between each transit (period)

def period_dif(found_transits):
    i=0
    
    output_name = []
    
    for i in range(len(found_transits)):
        dif = found_transits[i] - found_transits[i-1] #difference between neighbouring values
        output_name.append(dif)
        
    
    return(output_name)


"3.53 Function that finds residuals and creates a set of x data to plot residuals"
def clean_periods(periods, true_period):
    #periods: found periods i.e 13.1, 13.2, 13.0, 12.8
    #up: upper limit of expected value i.e 15
    #low: lower limit of expected value i.e 10
    #true period: period found folding the lightcurve i.e 13.170
    
    fin_y = []
    fin_x = []
    x2 = []
    i = 1
    
    up = true_period + 3
    low = true_period -3 
    
    array = np.array(periods) #change data type
    mask_period = (up>array) & (array>low) #upper and lower limit of expected periods, masking anomolies
    y = array[mask_period]
    fin_y = true_period - y
    
    for i in range (len(array)): #making x values for plotting
        x=i
        x2.append(x)
  
    
    x3 = np.array(x2)
    fin_x = x3[mask_period]
    
    return(fin_x, fin_y)


"3.54 Determining Standard Deviation in residual and Error in Period"
def std_period(x, y):
    
    n = len(y)
    sq=[]
    mean_n = np.mean(x)
    snn_ar = []
    
    for i in range (len(y)):
        y_sq = y[i]**2
        sq.append(y_sq)
    
    y_sq_arr = np.array(sq)
    total = np.sum(y_sq_arr)
    
    for i in range(len(x)):
        sum_sq = (x[i] - mean_n)**2
        snn_ar.append(sum_sq)
        
    snn_ar2 = np.array(snn_ar)
    Snn = np.sum(snn_ar2)
    
    err_res = np.sqrt(total/ (n-2)) #err residual = sqrt( sum(measured period - known period)^2 / n-2)
    err_period = err_res/ np.sqrt(Snn)
    print('Period Error [days]:', err_period)
    return(err_period)





"4) Other"
#Function investigate normalising methods, either SG or median method  and plot each lightcurve separately
#Option: used to determine which normalising method can be used option = 1 = SG method
def normtest(Flux, Time, Error, option, Wd_length, kernel_s):
    
    #remove nans
    missing=np.isfinite(Flux) 
    flux_output =Flux[missing]
    time_output= Time[missing]
    error_output=Error[missing]
  
    
    #If option = 1 run SG else run = Median method
    if option == 1: 
        
        interp_savgol = savgol_filter(flux_output, window_length=Wd_length, polyorder=3)
        #normalise
        flux_norm= flux_output/interp_savgol
        error_norm= error_output/interp_savgol
        norm= interp_savgol
        
    else:
        testfil = medfilt(flux_output, kernel_size = kernel_s)
        flux_norm = flux_output/testfil 
        error_norm = error_output/testfil
        norm = testfil
    
    
    #using standard deviation to remove outliers
    stdlc=np.std(flux_norm)
    stdmask=(flux_norm < 1+3*stdlc)
    
    #Remove outliers
    time_final = time_output[stdmask]
    flux_final = flux_norm[stdmask]
    error_final = error_norm[stdmask]
    
    plot1 = plt.figure(1)
    plt.plot(time_final, flux_final)
    plt.xlim(120,168)
    
    plot2 = plt.figure(2)
    plt.plot(time_final, flux_final)
    plt.xlim(170,260)
    
    plot3 = plt.figure(3)
    plt.plot(time_final, flux_final)
    plt.xlim(260,350)
    
    plot4 = plt.figure(4)
    plt.plot(time_final, flux_final)
    plt.xlim(350,445)
    
    plot5 = plt.figure(5)
    plt.plot(time_final, flux_final)
    plt.xlim(445,540)
    
    plot6 = plt.figure(6)
    plt.plot(time_final, flux_final)
    plt.xlim(540,630)
    
    plot7 = plt.figure(7)
    plt.plot(time_final, flux_final)
    plt.xlim(630,720)
    
    plot8 = plt.figure(8)
    plt.plot(time_final, flux_final)
    plt.xlim(720,630)
    
    plot9 = plt.figure(9)
    plt.plot(time_final, flux_final)
    plt.xlim(800,900)
    
    plot91 = plt.figure(91)
    plt.plot(time_final, flux_final)
    plt.xlim(910,1000)
    
    plot10 = plt.figure(10)
    plt.plot(time_final, flux_final)
    plt.xlim(1000,1100)
    
    plot11 = plt.figure(11)
    plt.plot(time_final, flux_final)
    plt.xlim(1175,1275)
    
    plot12 = plt.figure(12)
    plt.plot(time_final, flux_final)
    plt.xlim(1275,1290)
    
    plot13 = plt.figure(13)
    plt.plot(time_final, flux_final)
    plt.xlim(1300,1375)
    
    plot14 = plt.figure(14)
    plt.plot(time_final, flux_final)
    plt.xlim(1550,1600)
    
    return(time_final, flux_final, error_final)
             
         
