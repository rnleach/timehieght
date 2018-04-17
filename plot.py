import TimeSeries as ts
import Sounding

import datetime

import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import os.path as path
import os

from math import exp
from math import sin
from math import cos
from math import pi

################################################################################
#                        Global Configuration Items                            #
################################################################################
startdate = datetime.datetime(2017,9,1,0)
enddate = datetime.datetime(2017,9,5,0)

plot_width = 11
plot_height = 7

def val_hgts():
    return list(range(100,6000,150))

def lapse_hgts():
    return list(range(50,6100,150))

################################################################################
#                      End Global Configuration Items                          #
################################################################################
def rh(pair):
    t,dp = pair
    return 100*(exp((17.625*dp)/(243.04+dp))/exp((17.625*t)/(243.04+t)))

def spd_dir_to_u(pair):
    spd,direct = pair
    return -spd*sin(direct/180.0*pi)

def spd_dir_to_v(pair):
    spd,direct = pair
    return -spd*cos(direct/180.0*pi)

def interp_agl(sounding, tgt_height, attr):
    profile = zip(sounding.hgt, getattr(sounding,attr))
    profile = filter(lambda x: x[0] != -9999.0 and x[1] != -9999.0, profile)
    profile = map(lambda x: (x[0]-sounding.elevation, x[1]), profile)
    hgt, val = zip(*list(profile))

    return np.interp(tgt_height, hgt, val, -9999.0, -9999.0)

def local_lapse_rates_agl(sounding):
    hgt = lapse_hgts()
    tmp = [interp_agl(sounding,h,"temp") for h in hgt]

    lr = []
    asl = []

    for i in range(1,len(tmp)):
        h = hgt[i]
        t = tmp[i]
        h0 = hgt[i-1]
        t0 = tmp[i-1]
        lapse_rate = (t-t0)/(h-h0)*1000.0
        height = (h+h0)/2.0
        lr.append(lapse_rate)
        asl.append(height + sounding.elevation)

    return (np.array(lr), np.array(asl))

def sfc_lapse_rates_agl(sounding):
    hgt = lapse_hgts()
    tmp = [interp_agl(sounding,h,"temp") for h in hgt]

    h0 = hgt[0]
    t0 = tmp[0]

    lr = []
    asl = []

    for i in range(1,len(tmp)):
        h = hgt[i]
        t = tmp[i]
        lapse_rate = (t-t0)/(h-h0)*1000.0
        lr.append(lapse_rate)
        asl.append(h + sounding.elevation)

    return (np.array(lr), np.array(asl))

def omega_agl(sounding):
    hgt = val_hgts()
    omg = [interp_agl(sounding,h,"omega")*10.0 for h in hgt]
    asl = [h+sounding.elevation for h in hgt]

    return (np.array(omg), np.array(asl))

def wind_speed_agl(sounding):
    hgt = val_hgts()
    spd = [interp_agl(sounding,h,"windSpd")*1.15708 for h in hgt]
    asl = [h+sounding.elevation for h in hgt]

    return (np.array(spd), np.array(asl))

def u_wind_agl(sounding):
    hgt = val_hgts()
    spd = [interp_agl(sounding,h,"windSpd")*1.15708 for h in hgt]
    direct = [interp_agl(sounding,h,"windDir") for h in hgt]
    pairs = zip(spd,direct)
    u_vals = [spd_dir_to_u(p) for p in pairs]
    asl = [h+sounding.elevation for h in hgt]

    return (np.array(u_vals), np.array(asl))

def v_wind_agl(sounding):
    hgt = val_hgts()
    spd = [interp_agl(sounding,h,"windSpd")*1.15708 for h in hgt]
    direct = [interp_agl(sounding,h,"windDir") for h in hgt]
    pairs = zip(spd,direct)
    v_vals = [spd_dir_to_v(p) for p in pairs]
    asl = [h+sounding.elevation for h in hgt]

    return (np.array(v_vals), np.array(asl))
    

def rh_agl(sounding):
    hgt = val_hgts()
    tmp = [interp_agl(sounding,h,"temp") for h in hgt]
    dp = [interp_agl(sounding,h,"dewpoint") for h in hgt]
    pairs = zip(tmp,dp)
    rh_vals = [rh(p) for p in pairs]
    asl = [h+sounding.elevation for h in hgt]

    return (np.array(rh_vals), np.array(asl))

def get_profiles(target_dir):
    files = os.listdir(target_dir)

    profiles = [ts.TimeSeries.fromFile(path.join(target_dir,f)) for f in files]
    profiles = [item for sublist in profiles for item in sublist if item.time < enddate and item.time > startdate]
    p_dict = {}
    for p in profiles:
        key = p.time.strftime("%Y%m%d%H%M")
        if key not in p_dict.keys():
            p_dict[key] = p
        elif p.leadTime < p_dict[key].leadTime:
            p_dict[key] = p
    profiles = list(p_dict.values())
    profiles.sort(key=lambda x: x.time)

    return profiles

def make_lapse_rate_arrays(target_dir, lapse_rate_func):
    profiles = get_profiles(target_dir)

    vals = np.array([lapse_rate_func(p)[0] for p in profiles]).transpose()
    levels = local_lapse_rates_agl(profiles[0])[1]
    times = [x.time for x in profiles]
    base_time = times[0]
    dt = [(x.time-base_time).total_seconds()/3600 for x in profiles]

    return(times, dt, levels, vals)

def make_value_arrays(target_dir, values_func):
    profiles = get_profiles(target_dir)
    
    vals = np.array([values_func(p)[0] for p in profiles]).transpose()
    levels = omega_agl(profiles[0])[1]
    times = [x.time for x in profiles]
    base_time = times[0]
    dt = [(x.time-base_time).total_seconds()/3600 for x in profiles]

    return(times, dt, levels, vals)

def setup_y_axis(levels, ax):
    min_levels, max_levels = levels
    
    yticks = [i for i in xrange(0,20000,1000) if i < max_levels and i > min_levels]
    ylabels = [str(i / 1000) for i in yticks]

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylabel("km ASL")

def setup_x_axis(base_time, max_dt, ax):
    lst_shift = -7
    xticks = [i for i in range(0, int(max_dt)) if (i + lst_shift) % 24 == 0]
    xlabels = [base_time + datetime.timedelta(hours=(i+lst_shift)) for i in xticks]
    xlabels = [x.strftime("%m/%d") for x in xlabels]

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

def plot_lapse_rates(target_dir, f, ax, lapse_func, title):
    tms, dt, levels, vals = make_lapse_rate_arrays(target_dir, lapse_func)

    base_time = tms[0]

    X,Y = np.meshgrid(dt,levels)
    conts = [ -9.8, -7, 0]
    contsf = [v for v in np.arange(-15.0,25,0.25)]
    color_map = cm.get_cmap("gist_rainbow")

    # Filled contours
    CF = ax.contourf(X, Y, vals, contsf, cmap=color_map)
    CF.cmap.set_under('white')
    CF.cmap.set_over('white')
    cbar = f.colorbar(CF, ax=ax)

    # Lines
    CS = ax.contour(X, Y, vals,conts, colors='black', linestyles='solid', linewidths=0.5)
    cbar.add_lines(CS)
    cbar.set_ticks(conts)

    # Title and labels
    ax.set_title(title)
    cbar.ax.set_ylabel('C/km')

    setup_x_axis(base_time, max(dt), ax)
    setup_y_axis((min(levels),max(levels)), ax)
    
    # Show interesting times and levels
    plot_key_features(ax, (min(levels),max(levels)), base_time, max(dt))

def plot_values(target_dir, f, ax, val_func, conts, constsf, cmap_name, title, cbar_label):
    tms, dt, levels, vals = make_value_arrays(target_dir, val_func)

    base_time = tms[0]

    X,Y = np.meshgrid(dt,levels)
    color_map = cm.get_cmap(cmap_name)

    CF = ax.contourf(X, Y, vals, contsf, cmap=color_map)
    CF.cmap.set_under('white')
    CF.cmap.set_over('white')
    cbar = f.colorbar(CF, ax=ax)
    cbar.ax.set_ylabel(cbar_label)
    
    CS = ax.contour(X, Y, vals,conts, colors='black', linestyles='solid', linewidths=0.5)
    cbar.add_lines(CS)
    cbar.set_ticks(conts)
    
    ax.set_title(title)

    setup_x_axis(base_time, max(dt), ax)
    setup_y_axis((min(levels),max(levels)), ax)
    
    plot_key_features(ax, (min(levels),max(levels)), base_time, max(dt))

def plot_key_features(ax, elevations, base_time, max_dt):
    min_elevation, max_elevation = elevations
    
    # Big run on the evening of the third
    big_run_time = datetime.datetime(2017,9,4,2) # UTC time
    # The evening of the second, also a red flag warning
    evening_of_second = datetime.datetime(2017,9,3,0) # UTC time
    evening_of_27 = datetime.datetime(2017,8,28,0)
    evening_of_26 = datetime.datetime(2017,8,27,0)

    times_of_interest = [evening_of_second, big_run_time,evening_of_27,evening_of_26]
    for tm in times_of_interest:
        td = (tm - base_time).total_seconds() / 3600
        if td < 0 or td > max_dt:
            continue
        ax.plot([td,td],[min_elevation,max_elevation], linestyle='dashed', linewidth=2.0, color='black')

    # Seeley Lake
    seeley = 1218
    mtns = 2000

    elevations = [seeley, mtns]
    for e in elevations:
        if e < min_elevation or e > max_elevation:
            continue
        ax.plot([0,max_dt],[e,e], linestyle='dashed', linewidth=0.5, color='black')

if __name__ == "__main__":

    target_dirs = ["nam4km_mso", "nam_mso", "gfs3_mso"]
   
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    font = {
        'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 10,
        }

    matplotlib.rc('font', **font)

    for target_dir in target_dirs:

        plt.figure()

        f, axarr = plt.subplots(3, sharex=True)

        f.set_size_inches(plot_width,plot_height)
        f.set_dpi(300)

        # Lapse rate plots
        #plot_lapse_rates(target_dir, f, axarr[0], local_lapse_rates_agl, "Lapse rate")
        plot_lapse_rates(target_dir, f, axarr[0], sfc_lapse_rates_agl,"Surface to * average lapse rate")

        # Omega Plot
        conts = [ -20.0, -15, -10, -5, 0, 5, 10, 15, 20]
        contsf = [v for v in np.arange(-20.0,20,0.25)]
        plot_values(target_dir, f, axarr[1], omega_agl, conts, contsf, "RdBu",'Pressure vertical velocity','hPa/s')

        # RH plot
        conts = [ 0.0, 20.0, 40.0, 60.0, 80.0, 100.0]
        contsf = [v for v in np.arange(0.0,100.0,0.5)]
        plot_values(target_dir, f, axarr[2], rh_agl, conts, contsf, "RdYlGn", "Relative Humidity", "%")

        f.text(0.45,0.05,"Midnight MST (GMT-7)",horizontalalignment='center')

        # Done with first one
        f.savefig(target_dir+"_A.png")

        # Start the next one
        plt.figure()

        f, axarr = plt.subplots(3, sharex=True)

        f.set_size_inches(plot_width,plot_height)
        f.set_dpi(300)

        # U-Wind
        conts = [ -60.0, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60]
        contsf = [v for v in np.arange(-60.0,60,1.0)]
        plot_values(target_dir, f, axarr[0], u_wind_agl, conts, contsf, "RdBu_r"," U-wind speed","mph")

        # V-Wind
        plot_values(target_dir, f, axarr[1], v_wind_agl, conts, contsf, "RdBu_r","V-wind speed","mph")

        # Wind speed Plots
        conts = [ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        contsf = [v for v in np.arange(0.0,60.0,0.5)]
        plot_values(target_dir, f, axarr[2], wind_speed_agl, conts, contsf, "rainbow","Wind speed","mph")

        f.text(0.45,0.04,"Midnight MST (GMT-7)",horizontalalignment='center')

        # Done with first one
        f.savefig(target_dir+"_B.png")
        
        

        
