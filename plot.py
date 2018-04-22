from TimeSeries import TimeSeries
import Sounding

import datetime

import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import os.path as path
import os

from math import exp
from math import sqrt

################################################################################
#                        Global Configuration Items                            #
################################################################################
target_dirs = ["nam4km_mso", "nam_mso", "gfs3_mso"]

startdate = datetime.datetime(2017,9,1,0)
enddate = datetime.datetime(2017,9,5,0)

plot_width = 11
plot_height = 7

def val_hgts():
    return list(range(100,6000,150))

def lapse_hgts():
    return list(range(50,6100,150))

# Key features, times in UTC
evening_of_27 = datetime.datetime(2017,8,28,0) # large growth period
evening_of_26 = datetime.datetime(2017,8,27,0) # large growth period
evening_of_second = datetime.datetime(2017,9,3,0) # RFW, large growth period
big_run_time = datetime.datetime(2017,9,4,0) # RFW, largest growth burn period

times_of_interest = [
            evening_of_26,
            evening_of_27,
            evening_of_second, 
            big_run_time,
        ]

# Key features, elevations
seeley = 1218     # Seeley Lake
morrell_lo = 2377 # Morrell Lookout

key_elevations = [seeley, morrell_lo]

################################################################################
#                      End Global Configuration Items                          #
################################################################################
def rh(pair):
    t,dp = pair
    return 100*(exp((17.625*dp)/(243.04+dp))/exp((17.625*t)/(243.04+t)))

def make_profile(sounding, attr, hgt_func):
    profile = zip(sounding.hgt, getattr(sounding,attr))
    profile = filter(lambda x: x[0] != -9999.0 and x[1] != -9999.0, profile)
    profile = map(lambda x: (x[0]-sounding.elevation, x[1]), profile)
    snd_hgt, snd_val = zip(*list(profile))

    val = [np.interp(h, snd_hgt, snd_val, -9999.0,-9999.0) for h in hgt_func()]
    hgt = [h+sounding.elevation for h in hgt_func()]

    return (hgt,val)

def local_lapse_rates_agl(sounding):

    hgt, tmp = make_profile(sounding, "temp", lapse_hgts)

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
        asl.append(height)

    return (np.array(lr), np.array(asl))

def helicity(sounding):

    hgt, u_vals = make_profile(sounding, "uWind", lapse_hgts)
    _, v_vals = make_profile(sounding, "vWind", lapse_hgts)

    # Convert to m/s
    u_vals = list(map(lambda x: x * 0.5144, u_vals))
    v_vals = list(map(lambda x: x * 0.5144, v_vals))

    hy = []
    asl = []

    for i in range(1,min(len(v_vals),len(u_vals))):
        h = hgt[i]
        u = u_vals[i]
        v = v_vals[i]

        # h0 = hgt[i-1]
        v0 = v_vals[i-1]
        u0 = u_vals[i-1]

        # dz = h - h0
        du = (u-u0) 
        dv = (v-v0) 

        mu = (u+u0)/2.0
        mv = (v+v0)/2.0

        hlcty =  (mu * dv  - mv * du)
        
        hy.append(hlcty)
        asl.append(h)

    for i in range(len(hy)-1, -1,-1):
        hy[i] = sum(hy[:i])

    return (np.array(hy), np.array(asl))

def sfc_lapse_rates_agl(sounding):
    hgt, tmp = make_profile(sounding, "temp", lapse_hgts)

    h0 = hgt[0]
    t0 = tmp[0]

    lr = []
    asl = []

    for i in range(1,len(tmp)):
        h = hgt[i]
        t = tmp[i]
        lapse_rate = (t-t0)/(h-h0)*1000.0
        lr.append(lapse_rate)
        asl.append(h)

    return (np.array(lr), np.array(asl))

def sfc_to_level_bulk_shear(sounding):
    hgt, u_vals = make_profile(sounding, "uWind", lapse_hgts)
    _, v_vals = make_profile(sounding, "vWind", lapse_hgts)

    # h0 = hgt[0]
    u0 = u_vals[0]
    v0 = v_vals[0]

    bs = []
    asl = []

    for i in range(1,min(len(u_vals),len(v_vals))):
        h = hgt[i]
        u = u_vals[i]
        v = v_vals[i]

        u_shear = u - u0
        v_shear = v - v0

        shear_mag = sqrt(u_shear * u_shear + v_shear * v_shear)

        bs.append(shear_mag)
        asl.append(h)

    return (np.array(bs), np.array(asl))

def omega_agl(sounding):
    asl, omg = make_profile(sounding, "omega", val_hgts)
    omg = list(map(lambda x: x*10.0, omg))
    return (np.array(omg), np.array(asl))

def wind_speed_agl(sounding):
    asl, spd = make_profile(sounding, "windSpd", val_hgts)
    spd = list(map(lambda x: x*1.15708, spd))
    return (np.array(spd), np.array(asl))

def u_wind_agl(sounding):
    asl, u_vals = make_profile(sounding, "uWind", val_hgts)
    u_vals = list(map(lambda x: x*1.15708, u_vals))
    return (np.array(u_vals), np.array(asl))

def v_wind_agl(sounding):
    asl, v_vals = make_profile(sounding, "vWind", val_hgts)
    v_vals = list(map(lambda x: x*1.15708, v_vals))
    return (np.array(v_vals), np.array(asl))
    
def rh_agl(sounding):
    asl, tmp = make_profile(sounding, "temp", val_hgts)
    _, dp = make_profile(sounding, "dewpoint", val_hgts)
    pairs = zip(tmp,dp)
    rh_vals = [rh(p) for p in pairs]

    return (np.array(rh_vals), np.array(asl))

def get_profiles(target_dir):
    # Remember what directory we loaded last, and cache it.
    if "loadedDirectory" not in get_profiles.__dict__ or get_profiles.loadedDirectory != target_dir:
        get_profiles.loadedDirectory = target_dir
    
        # Load the files
        files = os.listdir(target_dir)

        paths = (path.join(target_dir,f) for f in files)
        time_series_list = map(TimeSeries.fromFile, paths)

        p_dict = {}
        for ts in time_series_list:
            for sounding in ts:
                if sounding.time > enddate or sounding.time < startdate:
                    continue
                key = sounding.time.strftime("%Y%m%d%H%M")
                if key not in p_dict.keys():
                    p_dict[key] = sounding
                elif sounding.leadTime < p_dict[key].leadTime:
                    p_dict[key] = sounding

        profiles_list = list(p_dict.values())
        profiles_list.sort(key=lambda x: x.time)

        get_profiles.profiles = profiles_list

    return get_profiles.profiles

def make_value_arrays(target_dir, values_func):
    profiles = get_profiles(target_dir)
    
    vals = np.array([values_func(p)[0] for p in profiles]).transpose()
    levels = values_func(profiles[0])[1]
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
    
    for tm in times_of_interest:
        td = (tm - base_time).total_seconds() / 3600
        if td < 0 or td > max_dt:
            continue
        ax.plot([td,td],[min_elevation,max_elevation], linestyle='dashed', linewidth=2.0, color='black')

    for e in key_elevations:
        if e < min_elevation or e > max_elevation:
            continue
        ax.plot([0,max_dt],[e,e], linestyle='dashed', linewidth=0.5, color='black')

################################################################################
#                           Declaration of plots                               #
################################################################################
plots = [
#       (function_to_plot,        contours,                              filled countours,                      color map name, title,                             colorbar label),
    [
        (sfc_lapse_rates_agl,     [ -9.8, -7, 0],                        [v for v in np.arange(-15.0,25,0.25)], "gist_rainbow", "Surface to * average lapse rate", "C/km"        ),
        (omega_agl,               [v for v in np.arange(-20.0, 21, 5)],  [v for v in np.arange(-20.0,20,0.25)], "RdBu",         "Pressure vertical velocity",      "hPa/s"       ),
        (rh_agl,                  [v for v in np.arange(0.0, 110, 20)],  [v for v in np.arange(0.0,100.0,0.5)], "RdYlGn",       "Relative Humidity",               "%"           ),
    ],
    [
        (u_wind_agl,              [v for v in np.arange(-60.0,60,10.0)], [v for v in np.arange(-60.0,60,1.0)], "RdBu_r",        "U-wind speed",                    "mph"         ),
        (v_wind_agl,              [v for v in np.arange(-60.0,60,10.0)], [v for v in np.arange(-60.0,60,1.0)], "RdBu_r",        "V-wind speed",                    "mph"         ),
        (wind_speed_agl,          [v for v in np.arange(0.0,60,10.0)],   [v for v in np.arange(0.0,60.0,0.5)], "rainbow",       "Wind Speed",                      "mph"         ),
    ],
    [
        (sfc_to_level_bulk_shear, [0.0, 10, 20, 30, 40, 50],             [v for v in np.arange(0.0,50,0.50)],   "rainbow",      "Sfc to * bulk shear",             "mph"         ),
        (local_lapse_rates_agl,   [ -9.8, -7, 0],                        [v for v in np.arange(-15.0,25,0.25)], "gist_rainbow", "Lapse Rate",                      "C/km"        ),
        (helicity,                [h for h in np.arange(-150, 150, 30)], [v for v in np.arange(-150, 150, 3)],  "gist_rainbow", "Helicity",                        "$m^2$/$s^2$" ),
    ],
]

if __name__ == "__main__":

    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    font = {
        'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 10,
        }

    matplotlib.rc('font', **font)

    for target_dir in target_dirs:

        for page_num, page in enumerate(plots):

            # Start the next one
            plt.figure()

            f, axarr = plt.subplots(len(page), sharex=True)

            f.set_size_inches(plot_width,plot_height)
            f.set_dpi(300)

            for row, plot in enumerate(page):

                func, conts, contsf, cmap_name, title, cbar_label = plot
                plot_values(target_dir, f, axarr[row], func, conts, contsf, cmap_name, title, cbar_label)            

            f.text(0.45,0.04,"Midnight MST (GMT-7)",horizontalalignment='center')

            # Save it.
            f.savefig(target_dir + "_" + str(page_num) + ".png")
        