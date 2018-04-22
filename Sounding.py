import re
import datetime as dt

import numpy as np

from math import sin
from math import cos
from math import pi

class Sounding:
    '''
    classdocs
    '''

    def __init__(self, rawText):
        '''
        Constructor
        '''
        
        self.rawText = rawText
        
        # Location Data
        self.stid      = None  # station identifier (or stnm, see below)
        self.stnm      = None  # station name, or id number
        self.time      = None  # valid time
        self.lat       = None  # latitude
        self.lon       = None  # longitude
        self.elevation = None  # Station elevation in feet
        self.leadTime  = None  # Lead time in hours
        
        # Sounding Parameters
        self.show = None       # Showalter index
        self.li   = None       # Lifted Index
        self.swet = None       # SWEAT Index
        self.kinx = None       # K index
        self.lcl  = None       # LCL in mb
        self.pwat = None       # precipitable water in inches
        self.totl = None       # Total Totals index
        self.cape = None       # CAPE J/Kg
        self.lclt = None       # Potential Temperature at LCL, kelvin
        self.cin  = None       # CIN J/Kg
        self.eql  = None       # Equilibrium level, mb
        self.lfc  = None       # Level of free convection, mb
        self.brch = None       # Bulk Richardson Number
        
        # Sounding Data
        self.pressure = None   # Pressure in mb
        self.temp     = None   # Temp in deg C
        self.wbt      = None   # Wet bulb temp in deg C
        self.dewpoint = None   # Dew point in deg C
        self.thetaE   = None   # equivalent potential temp in K
        self.windDir  = None   # wind direction
        self.windSpd  = None   # wind speed in knots
        self.uWind    = None   # west to east wind
        self.vWind    = None   # south to north wind
        self.omega    = None   # vertical velocity in microbars/sec
        self.cloud    = None   # Cloud fraction in percent
        self.hgt      = None   # Height in feet
        
        # Break the sounding into 3 sections
        regEx = re.compile(r"^STID(.|\n)*?^\s*$",  re.MULTILINE)
        match = regEx.search(self.rawText)
        statInfo = ""
        if match is not None:
            statInfo = match.group(0)

        regEx = re.compile(r"^STID(.|\n)*?^\s*$((.|\n)*?)^PRES",  re.MULTILINE)
        match = regEx.search(self.rawText)
        sndParams = ""
        if match is not None:
            sndParams = match.group(2)

        regEx = re.compile(r"^PRES(.|\n)*?^\s*$",  re.MULTILINE)
        match = regEx.search(self.rawText)
        sndData = ""
        if match is not None:
            sndData = match.group(0)
        
        # Get the station information
        regEx = re.compile(r"STID = ([A-Z]{3,4})?\s*(?=STNM)")
        match = regEx.search(statInfo)

        if match.group(1) is not None:
            self.stid = match.group(1)
        
        regEx = re.compile(r"STNM = ([0-9]+)\s+(?=TIME)")
        match = regEx.search(statInfo)
        self.stnm = match.group(1)
        
        if self.stid is None:
            self.stid = self.stnm
            
        regEx = re.compile(r'TIME = ([0-9]{2})([0-9]{2})([0-9]{2})/([0-9]{2})00')
        match = regEx.search(statInfo)
        year = int(match.group(1)) + 2000
        month = int(match.group(2))
        day = int(match.group(3))
        hour = int(match.group(4))
        self.time = dt.datetime(year, month, day, hour)
        
        regEx = re.compile(r"SLAT = (.*)\s+(?=SLON)")
        match = regEx.search(statInfo)
        self.lat = float(match.group(1))
        
        regEx = re.compile(r"SLON = (.*)\s+(?=SELV)")
        match = regEx.search(statInfo)
        self.lon = float(match.group(1))
        
        regEx = re.compile(r"SELV = (.*)\s+(?=STIM)")
        match = regEx.search(statInfo)
        self.elevation = float(match.group(1))
        
        regEx = re.compile(r"STIM = ([0-9]+)\s+")
        match = regEx.search(statInfo)
        self.leadTime = int(match.group(1))
        
        # Get the sounding parameters
        regEx = re.compile(r"SHOW = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.show = float(match.group(1))
            if self.show == -9999.0 or self.show == 9999.0:
                self.show = None
            
        regEx = re.compile(r"LIFT = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.li = float(match.group(1))
            if self.li == -9999.0 or self.li == 9999.0:
                self.li = None
            
        regEx = re.compile(r"SWET = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.swet = float(match.group(1))
            if self.swet == -9999.0 or self.swet == 9999.0:
                self.swet = None
            
        regEx = re.compile(r"KINX = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.kinx = float(match.group(1))
            if self.kinx == -9999.0 or self.kinx == 9999.0:
                self.kinx = None
            
        regEx = re.compile(r"LCLP = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.lcl = float(match.group(1))
            if self.lcl == -9999.0 or self.lcl == 9999.0:
                self.lcl = None
            
        regEx = re.compile(r"PWAT = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.pwat = float(match.group(1)) / 25.4
            
        regEx = re.compile(r"TOTL = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.totl = float(match.group(1))
            if self.totl == -9999.0 or self.totl == 9999.0:
                self.totl = None
            
        regEx = re.compile(r"CAPE = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.cape = float(match.group(1))
            if self.cape == -9999.0 or self.cape == 9999.0:
                self.cape = None
            
        regEx = re.compile(r"LCLT = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.lclt = float(match.group(1))
            if self.lclt == -9999.0 or self.lclt == 9999.0:
                self.lclt = None
            
        regEx = re.compile(r"CINS = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.cin = float(match.group(1))
            if self.cin == -9999.0 or self.cin == 9999.0:
                self.cin = None
            
        regEx = re.compile(r"EQLV = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.eql = float(match.group(1))
            if self.eql == -9999.0 or self.eql == 9999.0:
                self.eql = None
            
        regEx = re.compile(r"LFCT = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.lfc = float(match.group(1))
            if self.lfc == -9999.0 or self.lfc == 9999.0:
                self.lfc = None
            
        regEx = re.compile(r"BRCH = ([-+]?([0-9]*\.[0-9]+|[0-9]+))")
        match = regEx.search(sndParams)
        if match is not None:
            self.brch = float(match.group(1))
            if self.brch == -9999.0 or self.brch == 9999.0:
                self.brch = None

        # Parse the data section
        tokens = sndData.split()
        numCols = sum(1 for token in tokens if token.isalpha())
        headers = tokens[0:numCols]
        data = tokens[numCols:]

        # Make a mapping from the header to the column number 
        # and the list for each element    
        headerMapping = {}
        for i in range(len(headers)):
            headerMapping[headers[i]] = i
        
        # Make list for the present elements
        if "PRES" in headers:
            self.pressure = list(map(float, data[headerMapping["PRES"]::numCols]))
        if "TMPC" in headers:
            self.temp = list(map(float, data[headerMapping["TMPC"]::numCols]))
        if "TMWC" in headers:
            self.wbt = list(map(float, data[headerMapping["TMWC"]::numCols]))
        if "DWPC" in headers:
            self.dewpoint = list(map(float, data[headerMapping["DWPC"]::numCols]))
        if "THTE" in headers:
            self.thetaE = list(map(float, data[headerMapping["THTE"]::numCols]))
        if "DRCT" in headers:
            self.windDir = list(map(float, data[headerMapping["DRCT"]::numCols]))
        if "SKNT" in headers:
            self.windSpd = list(map(float, data[headerMapping["SKNT"]::numCols]))
        if "OMEG" in headers:
            self.omega = list(map(float, data[headerMapping["OMEG"]::numCols]))
        if "CFRL" in headers:
            self.cloud = list(map(float, data[headerMapping["CFRL"]::numCols]))
        if "HGHT" in headers:
            self.hgt = list(map(float, data[headerMapping["HGHT"]::numCols]))

        # Function to translate (spd,dir) into (u,v)
        def spd_dir_to_uv(pair):
            spd,direct = pair
            direct_rad = direct/180.0*pi
            return (-spd*sin(direct_rad), -spd*cos(direct_rad))

        if self.windSpd is not None and self.windDir is not None:
            sd_pairs = zip(self.windSpd, self.windDir)
            uv_pairs = map(spd_dir_to_uv, sd_pairs)
            self.uWind, self.vWind = zip(*list(uv_pairs))
        
        
    def __str__(self):
        toRet = "Sounding:"
        toRet += " STID = " + str(self.stid)
        toRet += " STNM = " + str(self.stnm)
        toRet += " Valid = " + str(self.time)
        toRet += " Lead time = " + str(self.leadTime)
        
        return toRet
