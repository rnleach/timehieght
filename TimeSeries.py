from Sounding import Sounding

class TimeSeries:
    '''
    classdocs
    '''

    def __init__(self, soundingList):
        '''
        Constructor
        '''
        self.data = soundingList
        
    @classmethod
    def fromFile(cls, path):
        
        # A list to store the soundings in
        soundingList = []
        
        # file handle for reading
        f = None
        
        try:
            f = open(path, 'r')
            text = f.read()
        
        except IOError as e:
            print("File Error.")
            raise e
        finally:
            if f and not f.closed:
                f.close()
                    
        soundings = text.split('STN YYMMDD/HHMM')[0]

        # Split at the station id section
        soundings = soundings.split("STID ")
            
        # Add the STID back to the front, and split out
        # the surface data at the end of the file
        for i in range(len(soundings)):
            if soundings[i].startswith("="):
                soundings[i] = "STID " + soundings[i]
            if soundings[i].find('STN YYMMDD/HHMM') > 0:
                shortList = soundings[i].split('STN YYMMDD/HHMM')
                soundings[i] = shortList[0].strip() + "\n"
                soundings.append('STN YYMMDD/HHMM' + shortList[1])
                    
        # Add the relevant sections to the sounding list
        for i in range(len(soundings)):
            if soundings[i].startswith("STID ="):
                snd = Sounding(soundings[i])
                soundingList.append(snd)
                    
        return TimeSeries(soundingList)
    
    def __iter__(self):
        '''This method and the next() method below make this iterable.'''
        self.counter = 0
        return self
    
    def __next__(self): # Python3
        '''Used for iteration in "for x in some_timeseries_data_object"'''
        self.counter = self.counter + 1
        if self.counter > len(self.data):
            raise StopIteration
        return self.data[self.counter - 1]

    def next(self): # Python2
        '''Used for iteration in "for x in some_timeseries_data_object"'''
        self.counter = self.counter + 1
        if self.counter > len(self.data):
            raise StopIteration
        return self.data[self.counter - 1]

    def __getitem__(self, sliced):
        return self.data[sliced]
                
if __name__ == "__main__":
    
    print("Testing TimeSeries: ")
    path = "nam4km_mso/17082600.nam4km_kmso.buf"
    ts = TimeSeries.fromFile(path)

    print("Number of items: %d" % sum(1 for _ in iter(ts)))
    print("Number 3 is %s" % str(ts[2]))
