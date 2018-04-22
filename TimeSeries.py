import re

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

        with open(path, 'r') as f:

            text = f.read()

            soundings = re.finditer(
                r'(^STID(.|\n)*?)(?=^(STID|STN YYMMDD/HHMM))', 
                text,  
                re.MULTILINE)

            soundings = map(lambda snd: snd.group(1), soundings)
            soundings = map(Sounding, soundings)

            return TimeSeries(soundings)
    
    def __iter__(self):
        return self.data.__iter__()

    def __getitem__(self, sliced):
        return self.data[sliced]
                
if __name__ == "__main__":
    
    print("Testing TimeSeries: ")
    path = "nam4km_mso/17082600.nam4km_kmso.buf"
    ts = TimeSeries.fromFile(path)

    print("Number of items: %d" % sum(1 for _ in iter(ts)))
    print("Number 3 is %s" % str(ts[2]))
