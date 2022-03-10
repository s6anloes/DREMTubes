# Conversion fron ascii data format to DREvent class

class DREvent:
  ''' Class that represent a Dual Readout event at TB 2021 @H8 '''
  
  def __init__(self):
    ''' Constructor '''
    self.EventNumber = 0
    self.NumOfPhysEv = 0
    self.NumOfPedeEv = 0
    self.NumOfSpilEv = 0
    self.TriggerMask = 0
    self.ADCs = {}        # Simple dict      key(channel) : value
    self.TDCs = {}        # Dict key:tuple   key(channel) : ( value, check )

  def headLine(self):
    """Write header in ascii data dump"""
    return "%8s, %8s, %8s, %8s, %10s, %6s, %6s" % (
           "evNum", "#phEv", "#peEv", "#spEv", "trigM",
           "#ADCs", "#TDCs", )
  
  def __str__(self):
    """Overload of str operator (for data dump)"""
    return "%8d, %8d, %8d, %8d, 0x%08x, %6d, %6d" % (
          self.EventNumber, self.NumOfPhysEv, self.NumOfPedeEv, self.NumOfSpilEv, self.TriggerMask,
          len(self.ADCs), len(self.TDCs) )

  def getAdcChannel(self, ch):
    """get the data value for ADC channel ch"""
    return self.ADCs[ch]

  def getTdcChannel(self, ch):
    """get the data value for TDC channel ch"""
    return self.TDCs[ch]  # Tuple: (value, check)

# Parse the evLine and return a DREvent object
def DRdecode(evLine):
  """Function that converts a raw data record (event) from
     ascii to object oriented representation: DREvent class"""
  # Create new DREvent
  e = DREvent()

  # Split sections
  strHeader  = evLine.split(":")[0]
  strPayload = evLine.split(":")[1]
  strPayload = strPayload.replace("TDC", ":TDC")
  strADCs    = strPayload.split(":")[0]
  strTDCs    = strPayload.split(":")[1]
	
  # Parse strHeader
  hList = strHeader.split()
  e.EventNumber = int( hList[2] )
  e.NumOfPhysEv = int( hList[5] )
  e.NumOfPedeEv = int( hList[6] )
  e.NumOfSpilEv = int( hList[7] )
  try:
    e.TriggerMask = int( hList[10], 16 )
  except ValueError:
    print('ERROR: INVALID TRIGGER MASK:', hList[10]) 
    e.TriggerMask = 0xFFFFFFFF 
   
  # Parse ADC 
  listADCs = strADCs.split()
  for i in range(len(listADCs)):
    if i%2 ==0:
      ch = int(listADCs[i]  , 10)
      val= int(listADCs[i+1], 16)
      e.ADCs[ch]=val 

  # Parse TDC  
  entries  = -1
  try: entries  =  int(strTDCs.split()[2], 10)
  except ValueError: 
    # In the 1st runs the TCD size was exadecimal, then it was changed in decimal
    print("WARNING: TCD size not in decimal format, trying exadecimal")
    try: entries  =  int(strTDCs.split()[2], 16)
    except ValueError: 
      print("WARNING: TCD size with unknown format.")
      pass
     

  strTDCs  = strTDCs[ strTDCs.find("val.s") + 6 : ]
  listTDCs = strTDCs.split()
  if entries > 0:
    for i in range(len(listTDCs)):
      if i%3 ==0:
        ch  = int( listTDCs[i+0], 10 )  # Channel
        ver = int( listTDCs[i+1], 10 )  # Varification number
        val = int( listTDCs[i+2], 10 )  # Value
        e.TDCs[ch] = ( val, ver) 
 
  return e


# Main for testing purpose
if __name__ == "__main__":
  import sys
  if len(sys.argv) < 2:
    print("Usage:", sys.argv[0], "filename [v=verbose]")
 
  verbose = False
  if len(sys.argv) == 3:
    verbose = True 
  for i, line in enumerate( open( sys.argv[1] ) ):
    ev = DRdecode(line)
    if verbose:
      if i%30 == 0:
        print(ev.headLine())
      print(ev)
    

  
