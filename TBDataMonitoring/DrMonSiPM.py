#!/usr/bin/env python
import DRSiPMEvent
import sys
import ROOT
import getopt
import glob
import os
import re

import DrMon
from DrMon import BLUBOLD, BOLD, BLU, RED, YELLOW, NOCOLOR
import importlib
importlib.reload(DrMon)

PathToData='./'

################################################################
# SIGNAL HANDLER ###############################################
################################################################
from signal import signal, SIGINT
stop=False
def handler(signal_rcv, frame):
  global stop
  stop=True


class DrMonSiPM(DrMon.DrMon):
  '''Data monitoring class for DualReadout 2021 Test Beam @H8  '''

  ##### DrMon method #######
  def __init__(self, fname, acqMode, maxEvts, sample):
    '''Constructor '''
    self.fname      = fname    # File name
    self.acqMode    = acqMode  # Acquisition mode of FERS boards
    self.maxEvts    = maxEvts  # Max number of events to process
    self.sample     = sample   # Sampling fraction
    self.evtDict    = {}       # Dictionary of events
    self.hDict      = {}       # Dictionary of histograms
    self.runNum = "0"
    # last number in filename should be run number
    tmp=re.findall('\d+',fname)
    if len(tmp) != 0:
      self.runNum = tmp[-1]
    


##### DrMon method #######
  def bookBoardHistos(self):
    '''Book histograms which were not booked during reading of file '''
    self.book1D( "boardID", self.hDict, 5, -0.5, 4.5, "BoardID", ymin=0)
    self.book1D( "numBoard", self.hDict, 5, 0.5, 5.5, "Num. of Boards per Event")

##### DrMon method #######
  def evtFill(self, evt):
    '''Fill the event dictionary for combining board information into full event'''
    key = str(evt.TriggerID)
    if key in self.evtDict:
      self.evtDict[key].append(evt)
    else:
      self.evtDict[key] = [evt]

##### DrMon method #######
  def checkEntryForOverUnderFlow(self, h, entry, trigID):
    '''Method for filling one entry into the histogram and checking for under-/overflow while filling'''
    #h.Fill(entry)
    goodentry = True
    nbins = h.GetNbinsX()
    entry_bin = h.FindBin(entry)
    hname = h.GetName()
    if entry_bin == 0:
      print(BOLD, YELLOW, hname+": Underflow", NOCOLOR)
      print(f"triggerID {trigID}; value: {entry}\n")
      goodentry=False
    elif entry_bin == nbins+1:
      print(BOLD, YELLOW, hname+": Overflow", NOCOLOR)
      print(f"triggerID {trigID}; value: {entry}\n")
      goodentry=False

    return goodentry

##### DrMon method #######
  def hSiPMFill(self, event):
    '''Fill the histograms '''
    
    # These are the histograms which can be filled with "board info" only
    # i.e. no need to combine boards into full events
    self.checkEntryForOverUnderFlow(self.hDict["boardID"], event.BoardID, event.TriggerID)
    self.hDict["boardID"].Fill(event.BoardID)
    self.hDict["triggerID"].Fill(event.TriggerID)
    self.hDict["triggerTimeStamp"].Fill(event.TriggerTimeStamp)

##### DrMon method #######
  def hFillSiPMEvent(self):
    '''Fill the histograms with combined event information (after looping and combining event information) '''
    
    lglist = []
    hglist = []

    maxlg = 0
    maxhg = 0

    badevtcounter = 0

    # one key is one triggerID
    for key in list(self.evtDict.keys()):
      
      # for counting trigID independent of number of boards
      uniqueTrigID = self.evtDict[key][0].TriggerID
        
      # Number of boards fired in this event
      # if numboards exceeds 5 the entire event is skipped for all histograms
      numboards = len(self.evtDict[key])
      dofill = self.checkEntryForOverUnderFlow(self.hDict["numBoard"], numboards, uniqueTrigID)
      if dofill: 
        self.hDict["numBoard"].Fill(numboards)
      else: 
        badevtcounter += 1
        del self.evtDict[key]
        continue

      self.hDict["uniqueTrigID"].Fill(uniqueTrigID)

      lgPhaSum = 0
      hgPhaSum = 0
      # one "evt" is the information of one board 
      for evt in self.evtDict[key]:
        self.hSiPMFill(evt)
        lgPhaSum += sum(evt.lgPha)
        hgPhaSum += sum(evt.hgPha)
        lglist.append(lgPhaSum)
        hglist.append(hgPhaSum)
        
      if lgPhaSum > maxlg: maxlg = lgPhaSum 
      if hgPhaSum > maxhg: maxhg = hgPhaSum

    # can only book these histos after the maximal value is known
    self.book1D("lgPhaSum", self.hDict, 4096, 0, maxlg+1, "lgPha Sum")
    self.book1D("hgPhaSum", self.hDict, 8192, 0, maxhg+1, "hgPha Sum")
    self.book1D("lgPhaSumZoom", self.hDict, 512, 2000, 5000, "lg Pha Sum Zoom")

    # need to loop over the list again to fill histos
    for lg, hg in zip(lglist, hglist):
      self.hDict["lgPhaSum"].Fill(lg)
      self.hDict["hgPhaSum"].Fill(hg)
      self.hDict["lgPhaSumZoom"].Fill(lg)

    print(f"Skipped {badevtcounter} events with more than 5 boards")


##### DrMon method #######
  def readSiPMFile(self, fname=None):
    '''Read raw ascii data from file, call the decoding function, fill the histograms'''

    if fname is None:
        fname = self.fname

    print("Read and parse. Type CTRL+C to interrupt")

    if self.acqMode not in [0,1]: self.acqMode=0
    
    maxtrigid = 0
    maxtrigtime = 0
   
    badboardcounter = 0
    
    with open(fname, "rb") as infile:
      # reading in file header (not used atm)
      data=infile.read(14)
      i = 0
      step = 1
      while data:
        data=infile.read(411)
        if len(data)!=411: 
            print("END OF FILE; "+str(len(data))+" bytes not read")
            print(str(badboardcounter)+" entries with event size not matching")
            continue
        ev, offset = DRSiPMEvent.DRdecode(data, self.acqMode)
        if offset:
          infile.seek(-411, 1)
          data=infile.read(411+offset)
          badboardcounter += 1
          continue
        
        # need to find max values for creating the histograms
        if ev.TriggerID > maxtrigid: maxtrigid = ev.TriggerID
        if ev.TriggerTimeStamp > maxtrigtime: maxtrigtime = ev.TriggerTimeStamp

        # counter for number of entries read
        i +=1
        if i>step*10: step*=10
        if i%step==0: print(i, "\t entries read")
        
        
        # add ev to evtDict for combining into full event information (1 event consists out of (up to) 5 boards which are treated as one event in this loop) 
        self.evtFill(ev)

        if i>=self.maxEvts: break
        
    
    self.book1D("triggerID", self.hDict, 50, 0, maxtrigid+1, "TriggerID", ymin=0)
    self.book1D("uniqueTrigID", self.hDict, 50, 0, maxtrigid+1, "UniqueTriggerID", ymin=0)
    self.book1D("triggerTimeStamp", self.hDict, 50, 0, maxtrigtime+1, "TriggerTimeStamp", ymin=0)
    self.hFillSiPMEvent()
        


  ##### DrMon method #######
  def DrawAll(self):
    '''Draw all event data'''
    if self.canNum != 9:
      self.createCanvas(9)
    self.canvas.cd(4); self.hDict['boardID'].Draw()#; chechOverUnderFlow(self.hDict['boardID'])
    self.canvas.cd(1); self.hDict['numBoard'].Draw()
    self.canvas.cd(6); self.hDict['triggerTimeStamp'].Draw()
    self.canvas.cd(5); self.hDict['triggerID'].Draw()
    self.canvas.cd(2); self.hDict['lgPhaSum'].Draw()
    self.canvas.cd(3); self.hDict['hgPhaSum'].Draw()
    self.canvas.cd(8); self.hDict['uniqueTrigID'].Draw()
    #if self.acqMode == 1:
      #self.canvas.cd(7); self.hDict[''].Draw()
    self.canvas.Update()

  ##### DrMon method #######
  def DrawBoard(self):
    '''Draw the distribution of BoardIDs and the num. of boards per event'''
    if self.canNum != 2:
      self.createCanvas(2)
    self.canvas.cd(1);  self.hDict["boardID"].Draw()
    self.canvas.cd(2);  self.hDict["numBoard"].Draw()
    self.canvas.Update()

  ##### DrMon method #######
  def DrawPHA(self):
    '''Draw lgPhaSum and hgPhaSum distribution'''
    if self.canNum != 2:
      self.createCanvas(2)
    self.canvas.cd(1);  self.hDict["lgPhaSum"].Draw()
    self.canvas.cd(2);  self.hDict["hgPhaSum"].Draw()
    self.canvas.Update()

  ##### DrMon method #######
  def DrawTrigger(self):
    '''Draw triggerID and triggerTimeStamp distribution'''
    if self.canNum != 3:
      self.createCanvas(3)
    self.canvas.cd(1);  self.hDict["triggerID"].Draw()
    self.canvas.cd(2);  self.hDict["uniqueTrigID"].Draw()
    self.canvas.cd(3);  self.hDict["triggerTimeStamp"].Draw()
    self.canvas.Update()
    
  ##### DrMon method #######
  def PrintHelp(self):
    '''Help: shows available commands and histograms'''
    print(BLU, "Available histograms:", NOCOLOR)
    for key in self.hDict:
      print(BLU, key, NOCOLOR)
    print(BLU, "\nenter histogram name to draw single histogram", NOCOLOR)
    print("\n--------------------------------------------------")
    print(BLU, "Available commands:", NOCOLOR)
    for key in self.CMDSHCUTS:
      print(BLU, key, NOCOLOR)
    print(BLU, "To exit type q", NOCOLOR)

  CMDSHCUTS  = {        # Mapping between command shortCuts and commands
      "all"     : DrawAll,
      "board"   : DrawBoard,
      "pha"     : DrawPHA,
      "trigger" : DrawTrigger,
      "help"    : PrintHelp,
    }   

##### DrMon method #######
  def commander(self):
    '''Interactive commander '''
    
    self.createCanvas(1)
    while True:
      # Command
      #self.dumpHelp()
      print("\n__________________________________________________")
      line = input( 'DrMonSiPM prompt --> ' + NOCOLOR)
      pars = line.split()
      if len(pars) == 0: continue
      cmd = pars[0]
      opt = ""
      if len(pars) > 1:
        opt = pars[1]

      if   cmd == "q": print("Bye"); sys.exit(0)
      elif cmd in self.CMDSHCUTS:  # SHORTCUTS
        DrMonSiPM.CMDSHCUTS[cmd](self)
      elif cmd.isdigit():          # SHORTCUTS WITH DIGITS
        hIdx = int(cmd)
        if hIdx < len(self.CMDSHCUTS): 
          self.CMDSHCUTS[ self.cmdShCutsV[hIdx] ]()
      else: 
        self.DrawSingleHisto(cmd, opt)


################################################################
# MAIN #########################################################
################################################################

def Usage():
  print("Read raw data from text file and create monitor histograms")
  print("Usage: DrMon.py [options]")
  print("   -f fname      Data file to analize (def=the latest data file)")
  print("   -a acqMode    Data acquisition mode [0=spectroscopy, 1=spec+timing] (def=0)")
  print("   -e maxEv      Maximum numver of events to be monitored (def=inf)")  
  print("   -s sample     Analyze only one event every 'sample'")
  print("   -r runNbr     Analyze run number runNbr")
  sys.exit(2)


def main(fname, acqMode, events, sample):
  print('Analyzing ', fname)
  drMon = DrMonSiPM(fname, acqMode, events, sample)
  drMon.bookBoardHistos()
  drMon.readSiPMFile()

  return drMon
  

if __name__ == "__main__":
  # Parse command line
  fname   = ""
  acqMode = 0
  events  = 999999999
  sample  = 1
  run     = 0
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hf:a:e:s:r:")
  except getopt.GetoptError as err:
    print(str(err))
    Usage()

  try:
    for o,a in opts:
      if    o == "-h": Usage()
      elif  o == "-f": fname    = a
      elif  o == "-a": acqMode  = int(a)
      elif  o == "-e": events   = int(a)
      elif  o == "-s": sample   = int(a)
      elif  o == "-r": run      = a
  except(AttributeError): pass

  if run > 0:
    fname = f"Run{run}_list.dat"
    fname = PathToData + fname
  elif len(fname) < 1:
    list_of_files = glob.glob( PathToData + 'Run*_list.dat' ) 
    fname = max(list_of_files, key=os.path.getctime)

  if not os.path.isfile(fname):
    print(RED, "[ERROR] File ", fname, " not found", NOCOLOR)
    sys.exit(404)

  # Install signal handler to interrupt file reading
  signal(SIGINT, handler)
  drMon = main(fname, acqMode, events, sample)
  drMon.commander()
