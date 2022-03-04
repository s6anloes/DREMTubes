#!/usr/bin/env python
import DRSiPMEvent
import sys
import ROOT
import getopt
import glob
import os
import re

PathToData='./'
BLUBOLD='\033[94m\033[1m'
BOLD   ='\033[1m'
BLU    ='\033[94m'
RED    ='\033[31m'
YELLOW ='\033[33m'
NOCOLOR='\033[0m'

# CONFIGURATION
NumAdcChannels = 96
NumTdcChannels = 48
ns_TdcCounts   = 0.139063
ns_mm          = 5.333333  
mm_ns          = 1./ns_mm

################################################################
# SIGNAL HANDLER ###############################################
################################################################
from signal import signal, SIGINT
stop=False
def handler(signal_rcv, frame):
  global stop
  stop=True


class DrMonSiPM:
  '''Data monitoring class for DualReadout 2021 Test Beam @H8  '''

  ##### DrMon method #######
  def __init__(self, fname, acqMode, maxEvts, sample, trigCut=0):
    '''Constructor '''
    self.fname      = fname    # File name
    self.acqMode    = acqMode
    self.maxEvts    = maxEvts  # Max number of events to process
    self.sample     = sample   # Sampling fraction
    self.evtDict    = {}       # Dictionary of events
    self.hDict      = {}       # Dictionary of histograms
    self.lastEv     = None     # Last DREvent object
    self.canvas     = None     # ROOT canvas
    self.canNum     = 0        # Number of pads in canvas
    self.numOfLines = 0        # Number of lines in the file
    self.lastLine   = 0        # Last line read
    self.runNum = "0"
    # last number in filename should be run number
    tmp=re.findall('\d+',fname)
    if len(tmp) != 0:
      self.runNum = tmp[-1]
    self.cmdShCuts  = {        # Mapping between command shortCuts and commands
      "all"     : self.DrawAll,
      "board"   : self.DrawBoard,
      "pha"     : self.DrawPHA,
      "trigger" : self.DrawTrigger,
      "help"    : self.PrintHelp,
    }


##### DrMon method #######
  def book1D(self, hname, bins, mi, ma, axTitle="", ymin=None):
    '''Utility to book histograms 1D '''
    h = ROOT.TH1I(hname, hname, bins, mi, ma) 
    h.GetXaxis().SetTitle(axTitle)
    if ymin is not None:
      h.SetMinimum(ymin)
    self.hDict[hname] = h
    return h

##### DrMon method #######
  def bookOthers(self):
    '''Book histograms which were not booked during reading of file '''
    self.book1D( "boardID", 5, -0.5, 4.5, "BoardID", ymin=0)
    self.book1D( "numBoard", 5, 0.5, 5.5, "Num. of Boards per Event")

##### DrMon method #######
  def evtFill(self, evt):
    '''Fill the event dictionary for combining board information into full event'''
    key = str(evt.TriggerID)
    if key in self.evtDict:
      self.evtDict[key].append(evt)
    else:
      self.evtDict[key] = [evt]

##### DrMon method #######
  def checkOverUnderFlow(self, h):
    '''Check and print out over- and underflow of histogram '''
    nbins = h.GetNbinsX()
    hname = h.GetName()
    uFlow = h.GetBinContent(0)
    oFlow = h.GetBinContent(nbins+1)
    if uFlow+oFlow > 0:
      print(BOLD)
      print(hname, ":")
      print(YELLOW)
      print("Underflow:", uFlow)
      print("Overflow :", oFlow)
      print(NOCOLOR)

##### DrMon method #######
  def checkOverUnderFlowFill(self, h, entry, trigID):
    '''Method for filling one entry into the histogram and cheking for under-/overflow '''
    h.Fill(entry)
    nbins = h.GetNbinsX()
    entry_bin = h.FindBin(entry)
    hname = h.GetName()
    if entry_bin == 0:
      print(BOLD, YELLOW, hname+": Underflow", NOCOLOR)
      print(f"triggerID {trigID}; value: {entry}\n")
    elif entry_bin == nbins+1:
      print(BOLD, YELLOW, hname+": Overflow", NOCOLOR)
      print(f"triggerID {trigID}; value: {entry}\n")

##### DrMon method #######
  def hFill(self, event):
    '''Fill the histograms '''
    
    # These are the histograms which can be filled with "board info" only
    # i.e. no need to combine boards into full events
    self.checkOverUnderFlowFill(self.hDict["boardID"], event.BoardID, event.TriggerID)
    self.hDict["triggerID"].Fill(event.TriggerID)
    self.hDict["triggerTimeStamp"].Fill(event.TriggerTimeStamp)

##### DrMon method #######
  def hFillEvent(self):
    '''Fill the histograms with combined event information (after looping and combining event information) '''
    
    lglist = []
    hglist = []

    maxlg = 0
    maxhg = 0

    # one "key" is one triggerID
    for key in self.evtDict:
      lgPhaSum = 0
      hgPhaSum = 0
      
      # one "evt" is the information of one board 
      for evt in self.evtDict[key]:
        self.hFill(evt)
        lgPhaSum += sum(evt.lgPha)
        hgPhaSum += sum(evt.hgPha)
        lglist.append(lgPhaSum)
        hglist.append(hgPhaSum)
        
      if lgPhaSum > maxlg: maxlg = lgPhaSum 
      if hgPhaSum > maxhg: maxhg = hgPhaSum

      # for counting trigID independent of number of boards
      uniqueTrigID = self.evtDict[key][0].TriggerID
      self.hDict["uniqueTrigID"].Fill(uniqueTrigID)
        
      # Number of boards fired in this event
      numboards = len(self.evtDict[key])
      self.checkOverUnderFlowFill(self.hDict["numBoard"], numboards, uniqueTrigID)

    # can only book these histos after the maximal value is known
    self.book1D("lgPhaSum", 4096, 0, maxlg, "lgPha Sum")
    self.book1D("hgPhaSum", 8192, 0, maxhg, "hgPha Sum")
    self.book1D("lgPhaSumZoom", 512, 2000, 5000, "lg Pha Sum Zoom")

    # need to loop over the list again to fill histos
    for lg, hg in zip(lglist, hglist):
      self.hDict["lgPhaSum"].Fill(lg)
      self.hDict["hgPhaSum"].Fill(hg)
      self.hDict["lgPhaSumZoom"].Fill(lg)


##### DrMon method #######
  def readFile(self, offset=0):
    '''Read raw ascii data from file, call the decoding function, fill the histograms'''
    print("Read and parse. Type CTRL+C to interrupt")

    if self.acqMode not in [0,1]: self.acqMode=0
    
    maxtrigid = 0
    maxtrigtime = 0
   
    badboardcounter = 0
    
    with open(self.fname, "rb") as infile:
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
        
    
    self.book1D("triggerID", 50, 0, maxtrigid, "TriggerID", ymin=0)
    self.book1D("uniqueTrigID", 50, 0, maxtrigid, "UniqueTriggerID", ymin=0)
    self.book1D("triggerTimeStamp", 50, 0, maxtrigtime, "TriggerTimeStamp", ymin=0)
    self.hFillEvent()
        

  ##### DrMon method #######
  def createCanvas(self, dim):
    '''Create a canvas'''
    s="Run" + self.runNum + ": "
    if   dim == 4: 
      self.canvas = ROOT.TCanvas('c4', s+'IDEA-DR4', 0, 0, 800, 800)
      self.canvas.Divide(2,2)
    elif dim == 6: 
      self.canvas = ROOT.TCanvas('c6', s+'IDEA-DR6', 0, 0, 1000, 600)
      self.canvas.Divide(3,2)
    elif dim == 2: 
      self.canvas = ROOT.TCanvas('c2', s+'IDEA-DR2', 0, 0, 1000, 500)
      self.canvas.Divide(2,1)
    elif dim == 3: 
      self.canvas = ROOT.TCanvas('c3', s+'IDEA-DR3', 0, 0, 1500, 500)
      self.canvas.Divide(3,1)
    elif dim == 9: 
      self.canvas = ROOT.TCanvas('c9', s+'IDEA-DR9', 0, 0, 900, 900)
      self.canvas.Divide(3,3)
    else: 
      self.canvas = ROOT.TCanvas('c1', s+'IDEA-DR', 0, 0, 800, 600)
    self.canNum = dim

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
  def DrawSingleHisto(self, cmd, opt):
    '''Draw single histogram'''

    if self.canNum != 1:
      self.createCanvas(1)
        
    h = None
    hname = cmd
    if hname in self.hDict:
      h = self.hDict[hname]
    else:
      print(RED, BOLD, 'Unknown histogram', NOCOLOR)
      return
    if opt == "same":
      h.SetFillColor( h.GetFillColor() + 3 )
    h.Draw(opt)
    self.checkOverUnderFlow(h)
    """
    nBins=h.GetNbinsX()
    uFlow = h.GetBinContent(0)
    oFlow = h.GetBinContent(nBins+1)
    if uFlow+oFlow > 0:
      print(BOLD)
      print("Underflow:", uFlow)
      print("Overflow :", oFlow)
      print(NOCOLOR)
    """
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
    for key in self.cmdShCuts:
      print(BLU, key, NOCOLOR)
    print(BLU, "To exit type q", NOCOLOR)

    



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
      elif cmd in self.cmdShCuts:  # SHORTCUTS
        self.cmdShCuts[cmd]()
      elif cmd.isdigit():          # SHORTCUTS WITH DIGITS
        hIdx = int(cmd)
        if hIdx < len(self.cmdShCuts): 
          self.cmdShCuts[ self.cmdShCutsV[hIdx] ]()
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
  print("   -t trigCut    Trigger cut [5=phys, 6=pede]")
  sys.exit(2)

# Parse command line
fname  = ""
events = 999999999
sample = 1
run    = 0
trigCut= 0
try:
  opts, args = getopt.getopt(sys.argv[1:], "hf:a:e:s:r:t:")
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
    elif  o == "-t": trigCut  = int(a)
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

print('Analyzing ', fname)
drMon = DrMonSiPM(fname, acqMode, events, sample, trigCut)
drMon.bookOthers()
drMon.readFile()
drMon.commander()
