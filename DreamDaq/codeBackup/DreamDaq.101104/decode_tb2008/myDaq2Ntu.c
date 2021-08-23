/*************************************************************************

        myDaq2Ntu.c
        -----------

        myDaq --> Ntuple conversion program

        Version 0.1,      A.Cardini C. Deplano 29/11/2001

                          E. Meoni L. La Rotonda 25/07/2006
*********************************************************************/

#include "../myRawFile.h"
// DRS
#include "./MCommon/MagicTypes.h"
#include "./MCommon/MagicRawData.h"
#include <sys/times.h>
#include <TH1F.h>
#include <TROOT.h>
#include <fstream>
// end DRS
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>     
#include <TH1.h> 
#include <TMath.h>

#define TOTFADC 1
#define TOTADC  6
#define TOTSCA  1
#define TOTTDC  1
#define TOTKLOETDC 1
#define TOTTH03  1
#define TOTADCN 2

TTree *tree;
/* Global Variables for NTUPLE (not all used in the ntuple at the moment) */
bool phys;
bool drs;
bool drs_adc;
bool drs_v;
bool drs_leak;
int Nrunnumber;           // Number of the Run  
int BegTimeEvs;             // Time of first event in Seconds
int BegTimeEvu;             // Time of first event in MicroSeconds
int TimeEvs;                // Time of event in Seconds
int TimeEvu;                // Time of event in MicroSeconds
unsigned int Nevtda =0 ;             // Number of data event
unsigned int NSCA;              // Total number of counts in Scaler 
unsigned int CHSCA[16];         // Channels in the scaler (the channel of the count) 
unsigned int  COUNTSCA[16];     //  Counts in the Scaler  (the value of the count)


unsigned int NhitTH03;          // Total number of measurements in the TH03 (always 3) 
unsigned int CHTH03[3];         // Channels in the TH03 (the channel of the measurement) 
                                // (CHTH03[i] is the channel of the i-th measurement) 

UShort_t lowth03;     
UShort_t highth03;     
UShort_t RAWTH03[3];     
                         
float DATATH03[3];              // Measurements in the TH03 (the value of the mesurement)
                                // (DATATH03[i] is the value of the i-th measurement) 


unsigned int NhitOSC;           // Total number of points in the Oscilloscope 
                                // (=1128 if there are 4 channels each with 282 points )
unsigned int CHOSC[1128];       // Channels in the oscilloscope
                                // CHOSC[i] is the channel of the i-th point )  
Char_t DATAOSC[1128];              // Values of the poit in the oscilloscope
// (DATAOSC[i] is the value of the i-th point)
unsigned int NUMOSC[1128];      // Number of the point in the channel (max value =282)
                                
int POSOSC[4];                  // Position for Oscilloscope (as in myOscDataFile.h) 
unsigned int TDOSC[4];          // Time Diff for Oscilloscope (as in myOscDataFile.h)
unsigned int SCALEOSC[4];       // Scale of the Oscilloscope (as in myOscDataFile.h)  
unsigned int CHFOSC[4];         // Oscilloscope Channel Flag (0= channel enabled, 1 =channel abled)  
unsigned int PTSOSC;            // Oscilloscope pts (as in myOscDataFile.h , at the moment =282)  
unsigned int SROSC;             // Sampling rate of the Oscilloscope (as in myOscDataFile.h)
float PED_MEAN_OSC[4];
float PED_RMS_OSC[4];

unsigned int NhitADCN0;         // Hits in ADCN0 (the total number of couts in ADCNO, MODULE 0 of V792AC)  
unsigned int CHADCN0[32];       // Channel in ADCN0 (the channel of cout)
                                // (CHADCN0[i] is the channel of the i-th count )
unsigned int CHARGEADCN0[32] ;  // Charge in ADCN0  (the value of the count)
                                // (CHARGEADCN0[i] is the  value of the i-th cout )
UChar_t OVERADCN0[32] ; 
UChar_t UNDERADCN0[32] ; 


float   PED_MEAN_ADCN0[32];
float   PED_RMS_ADCN0[32];

unsigned int NhitADCN1;         // Hits in ADCN1 (the total number of couts in ADCNO, MODULE1  of V792AC)
unsigned int CHADCN1[32];       // Channel in ADCN1 (the channel of the count)
                                // (CHADCN1[i] is the channel of the i-th count )
unsigned int CHARGEADCN1[32] ;  // Charge in ADCN1  (the value of the count)
                                // (CHARGEADCN1[i] is the  value of the i-th cout )

UChar_t OVERADCN1[32] ; 
UChar_t UNDERADCN1[32] ; 
 
float   PED_MEAN_ADCN1[32];
float   PED_RMS_ADCN1[32];

//unsigned int NhitADC0;          // Hits in ADC0  (the total number of couts in ADCO, MODULE 0 of ADC L1182 )
//unsigned int CHADC0[8];         // Channel in ADC0  (the channel of the count)
                                //(CHADC0[i] is the channel of the i-th count )
//unsigned int CHARGEADC0[8] ;    // Charge in ADC0 (the value of the count)
                                //(CHARGEADC0[i] is the  value of the i-th cout )

//float   PED_MEAN_ADC0[8];
//float   PED_RMS_ADC0[8];
                                                                                                               
unsigned int NhitTDC;           // Hits in TDC  (the total number of couts in TDC, MODULE 1176 )
                                // Remember TDC have 16 channel multihit 
                                //(at the moment used only 8 channels)  
unsigned int CHTDC[64];         // Channel in TDC (the channel of the count)
                                // (CHTDC[i] is the channel of the i-th count)
unsigned int VALIDTDC[64];      //  0=good 1=corrupt
unsigned int EDGETDC[64];       // Edge in TDC (0=falling 1=rising)
unsigned int COUNTTDC[64];      // Count in TDC (the value of the count)
                                //(COUNTTDC[i] is the  value of the i-th cout) 
float X_UP[16];
float Y_UP[16];                                                                        
float X_DW[16];
float Y_DW[16];

unsigned int N_X_UP;
unsigned int N_Y_UP;
unsigned int N_X_DW;
unsigned int N_Y_DW;

UInt_t eventnumDRS=0;
UInt_t TimeSecDRS;
UInt_t TimeSubSecDRS;
unsigned int  drs_flag; //0-is the DRS event is not matched; 1-is the DRS event is matched
UShort_t drs1a[10][1024]; // DRS ADC counts
UShort_t drs2a[10][1024]; // ...
UShort_t drs3a[10][1024]; //      
UShort_t drs4a[10][1024]; //
UShort_t drs5a[10][1024]; //
UShort_t drs6a[10][1024]; //
UShort_t drs7a[10][1024]; //
UShort_t drs8a[10][1024]; //

Short_t drs1v[10][1024]; // DRS ADC counts
Short_t drs2v[10][1024]; // ...
Short_t drs3v[10][1024]; //
Short_t drs4v[10][1024]; //
Short_t drs5v[10][1024]; //
Short_t drs6v[10][1024]; //
Short_t drs7v[10][1024]; //
Short_t drs8v[10][1024]; //



TH1F *hmeanped[80];//istograms with the mean pedestal computed using the second part of the file
TH1F *hrmsped[80];//istograms with the rms pedestal computed using the second part of the file
TH1F *hmeanped1[80];//istograms with the mean pedestal computed using 
                    //the first part of the file (not stored)
TH1F *hrmsped1[80];//istograms with the rms pedestal computed using 
//the first part of the file (not stored)

float ped_drs[80][1024]; // mean pedestal computed using the first part of the file
float rms_drs[80][1024]; // rms pedestal computed using the first part of the file

TH1F *hmeanpedv[80];//istograms with the mean pedestal computed using the second part of the file
TH1F *hrmspedv[80];//istograms with the rms pedestal computed using the second part of the file
TH1F *hmeanped1v[80];//istograms with the mean pedestal computed using
//the first part of the file (not stored)
TH1F *hrmsped1v[80];//istograms with the rms pedestal computed using
//the first part of the file (not stored)

float ped_drsv[80][1024]; // mean pedestal computed using the first part of the file
float rms_drsv[80][1024]; // rms pedestal computed using the first part of the file

float Adc[80][1024][20];
float Intercept[80][1024][20];
float Slope[80][1024][20];

/*
// Variables for Other Modules ADC L1182 present in testbeam 2006  
 
unsigned int NhitADC1;        // Hits in ADC1
unsigned int CHADC1[8];       // Channel in ADC1
unsigned int CHARGEADC1[8];   // Charge in ADC1

unsigned int NhitADC2;        // Hits in ADC2
unsigned int CHADC2[8];       // Channel in ADC2
unsigned int CHARGEADC2[8];   // Charge in ADC2

unsigned int NhitADC3;        // Hits in ADC3
unsigned int CHADC3[8];       // Channel in ADC3
unsigned int CHARGEADC3[8];   // Charge in ADC3
                                                                                
unsigned int NhitADC4;        // Hits in ADC4
unsigned int CHADC4[8];       // Channel in ADC4
unsigned int CHARGEADC4[8];   // Charge in ADC4
                                                                                
unsigned int NhitADC5;        // Hits in ADC5
unsigned int CHADC5[8];       // Channel in ADC5
unsigned int CHARGEADC5[8];   // Charge in ADC5

// Variables for the FLASH ADC present in testbeam 2006                                                                               
unsigned int NhitFADC;        // Hits in FADC
unsigned int CHFADC[128];     // Channel in FADC
float VALUEFADC[128];         // Value in FADC
unsigned int NUMFADC[128];   // Number of sampling in FADC                                                                     

// Variables for KLOE TDC 
int NhitKTDC;               // Hits in KLOE TDC
int CHKTDC[64];             // Channel in KLOE TDC
int OVERKTDC[64];           // Overflow in KLOE TDC
int EDGEKTDC[64];           // Edge in KLOE TDC
float COUNTKTDC[64];        // Count in KLOE TDC
*/

/************************************************************************/
void ntbook()
/************************************************************************/
{                                                                                                                             
  //  TTree::SetMaxTreeSize(1000*Long64_t(2000000000));
  // Declare ntuple identifier and title
  tree = new TTree("DREAM","Title Dream 2008");

                                                                                                 
  // Describe data to be stored in the Ntuple (see comments before) 
 
  tree->Branch("Nrunnumber",&Nrunnumber,"Nrunnumber/i");
  tree->Branch("TimeEvs",&TimeEvs,"TimeEvs/i");
  tree->Branch("TimeEvu",&TimeEvu,"TimeEvu/i");
  tree->Branch("Nevtda",&Nevtda,"Nevtda/i");
 
  tree->Branch("HIGHREFTH03",&highth03,"HIGHREFTH03/s");
  tree->Branch("LOWREFTH03",&lowth03,"LOWREFTH03/s");
  tree->Branch("RAWTH03",RAWTH03,"RAWTH03[3]/s");
  tree->Branch("DATATH03",DATATH03,"DATATH03[3]/F");

  if (!drs) {
    tree->Branch("DATAOSC",DATAOSC,"DATAOSC[1128]/B");
    tree->Branch("CHFOSC",CHFOSC,"CHFOSC[4]/i");
    tree->Branch("POSOSC",POSOSC,"POSOSC[4]/I");
    tree->Branch("SCALEOSC",SCALEOSC,"SCALEOSC[4]/i");
    tree->Branch("TDOSC",TDOSC,"TDOSC[4]/i");
    tree->Branch("PTSOSC",&PTSOSC,"PTSOSC/i");
    tree->Branch("SROSC",&SROSC,"SROSC/i");
    }
    tree->Branch("CHADCN0",CHADCN0,"CHADCN0[32]/i");
    tree->Branch("CHARGEADCN0",CHARGEADCN0,"CHARGEADCN0[32]/i");
    tree->Branch("OVERADCN0",OVERADCN0,"OVERADCN0/b");
    tree->Branch("UNDERADCN0",UNDERADCN0,"UNDERADCN0/b");

   if (!drs){
    tree->Branch("CHADCN1",CHADCN1,"CHADCN1[32]/i");
    tree->Branch("CHARGEADCN1",CHARGEADCN1,"CHARGEADCN1[32]/i");
    tree->Branch("OVERADCN1",OVERADCN1,"OVERADCN1/b");
    tree->Branch("UNDERADCN1",UNDERADCN1,"UNDERADCN1/b");
  }
  tree->Branch("NhitTDC",&NhitTDC,"NhitTDC/i"); 
  tree->Branch("CHTDC",CHTDC,"CHTDC[64]/i");
  tree->Branch("COUNTTDC",COUNTTDC,"COUNTTDC[64]/i");
  tree->Branch("VALIDTDC",VALIDTDC,"VALIDTDC[64]/i");
  tree->Branch("EDGETDC",EDGETDC,"EDGETDC[64]/i");
  
  tree->Branch("X",X_DW,"X[16]/F");
  tree->Branch("N_X",&N_X_DW,"N_X/i");
  tree->Branch("Y",Y_DW,"Y[16]/F");
  tree->Branch("N_Y",&N_Y_DW,"N_Y/i");

  tree->Branch("OPT_X",X_UP,"OPT_X[16]/F");
  tree->Branch("N_OPT_X",&N_X_UP,"N_OPT_X/i");
  tree->Branch("OPT_Y",Y_UP,"OPT_Y[16]/F");
  tree->Branch("N_OPT_Y",&N_Y_UP,"N_OPT_Y/i");
  if (!drs) { 
//TB 2007
//  tree->Branch("CHADC0",CHADC0,"CHADC0[8]/i");
//  tree->Branch("CHARGEADC0",CHARGEADC0,"CHARGEDC0[8]/i");
  } 
  tree->Branch("CHSCA",CHSCA,"CHSCA[16]/i");
  tree->Branch("COUNTSCA",COUNTSCA,"COUNTSCA[16]/i");
  //DRS
  if (drs) {
    tree->Branch("eventnumDRS",&eventnumDRS,"eventnumDRS/i");
    tree->Branch("TimeSecDRS",&TimeSecDRS,"TimeSecDRS/i");
    tree->Branch("TimeSubSecDRS",&TimeSubSecDRS,"TimeSubSecDRS/i");
    if (drs_adc) {
      //      tree->Branch("drs_flag",&drs_flag,"drs_flag/s"); 
      tree->Branch("drs1a",drs1a,"drs1a[10][1024]/s");
      tree->Branch("drs2a",drs2a,"drs2a[10][1024]/s");
      tree->Branch("drs3a",drs3a,"drs3a[10][1024]/s");
      tree->Branch("drs4a",drs4a,"drs4a[10][1024]/s");
      if (drs_leak) {
        tree->Branch("drs5a",drs5a,"drs5a[10][1024]/s");
      // tree->Branch("drs6a",drs6a,"drs6a[10][1024]/s");
      // tree->Branch("drs7a",drs7a,"drs7a[10][1024]/s");
      //  tree->Branch("drs8a",drs8a,"drs8a[10][1024]/s");
        }
    }
    if (drs_v) {
      tree->Branch("drs1v",drs1v,"drs1v[10][1024]/S");
      tree->Branch("drs2v",drs2v,"drs2v[10][1024]/S");
      tree->Branch("drs3v",drs3v,"drs3v[10][1024]/S");
      tree->Branch("drs4v",drs4v,"drs4v[10][1024]/S");
       if (drs_leak){
        tree->Branch("drs5v",drs5v,"drs5v[10][1024]/S");
      // tree->Branch("drs6v",drs6v,"drs6v[10][1024]/S");
      //  tree->Branch("drs7v",drs7v,"drs7v[10][1024]/S");
      //  tree->Branch("drs8v",drs8v,"drs8v[10][1024]/S");
        }
    }
    int k=0;
    for (int i=1;i<9;i++){
      for (int j=0;j<10;j++){
        if (drs_adc){
	  hmeanped[k] = new TH1F(Form("hmeanped_drs%01da_%01d",i,j),Form("hmeanped_drs%01da_%01d",i,j),1024,0,1024);
	  hrmsped[k] =  new TH1F(Form("hrmsped_drs%01da_%01d",i,j),Form("hrmsped_drs%01da_%01d",i,j),1024,0,1024);
        }
        if (drs_v){
	  hmeanpedv[k] = new TH1F(Form("hmeanped_drs%01dv_%01d",i,j),Form("hmeanped_drs%01dv_%01d",i,j),1024,0,1024);
	  hrmspedv[k] =  new TH1F(Form("hrmsped_drs%01dv_%01d",i,j),Form("hrmsped_drs%01dv_%01d",i,j),1024,0,1024);
        }
	if (!phys) {
          if (drs_adc){
	    hmeanped1[k] = new TH1F(Form("hmeanped_drs%01da_%01d_1",i,j),Form("hmeanped_drs%01da_%01d_1",i,j),1024,0,1024);
	    hrmsped1[k] =  new TH1F(Form("hrmsped_drs%01da_%01d_1",i,j),Form("hrmsped_drs%01da_%01d_1",i,j),1024,0,1024);
          }
          if (drs_v){
	    hmeanped1v[k] = new TH1F(Form("hmeanped_drs%01dv_%01d_1",i,j),Form("hmeanped_drs%01dv_%01d_1",i,j),1024,0,1024);
	    hrmsped1v[k] =  new TH1F(Form("hrmsped_drs%01dv_%01d_1",i,j),Form("hrmsped_drs%01dv_%01d_1",i,j),1024,0,1024);
          }
	}
	k++;
      }
    }
  }
  /*
    tree->Branch("CHADC1",CHADC1,"CHADC1[8]/i");
    tree->Branch("CHARGEADC1",CHARGEADC1,"CHARGEDC1[8]/i");
    tree->Branch("CHADC2",CHADC2,"CHADC2[8]/i");
    tree->Branch("CHARGEADC2",CHARGEADC2,"CHARGEDC2[8]/i");
    tree->Branch("CHADC3",CHADC3,"CHADC3[8]/i");
    tree->Branch("CHARGEADC3",CHARGEADC3,"CHARGEDC3[8]/i");
    tree->Branch("CHADC4",CHADC4,"CHADC4[8]/i");
    tree->Branch("CHARGEADC4",CHARGEADC4,"CHARGEDC4[8]/i");
    tree->Branch("CHADC5",CHADC5,"CHADC5[8]/i");
    tree->Branch("CHARGEADC5",CHARGEADC5,"CHARGEDC5[8]/i");

    tree->Branch("CHFADC",CHFADC,"CHFADC[128]/i");
    tree->Branch("VALUEFADC",VALUEFADC,"VALUEFADC[128]/F");
    tree->Branch("NUMFADC",NUMFADC,"NUMFADC[128]/i");
  */
  if (phys) {
     if (!drs) {
// UNtil 2007
//      tree->Branch("PED_MEAN_ADC0",PED_MEAN_ADC0,"PED_MEAN_ADC0[8]/F");
//      tree->Branch("PED_RMS_ADC0",PED_RMS_ADC0,"PED_RMS_ADC0[8]/F");
      }
      tree->Branch("PED_MEAN_ADCN0",PED_MEAN_ADCN0,"PED_MEAN_ADCN0[32]/F");
      tree->Branch("PED_RMS_ADCN0",PED_RMS_ADCN0,"PED_RMS_ADCN0[32]/F");
     if (!drs) {
      tree->Branch("PED_MEAN_ADCN1",PED_MEAN_ADCN1,"PED_MEAN_ADCN1[32]/F");
      tree->Branch("PED_RMS_ADCN1",PED_RMS_ADCN1,"PED_RMS_ADCN1[32]/F");
      tree->Branch("PED_MEAN_OSC",PED_MEAN_OSC,"PED_MEAN_OSC[4]/F");
      tree->Branch("PED_RMS_OSC",PED_RMS_OSC,"PED_RMS_OSC[4]/F");
    }

    // Fill DRS ped istogram in the ntuple
    if (drs) {
      int k=0;
      for (int l=1;l<9;l++) {
	for (int i=0;i<10;i++){
	  for (int j=0;j<1024;j++){
            if (drs_adc){  
	      hmeanped[k]->SetBinContent(j+1,ped_drs[k][j]);
	      hrmsped[k]->SetBinContent(j+1,rms_drs[k][j]);
            }
            if (drs_v){
	      hmeanpedv[k]->SetBinContent(j+1,ped_drsv[k][j]);
	      hrmspedv[k]->SetBinContent(j+1,rms_drsv[k][j]);
            }

	  }
          if (drs_adc){
             if (!drs_leak){
	    if (k<40){ 
	      hmeanped[k]->Write();
	      hrmsped[k]->Write();
	    }
            }
            if (drs_leak) {
             if (k<50) {
             hmeanped[k]->Write();
              hrmsped[k]->Write();
             }
            }
          }
          if (drs_v){
            if (!drs_leak){
	    if (k<40){
	      hmeanpedv[k]->Write();
	      hrmspedv[k]->Write();
	    }
            }
           if (drs_leak) {
              if (k<50){
              hmeanpedv[k]->Write();
              hrmspedv[k]->Write();
              }
           } 
          }
           
	  k++;
	}
      }
    }


  } 
}

/***/
int ReadCalibCoeff(char* drs_cal_coeff_name){
  ifstream ifile;
  int dummy;
  //  char drs_cal_coeff_name[256];

  //  sprintf(drs_cal_coeff_name, "/home/dreamtest/working/drs/NewCalibrationCoeff.txt");
  ifile.open(drs_cal_coeff_name);
  if (ifile.fail()){
    fprintf(stderr, "Cannot open DRS Calibration Coefficient data file %s\n", drs_cal_coeff_name);
    return -1;
  }else{

    for (int iPixel=1;iPixel<=80;iPixel++){
      for (int iCell=0;iCell<1024;iCell++){
	ifile>>dummy>>Adc[iPixel-1][iCell][0];
	for (int iLevel=0;iLevel<20;iLevel++){
	  ifile>>Intercept[iPixel-1][iCell][iLevel]>>Slope[iPixel-1][iCell][iLevel]>>Adc[iPixel-1][iCell][iLevel+1];
	}
      }
    }
  }
  return 0;
}
/***/

/************************************************************************/
int main(int argc, char **argv)
/************************************************************************/
{
  int rc;
  unsigned int j;
  int hits;
  unsigned int buf[10000];
  unsigned int *addr;

  myTRIGNUM trignumData;
  mySCA scaData;
  myADCN adcnData;
  myTH03 th03Data;
  myTEKOSC tekoscData;
  myTDC tdcData;

/*  myADC adcData;
  myFADC fadcData;
    myKTDC ktdcData;
  */

  char datadir[256];
  char datafilename[256];
  char datafilenamedrs[256];
  char ntdir[256];
  char ntfilename[256];
  char drs_cal_coeff_name[256];

  int isFirstEvent = 1;
  //  int istat;
  int carry;

  if(argc<3){
    std::cout << "Usage: " << argv[0] << " runnumber datatype" 
	      << " (pedestal or data)" << std::endl;
    std::cout << "PLEASE set correctly the following environmental variables"<< std::endl;
   
    std::cout << "Input binary files are searched in:"<< std::endl;
    std::cout << "\t $DATADIR/data"<< std::endl;
    std::cout << "\t $DATADIR/pedestal"<< std::endl;
    std::cout << "Input root files for pedestal math are searched in:"<< std::endl;
    std::cout << "\t $NTUPLEDIR/"<< std::endl;
    std::cout << "Output files will be placed in:"<< std::endl;
    std::cout << "\t $NTUPLEDIR/"<< std::endl;
    exit(1);
  }
  
  if(!strcmp(argv[2],"pedestal")){
    phys=false;
  }else if(!strcmp(argv[2],"data")){
    phys=true;
  }else{
    std::cout << "Unknown data type: " << argv[1] << std::endl;
    exit(1);
  }
  drs=true;
  drs_adc=false;
  drs_v=true;
  drs_leak=false;
 
  if (phys) {
    if (argc < 4) {
      drs=true;
    } else if (argc ==4) {
      if (!strcmp(argv[3],"drs")){
	drs=true;
      }else if (!strcmp(argv[3],"nodrs")){
	drs=false;
        drs_adc=false;
        drs_v=false;
      } else if (!strcmp(argv[3],"drs_leak")){
       drs=true;
       drs_leak=true;
      }
    }
    else if (argc ==5) {
      if (!strcmp(argv[3],"drs")){
	drs=true;
      }else if (!strcmp(argv[3],"nodrs")){
	drs=false;
        drs_adc=false;
        drs_v=false;
      }else if (!strcmp(argv[3],"drs_leak")){
       drs=true;
       drs_leak=true;
      }

      if(!strcmp(argv[4],"adc")){
	drs_adc=true;
	drs_v=false;
      } else if (!strcmp(argv[4],"volt")) {
	drs_adc=false;
	drs_v=true;
      }
    }
  }

  if (!phys) {

    if (argc < 4) {

      drs=true;
    } else if (argc ==4) {
      if(!strcmp(argv[3],"all")){
	drs=true;
      }else if (!strcmp(argv[3],"drs")){       
	drs=true;
      } else if (!strcmp(argv[3],"nodrs")){
	drs=false;
      }else if (!strcmp(argv[3],"drs_leak")){
        drs=true;
        drs_leak=true;
      }

    }else if (argc ==5) {
      if (!strcmp(argv[3],"drs")){
	drs=true;
      } else if (!strcmp(argv[3],"nodrs")){
	drs=false;
      } else if (!strcmp(argv[3],"drs_leak")){
        drs=true;
        drs_leak=true;
      }
      if(!strcmp(argv[4],"adc")){
	drs_adc=true;
	drs_v=false;
      } else if (!strcmp(argv[4],"volt")) {
	drs_adc=false;
	drs_v=true;
      }

    } 
  } 
      


  // Create Input/Output filenames
  if (getenv("DATADIR") == NULL){
    if(phys)
      sprintf(datadir, "/home/dreamtest/working");
    else
      sprintf(datadir, "/home/dreamtest/working");
  }else
    sprintf(datadir, "%s", getenv("DATADIR"));
  
  if (getenv("NTUPLEDIR") == NULL)
    sprintf(ntdir, "/home/dreamtest/working/ntuple");
  else
    sprintf(ntdir, "%s", getenv("NTUPLEDIR"));
   
  
 
  if(phys){
    if (!drs){
      sprintf(datafilename, "%s/data/datafile_run%s.dat", datadir, argv[1]);
      sprintf(ntfilename, "%s/datafile_ntup_run%s.root", ntdir, argv[1]);
    } else if (drs) {
      sprintf(datafilename, "%s/data/datafile_run%s.dat", datadir, argv[1]);
      sprintf(ntfilename, "%s/datafile_ntup_drs_run%s.root", ntdir, argv[1]);
    }
  }else{
    if (!drs){
      sprintf(datafilename, "%s/pedestal/pedestal_run%s.dat", datadir, argv[1]);
      sprintf(ntfilename, "%s/pedestal_ntup_run%s.root", ntdir, argv[1]);
    }else if (drs){
      sprintf(datafilename, "%s/pedestal/pedestal_run%s.dat", datadir, argv[1]);
      sprintf(ntfilename, "%s/pedestal_ntup_drs_run%s.root", ntdir, argv[1]);
    } 

  }
  //DRS
  if (drs){
    sprintf(datafilenamedrs, "%s/drs/drs_datafile_run%s.dat", datadir, argv[1]);
    if (drs_v) {
      sprintf(drs_cal_coeff_name, "%s/drs/NewCalibrationCoeff.txt",datadir);
    }
  }
  //

  if(phys){
    // Read pedestal information for the different modules from Pedestal Ntuple
    char pedntfilename[256],channel[50],charge_modul[50];
    float mean_adcn[2][32],rms_adcn[2][32],mean_adc[8],rms_adc[8],mean_osc[4],rms_osc[4];

    for (unsigned int i=0;i<8;i++){
      mean_adc[i]=0;
      rms_adc[i]=0;
      if (i<4) {
	mean_osc[i]=0;
	rms_osc[i]=0;
      } 
      if( i<2) {
	for (unsigned int j=0;j<32;j++){
	  mean_adcn[i][j]=0;
	  rms_adcn[i][j]=0;
	}
      }
    }//for i=0,8
    //    if (drs) { 
    for (unsigned int i=0;i<80;i++){   
      for (unsigned int j=0;j<1024;j++){
	ped_drs[i][j]=0;
	rms_drs[i][j]=0;
	ped_drsv[i][j]=0;
	rms_drsv[i][j]=0;

      }
    }
    //    }
    if (!drs) {  
    sprintf(pedntfilename, "%s/pedestal_ntup_run%s.root", ntdir, argv[1]);
    } else if(drs){
    sprintf(pedntfilename, "%s/pedestal_ntup_drs_run%s.root", ntdir, argv[1]);
    }
     
    TFile *fileped = new TFile(pedntfilename,"READ");
    if (!(fileped->IsOpen())){
      printf("myDaq2Ntu: cannot open pedestal ntuple file %s\n",pedntfilename);
    }//file open
    if(fileped->IsOpen()) {
//      if (!drs) {
	TTree *t=(TTree *)gDirectory->Get("DREAM");
      
	if(t!=NULL && t->GetEntries()>0){
	
	  TH1F *ht=new TH1F("ht","ADC",4096,-0.5,4095.5); 
	  for (int j=0;j<2;j++){
            if ( ((j==0)&&(drs))|| (!drs) ){
	    for (int i=0;i<32;i++) {
	      sprintf(charge_modul,"CHARGEADCN%d",j);
	      sprintf(channel,"CHADCN%d==%d",j,i);
	      if(t->Project("ht",charge_modul,channel)){
		mean_adcn[j][i]= ht->GetMean();
		rms_adcn[j][i] = ht->GetRMS();
	      }else{
		mean_adcn[j][i]=0;
		rms_adcn[j][i]=0;
	      }
	      ht->Reset();
	    }
            }// !drs
	  }
	
	  for (int i=0;i<32;i++) {
	    PED_MEAN_ADCN0[i]=mean_adcn[0][i];
	    PED_RMS_ADCN0[i]=rms_adcn[0][i];
            if (!drs) {
	    PED_MEAN_ADCN1[i]=mean_adcn[1][i];
	    PED_RMS_ADCN1[i]=rms_adcn[1][i];
            }
	  }
          if (!drs) {	
/*	  int j=0;
	  for (int i=0;i<8;i++) {
	    sprintf(charge_modul,"CHARGEADC%d",j);
	    sprintf(channel,"CHADC%d==%d",j,i);
	    if(t->Project("ht",charge_modul,channel)){
	      mean_adc[i]= ht->GetMean();
	      rms_adc[i] = ht->GetRMS();
	    }else{
	      mean_adc[i]=0;
	      rms_adc[i]=0;
	    }
	    ht->Reset();
	  }

	  for (int i=0;i<8;i++) {
	    PED_MEAN_ADC0[i]=mean_adc[i];
	    PED_RMS_ADC0[i]=rms_adc[i];
	  }
*/
	  delete ht;


	  UInt_t ptsosc,chfosc[4];
	  Char_t dataosc[1128];
	  t->SetBranchAddress("PTSOSC",&ptsosc);
	  t->SetBranchAddress("CHFOSC",chfosc);
	  t->SetBranchAddress("DATAOSC",dataosc);
	  t->GetEntry(0);
	
	  unsigned int k[4];
	  for (int i=0;i<4;i++) {
	    if (chfosc[i]==1)
	      k[i]=i;
	    else
	      k[i]=0;
	  }
	
	  char histo[20],name[20];
	  TH1F *htosc[4];
	  for(unsigned int j=0;j<4;j++) {
	    if (chfosc[j]==1){
	      sprintf(histo,"ht%d",j);
	      sprintf(name,"CHFOSC%d",j);
	      htosc[j]=new TH1F(histo,name,256,-128.5,127.5);
	    }
	  } 
	
	  for (unsigned int j=0;j<t->GetEntriesFast();j++){
	    t->GetEntry(j);
	    for (int i=0;i<4;i++) {
	      if (chfosc[i]==1) {
		for (unsigned l=0;l<ptsosc;l++) {
		  htosc[i]->Fill(dataosc[(k[i]*ptsosc)+l]);
		} 
	      } 
	    }
	  }
	
	  for (int i=0;i<4;i++){
	    if (chfosc[i]==1){
	      PED_MEAN_OSC[i]=htosc[i]->GetMean();
	      PED_RMS_OSC[i]=htosc[i]->GetRMS();
	      delete htosc[i];
	    }else{
	      mean_osc[i]=0;
	      rms_osc[i]=0;
	      PED_MEAN_OSC[i]=0;
	      PED_RMS_OSC[i]=0;

	    }
	  }
               }   //(!drs)
	}//entries>0
//      } //(!drs)
      //DRS PED
      if (drs){
	int icha=0;
	for (int icha1=1;icha1<9;icha1++){
	  for (int icha2=0;icha2<10;icha2++){
	    for(int isam=0;isam<1024;isam++){
	      char cmean[30],crms[30];
              if (drs_adc){  
		sprintf(cmean,"hmeanped_drs%01da_%01d",icha1,icha2);
		sprintf(crms,"hrmsped_drs%01da_%01d",icha1,icha2);
		TH1F *hmean = (TH1F*)gDirectory->Get(cmean);
		TH1F *hrms = (TH1F*)gDirectory->Get(crms);
		if (hmean!=0){ 
		  ped_drs[icha][isam]=hmean->GetBinContent(isam+1);
		}
		if (hrms!=0){
		  rms_drs[icha][isam]=hrms->GetBinContent(isam+1);
		}
              }
              if (drs_v){
		sprintf(cmean,"hmeanped_drs%01dv_%01d",icha1,icha2);
		sprintf(crms,"hrmsped_drs%01dv_%01d",icha1,icha2);
		TH1F *hmeanv = (TH1F*)gDirectory->Get(cmean);
		TH1F *hrmsv = (TH1F*)gDirectory->Get(crms);
		if (hmeanv!=0){
		  ped_drsv[icha][isam]=hmeanv->GetBinContent(isam+1);
		}
		if (hrmsv!=0){
		  rms_drsv[icha][isam]=hrmsv->GetBinContent(isam+1);
		}
              }
	    }
	    icha++;
	  } 
	}
      }
      fileped->Close(); 
    }//file open
  }//if phys
 
   // Open Data file
  rc = RawFileOpen(datafilename);
  if (rc) {
    printf("myDaq2Ntu: cannot open data file %s\n", argv[1]);
    return -1;
  }
 
  TFile *file = new TFile(ntfilename,"RECREATE");

  if (!file->IsOpen()) {
    printf("myDaq2Ntu: cannot open ntuple file %s\n", ntfilename);
    return -1;
  }
  if (drs) {
    file->SetCompressionLevel(1);
  }else{
    file->SetCompressionLevel(4);
  }

   
  //DRS- Open and read the file
  FILE *fptr=NULL;
  unsigned int  iChannel;
  unsigned int evt_tot=0;
    
  if (drs) {
    if ((fptr = fopen(datafilenamedrs, "r")) == NULL) {
      printf("myDaq2Ntu: DRS: Could not open DRS file: %s\n",datafilenamedrs);
      return -1;
    }
  }
  run_header_str    *RunHeader = NULL;
  event_header_str   EventHeader;
  channel_header_str ChannelHeader;
  U16 *ChannelData =0;
  int NumSamplesInChannel=0;

  if (drs) {  
    if (!(RunHeader = fReadRunHeader(fptr))) {
      printf("myDaq2Ntu: DRS: Cannot handle format version: %d\n", RunHeader->FormatVersion);
      exit(RunHeader->FormatVersion);
    }
    PrintRunHeader(RunHeader, stdout);

    NumSamplesInChannel = RunHeader->NumPixInChannel * RunHeader->NumSamplesPerPixel;
    ChannelData = (U16 *)malloc(NumSamplesInChannel * NUM_BYTES_PER_SAMPLE);

    U32 nRun = RunHeader->RunNumber;
    U32 nEvts = RunHeader->NumEvents;
    evt_tot =nEvts; 
    std::cout<<"myDaq2Ntu: DRS : Run Number = "<<nRun<<std::endl;
    std::cout<<"myDaq2Ntu: DRS : Number of events = "<<nEvts<<std::endl;
  }
  //
  // Book Ntuple
  ntbook();

  // Some DRS variables
  unsigned int Counter;
  unsigned int Triggnumb;

  int evt_DRS=0; //number of the matched DRS events
  int evt_DRS_fp=0; //number of the matched DRS events in the first part of the file
  int GoodEvts[80][1024]; //number of the matched DRS events in the second part of the file
  int GoodEvtsv[80][1024];
  double    muv[80][1024]; // ped mean computed using the first part of the file
  double    sigv[80][1024]; // ped rms computed using the first part of the file
  double    muvv[80][1024]; 
  double    sigvv[80][1024];
  bool check_drs_half=true;
  bool check_drs_halfv=true;

  for (int i=0;i<80;i++) {
    for (int j=0;j<1024;j++) {
      GoodEvts[i][j]=0;
      GoodEvtsv[i][j]=0;  
    }
  }
  if ((drs)&&(drs_v))  {
    int r=ReadCalibCoeff(drs_cal_coeff_name);
    if (r==-1) {return r;}
  }
  // Main Loop on all events
  for (;;) {
    Counter=0;
    Triggnumb=0;
    // Read Event and write it in buffer 'buf'
    rc = RawFileReadEventData(buf);
    if (rc == RAWDATAEOF) { 
      printf("Found EOF at event %d\n", GetEventNumber()); 
      break; }
    if (rc == RAWDATAUNEXPECTEDEOF) { printf("Unexpected EOF at event %d\n", GetEventNumber()); break; }
    rc = 0;  // 
    Nrunnumber=GetRunNumber();
    //    int Nevttot=GetTotEvts();
    // Write in the Ntuple  the Event Number and Time (seconds and uSec) w.r.t first event.
    if (isFirstEvent) {
      BegTimeEvs = GetEventTimes();
      BegTimeEvu = GetEventTimeu();
      isFirstEvent =0;
    }
    Nevtda   = GetEventNumber();
//    std::cout<<" Nevt "<<Nevtda<<std::endl;
    TimeEvu  = GetEventTimeu() - BegTimeEvu;
    if (TimeEvu<0.) {
      double TimeEvu_t= 1e6 + GetEventTimeu() - BegTimeEvu; 
      TimeEvu = (int) TimeEvu_t; 
      carry = 1; 
    }
    else 
      carry = 0;
    
    TimeEvs  = GetEventTimes() - BegTimeEvs - carry;


    //    printf("while fillng time info: Event %d TimeEvs %d TimeEvu %d\n", Nevtda, TimeEvs, TimeEvu);


    //  Decode module and compute the variables of the NTUPLE
    //

// FILL EVENTNUMER
  if (drs){
   hits = 0;
    addr = SubEventSeek(0x100023, buf);
    rc += DecodeTriggerNumber(addr, &trignumData);
      hits=trignumData.index;
      Counter=trignumData.counter;
      Triggnumb=trignumData.triggnumb;
/*
     if (Counter!=Triggnumb) {
      std::cout<<" myDaq2Ntu : Counter and Trigger number different! Counter="<<Counter<<" Triggnumb="<<Triggnumb<<std::endl;  return -1;}
     if (Nevtda!= Counter) {
     std::cout<<" myDaq2Ntu : Nevtda(GetEventNumber) and Trigger number different! Nevtda="<<Nevtda<<" Couner="<<Counter<<std::endl;  return -1;}
*/     if ((Nevtda+1)!= Triggnumb) {
 std::cout<<" myDaq2Ntu : Nevtda(GetEventNumber) and Trigger number different! Nevtda="<<Nevtda<<" Triggnumb="<<Triggnumb<<std::endl;  return -1;}

// std::cout<<" myDaq2Ntu NEW EVENT : Nevtda+1= "<<Nevtda+1<< " Triggnumb= "<<Triggnumb<<" Counter= "<<Counter<<std::endl;
 }
   // FILL SCALER 260
    for (int i=0;i<16;i++) {
      CHSCA[i] = 100;
    }

    hits = 0;
    addr = SubEventSeek(0x200003, buf);
    rc += DecodeV260(addr, &scaData);
    for (j=0; j<scaData.index; j++) {
      CHSCA[hits]  = scaData.channel[j];
      COUNTSCA[hits] = scaData.counts[j];
      hits++;
    }
    NSCA = hits;

    // FILL TH03
    hits = 0;
    addr =  SubEventSeek(0x0000ffff, buf);
    rc += DecodeTH03(addr, &th03Data);

    for (j=0; j<th03Data.index; j++) {
      CHTH03[hits]  = th03Data.channel[j];
      RAWTH03[hits]= th03Data.rawdata[j];
      DATATH03[hits] = th03Data.data[j];
      hits++;
    }
    lowth03=th03Data.lowref;
    highth03=th03Data.highref;
    NhitTH03 = hits;

//    if (!drs) {
      // FILL ADC NEW 0
      for (int i=0;i<32;i++){     
	CHADCN0[i] = 100;
	CHADCN1[i] = 100; 
      }
 
      hits = 0;
      addr =  SubEventSeek(0x04000005, buf);
      rc += DecodeV792AC(addr, &adcnData);
      
      for (j=0; j<adcnData.index; j++) {
	CHADCN0[hits]  = adcnData.channel[j];
	CHARGEADCN0[hits] = adcnData.data[j];
	OVERADCN0[hits] = adcnData.ov[j];
	UNDERADCN0[hits] = adcnData.un[j];
	hits++;
      }
      NhitADCN0 = hits;
     
      // FILL ADC NEW 1
     if (!drs){ 
     hits = 0;
      addr =  SubEventSeek(0x06000005, buf);
      rc += DecodeV792AC(addr, &adcnData);
                                                                                                 
      for (j=0; j<adcnData.index; j++) {
	CHADCN1[hits]  = adcnData.channel[j];
	CHARGEADCN1[hits] = adcnData.data[j];
	OVERADCN1[hits] = adcnData.ov[j];
	UNDERADCN1[hits] = adcnData.un[j];
	hits++;
      }
      NhitADCN1 = hits;
                                                                                                 
                                                                                           


      // FILL TEKOSC

      hits = 0;
      addr =  SubEventSeek(0x0000fafa, buf);
      rc += DecodeTEKOSC(addr, &tekoscData);
      unsigned int kp;
      for (j=0; j<tekoscData.index; j++) {
	CHOSC[hits]  = tekoscData.channel[j];
	DATAOSC[hits] = (Char_t)(tekoscData.data[j]/256);
	NUMOSC[hits] =tekoscData.num[j];
      
	hits++;
      }
      NhitOSC = hits;

      for (kp=0;kp<4;kp++) {
	SCALEOSC[kp]=tekoscData.scale[kp];
	POSOSC[kp]=tekoscData.position[kp];
	TDOSC[kp]=tekoscData.tdiff[kp];
	CHFOSC[kp]=tekoscData.chfla[kp];
      }
      PTSOSC=tekoscData.pts;
      SROSC=tekoscData.samplerate;
    }
//Until tb 2007 
/*    //FILL ADC LECROY L1182 0
    for (int i=0;i<8;i++) {
      CHADC0[i]  = 100;
    } 

    hits = 0;
    addr = SubEventSeek(0x10002, buf);
    rc += DecodeL1182(addr, &adcData);
    for (j=0; j<adcData.index; j++) {
      CHADC0[hits]  = adcData.channel[j];
      CHARGEADC0[hits]  = adcData.charge[j];
      hits++;
    }
    NhitADC0 = hits;
*/

    /*    //FILL ADC LECROY L1182 1
	  hits = 0;
	  addr = SubEventSeek(0x50002, buf);
	  rc += DecodeL1182(addr, &adcData);
	  for (j=0; j<adcData.index; j++) {
	  CHADC1[hits]  = adcData.channel[j];
	  CHARGEADC1[hits]  = adcData.charge[j];
	  hits++;
	  }
	  NhitADC1 = hits;
	  //FILL ADC LECROY L1182 2
	  hits = 0;
	  addr = SubEventSeek(0x40002, buf);
	  rc += DecodeL1182(addr, &adcData);
	  for (j=0; j<adcData.index; j++) {
	  CHADC2[hits]  = adcData.channel[j];
	  CHARGEADC2[hits]  = adcData.charge[j];
	  hits++;
	  }
	  NhitADC2 = hits;
	  //FILL ADC LECROY L1182 3
	  hits = 0;
	  addr = SubEventSeek(0x30002, buf);
	  rc += DecodeL1182(addr, &adcData);
	  for (j=0; j<adcData.index; j++) {
	  CHADC3[hits]  = adcData.channel[j];
	  CHARGEADC3[hits]  = adcData.charge[j];
	  hits++;
	  }
	  NhitADC3 = hits;
	  //FILL ADC LECROY L1182 4
	  hits = 0;
	  addr = SubEventSeek(0x20002, buf);
	  rc += DecodeL1182(addr, &adcData);
	  for (j=0; j<adcData.index; j++) {
	  CHADC4[hits]  = adcData.channel[j];
	  CHARGEADC4[hits]  = adcData.charge[j];
	  hits++;
	  }
	  NhitADC4 = hits;
	  //FILL ADC LECROY L1182 5
	  hits = 0;
	  addr = SubEventSeek(0x10002, buf);
	  rc += DecodeL1182(addr, &adcData);
	  for (j=0; j<adcData.index; j++) {
	  CHADC5[hits]  = adcData.channel[j];
	  CHARGEADC5[hits]  = adcData.charge[j];
	  hits++;
	  }
	  NhitADC5 = hits;
    */
    //FILL FADC SIS3320
    /*
      hits = 0;
      addr = SubEventSeek(0x20000013, buf);
      rc += DecodeSIS3320(addr, &fadcData);
      for (j=0; j<fadcData.index; j++) {
      CHFADC[hits]  = fadcData.channel[j];
      VALUEFADC[hits]  = fadcData.value[j];
      NUMFADC[hits]  = fadcData.num[j];
      hits++;
      }
      NhitFADC = hits;
    */
    //FILL TDC L1176
    for (int i=0;i<64;i++) {
      CHTDC[i] = 100;
    }

    hits = 0;
    addr = SubEventSeek(0x300020, buf);
    rc += DecodeL1176(addr, &tdcData);
//     printf("TDC TDC index= %d",tdcData.index);
    for (j=0; j<tdcData.index; j++) { 
      COUNTTDC[hits]  = tdcData.data[j];
      EDGETDC[hits]  = tdcData.edge[j];
      CHTDC[hits]  = tdcData.channel[j];
      VALIDTDC[hits]  = tdcData.valid[j];
//      printf("TDC Nevtda=%d hits= %d chtdc =%d valid=%d count=%d \n",Nevtda, hits,CHTDC[hits] , VALIDTDC[hits], COUNTTDC[hits]);
      hits++;
    }
    NhitTDC = hits;
    //    printf("TDC NUOVO EVENTO NHIT=%d\n",NhitTDC);

    // Convert time information of TDC L1176 in space coordinate
    unsigned int x,y,c[9];
    float count[9][4];
/*  CALIBRATION OF THE DWCs for TB 2007
   float hSlope_up=0.173;
    float vSlope_up=0.176;
    float hOffset_up=0.3;
    float vOffset_up=-1.3;

    float hSlope_dw=0.172;
    float vSlope_dw=0.174;
    float hOffset_dw=-1.2;
    float vOffset_dw=-0.9;
*/
//CALIBRATION OF THE DWCs for TB 2008
//    float hSlope_up=-0.183;
    float hSlope_up=0.183;
    float vSlope_up=0.183;
    float hOffset_up=0.06;
    float vOffset_up=0.79;

//    float hSlope_dw=-0.185;
    float hSlope_dw=0.185;
    float vSlope_dw=0.182;
    float hOffset_dw=0.92;
    float vOffset_dw=0.54;
 


//2007    for (j=0;j<9;j++) c[j]=0;
       for (j=0;j<8;j++) c[j]=0;
    for (j=0; j<NhitTDC; j++) {
      if (VALIDTDC[j]==0) {
//2007	for (unsigned int i=1;i<9;i++) {
         for (unsigned int i=0;i<8;i++) {
	  if (CHTDC[j]==i) {
	    count[i][c[i]]=COUNTTDC[j];
	    //           printf("TDC convers channel=%d c=%d count=%f\n",i,c[i],count[i][c[i]]);
	    c[i]++;
	  } 
	}
      } 
    }
    x=0;
//2007    for (unsigned int i1=0;i1<c[1] && x<=16;i1++)  {
//2007      for (unsigned int i2=0;i2<c[2] && x<=16;i2++) {
       for (unsigned int i1=0;i1<c[0] && x<=16;i1++)  {
      for (unsigned int i2=0;i2<c[1] && x<=16;i2++) {

//2007	X_UP[x]=((count[2][i2]-count[1][i1])*hSlope_up)+hOffset_up;
         X_UP[x]=((count[1][i2]-count[0][i1])*hSlope_up)+hOffset_up;
	//        printf("TDC x=%d i1=%d i2=%d X_UP=%f count1=%f count2=%f\n",x,i1,i2,X_UP[x],count[1][i1],count[2][i2]); 
         
	x++;
      }
    }
    N_X_UP=x;
    //printf("TDC N=%d \n", N_X_UP);
    y=0;
//2007    for (unsigned int i3=0;i3<c[3] && y<=16;i3++)  {
//2007      for (unsigned int i4=0;i4<c[4] && y<=16;i4++) {
       for (unsigned int i3=0;i3<c[2] && y<=16;i3++)  {
             for (unsigned int i4=0;i4<c[3] && y<=16;i4++) {

	//       Y_UP[y]=((count[3][i3]-count[4][i4])*vSlope_up)+vOffset_up;
//2007	Y_UP[y]=((count[4][i4]-count[3][i3])*vSlope_up)+vOffset_up; 
        Y_UP[y]=((count[3][i4]-count[2][i3])*vSlope_up)+vOffset_up;
	//       printf("TDC y=%d i3=%d i4=%d  Y_UP=%f count3=%f count4=%f\n",y,i3,i4,Y_UP[y],count[3][i3],count[4][i4]);
	y++;
      }
    }
    N_Y_UP=y;
    //printf("TDC N=%d \n", N_Y_UP);

    x=0;
//2007    for (unsigned int i5=0;i5<c[5] && x<=16;i5++)  {
//2007      for (unsigned int i6=0;i6<c[6] && x<=16;i6++) {
        for (unsigned int i5=0;i5<c[4] && x<=16;i5++)  {
        for (unsigned int i6=0;i6<c[5] && x<=16;i6++) {

//2007	X_DW[x]=((count[6][i6]-count[5][i5])*hSlope_dw)+hOffset_dw;
        X_DW[x]=((count[5][i6]-count[4][i5])*hSlope_dw)+hOffset_dw;
	//       printf("TDC x=%d i5=%d i6=%d X_DW=%f cont6=%f count5=%f\n",x,i5,i6,X_DW[x],count[6][i6],count[5][i5]);
	x++;
      }
    }
    N_X_DW=x;
    //printf("TDC N=%d \n", N_X_DW); 
    y=0;
//2007    for (unsigned int i7=0;i7<c[7] && y<=16;i7++)  {
//2007      for (unsigned int i8=0;i8<c[8] && y<=16;i8++) {
           for (unsigned int i7=0;i7<c[6] && y<=16;i7++)  {
        for (unsigned int i8=0;i8<c[7] && y<=16;i8++) {

	//       Y_DW[y]=((count[7][i7]-count[8][i8])*vSlope_dw)+vOffset_dw;
//2007	Y_DW[y]=((count[8][i8]-count[7][i7])*vSlope_dw)+vOffset_dw;
       Y_DW[y]=((count[7][i8]-count[6][i7])*vSlope_dw)+vOffset_dw;
	//        printf("TDC y=%d i7=%d i8=%d Y_DW=%f count7=%f count8=%f\n",y,i7,i8,Y_DW[y],count[7][i7],count[8][i8]); 
	y++;
      }
    }
    N_Y_DW=y;
    //printf("TDC N=%d \n", N_Y_DW);
    /*
    //FILL KLOE TDC

    hits = 0;
    addr = SubEventSeek(0x7C000012, buf);
    rc += DecodeKLOETDC(addr, &ktdcData);
    for (j=0; j<ktdcData.index; j++) {
    COUNTKTDC[hits]  = ktdcData.data[j];
    EDGEKTDC[hits]  = ktdcData.edge[j];
    CHKTDC[hits]  = ktdcData.channel[j];
    OVERKTDC[hits]  = ktdcData.over[j];

    hits++;
    }

    NhitKTDC = hits;
    */




    // FILL DRS
    if (drs) {
      // if no match between event of pedestal or phys. data file and event of drs file 
      // drs variables put equal to zero
      for (int i=0;i<10;i++) {
	for (int j=0;j<1024;j++){
	  drs1a[i][j]=0;
	  drs2a[i][j]=0;
	  drs3a[i][j]=0;
	  drs4a[i][j]=0;
	  drs5a[i][j]=0;
	  drs6a[i][j]=0;
	  drs7a[i][j]=0;
	  drs8a[i][j]=0;
          drs1v[i][j]=0;
          drs2v[i][j]=0;
          drs3v[i][j]=0;
          drs4v[i][j]=0;
          drs5v[i][j]=0;
          drs6v[i][j]=0;
          drs7v[i][j]=0;
          drs8v[i][j]=0;
	}
      }

      drs_flag=0;
      bool check_drs=true;
      while (check_drs) {
	//  if(fReadEventHeader(&EventHeader, fptr) == EVENT_HEADER_SIZE)
	//      PrintEventHeader(&EventHeader, stdout);
	//     fReadEventHeader(&EventHeader, fptr);
	if(fReadEventHeader(&EventHeader, fptr) != EVENT_HEADER_SIZE){
	  if (eventnumDRS == evt_tot) {
	    check_drs=false;
	  }else{
	    std::cout<<" myDaq2Ntu : DRS :ERROR Wronf header size..."<<std::endl;
	    return -1;}
	}

	eventnumDRS = EventHeader.EvtNumber;
	TimeSecDRS = EventHeader.TimeSec;
	TimeSubSecDRS = EventHeader.TimeSubSec;
//	std::cout<<"myDaq2Ntu: evt. number  (Nevtda+1)="<<Nevtda+1<<"  DRS evt. number= "<<(eventnumDRS)<<std::endl;
        unsigned int data=0;
        double datad=0;
        float datavo=0;
        int datav=0;
        double datadv=0;
	for (iChannel=0; iChannel<RunHeader->NumChannels; iChannel++) {

	  fReadChannelHeader(&ChannelHeader, fptr);
	  //    std::cout<<"myDaq2Ntu: DRS: === NumChannels:"<<RunHeader->NumChannels <<" Channel Number:"<<ChannelHeader.ChannelNumber<<std::endl;
	  fread(ChannelData, NumSamplesInChannel, NUM_BYTES_PER_SAMPLE, fptr);
	  // check if the match is ok
//          std::cout<<" Nevtda "<<Nevtda <<" eventnumDRS "<<eventnumDRS<<std::endl; 
	  if (Nevtda+1==eventnumDRS) {
//          if (Triggnumb==eventnumDRS) {
	    drs_flag=1;
	    for (unsigned int i=0; i<RunHeader->NumPixInChannel; i++) {
	      //fprintf(stdout, "Pixel %3d :",i);
	      //        std::cout<<"myDaq2Ntu: i "<<i<<" imax= "<<RunHeader->NumPixInChannel<<std::endl;
	      for (unsigned int j=0; j<RunHeader->NumSamplesPerPixel; j++) {
		//fprintf(stdout, " %5d", ChannelData[i*RunHeader->NumSamplesPerPixel + j]);
		data = (unsigned int) ChannelData[i*RunHeader->NumSamplesPerPixel+j];
		int ncha=(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i;   
		datad=(double) data/1000.;
		// Make conversion in 0.1milliVolt
		datav=0;
                datavo=0;
		datadv=0;  
                if (drs_v) {
		  for(int iLevel=0;iLevel<20;iLevel++){
		    if (data<Adc[ncha][j][iLevel+1]){
		      datavo = (int) 10*(data*Slope[ncha][j][iLevel]+Intercept[ncha][j][iLevel]);
		      break;
		    }
                  }
                  if (phys) {
		    if (ncha!=0) { 
		      datav= (int) (datavo - ped_drsv[ncha][j]);
		    } else if (ncha==0) {
		      datav= (int) data;
		    }
                  }else {
		    if (ncha!=0) {
		      datav= (int) datavo;
		    } else if (ncha==0) {
		      datav= (int) data;
		    }

                  }
		  datadv=(double) datav/1000.;   
		} 
		// fill histos
		if (!phys) {
		  // in case of pedestal file, fill the istograms with the first part of the events
		  if (eventnumDRS<=(unsigned int)(evt_tot/2)){
                    if (drs_adc) {
		      hmeanped1[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1,datad);
		      hrmsped1[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1,datad*datad);
                    }
                    if (drs_v){
		      hmeanped1v[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1,datadv);
		      hrmsped1v[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1,datadv*datadv);
                    }
		  }else{
		    //in case of pedestal file, fill the istograms with the second part of the events,
		    // using only GoodEvents (evts with mu-3sig<data<mu+3sig , mu and sig computed with the 
		    //first part of the events
		    //          float mu=hmeanped1[ncha]->GetBinContent(j+1);
		    //          float sig=hrmsped1[ncha]->GetBinContent(j+1);
                    if (drs_adc) {
		      if ( datad>(muv[ncha][j]-3*sigv[ncha][j]) && datad<(muv[ncha][j]+3*sigv[ncha][j]) && (!check_drs_half)) {
			GoodEvts[ncha][j]++;
			hmeanped[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1, datad);
			hrmsped[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1, datad*datad);
		      }
                    }
                    if (drs_v){
		      if ( datadv>(muvv[ncha][j]-3*sigvv[ncha][j]) && datadv<(muvv[ncha][j]+3*sigvv[ncha][j]) && (!check_drs_halfv)) {
			GoodEvtsv[ncha][j]++;
			hmeanpedv[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1, datadv);
			hrmspedv[(ChannelHeader.ChannelNumber-1)*RunHeader->NumPixInChannel+i]->AddBinContent(j+1, datadv*datadv);
		      }
		    }   
		  }

		}
		// fill rootple
		switch ((int)ChannelHeader.ChannelNumber)
		  {
		  case 1:
		    drs1a[i][j]=data;
                    drs1v[i][j]=datav;
		    break;
		  case 2:
		    drs2a[i][j]=data;
                    drs2v[i][j]=datav;
		    break;
		  case 3:
		    drs3a[i][j]=data;
                    drs3v[i][j]=datav; 
		    break;
		  case 4:
		    drs4a[i][j]=data;
                    drs4v[i][j]=datav;
		    break;
		  case 5:
		    drs5a[i][j]=data;
                    drs5v[i][j]=datav;
		    break;
		  case 6:
		    drs6a[i][j]=data;
                    drs6v[i][j]=datav;
		    break;
		  case 7:
		    drs7a[i][j]=data;
                    drs7v[i][j]=datav; 
		    break;
		  case 8:
		    drs8a[i][j]=data;
                    drs8v[i][j]=datav;
		    break;
		  default:
		    std::cout<<"myDaq2Ntu : DRS: unknown channel number :"<<ChannelHeader.ChannelNumber<<" ... I quit"<<std::endl;
		    return 1;
		    break;
		  }//swich
          
	      }//for
	      //fprintf(stdout, "\n");
	    }//for
	  }//if (of the match)
	}//for
	// Exit on the loop on DRS events if you are at the end of the file or if the match is Ok
	if (eventnumDRS==evt_tot) check_drs=false;  
	if (drs_flag==1) check_drs=false;

      }//while
      // Count DRS event matched
      if (drs_flag==1){
	evt_DRS++;
      }
      if (!phys) {

	//Count DRS event matched in the first part of the file
	if (drs_flag==1){
	  if (eventnumDRS<=(unsigned int)(evt_tot/2)) {
	    evt_DRS_fp++;} 
	}

	// Compute ped mu and ped rms using the first part of the file
        if (drs_adc) {
	  if ((eventnumDRS>=(unsigned int)(evt_tot/2))&&(check_drs_half)) {
	    check_drs_half=false;
	    for(int i=0;i<80;i++){
	      hmeanped1[i]->Scale(1./(evt_DRS_fp));
	      hrmsped1[i]->Scale(1./(evt_DRS_fp));
	      for(int j=0;j<1024;j++){
		double mu=hmeanped1[i]->GetBinContent(j+1);
		double sig2=hrmsped1[i]->GetBinContent(j+1);
		double diff2=sig2-mu*mu;
		if (diff2<0)  diff2=0.000001;
		hrmsped1[i]->SetBinContent(j+1,1000*TMath::Sqrt(diff2));
		hmeanped1[i]->SetBinContent(j+1,1000*mu);
		muv[i][j]=mu;
		sigv[i][j]=sqrt(diff2);
	      }

	    }
	  }
        }
	// Compute ped mu and ped rms using the first part of the file (millVolt)
	if (drs_v) {
	  if ((eventnumDRS>=(unsigned int)(evt_tot/2))&&(check_drs_halfv)) {
	    check_drs_halfv=false;
	    for(int i=0;i<80;i++){
	      hmeanped1v[i]->Scale(1./(evt_DRS_fp));
	      hrmsped1v[i]->Scale(1./(evt_DRS_fp));
	      for(int j=0;j<1024;j++){
		double muv=hmeanped1v[i]->GetBinContent(j+1);
		double sig2v=hrmsped1v[i]->GetBinContent(j+1);
		double diff2=sig2v-muv*muv;
		if (diff2<0) diff2=0.000001; 
		hrmsped1v[i]->SetBinContent(j+1,1000*TMath::Sqrt(diff2));
		hmeanped1v[i]->SetBinContent(j+1,1000*muv);
		muvv[i][j]=muv;
		sigvv[i][j]=sqrt(diff2);
	      }

	    }
	  }
        }
      } // if !phys
    } // if drs=true
    //end DRS
    tree->Fill();
  } // End of loop forever

  if (drs) {
    // For pedestal file fill istograms with mean and rms using the second part of the file
    if (!phys) {
      if (drs_adc) {
	for(int i=0;i<80;i++){
	  for(int j=0;j<1024;j++){
	    double mu,sig2,diff2;
	    if (GoodEvts[i][j]!=0) {
	      mu=hmeanped[i]->GetBinContent(j+1)/GoodEvts[i][j];
	      sig2=hrmsped[i]->GetBinContent(j+1)/GoodEvts[i][j];
	      diff2=sig2-mu*mu;
	      if (diff2<0) diff2=0.000001; 
            }else{
	      mu=0;
	      sig2=0;
	      diff2=0;
	    }
	    hmeanped[i]->SetBinContent(j+1,1000*mu);
	    hrmsped[i]->SetBinContent(j+1,1000*TMath::Sqrt(diff2));
  
	  }
          if (!drs_leak) {
	  if (i<40){
	    hmeanped[i]->Write();
	    hrmsped[i]->Write();
	  }
          }
          if (drs_leak) {
           if (i<50){
            hmeanped[i]->Write();
            hrmsped[i]->Write();
            }
          }

	}
      } 
      //Pedestal milliVolt
      if (drs_v) { 
	for(int i=0;i<80;i++){
	  for(int j=0;j<1024;j++){
	    double muv,sig2v,diff2;
	    if (GoodEvtsv[i][j]!=0) {
	      muv=hmeanpedv[i]->GetBinContent(j+1)/GoodEvtsv[i][j];
	      sig2v=hrmspedv[i]->GetBinContent(j+1)/(GoodEvtsv[i][j]);
	      diff2=sig2v-muv*muv;
              if (diff2<0) diff2=0.000001; 
            } else {
	      muv=0;
	      sig2v=0;
	      diff2=0; 
	    }
            hmeanpedv[i]->SetBinContent(j+1,1000*muv);
            hrmspedv[i]->SetBinContent(j+1,1000*TMath::Sqrt(diff2));
	  }
          if (!drs_leak){
	  if (i<40) {
	    hmeanpedv[i]->Write();
	    hrmspedv[i]->Write();
	  }
          }
         if(drs_leak) {
            if (i<50) {
            hmeanpedv[i]->Write();
            hrmspedv[i]->Write();
            }
          } 
	}

      }
    } 
    //DRS Close file
    fclose(fptr);
  } //if drs
  // Close Data file
  RawFileClose();
  // Write and close Ntuple
  //  tree->Write();
  tree->Write("",TObject::kWriteDelete);
  file->Close();

  return 0;

}
