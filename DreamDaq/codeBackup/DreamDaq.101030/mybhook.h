//NTUPLE VAR
#ifndef _MYHBOOK_
#define _MYHBOOK_
/* Global Variables for HBOOK */
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <string>

TFile *hfile;
char datadir[256];
char datafilename[256];
char ntdir[256];
char peddir[256]; 
char ped_run_name[256];
char ntfilename[256];
int evt_max; 
int ped_run_number;
float pmean[66];
int Nrunnumber;
int BegTimeEvs;             // Time of first event in Seconds
int BegTimeEvu;             // Time of first event in MicroSeconds
int TimeEvs;                // Time of event in Seconds
int TimeEvu;                // Time of event in MicroSeconds
int Nevtda;                // Number of data event

unsigned int NSCA;                   //
unsigned int CHSCA[16];              //
unsigned int  COUNTSCA[16];        //

unsigned int NhitTH03;
unsigned int CHTH03[3];
float DATATH03[3];
                                                                                
unsigned int NhitOSC;
unsigned int CHOSC[1128];
unsigned int NUMOSC[1128];
int POSOSC[4];
unsigned int TDOSC[4];
unsigned int SCALEOSC[4];
unsigned int CHFOSC[4];
unsigned int PTSOSC;
unsigned int MAPOSC;
unsigned int SROSC;
    
float sum1[1128];
float sum2[1128];
float sum3[1128];
float sum4[1128];
float sum1v[1128];
float sum2v[1128];
float sum3v[1128];
float sum4v[1128];
float dato1[1128];
float dato2[1128];
float dato3[1128];
float dato4[1128];
float osc_int[4];
float osc_intf[4];
unsigned int num1;
unsigned int num2;
unsigned int num3;
unsigned int num4;
unsigned int counto;

                                                                            
unsigned int NhitADCN0;               // Hits in ADCN0
unsigned int CHADCN0[32];            // Channel in ADCN0
unsigned int CHARGEADCN0[32] ;     // Charge in ADCN0
                                                                                
unsigned int NhitADCN1;               // Hits in ADCN1
unsigned int CHADCN1[32];            // Channel in ADCN1
unsigned int CHARGEADCN1[32] ;     // Charge in ADCN1
                                                                                

unsigned int NhitADCN2;               // Hits in ADCN1
unsigned int CHADCN2[32];            // Channel in ADCN1
unsigned int CHARGEADCN2[32] ;     // Charge in ADCN1


unsigned int NhitADC0;               // Hits in ADC0
unsigned int CHADC0[8];            // Channel in ADC0
unsigned int CHARGEADC0[8] ;     // Charge in ADC0
                                                                                
unsigned int NhitADC1;               // Hits in ADC1
unsigned int CHADC1[8];            // Channel in ADC1
unsigned int CHARGEADC1[8];      // Charge in ADC1
                                                                                
unsigned int NhitADC2;               // Hits in ADC2
unsigned int CHADC2[8];            // Channel in ADC2
unsigned int CHARGEADC2[8];      // Charge in ADC2
                                                                                
unsigned int NhitADC3;               // Hits in ADC3
unsigned int CHADC3[8];            // Channel in ADC3
unsigned int CHARGEADC3[8];      // Charge in ADC3
                                                                                
unsigned int NhitADC4;               // Hits in ADC4
unsigned int CHADC4[8];            // Channel in ADC4
unsigned int CHARGEADC4[8];      // Charge in ADC4
                                                                                
unsigned int NhitADC5;               // Hits in ADC5
unsigned int CHADC5[8];            // Channel in ADC5
unsigned int CHARGEADC5[8];      // Charge in ADC5
                                                                                
unsigned int NhitFADC;               // Hits in FADC
unsigned int CHFADC[128];            // Channel in FADC
float VALUEFADC[128];       // Value in FADC
unsigned int NUMFADC[128];           // Number of sampling in FADC                                                                                 

                                                                                
                                                                                
unsigned int NhitTDC;                // Hits in TDC
unsigned int CHTDC[64];              // Channel in TDC
unsigned int VALIDTDC[64];           //  0=good 1=corrupt
unsigned int EDGETDC[64];            // Edge in TDC (0=falling 1=rising)
unsigned int COUNTTDC[64];         // Count in TDC
float X_UP[16];
float Y_UP[16];
float X_DW[16];
float Y_DW[16];

unsigned int N_X_UP;
unsigned int N_Y_UP;
unsigned int N_X_DW;
unsigned int N_Y_DW;



/*
int NhitKTDC;               // Hits in KLOE TDC
int CHKTDC[64];             // Channel in KLOE TDC
int OVERKTDC[64];           // Overflow in KLOE TDC
int EDGEKTDC[64];           // Edge in KLOE TDC
float COUNTKTDC[64];        // Count in KLOE TDC
*/
TH1F* hrn;
TH1F* hprn;
TH1F* hadc[48];
TH1F* hadcp[48];
TH1F* hadc0n[32];
TH1F* hadc1n[32];
TH1F* hadc0np[32];
TH1F* hadc1np[32];
TH1F* hadc0ovmap;
TH1F* hadc1ovmap;
TH1F* hs;
TH1F* hsp;
TH1F* hq;
TH1F* hqp;
TH1F* ha;
TH1F* hap;
TH1F* hfadc[8];
TH1F* htdc[16];
TH1F* htemp[3];
TH1F* htemptrend[3];
TH1F* hosc[4];
TH1F* hunder[4];
TH2F* hfadcsetup1[8];
TH2F* hfadcsetup2[8];

TH1F* hosc2[100][4];
TH1F* hosc_mean[4];
TH1F* hosc_meanv[4];

TH1F* htdccorr;
TH1F* htdcupud;
TH1F* htdcuplr;
TH1F* htdcdwud;
TH1F* htdcdwlr;
TH1F* h_x_dw;
TH1F* h_y_dw;
TH2F* h_xy_dw;
TH1F* h_x_up;
TH1F* h_y_up;
TH2F* h_xy_up;
TH1F* hosc_int_mean_v[4];
/*TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
TH1F* ;
 */
#endif
// 
  
