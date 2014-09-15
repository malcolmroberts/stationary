/* ========================================================================== */
/*                                                                            */
/*   PeriodicStochasticEquilibrium_link_new.cpp                                   */
/*   (c) 2014 Schaff                                                          */
/*       corresponding author: frederik.schaff@fernuni-hagen.de               */
/*   08.09.2014 13:55:48                                                      */
/*                                                                            */
/*                                                                            */
/*   Usage: Save Data with Windows Machine.                                   */
/*                                                                            */
/* ========================================================================== */
/* ========================================================================== */
//Declaration of Variables. //exact the same labels as in the functions. Makes things easier.

#define LogNow true

char PLOGTXT[500];



//for saving data
#include <iostream>
#include <fstream>
using namespace std;

void WriteToLog( const std::string &text, int cur_time );
void print_MSER_batch(double MSER_data[], int Tail_end, int number_bins, int batch_size);
void print_data ();
char* statistics_time(char LabelName[32],bool start);
char* statistics_time_complete(char LabelName[32],bool start);
char* statistics_time_local(char LabelName[32],bool start);



#include "PeriodicStatisticEquilibrium_testProcedure.cpp"            //.......8



 
int CurrentRun = 0;

#include <windows.h>
void init_Test(){
CurrentRun++;
char buffer[32]; // The filename buffer.
    // Put "file" then k then ".txt" in to filename.
snprintf(buffer, sizeof(char) * 32, "TestData_%i", CurrentRun);
CreateDirectory(buffer,NULL);
}; 
 
 
                                                  
int LastLogTime = 0;                                                                               
void WriteToLog( const std::string &text, int cur_time ){
    if(LogNow){
    
        //Filename will include the current seed.
    char buffer[32]; // The filename buffer.
    // Put "file" then k then ".txt" in to filename.
    snprintf(buffer, sizeof(char) * 32, "TestData_%i/DebugLog.txt", CurrentRun);

    
      if (LastLogTime!= cur_time){
        LastLogTime = cur_time;
        char txt2[32];
        snprintf(txt2, sizeof(char) * 32, "Logging for time: %i", cur_time);
        std::ofstream log_file(
          buffer, std::ios_base::out | std::ios_base::app );
      log_file << txt2 << std::endl;
      }
      std::ofstream log_file(
          buffer, std::ios_base::out | std::ios_base::app );
      log_file << text << std::endl;
   }
}

/* Bunch of global variables */
#define NumberOfBins SPE_MaximumNumberOfSteps+1
double Residual_Data[NumberOfBins], Transform_Data[NumberOfBins], Transform_StdDev_Data[NumberOfBins],    Time_Series_Data[NumberOfBins], ACF_Tail_Data[NumberOfBins], ACF_Complete_Data[NumberOfBins];
double Gamma_Transform, MSE_Transform, MSE_Simple_Mean, Tail_Mean;


void print_data (){
statistics_time("",true);
/*Save data to harddisk */
save_Data(Residual_Data, Transform_Data, Transform_StdDev_Data,         Time_Series_Data,  ACF_Tail_Data,  ACF_Complete_Data, NumberOfBins, Gamma_Transform, MSE_Transform, MSE_Simple_Mean, Tail_Mean);
WriteToLog(statistics_time("Time needed to calc data to save",false),SPE_TempTailEnd);
statistics_time("",true);
  char txtA[500],txtB[500];
  char buffer[64]; // The filename buffer.
  ofstream outputfile;
  
  /*Autocorrelation */
  
    // Put "file" then k then ".txt" in to filename.
    snprintf(buffer, sizeof(char) * 64, "TestData_%i/Autocorrdata.txt", CurrentRun);
    
    outputfile.open(buffer);
    
    // The Column header
    sprintf(txtA,"#Time\t#ACF_Complete\t#ACF_Tail");
    outputfile << txtA << endl;

    // The Data
    for (int i = 0; i<NumberOfBins; i++){
      sprintf(txtB,"%i\t%.15g\t%.15g",i, ACF_Complete_Data[i],ACF_Tail_Data[i]);
      outputfile << txtB << endl; 
    }
    outputfile.close();
    
    
  /*Data, Transform and Residual */
  
    // Put "file" then k then ".txt" in to filename.
    snprintf(buffer, sizeof(char) * 64, "TestData_%i/Data_and_Transform.txt", CurrentRun);
    
    outputfile.open(buffer);
    
    // The Column header
    sprintf(txtA,"#Time\t#Data\t#Transform\t#StdDev_Transform\t#Residual");
    outputfile << txtA << endl;


    // The Data
    for (int i = 0; i<NumberOfBins; i++){
      sprintf(txtB,"%i\t%.15g\t%.15g\t%.15g\t%.15g",i, Time_Series_Data[i],Transform_Data[i], Transform_StdDev_Data[i],Residual_Data[i]);
      outputfile << txtB << endl; 
    }
    outputfile.close();
    
  /*Statistics - Pseudo time var.*/
  
    // Put "file" then k then ".txt" in to filename.
    snprintf(buffer, sizeof(char) * 64, "TestData_%i/Statistics.txt", CurrentRun);
    
    outputfile.open(buffer);
    
    // The Column header
    sprintf(txtA,"%-16s%-16s%-16s%-16s","#Gamma_Transform","#MSE_Transform","#MSE_simple_mean","#Tail_Mean");
    outputfile << txtA << endl;

             
    // The Data
    sprintf(txtB,"%-16g%-16g%-16g%-16g",Gamma_Transform, MSE_Transform, MSE_Simple_Mean, Tail_Mean);
      outputfile << txtB << endl; 
    
    sprintf(txtA,"%-16s%-16s%-16s%-16s","#CycleLength","#Phase_Start","#Wilcoxon","#Wilcoxon_zValue");
    outputfile << txtA << endl;
    
     sprintf(txtB,"%-16i%-16i%-16g%-16g",SPE_TempCycleLength,SPE_AvgCycle_Phase_start, SPE_Wilcoxon,SPE_Wilcoxon_zVal);
      outputfile << txtB << endl; 
    
    outputfile.close();        
    
WriteToLog(statistics_time("Time needed to save to HardDisk",false),SPE_TempTailEnd);

}




int statistics_startTime,statistics_startTime_complete,statistics_startTime_local;
char TimeNeeded[500];
  char* statistics_time(char LabelName[32],bool start){
    if (start){
      statistics_startTime = clock();
    } else {
      sprintf(TimeNeeded,"%s took %.15g seconds. Time called was %i.",LabelName,(double) (clock()-statistics_startTime)/CLOCKS_PER_SEC,SPE_TempTailEnd);

      return TimeNeeded;
    }
  }
  char* statistics_time_complete(char LabelName[32],bool start){
    if (start){
      statistics_startTime_complete = clock();
    } else {
      sprintf(TimeNeeded,"%s took %.15g seconds. Time called was %i.",LabelName,(double) (clock()-statistics_startTime_complete)/CLOCKS_PER_SEC,SPE_TempTailEnd);

      return TimeNeeded;
    }
  }
  
  //Does not work...

  char* statistics_time_local(char LabelName[64],bool start){
    if (start){
      statistics_startTime_local = clock();
      sprintf(TimeNeeded,"%s",LabelName);
    } else {
      sprintf(TimeNeeded,"%s took %.15g seconds. Time called was %i.",LabelName,(double) (clock()-statistics_startTime_local)/CLOCKS_PER_SEC,SPE_TempTailEnd);

      return TimeNeeded;
    }
  }

  //int MSER_count = 0;
  
  void print_MSER_batch(double MSER_data[], int Tail_end, int number_bins, int batch_size){
  //MSER_count++;
  char txtA[500],txtB[500];
  char buffer[64]; // The filename buffer.
  ofstream outputfile;
  
  /*Autocorrelation */
  
    // Put "file" then k then ".txt" in to filename.
    snprintf(buffer, sizeof(char) * 64, "TestData_%i/MSER_Data.txt", CurrentRun);
    
    outputfile.open(buffer);
    
    // The Column header
    sprintf(txtA,"#Time\t#MSER");
    outputfile << txtA << endl;

    // The Data
    double temp_MSER;
    int Tail_start = Tail_end-batch_size*(number_bins-1)+1;
    for (int i = 0, batch_lag = 0, batch_number=0; i<NumberOfBins; i++, batch_lag++){
      if (i<Tail_start || i >Tail_end){
        temp_MSER = SPE_NotANumber;
      } else {
        if (i == Tail_start){
          batch_lag = 0;
        }
        if (batch_lag == batch_size){
          batch_lag = 0;
          batch_number++;
        }
        temp_MSER = MSER_data[batch_number];
      }
      
      sprintf(txtB,"%i\t%.15g",i, temp_MSER);
      outputfile << txtB << endl; 
    }
    outputfile.close();
  
  
  
  }
