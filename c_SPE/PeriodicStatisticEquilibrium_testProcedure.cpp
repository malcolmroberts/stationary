/* ========================================================================== */
/*                                                                            */
/*   PeriodicStochasticEquilibrium_coreAlgorithms.cpp                         */
/*   (c) 2014 Klabunde / Roberts / Schaff                                     */
/*   05.09.2014 14:01:34                                                      */
/*                                                                            */
/*   Test procedure. Whereas the coreAlgorithms are pure stand alone, the     */
/*   test procedure makes use of global parameters defined in the CONFIG.     */
/*                                                                            */
/*   Usage on the run: Each time the simulated time series is updated, also   */
/*   update the test by calling:                                              */
/*         bool OnTheRun_Test(double Data_value,int data_time);               */
/*   The test returns "true" when a periodic component has been found, this   */
/*   should be taken as stop-condition for the simulation-run.                */
/*                                                                            */
/*   The function: void save_Data(double Residual_Data[], double 
        Transform_Data[], double Transform_StdDev_Data[], double 
        Time_Series_Data[], double ACF_Tail_Data[], double ACF_Complete_Data[], 
        int NumberOfBins, double & Gamma_Transform, double & MSE_Transform, 
        double & MSE_Simple_Mean, double & Tail_Mean);                                                                         
/*   Can be used to save the results. Not that it may only be used correctly, */
/*   once the stopcondition is met.                                           */
/*                                                                            */
/*                                                                            */
/* ========================================================================== */
/* ========================================================================== */

/* Include Config */

#include "PeriodicStatisticEquilibrium_CONFIG.txt"


/* Declarations */
void init_TestVariables();
int calc_WindowSize(int Data_start, int Data_end);
bool preTest_onTheRun(double Data_input[], int Data_end);
bool OnTheRun_Test (double Data_value, int Data_time);
void save_Data(double Residual_Data[], double Transform_Data[], double Transform_StdDev_Data[], double Time_Series_Data[], double ACF_Tail_Data[], double ACF_Complete_Data[], int NumberOfBins, double & Gamma_Transform, double & MSE_Transform, double & MSE_Simple_Mean, double & Tail_Mean);
int resize_TailSimple(double Data_input[],int Data_end, int batch_size);
int resize_TailFinal(double Data_input[],int Tail_start,int Data_end, double AvgCycle_Input[], int Cyclelength, int Phase_start);
int find_optCyclelength(double Data_input[],int Tail_start, int Data_end, int maxCycleLength, double ACF_temp[]);
bool check_SampleSize(int Cyclelength, int Tail_start, int Data_end, double AvgCycle_Input[], int Phase_start);
 int find_optAvgCycle(double ACF_input[],int ACF_bins, int ACF_dominantCycle,double Data_input[], int Tail_start, int Data_end);
bool Tail_is_stationary_Res (double Data_Input[], int Tail_start, int Tail_end, double AvgCycle_Input[], int Cyclelength, int Phase_start);
bool Tail_is_stationary (double Data_Input[], int Tail_start, int Tail_end);
bool Periodicity_is_significant (double Data_Input[], int Tail_start, int Tail_end, double AvgCycle_Input[], int Phase_start, int Cyclelength);
void calc_opt_Cycle ();
  /* Global Variables */
  
  //Used by the preTest_onTheRun and *_offTheRun
  int SPE_CurrentTestWindowSize, SPE_preTestStart, SPE_NextTestTime, SPE_TempTailStart, SPE_TempCycleLength, SPE_TempTailEnd, SPE_AvgCycle_Phase_start, SPE_MaxCycleLength;
  double SPE_CurrentMeanDelta, SPE_CurrentVarDelta, SPE_FrameASum, SPE_FrameBSum, SPE_Wilcoxon, SPE_Wilcoxon_zVal, SPE_Gamma, SPE_relativeMSE;
  double SPE_Data[SPE_MaximumNumberOfSteps+1],SPE_AvgCycle[SPE_MaximumNumberOfSteps+1],SPE_AvgCycle_StdDev[SPE_MaximumNumberOfSteps+1], SPE_ACF[SPE_MaximumNumberOfSteps+1], SPE_MSER[SPE_MaximumNumberOfSteps+1],SPE_Residual[SPE_MaximumNumberOfSteps+1];
  bool SPE_PreTestDue, SPE_TestOK;
  void init_TestVariables(){
/* Initialises global variables */
    SPE_PreTestDue = false;
    SPE_preTestStart = 1;
    SPE_CurrentTestWindowSize = 0;
    SPE_NextTestTime = SPE_MinimumNumberOfSteps;
    SPE_MaxCycleLength = (int) floor(SPE_MaximumNumberOfSteps/SPE_CONFIG_MinimumNumberOfCycles);
  }

/* Include CORE Algorithms */

#include "PeriodicStatisticEquilibrium_coreAlgorithms.cpp"


/* Algorithms */

/* 07.09.2014 14:33:07 Simpler Version as before, because the Tail-Refinement 
will take care of the "correct" truncation point anyway and the sole reason of performing a pre-test is to provide a first time for the analysis. */
                                                                                                              
int calc_WindowSize(int Data_start, int Data_end){
/* Calculates the confirmation window size. The ratio of tail to transient as 
given by SPE_MinimumTailRatio is considered, e.g. the window size will only be 
calculated such that the tail under scrutiny is as long as the inverse of this 
ratio. (The inverse condition) This will allow a faster initial detection. 

  
  /* The Tail is split in two parts, the ConfirmationWindows A and B plus some 
  residual, for we only consider log2 sizes. The reason is that this allows for 
  simple updating plus a doubling of window sizes will never decrease the 
  statistical homogeneity of the data in the window, given the data is periodic 
  stationary (??)*/
  
  int ConfirmationWindowSize = 2; //This is the minimum size.
  while (ConfirmationWindowSize*2 <= (Data_end - Data_start +1)){
    ConfirmationWindowSize *= 2;
  }
  return ConfirmationWindowSize;
}

bool preTest_onTheRun(double Data_input[], int Data_end){
/* The preTest statistic is simply the relative difference of the mean and variance of the two confirmation windows. */
  
  int oldStart = SPE_preTestStart;
  SPE_CurrentTestWindowSize = calc_WindowSize(1,Data_end);
  SPE_preTestStart = Data_end-SPE_CurrentTestWindowSize+1;
  int FrameBStart = SPE_preTestStart + SPE_CurrentTestWindowSize/2; 
  
  /* If the TestWindowSize didn't change we can update the mean very fast */
  if (oldStart == SPE_preTestStart-1){
    SPE_FrameASum = SPE_FrameASum - Data_input[oldStart] + Data_input[FrameBStart-1];
    SPE_FrameBSum = SPE_FrameBSum - Data_input[FrameBStart-1] + Data_input[Data_end];
  } else {
    SPE_FrameASum=0, SPE_FrameBSum=0;
    for (int i = 0; i < SPE_CurrentTestWindowSize/2; i++){
      SPE_FrameASum+=Data_input[SPE_preTestStart+i];
      SPE_FrameBSum+=Data_input[FrameBStart+i];
    }
  }
  
  /* Next the relative difference from frames A to B (centered) is calculated. 
  This is the difference from the frames (symetric) to the window mean */
  SPE_CurrentMeanDelta =  abs(SPE_FrameASum-SPE_FrameBSum)/((SPE_FrameASum+SPE_FrameBSum)/2);                                          

  /* If the FastPreTest is used, the Difference in Variance is only calculated if the Difference in Mean is below the tolerance level. */
  
  if (!SPE_FAST_VAR || SPE_CurrentMeanDelta <= SPE_MaximumDeltaMeanTolerance){
    double A_mean = SPE_FrameASum/SPE_CurrentTestWindowSize*2;
    double B_mean = SPE_FrameBSum/SPE_CurrentTestWindowSize*2;
    double A_var=0, B_var=0;
  
    for (int i = 0; i < SPE_CurrentTestWindowSize/2; i++){
      A_var+=pow(Data_input[SPE_preTestStart+i]-A_mean,2);
      B_var+=pow(Data_input[FrameBStart+i]-B_mean,2);
    }
    SPE_CurrentVarDelta = abs(A_var-B_var)/((A_var+B_var)/2); 
  } else {
    SPE_CurrentVarDelta = 2*SPE_MaximumDeltaVarTolerance;    //not necessary...
  }
    //sprintf(PLOGTXT,"\nSPE_CurrentMeanDelta: %g\nSPE_CurrentVarDelta: %g",SPE_CurrentMeanDelta,SPE_CurrentVarDelta);
    //plog(PLOGTXT);
  



  /*11.09.2014 07:48:49
    The pre-test is advanced in a similar way as before due to detecting a 
    possible SPE to early. 
    Considerations: In face of periodic components, to small batch sizes will 
    cut in between cycles and thus variance among the batches will be too big.
    Therefore, as a rough guess, we will define as many test-batches as [one 
    half] the number of Cycles necessary (m) for SPE detection. We will 
    then aportion the test-window in m equal sized neighbouring test-batches 
    (starting from the right) and check if the means are within a given 
    tolerance level of the window-mean. If not, we will decrease the batch-size 
    until either the tolerance level is met or the batch-size was decrease by 
    one-third, where we will stop. Only the man, not the variance, is checked.                                                  
  */
  
  int NoOfTestBatches = SPE_CONFIG_MinimumNumberOfCycles;
  int SizeOfTestBatches = (int) floor(SPE_CurrentTestWindowSize/ NoOfTestBatches);
  int minSizeOfTestBatches = (int) ceil(SizeOfTestBatches*2/3);
  //We can use the unnormalised mean
  double MeanOfTestBatches[NoOfTestBatches];
    //Initialise
  for (int i = 0; i <NoOfTestBatches; i++){
    MeanOfTestBatches[i]=0.0;
  }
    //Calculate the "starting" mean
  for (int i =  Data_end; i > Data_end - SizeOfTestBatches; i --){
    for (int j = 0; j < NoOfTestBatches; j++){
      MeanOfTestBatches[j]+=Data_input[i-j*SizeOfTestBatches];
    }  
  }
    double ReferenceMean;
  
    bool BatchTestOk = false, BatchTestAtEnd = false;
    int count, count2;
    while (!BatchTestOk && !BatchTestAtEnd){
      ReferenceMean=0.0;
      for (count = 0; count < NoOfTestBatches; count++){
        ReferenceMean+=MeanOfTestBatches[count]/NoOfTestBatches; 
      }
      for (count = 0; count < NoOfTestBatches; count++){
        if ( abs(1-MeanOfTestBatches[count]/ReferenceMean)<SPE_MaximumDeltaMeanTolerance  ) {
          BatchTestOk = true;
        } else {
          BatchTestOk = false;
          break;
        }
      }
      

      if (SizeOfTestBatches-1 < minSizeOfTestBatches){
        BatchTestAtEnd = true;
      }
      
        //Resize batches and change info. The batches are shifted and decreased.
      if (!BatchTestOk && !BatchTestAtEnd){

        
        for (count = 0; count < NoOfTestBatches; count++){
          //first, take away data on left side
          for (count2 = count; count2 >=0; count2--){
            MeanOfTestBatches[count]-= Data_input[Data_end-(count+1)*SizeOfTestBatches+1-count2];
          }
          //next, add data on right size 
          for (count2 = count; count2 >0; count2--){
            MeanOfTestBatches[count]+= Data_input[Data_end-(count)*SizeOfTestBatches+count2];           
          }          
        }
        SizeOfTestBatches--;
      }
  
    } 
    

  return (SPE_CurrentMeanDelta<=SPE_MaximumDeltaMeanTolerance && SPE_CurrentVarDelta <=SPE_MaximumDeltaVarTolerance && BatchTestOk);

}




bool OnTheRun_Test (double Data_value, int Data_time){
/* This is the top-level procedure for the on-the run test. It will return true when a stochastic periodic equilibrium is found and false else. It needs consecutive updating at each time a new point of data is simulated.*/

  SPE_TestOK = false;
  
  if (Data_time <= SPE_MaximumNumberOfSteps){
  /* First, we update the data */
   
    if (Data_time != SPE_TempTailEnd+1 && Data_time != 1){
    plog("\nThe data is not being updated continuously. Provide a continuous time_scale for the test to work correctly.");
    }
  SPE_TempTailEnd = Data_time; //Save the end point.
   /* Initilise at start */
   if (SPE_TempTailEnd == 1){
    init_TestVariables();
    statistics_time_complete("",true);
   }
  
  SPE_Data[SPE_TempTailEnd] = Data_value;


  bool SampleSizeOk = false;

  /* Check if a new test is warranted */
    if (SPE_TempTailEnd >= SPE_NextTestTime){
  
      /* As long as the PreTest did not find an initial truncation point, it 
      will be performed each time the test is called. */
      if (!SPE_PreTestDue) {
        //WriteToLog("Pretest started",Data_time);
          SPE_PreTestDue=preTest_onTheRun(SPE_Data, SPE_TempTailEnd);
          SPE_TempTailStart=SPE_preTestStart;             
      }
      
      //12.09.2014 13:07:00 Added: After pre-test is O.K. a first refinement based on the MSER-5 is made
      if (SPE_PreTestDue && SPE_TempTailStart==SPE_preTestStart){
        if (SPE_FixedTruncation == 0){
        SPE_TempTailStart=resize_TailSimple(SPE_Data, SPE_TempTailEnd, 5); 
        } else {
          SPE_TempTailStart=SPE_FixedTruncation;
        }   
      }
      
      
      /* Whenever a test is indicated (first time after Pretest OK) a first 
      rough analysis is performed to refine the Tail_start and check if the 
      sample_size is big enough. */  
      if (SPE_PreTestDue) {
        WriteToLog("Pretest due",Data_time);
        /* The very first iteration is to calculate the dominant cycle of the 
        tail as is and resize the tail on the basis of a batch-mean statistic 
        with batch-size equal to the dominant cycle length*/
        
        // SPE_TempTailStart = resize_Tail(SPE_Data, SPE_TempTailStart, SPE_TempTailEnd);
        
        //11.09.2014 09:51:46 TEST vs Batch Size = CycleSize, if any
        
        /* We analyse the tail and find the optimal cycle (if significant, else mean aka pseudo-one-cyle)*/
       calc_opt_Cycle();
        
        /* We resize the tail */
        SPE_TempTailStart = resize_TailFinal(SPE_Data, SPE_TempTailStart,SPE_TempTailEnd, SPE_AvgCycle, SPE_TempCycleLength, SPE_AvgCycle_Phase_start); 
        sprintf(PLOGTXT,"\tTempTail Start refined to %i",SPE_TempTailStart);
        WriteToLog(PLOGTXT,Data_time);
        /* And once again calc. the optimal cycle based on the resisement */
        calc_opt_Cycle (); 
        
        
        /* and check, wether the Sample is big enough to accomodate the minimum 
        number of cycles. If not, the new next testtime is set.*/
        SampleSizeOk = check_SampleSize(SPE_TempCycleLength,SPE_TempTailStart,SPE_TempTailEnd,SPE_AvgCycle,SPE_AvgCycle_Phase_start);
      }
      
     /* If the Sample is big enough (implying a short-enough cycle exists), the 
     average cycle is calculated and the tail is refined one more time. Only now 
     do we not allow to long cycles. */
     if (SampleSizeOk){
       WriteToLog("SampleSize is now Ok.",Data_time);     
       SPE_TestOK = true; //StopCondition.
     }
  
  
   }
  
  }  else {
  plog("\n The maximum Number of Datapoints defined is reached. Rescale Test.");
  }
  
  if (SPE_TestOK || SPE_TempTailEnd == SPE_MaximumNumberOfSteps) {
  WriteToLog(statistics_time_complete("\n::::::::::::::::::\n\t\tThe complete TEST ",false),SPE_TempTailEnd);
  }
  
  return SPE_TestOK;
    

}

void save_Data(double Residual_Data[], double Transform_Data[], double Transform_StdDev_Data[], double Time_Series_Data[], double ACF_Tail_Data[], double ACF_Complete_Data[], int NumberOfBins, double & Gamma_Transform, double & MSE_Transform, double & MSE_Simple_Mean, double & Tail_Mean){
/* This procedure calculates the transform and other stuff and simply saves the 
data to the pointed variables and arrays. Note, that the arrays need all to be 
of the same size "NumberOfBins". */
     WriteToLog("\t:::\tNow in save_Data()",SPE_TempTailEnd);
/*In the Arrays holding the data, the first bin will always be SPE_NaN and all 
non-used bins the same.*/

  /* First we calculate and save the Autocorrelation Info */

         statistics_time_local("calc_ACF_Complete_time",true);
  calc_ACF(ACF_Complete_Data, SPE_Data, 1, SPE_TempTailEnd, 0,SPE_MaxCycleLength);
         WriteToLog(statistics_time_local("calc_ACF_Complete_time",false),SPE_TempTailEnd);
        statistics_time_local("calc_ACF_Tail_time",true);
  double temp_ACF_Tail[SPE_TempTailEnd-SPE_TempTailStart+1];
  calc_ACF(temp_ACF_Tail, SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, 0,min((SPE_TempTailEnd-SPE_TempTailStart+1)/2,SPE_MaxCycleLength));
        WriteToLog(statistics_time_local("calc_ACF_Tail_time",false),SPE_TempTailEnd);
  //AC goes from 0 on.
  statistics_time_local("save ACF Data to arrays",true);
  for (int i = 0, tbin = 0; i < NumberOfBins; i++){
    if (i >= SPE_TempTailEnd){
      ACF_Complete_Data[i]=SPE_NotANumber;
      ACF_Tail_Data[i]=SPE_NotANumber;
    } else if (i < SPE_TempTailStart){
      ACF_Tail_Data[i]=SPE_NotANumber;
    } else {
      ACF_Tail_Data[i]=temp_ACF_Tail[tbin];
      tbin++;
    }
  }
  WriteToLog(statistics_time_local("save ACF Data to arrays",false),SPE_TempTailEnd);
       WriteToLog("\t\t\t ... saved the ACF..",SPE_TempTailEnd);

  /* Next the rest of the Data */
    statistics_time_local("save other Data to arrays",true);
  if (SPE_TestOK){

    for (int bin = 0, lag = Phase_Info(SPE_TempTailStart,SPE_TempCycleLength,0); bin < NumberOfBins; bin++, lag++){
      if (lag == SPE_TempCycleLength){
        lag = 0;
      }
      
      if (bin == 0 || bin >SPE_TempTailEnd){
        Time_Series_Data[bin]=SPE_NotANumber;
        Residual_Data[bin]=SPE_NotANumber;
        Transform_Data[bin]=SPE_NotANumber;
        Transform_StdDev_Data[bin]=SPE_NotANumber;
    
      } else {
        Time_Series_Data[bin]=SPE_Data[bin];
        Residual_Data[bin]=SPE_Data[bin]-SPE_AvgCycle[lag];
        
        if (bin < SPE_TempTailStart){
          Transform_Data[bin]=SPE_NotANumber;
          Transform_StdDev_Data[bin]=SPE_NotANumber;
        } else {
          Transform_Data[bin]=SPE_AvgCycle[lag];
          Transform_StdDev_Data[bin]= SPE_AvgCycle_StdDev[lag];
        }
      }
    }
  
   } else {
    for (int bin = 0; bin < NumberOfBins; bin++){
    if (bin == 0 || bin >SPE_TempTailEnd){
        Time_Series_Data[bin]=SPE_NotANumber;
        Residual_Data[bin]=SPE_NotANumber;
        Transform_Data[bin]=SPE_NotANumber;
        Transform_StdDev_Data[bin]=SPE_NotANumber;
      } else {
        Time_Series_Data[bin]=SPE_Data[bin];
        Residual_Data[bin]=SPE_Data[bin];
        Transform_Data[bin]=0;
        Transform_StdDev_Data[bin]=0;
      }
    }
   
   } 
   WriteToLog(statistics_time_local("save other Data to arrays",false),SPE_TempTailEnd);
           WriteToLog("\t\t\t ... saved the Data..",SPE_TempTailEnd);
  /* Finally, the basic statistics are calculated and saved. */
  if (SPE_TempCycleLength>0){
  Gamma_Transform = calc_gamma_corr(SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, SPE_AvgCycle, SPE_TempTailStart, SPE_TempCycleLength);
  MSE_Transform = calc_MSE_AvgCycle(SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, SPE_AvgCycle, SPE_TempTailStart, SPE_TempCycleLength);
  Tail_Mean=0;
  for (int bin = 0; bin < SPE_TempCycleLength; bin++){
    Tail_Mean+=SPE_AvgCycle[bin];
  }
  Tail_Mean /=  SPE_TempCycleLength;
  MSE_Simple_Mean = calc_MSE_AvgCycle(SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, & Tail_Mean, SPE_TempTailStart, 1);
  } else {
  Gamma_Transform = -1;
  MSE_Transform = -1;  
  Tail_Mean = -1;
  MSE_Simple_Mean = -1;  
  }
  
  if (!SPE_TestOK) {

    WriteToLog("!!!!!!!!!\nTest not OK but at end.\n!!!!!!!",SPE_TempTailEnd);
  }

}
 

int resize_TailSimple(double Data_input[], int Data_end, int batch_size){
/* This procedure resizes the tail. An initial starting point SPE_TempTailStart 
is used to calculate the dominant cycle (if none, 1-cycle), calculate the batch 
means statistic with a batch size according to the dominant cycle and calculate 
the MSER-batch statistic. The minimum MSER and on tie the earliest one is 
selected as starting point for further investigations. */

  /* a) Calculate maximal Cyclelength such that at least two cycles are covered. */
  WriteToLog("\t:::\tNow in resize_TailSimple()",SPE_TempTailEnd);
  //int temp_maxCycleLength = (int) floor(Data_end/SPE_CONFIG_MinimumNumberOfCycles); //E.g. minimum according number of batches in complete data.
  
  /* b) And recieve the optimal cyclelength based on this */
    //double temp_ACF[Data_end-Tail_start+1];
  //int temp_dominant_cyclelength = find_optCyclelength(Data_input,Tail_start, Data_end,temp_maxCycleLength,temp_ACF );
                                                    
  /* c) This is used to calculate the batch-averages for the complete data, starting from the end and possibly neglecting partial batches at the beginning. */
  
  //09.09.2014 13:28:09 Use simple MSER-5 Test on raw data.
  
  int temp_truncPoint = 1;
  
  if (SPE_FixedTruncation == 0){

  
    sprintf(PLOGTXT,"\t\t\tBatch-%i MSER Stats with batches are calculated",batch_size);
    WriteToLog(PLOGTXT,SPE_TempTailEnd);
    
  int batch_number = (int) floor(Data_end/batch_size);
  double Batch_Mean[batch_number]; //Initialised to 0.
  for (int i = 0; i<batch_number; i++){
    Batch_Mean[i]=0.0;                                 
  }
   /* Unnormalised, e.g Sums */
  for (int i = batch_number, bin, data=Data_end; i> 0; i--){
    for ( bin = batch_size; bin >0; bin--){
      Batch_Mean[i]+=Data_input[data];
      data--;                                       
    }
  }
      WriteToLog("\t\t\tBatch Means have been calculated",SPE_TempTailEnd);
 
 /* d) We derive the MSER-Batch Statistic */
 double Batch_MSER[batch_number+1];

 calc_MSER_complete(Batch_Mean, batch_number, Batch_MSER);     //Check speed
      WriteToLog("\t\t\tA MSER Stats due",SPE_TempTailEnd);
 
 //print_MSER_batch(Batch_MSER, Data_end, batch_number+1, batch_size);
 
 
 /* e) and search for the minimum MSER which is most to the left, neglecting a number of batches equal to half the number of batches inside the pre-tail  */
 //Ignore first 10% of  complete series
 int ignore_batches = (int) floor(batch_number/10);
 double MSER_min =  Batch_MSER[1];
 int MSER_min_batch = 1;
 for (int batch = batch_number-ignore_batches; batch>=1;batch--){
                  //-1 because arrays index from 0 to n-1.
  if (Batch_MSER[batch]<=MSER_min){
    MSER_min_batch = batch;
    MSER_min =Batch_MSER[batch];   
  }  
 }
 
 /* f) Now we check for the first point left of the minMSER position which is outside the sensitivity of our MSE at minMSER */
 double thresholdMSE = MSER_min*(batch_number-MSER_min_batch+1) * (1+SPE_MSESensitivity);
 int MSE_min_batch = MSER_min_batch;
 for (int batch = MSER_min_batch; batch >= 1;batch--){
  if (Batch_MSER[batch]*(batch_number-batch+1) <= thresholdMSE){
    MSE_min_batch = batch;
  }
 }
 sprintf(PLOGTXT,"The new sensitivity check refined the truncation point from batch %i with MSER %.15g to batch %i. with MSER %.15g",MSER_min_batch,Batch_MSER[MSER_min_batch],MSE_min_batch,Batch_MSER[MSE_min_batch]);
 WriteToLog(PLOGTXT,Data_end);
 /*We do now have the new truncation point as the first bin in the batch received*/
 temp_truncPoint = Data_end - (batch_number-MSE_min_batch-1)*batch_size +1;
sprintf(PLOGTXT,"The new start of the time-series is: %i",temp_truncPoint);
WriteToLog(PLOGTXT,Data_end);
 } else {
   temp_truncPoint = SPE_FixedTruncation;
 }
 
 
 return temp_truncPoint;

}

int resize_TailFinal(double Data_input[],int Tail_start,int Data_end, double AvgCycle_Input[], int Cyclelength, int Phase_start){
/* resizes the tail based on the residual of the data and the transform of the 
avg cycle, the latter beeing calculated first (for the tail). */
     WriteToLog("\t:::\tNow in resize_TailFinal()",SPE_TempTailEnd);
  /* First, the current avg cycle is (re)calculated */
     
     int new_Start = Tail_start;
     
     if (SPE_FixedTruncation == 0){
 
  /* Next the MSER-R statistic is calculated. For now the whole range is covered, 
  but later only the part considered (e.g. without the last 
  SPE_CONFIG_MinimumNumberOfCycles cycles) should be calculated for performance 
  reasons.) */
                                //for simplicity index=time. 0 is not in use.
  double MSER_statistic[Data_end+1];

  calc_MSER_completeRes(Data_input, Data_end, AvgCycle_Input, Cyclelength, Phase_start, MSER_statistic);
  print_MSER_batch(MSER_statistic, Data_end, Data_end+1, 1);
  
  /* Search first minimum in valid interval of MSER-R Statistic */

  double min_MSER = MSER_statistic[Tail_start];
  for (int i = 1;i<Data_end-Cyclelength+1;i++){
    if (MSER_statistic[i]<min_MSER){          
      min_MSER=MSER_statistic[i];
      new_Start = i;
    }
  }
  sprintf(PLOGTXT,"\t\t The cycle was of length %i and the min MSER was found at %i",Cyclelength,new_Start);
  WriteToLog(PLOGTXT,Data_end);
  
   /* f) Now we check for the first point left of the minMSER position which is outside the sensitivity of our MSE at minMSER */
 double thresholdMSE = min_MSER*(Data_end-new_Start+1) * (1+SPE_MSESensitivity);
 int MSE_min_pos = new_Start;
 for (int i = MSE_min_pos; i >= 1;i--){
  if (MSER_statistic[i]*(Data_end-i+1) <= thresholdMSE){
    MSE_min_pos = i;
  }
 }
 sprintf(PLOGTXT,"The new sensitivity check refined the truncation point from %i with MSER %.15g to %i. with MSER %.15g",new_Start,MSER_statistic[new_Start],MSE_min_pos,MSER_statistic[MSE_min_pos]);
 WriteToLog(PLOGTXT,Data_end);
 new_Start = MSE_min_pos; 
  
  } else {
    new_Start=SPE_FixedTruncation;
  }
  
  return new_Start;

}

int find_optCyclelength(double Data_input[],int Tail_start, int Data_end, int maxCycleLength, double ACF_temp[]){
     WriteToLog("\t:::\tNow in find_optCyclelength()",SPE_TempTailEnd);
/* a) Calculate ACF */
          statistics_time_local("calc_ACF",true);
          int ACF_bins = min((Data_end-Tail_start+1)/2,maxCycleLength);
  calc_ACF(ACF_temp, Data_input, Tail_start, Data_end, 0,maxCycleLength);
          WriteToLog(statistics_time_local("\t\tcalc_ACF ",false),Data_end);
  /* b) Calculate current dominant cycle length based on ACF and maxlength*/   
                  statistics_time_local("calc_dominant_cycle_length",true);
  int temp_dominant_cyclelength = calc_dominant_cycle_length(ACF_temp, ACF_bins, SPE_AutocorrelationSensitivity, SPE_AutocorrelationPeakSensitivity, 0, maxCycleLength);
  WriteToLog(statistics_time_local("calc_dominant_cycle_length",false),Data_end);
    sprintf(PLOGTXT,"\t:::\t opt cycle length calculated as %i with TailStart %i",temp_dominant_cyclelength, Tail_start);
       WriteToLog(PLOGTXT,SPE_TempTailEnd);
  return temp_dominant_cyclelength;
}


bool check_SampleSize(int Cyclelength, int Tail_start, int Data_end, double AvgCycle_Input[], int Phase_start){
/* This test procedure checks if an SPE has been found or else updates the next 
time for testing. */
       WriteToLog("\t:::\tNow in check_SampleSize()",Data_end);
  
  bool temp_OK = true;
  
  
  
  
  /*a) Is the sample big enough? If not, reframe next stopping point.*/
  
    //a i) Absolut Sample Size 
  if ((Data_end - Tail_start +1) < SPE_MinimumSampleSize){
    temp_OK = false;
    SPE_NextTestTime = (int) min(SPE_MaximumNumberOfSteps, Tail_start+SPE_MinimumSampleSize*(1+SPE_TestTimeBonus));                                            
  }
  
  // b) THe ratio of the tail to the transient-part
  if ( (Data_end-Tail_start+1)/Tail_start < SPE_MinimumTailRatio){
    temp_OK = false;
    SPE_NextTestTime = (int) min(SPE_MaximumNumberOfSteps, (1+SPE_MinimumTailRatio)*Tail_start*(1+SPE_TestTimeBonus));
  }
  
     
       
    //c if There is a cyclic component but we do not yet have enough cycles 
    
    if (Cyclelength >0 && (Data_end - Tail_start +1)/Cyclelength < SPE_CONFIG_MinimumNumberOfCycles) {
    temp_OK = false;
    SPE_NextTestTime = (int) max(SPE_NextTestTime, Tail_start+SPE_CONFIG_MinimumNumberOfCycles*Cyclelength*(1+SPE_TestTimeBonus));
    SPE_NextTestTime = min (SPE_NextTestTime, SPE_MaximumNumberOfSteps);
  }
  
  //d Wilcoxon_signed_rank test if the tail is stationary. If yes, o.k., if not, go on.
  if (temp_OK) {
    temp_OK = Tail_is_stationary_Res (SPE_Data, Tail_start, Data_end,  AvgCycle_Input, Cyclelength, Phase_start );
  }
  
  
  /* MInimum increase in test-time */
  if (!temp_OK){
  SPE_NextTestTime = (int) max(SPE_NextTestTime,Data_end+floor((Data_end-Tail_start+1)/4) );    //ToDo Parameterise
  SPE_NextTestTime = min (SPE_NextTestTime, SPE_MaximumNumberOfSteps);
  }
  
  
  return temp_OK;

}

int find_optAvgCycle(double ACF_input[],int ACF_bins, int ACF_dominantCycle,double Data_input[], int Tail_start, int Data_end){
/* Finds the optimal avg Cycle. */
       WriteToLog("\t:::\tNow in find_optAvgCycle()",SPE_TempTailEnd);
  /* We define the maximum testable cyclelength although the SampleSize has been 
  tested already, because we want to try cycles longer than the dominant cycle 
  s.t. ACF analysis. */
                                                                  //At this stage we do ensure that no bigger cycles are found.
  int maxTestedCyclelength = (int) min(SPE_MaxCycleLength,(floor (Data_end-Tail_start+1)/2));//SPE_CONFIG_MinimumNumberOfCycles); 
  
  int MinCycleLength, MaxCycleLength;
  
  //This procedure calculates the left and right border of the test interval.
  calc_cyclelength_interval(ACF_input, ACF_bins, ACF_dominantCycle, SPE_PeriodicityACtolerance, SPE_AutocorrelationSensitivity, 0, maxTestedCyclelength, MinCycleLength, MaxCycleLength);

  /* The following procedure calculates all possible cycles in the range given 
  above for the tail given and saves them to a temporary array. It also 
  calculates the MSE for each single one (Transform vs Data for Tail) and saves 
  the statistics. 
  
  Considerations: The rms of all different transforms/avg cycles is the same, 
  for the transforms are energy preserving per se. The gamma correlation will be 
  different and it might even be the case that the cycle with the highest gamma 
  differs from the cycle with the lowest MSE. However, the MSE is a better 
  predictor / more important. The gamma is only a rank correlation whereas the 
  MSE considers the quantitative difference. However, deeper (theoretical) 
  investigation might be of interest.
  */
  
  /* At the moment the Gamma and RMS are also calculated. In a later, efficient 
  version, this should be deleted. Also the optimal avg cycle is calculated two 
  times atm.*/
  
  
  double cycle_statistics[MaxCycleLength-MinCycleLength+1][2];  //0: Cyclelength 1: MSE
  double temp_AvgCycle[MaxCycleLength];
  int Phase_start; //holds info of tail start.
  //double rms_data =  calc_rms(Data_input, Tail_start, Data_end); 
  
  for (int length = MinCycleLength, i = 0; length <= MaxCycleLength; length++, i++){
    //cyclelength
    cycle_statistics[i][0]=length;
    
    calc_AvgCycle (Data_input, Tail_start, Data_end, length, temp_AvgCycle, Phase_start);
    
    //MSE 
    cycle_statistics[i][1]= calc_MSE_AvgCycle(Data_input, Tail_start, Data_end, temp_AvgCycle, Phase_start, length);
    
    //RMS delta -> should be the same for all transforms! (Test)
    //cycle_statistics[i][2]=calc_rms_Transform(temp_AvgCycle, length, Phase_start, Tail_start, Data_end);
    //cycle_statistics[i][2]=(rms_data-cycle_statistics[i][2])/((rms_data+cycle_statistics[i][2])/2);
    
    //Gamma corr
    //cycle_statistics[i][3]=calc_gamma_corr(Data_input, Tail_start, Data_end, temp_AvgCycle, Phase_start, length);
    
  }
  
  /* Select the Cycle with the lowest MSE on tie the shortest one*/
  
  int optimal_cyclelength= (int) cycle_statistics[0][0];
  double min_MSE = cycle_statistics[0][1];
  for (int i = 0; i < MaxCycleLength-MinCycleLength+1; i++){
    if (cycle_statistics[i][1]<min_MSE){
      optimal_cyclelength=(int)cycle_statistics[i][0];
      min_MSE = cycle_statistics[i][1];
    }
  }
  
  sprintf(PLOGTXT,"The optimal cycle length detected was: %i with MSE of %.15g",optimal_cyclelength,min_MSE);
  WriteToLog(PLOGTXT,Data_end);
  
  return optimal_cyclelength;
   

}

bool Tail_is_stationary_Res (double Data_Input[], int Tail_start, int Tail_end, double AvgCycle_Input[], int Cyclelength, int Phase_start){
         sprintf(PLOGTXT,"\tTail_is_stationary RESIDUAL is tested with Tail_start %i",Tail_start);
         WriteToLog(PLOGTXT,Tail_end); 
 int signed_ranks;
 
 calc_Residual(Data_Input, Tail_end, AvgCycle_Input, Cyclelength, Phase_start, SPE_Residual);
 
 SPE_Wilcoxon =calc_Wilcoxon_signed_rank(SPE_Residual, Tail_start, Tail_end, &signed_ranks);
 double ranks = (double)signed_ranks; 
 SPE_Wilcoxon_zVal =   (SPE_Wilcoxon-0.5)/ ( sqrt( (ranks*(ranks+1)*(2*ranks+1))/6   ) );
 
 sprintf(PLOGTXT,"!!!!:\tChecking for stationarity at time %i.Wilcoxon is %15.f and z-level is %15.f. Ranks were %i",Tail_end,SPE_Wilcoxon,SPE_Wilcoxon_zVal,signed_ranks);
 WriteToLog(PLOGTXT,Tail_end);
 
 return  SPE_Wilcoxon_zVal>SPE_Wilcoxon_Critical_zValue;

}

bool Tail_is_stationary (double Data_Input[], int Tail_start, int Tail_end){
         sprintf(PLOGTXT,"\tTail_is_stationary is tested with Tail_start %i",Tail_start);
         WriteToLog(PLOGTXT,Tail_end); 
 int signed_ranks;
 SPE_Wilcoxon =calc_Wilcoxon_signed_rank(Data_Input, Tail_start, Tail_end, &signed_ranks);
 double ranks = (double)signed_ranks; 
 SPE_Wilcoxon_zVal =   (SPE_Wilcoxon-0.5)/ ( sqrt( (ranks*(ranks+1)*(2*ranks+1))/6   ) );
 
 sprintf(PLOGTXT,"!!!!:\tChecking for stationarity at time %i.Wilcoxon is %15.f and z-level is %15.f. Ranks were %i",Tail_end,SPE_Wilcoxon,SPE_Wilcoxon_zVal,signed_ranks);
 WriteToLog(PLOGTXT,Tail_end);
 
 return  SPE_Wilcoxon_zVal>SPE_Wilcoxon_Critical_zValue;

}

bool Periodicity_is_significant (double Data_Input[], int Tail_start, int Tail_end, double AvgCycle_Input[], int Phase_start, int Cyclelength){
         WriteToLog("\tPeriodicity_is_significant is tested",Tail_end);
         

SPE_Gamma = calc_gamma_corr(Data_Input, Tail_start, Tail_end, AvgCycle_Input, Phase_start, Cyclelength);
          /* We need to check if the gamma is critical, also. */
          
          
          WriteToLog("\tcalc_MSE_AvgCycle_vs_Mean is tested",Tail_end);
 SPE_relativeMSE = calc_MSE_AvgCycle_vs_Mean(Data_Input, Tail_start, Tail_end, AvgCycle_Input, Phase_start, Cyclelength);
          sprintf(PLOGTXT,"\t 1-BestCycleMSE/simpleMeanMSE is: %g and Gamma Correlation is %g",SPE_relativeMSE,SPE_Gamma);
          WriteToLog(PLOGTXT,Tail_end);                          
return (SPE_relativeMSE > SPE_Critical_relMSE); 
}



void calc_opt_Cycle (){
/* This procedure calculates the optimal cycle. */
   WriteToLog("::::calc_opt_Cycle() called.",SPE_TempTailEnd);
  /* a) Conduct an ACF analysis for to find the prior of the dominant cycle */

  SPE_TempCycleLength = find_optCyclelength(SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, SPE_MaxCycleLength, SPE_ACF);
        sprintf(PLOGTXT,"\ta) Prior to Cycle length via ACF analyis calculated as: %i with tail start at %i",SPE_TempCycleLength,SPE_TempTailStart);
        WriteToLog(PLOGTXT,SPE_TempTailEnd);
        
  /* b) Calculate optimal avg cycle based on prior */
  SPE_TempCycleLength = find_optAvgCycle(SPE_ACF, SPE_TempTailEnd-SPE_TempTailStart+1, SPE_TempCycleLength, SPE_Data, SPE_TempTailStart, SPE_TempTailEnd);
       sprintf(PLOGTXT,"\tb) Refined (optimal) Cycle length calculated as: %i with tail start at %i",SPE_TempCycleLength,SPE_TempTailStart);
       WriteToLog(PLOGTXT,SPE_TempTailEnd); 
       
       /* sub b) And save it (atm recalc)*/
       
        calc_AvgCycle(SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, SPE_TempCycleLength, SPE_AvgCycle, SPE_AvgCycle_Phase_start);
        calc_AvgCycle_StdDev (SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, SPE_TempCycleLength, SPE_AvgCycle, SPE_AvgCycle_Phase_start, SPE_AvgCycle_StdDev);
  /* c) Check if the cycle is a significant improvement vs. simple mean */
        if (SPE_TempCycleLength > 1){
           if (!Periodicity_is_significant(SPE_Data, SPE_TempTailStart, SPE_TempTailEnd, SPE_AvgCycle, SPE_AvgCycle_Phase_start, SPE_TempCycleLength)){
            //Change Cycle to pseudo one-cycle
            for (int i = 1; i< SPE_TempCycleLength; i++){
              SPE_AvgCycle[0]+= SPE_AvgCycle[i];
            }
            SPE_AvgCycle[0]/= SPE_TempCycleLength;
            SPE_TempCycleLength = 1;
            SPE_AvgCycle_Phase_start = SPE_TempTailStart;
          }
        }
}

