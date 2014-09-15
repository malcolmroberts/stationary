/* ========================================================================== */
/*                                                                            */
/*   PeriodicStochasticEquilibrium_coreAlgorithms.cpp                         */
/*   (c) 2014 Klabunde / Roberts / Schaff                                     */
/*   05.09.2014 14:01:34                                                      */
/*                                                                            */
/*   Core algorithms for the project                                          */
/*                                                                            */
/*                                                                            */
/* ========================================================================== */
/* ========================================================================== */
/* Declarations */

void calc_ACF(double ACF_output[], double ACF_input[],int ACF_input_start, int ACF_input_end, int Mode, int MaxNumberOfLags);
int Phase_Info (int Phase_start, int Cyclelength, int Info_time);
double calc_MSE_AvgCycle(double Data_input[], int Tail_start, int Tail_end, double AvgCycle[], int Phase_start, int Cyclelength);
double calc_MSE_AvgCycle_vs_Mean(double Data_input[], int Tail_start, int Tail_end, double AvgCycle[], int Phase_start, int Cyclelength);
double calc_gamma_corr(double Data_input[], int Tail_start, int Tail_end, double AvgCycle[], int Phase_start, int Cyclelength);
int calc_dominant_cycle_length(double ACF_input[], int ACF_bins, double ACF_Absolute_Sensitivity, double ACF_Relative_Sensitivity, int Mode, int maxCycleLength);   
void calc_cyclelength_interval(double ACF_input[], int ACF_bins, int ACF_dominantCycle, double ACF_RelativeNeighbourhood, double ACF_Absolute_Sensitivity, int Mode, int maxCycleLength, int & MinCyclelength, int & MaxCyclelength);
void calc_AvgCycle (double Data_input[], int Tail_start, int Tail_end, int Cyclelength, double AvgCycle_output[], int & Phase_start);
void calc_AvgCycle_StdDev (double Data_input[], int Tail_start, int Tail_end, int Cyclelength, double AvgCycle[], int Phase_start, double AvgCycleStdDev_output[]);
double calc_rms(double Input_Data[], int Input_start, int Input_end);
double calc_rms_Transform(double AvgCycle_Data[], int Cyclelength, int Phase_start, int Transform_start, int Transform_end);
void calc_Residual(double Input_Data[], int Data_end, double Input_AvgCyle[], int Cyclelength, int Phase_start, double Output_Residual[]);
void calc_MSER_completeRes(double Input_Data[], int Data_end, double Input_AvgCyle[], int Cyclelength, int Phase_start, double Output_Statistic[]);
void calc_MSER_complete(double Input_Data[], int Data_end, double Output_Statistic[]);
double calc_Wilcoxon_signed_rank (double Input_Data[],int Data_start, int Data_end, int *signed_number);
/* Algorithms */
                                                                                                              
void calc_ACF(double ACF_output[], double ACF_input[],int ACF_input_start, int ACF_input_end, int Mode, int MaxNumberOfLags){
/* Calculate the autocorrelation.  
in: ACF_input[], ACF_input_start, ACF_input_end
out: ACF_output[]

Only calculates the ACF for the given number of Lags. All rest is pur equal to zero.
*/
  //WriteToLog("\t:CORE:\tNow in calc_ACF()",SPE_TempTailEnd);
  //Number of bins
  int  ACF_bins =  ACF_input_end-ACF_input_start +1;
  int ACF_bins_used = min(ACF_input_end-ACF_input_start +1,MaxNumberOfLags);
  
  //Temporary variables
  double temp_mean_base,temp_mean_lagged, temp_var_base,temp_var_lagged;
  double temp_autocovar,temp_count; //    
    
  if (Mode == 0){    
    /* MODE 0 Calculation of the statistical autocorrelation function (ACF), see  
    mathworld.wolfram.com/StatisticalCorrelation.html for the correlation 
    function. The ACF is the correlation of the time-series (A) with a lagged 
    version of itself (B), starting with no lag (always 1) up to n-1 lags.
    for A truncated at the end, such that it has the same length as B. The 
    overlap of the series is thus always N-lag. Therefore it is, in addition, 
    only useful to consider the first half of the time-series for otherwise the
    part neglected in the study (lag) is bigger than half of the complete data.
    */
    
    /* The special cases, where the variance of one or both input series (laged 
    and non-lagged time-series) are 0 are treated as follows: if only one series 
    has a zero variance, the ACF is -1.5, if both are zero (a perfect 
    equilibrium for both), it is +1.5. */  
    
    
    //for each lag from 0 to N-1 we will test.
      for (int lag = 0; lag < ACF_bins_used; lag++){
      
      temp_mean_base = 0;
      temp_mean_lagged = 0;
      temp_var_base = 0;
      temp_var_lagged = 0;
      temp_autocovar = 0; //
      temp_count = 0;
      
      
        //First, calculate means:
        //base
        for (int bin = ACF_input_start; bin <= ACF_input_end-lag; bin++){
        temp_mean_base += ACF_input[bin];
        temp_mean_lagged += ACF_input[bin+lag]; 
        temp_count++;     
        }
        temp_mean_base/=temp_count;
        temp_mean_lagged/=temp_count;
        
        //Then variances (unnormalised e.g. only sum) & Autocovar:       
        //11.09.2014 19:51:06 Here was the error...
        for (int bin = ACF_input_start; bin < ACF_input_end-lag; bin++){
          temp_var_base += pow(ACF_input[bin]-temp_mean_base,2);
          temp_var_lagged += pow(ACF_input[bin+lag]-temp_mean_lagged,2);
          temp_autocovar += (ACF_input[bin]-temp_mean_base)*(ACF_input[bin+lag]-temp_mean_lagged);
        }
    
        
        //Finally, calculate autocorrelation
        
        if (temp_var_base != 0 && temp_var_lagged != 0){
          ACF_output[lag]= temp_autocovar / (sqrt(temp_var_base) * sqrt(temp_var_lagged));
        } else if (temp_var_base == 0 && temp_var_lagged == 0){
          ACF_output[lag]=+1.5; //not defined
        } else {
          ACF_output[lag]=-1.5; //the same but different / none at all
        }
      }
      
      for (int lag = ACF_bins_used+1; lag < ACF_bins; lag++){
        ACF_output[lag]=0.0;
      }                              
      
  } else {
  
  //ToDo
  
  
  }


}


int Phase_Info (int Phase_start, int Cyclelength, int Info_time){
/* The Function Calculates the Phase Info for the Info_time */
  ////WriteToLog("\t:CORE:\tNow in Phase_Info()",SPE_TempTailEnd);  
 int phase_start = Phase_start%Cyclelength; //lowest time with 0 info.
 int info_time = Info_time%Cyclelength;
 if (info_time < phase_start){
  info_time += Cyclelength;
 }
 /* Turn the clock */
 int phase_info = 0;
  while (info_time != phase_start){
  phase_info++;
  phase_start++;
  if (phase_info == Cyclelength){
    phase_info = 0;
  }  
  } 

/*  
  for (int lag = 0, time = ; time <= Info_time; time++;lag++){
  if ()
 }



  while (Info_time <= Phase_start){
    Info_time += Cyclelength; 
  }
  
  while (Phase_start+Cyclelength <= Info_time){
    Phase_start += Cyclelength;
  }                        
  
  
  for (phase_info = 0; phase_info < Cyclelength; phase_info++){
    if (Phase_start + phase_info == Info_time ){
    break;
    }
  } 
 */
 return phase_info;   
}
    
double calc_MSE_AvgCycle(double Data_input[], int Tail_start, int Tail_end, double AvgCycle[], int Phase_start, int Cyclelength){
/* Calculates the mean square error of prediction */
  //WriteToLog("\t:CORE:\tNow in calc_MSE_AvgCycle()",SPE_TempTailEnd);
  double MSE_pred = 0;
  
  for (int bin = Tail_start, lag = Phase_Info(Phase_start, Cyclelength,Tail_start); bin <= Tail_end; bin++,lag++){
    if (lag == Cyclelength){
      lag = 0;
    }
    MSE_pred += pow(Data_input[bin]-AvgCycle[lag],2);
  }
  MSE_pred /= Tail_end - Tail_start +1;
  
  return MSE_pred;

} 

double calc_MSE_AvgCycle_vs_Mean(double Data_input[], int Tail_start, int Tail_end, double AvgCycle[], int Phase_start, int Cyclelength){
/* Calculates the mean square error of prediction */
  //WriteToLog("\t:CORE:\tNow in calc_MSE_AvgCycle()",SPE_TempTailEnd);
  double MSE_mean = 0, MSE_cycle = 0;
  
  double Tail_mean = 0;
  for (int bin = 0; bin < Cyclelength; bin++){
    Tail_mean += AvgCycle[bin];
  }
  Tail_mean /= Cyclelength;
  
  for (int bin = Tail_start, lag = Phase_Info(Phase_start, Cyclelength,Tail_start); bin <= Tail_end; bin++,lag++){
    if (lag == Cyclelength){
      lag = 0;
    }
    MSE_cycle += pow(Data_input[bin]-AvgCycle[lag],2);
    MSE_mean += pow(Data_input[bin]-Tail_mean,2); 
  }
  MSE_cycle /= Tail_end - Tail_start +1;
  MSE_mean /= Tail_end - Tail_start +1;
  
  return 1-(MSE_cycle+0.0000001)/(MSE_mean+0.0000001);  //to cope with perfect equilibria.

}


double calc_gamma_corr(double Data_input[], int Tail_start, int Tail_end, double AvgCycle[], int Phase_start, int Cyclelength){
/*Calculates the simple gamma correlation function for a given time-series and the known average cycle plus the Phase-info at the start of the tail.*/
  //WriteToLog("\t:CORE:\tNow in calc_gamma_corr()",SPE_TempTailEnd);

/*Is only usefull if the Cyclelength is at least 2, obviously*/

  //temporary variables

  //The epsilon-tolerance is equal to the min of the std.dev. of the bin and the difference between two bins/2.
  int concordant=0, discordant=0;
  double x0,x1,y0,y1,gamma;
  
  for (int bin = Tail_start, lag = Phase_Info(Phase_start, Cyclelength,Tail_start), nextlag = lag+1; bin < Tail_end; bin++, lag++, nextlag++){
  
    
    if (lag == Cyclelength){
      lag = 0;
    }
    if (nextlag == Cyclelength){
      nextlag = 0;
    }
  
    x0 = Data_input[bin];
    x1 = Data_input[bin+1];
    y0 = AvgCycle[lag];
    y1 = AvgCycle[nextlag];
    
    
    if ( (x0-x1 > 0 && y0-y1 >0) || (x0-x1 < 0 && y0-y1 <0) ) {
      concordant++;
    } else {
      discordant++;
    }   
            
  }
  
  
  if (concordant+discordant != 0){
      gamma = (double)(concordant - discordant) / (double)(concordant + discordant) ;
  } else{
      gamma = 0.0;
  }

  return gamma;
}



int calc_dominant_cycle_length(double ACF_input[], int ACF_bins, double ACF_Absolute_Sensitivity, double ACF_Relative_Sensitivity, int Mode, int maxCycleLength){
/* Search for maximum corelation (with sensitivity) -> this is the dominant cycle length. Period. Output: Dominant cycle length*/
  //WriteToLog("\t:CORE:\tNow in calc_dominant_cycle_length()",SPE_TempTailEnd);

  /* IMPORTANT NOTE: The "ACF_bins" parameter should reflect that only "
  significant" ACF is considered (on whatever basis). This is not very important 
  if the normalised ACF MODE 1 (see calc_ACF) is used, for the ACF is then 
  "corrected" by a f(lag) dimnishing later ACF. In the MODE 0 ACF Version this 
  is not the case and only the first half of the ACF should be considered. In 
  order to make certain the correct MODE is selected (currently only MODE 0 
  implemented) an additional argument is provided. If MODE 0 is selected, only 
  the first half of the bins will be considered. */   
  
  /* Second Note: The ACF MODE 0 as implemented allows for ACF of 1.5 if a 
  perfect equilibrium is found. This is also recognised for there will be no 
  zero-crossing in this case.   */
  
  /* Third Note: If we end up in a condition where we want to find the optimal 
  Cycle GIVEN a specific minimum number of cycles, we need to define the maximum 
  cycle length. This should be kept in mind. In practise that could be because 
  we are interested in a specific time-frame and therefore don't bother to much 
  about very long cycles, ALTHOUGH ignoring these long cycles will ultimately   
  increase the residual of our transform - data and thus decrease precision. I 
  added the argument maxCycleLength to allow taking control of this.*/  
  
  /* Fourth Note (09.09.2014 09:29:33): We also need to check that, 
  for the complete ACF, there are several (at least two) cycles. To do so, 
  we truncate the ACF at the last zero_switch.
  */

  double temp_cor_max=0;
  bool first_switch = false;
  bool last_switch = false; 
  int last_switch_position; 
  int first_switch_position;

  int temp_dominant_cycle_length = 1; //No cycle
  double temp_dominant_cycle_AC = 0; 
                          
  int last_bin_checked = 0;
  
  if (Mode == 0){
      /* First, search for last zero crossing from top, if any. This would be 
      after the last complete cycle */
      for (int bin = ACF_bins-1; bin > 0; bin--){
        if (ACF_input[bin]<0 && ACF_input[bin-1]>0){
          last_switch_position = bin;
          last_switch = true;
          break;
        }
       } 
                            //Only the second half of the ACF is considered. This is already taken care of in the callee function for only the relevant ACF is calculated to start with.
     last_bin_checked = last_switch_position;
     
                            //Also the callee takes care that the number of bins under scrutiny is not bigger than the maximal cycle length. 
     
     
    /*The procedure is simple. First (a), we check for when the first zero- 
    crossing occurs for the dominant mode must be afterwards. Then we (b) find 
    the maximum AC. This is our prior assumption of the dominant mode. #1
    
    #1: I think it is not problematic that in case of a limited maxCycleLength 
    we might consider a mode as dominant that is, if we had checked the next 
    view bins, not a local maximum. The reason is (c) , see below.
    */ 
    
    if (last_switch){
      //WriteToLog("\t:---:\tcalc_dominant_cycle_length()\tLast Switch Ok",SPE_TempTailEnd);
        
      for (int bin = 1; bin <= last_bin_checked; bin++){
        if (ACF_input[bin-1]>0 && ACF_input[bin]<0 && !first_switch){
          first_switch_position = bin;
          first_switch = true;   
        }    
        
        if (first_switch){
          if (ACF_input[bin]>temp_cor_max && ACF_input[bin]>=ACF_Absolute_Sensitivity){
          temp_cor_max = ACF_input[bin];
          temp_dominant_cycle_length = bin;             
          temp_dominant_cycle_AC = temp_cor_max;
          }    
        }
      }
       if (!first_switch){                
       //WriteToLog("\t:---:\tcalc_dominant_cycle_length()\t No first Switch found",SPE_TempTailEnd);
       }
      
      
      /* (c) Now that we know a dominant cycle might exist (if there was a 
      first_switch AND the ACF was above the absolute significance level), we will 
      check if there is a fast/shorter Cycle that is not significantly different 
      from the maximum and consider this one as our cycle. The first local maximum 
      within range will be chosen. This is unproblematic for we will later also 
      define a range of potential cyclelength around this local maximum. (see 
      int calc_optimal_cyclelength() )*/ 
      
      if (first_switch && temp_cor_max >= ACF_Absolute_Sensitivity && first_switch_position != last_switch_position) {
            //WriteToLog("\t:---:\tcalc_dominant_cycle_length()\tFirst Switch Ok",SPE_TempTailEnd);
      bool inRange = false;
      temp_cor_max = 0.0;
      double Threshold_AC =  max (ACF_Absolute_Sensitivity, temp_dominant_cycle_AC* (1-ACF_Relative_Sensitivity));
      
  
        for (int bin = first_switch_position; bin <= temp_dominant_cycle_length; bin++){  
          if (ACF_input[bin]>Threshold_AC) {
            inRange = true;
          }
          if (inRange){   //region near local maximum entered
            if (ACF_input[bin] > temp_cor_max){
              temp_cor_max = ACF_input[bin];
              temp_dominant_cycle_length = bin;
            } else {
              break; //region near local maximum left, the dominant cycle is found.
            }
          }
        }
      } else {
        temp_dominant_cycle_length = 1;   //Pseudo one cycle.
      }
    
    /* (d) For the case that no dominant cycle has been found, we need to check if 
    perhaps we have a perfect equilibrium in which case the ACF of MODE 0 was 
    defined +1.5 for lag 1 and the (pseudo)one-cycle is perfect. If this is not 
    the case, we will instead give back 0. */
    
    if (temp_dominant_cycle_length == 0 && ACF_input[1]==1.5){
      temp_dominant_cycle_length = 1;
      //WriteToLog("\t:!!!!!!:\tcalc_dominant_cycle_length()\tPERFECT EQ FOUND.",SPE_TempTailEnd);
    } 
  }
  
  //Mode 1
  } else {
  
  //ToDo!
  }
  
  
  
  //We do now know the dominant cycle length.
  //WriteToLog("\t:---:\tcalc_dominant_cycle_length()\tAt end",SPE_TempTailEnd);
  return temp_dominant_cycle_length;
}

                                              

void calc_cyclelength_interval(double ACF_input[], int ACF_bins, int ACF_dominantCycle, double ACF_RelativeNeighbourhood, double ACF_Absolute_Sensitivity, int Mode, int maxCycleLength, int & MinCyclelength, int & MaxCyclelength){
/* This procedure calculates the minimum and maximum cyclelength for 
consideration as follows: Using the Information of the dominant Cycle length and 
the ACF, the algorithm searches to the left (right) of the dominant cycle and 
changes Min(Max)Cyclelength to the according lag (=cyclelength) until it gets +
out of range. The Range is defined by the ACF at the lag=ACF_dominantCycle*(1-
ACF_RelativeNeighbourhood) or ACF_Absolute_Sensitivity, whatever is lower.;
*/
  //WriteToLog("\t:CORE:\tNow in calc_cyclelength_interval()",SPE_TempTailEnd);
/* The mode is again considered, for with MODE 0 we would only want to check the 
first half of the ACF function.
In Addition, the maximumCycleLength is also considered by necessity (if any).*/
    
  double local_Threshold = max(ACF_input[ACF_dominantCycle]*(1-ACF_RelativeNeighbourhood),ACF_Absolute_Sensitivity);
  MinCyclelength = ACF_dominantCycle;
  MaxCyclelength = ACF_dominantCycle;
    //Search left                      
  for (int bin = ACF_dominantCycle; bin >0; bin--){
    if (ACF_input[bin]>local_Threshold){
      MinCyclelength = bin;
    } else {
      break;
    }
  }
  
  
  int last_bin_checked = maxCycleLength;
  
  if (Mode == 0){
    last_bin_checked = min (maxCycleLength , (int) (floor(ACF_bins/2) ) );
  }
  
    //Search right
  for (int bin = ACF_dominantCycle; bin <= last_bin_checked; bin++){
    if (ACF_input[bin]>local_Threshold){
      MaxCyclelength = bin;
    } else {
      break;
    }
  }
 
}

void calc_AvgCycle (double Data_input[], int Tail_start, int Tail_end, int Cyclelength, double AvgCycle_output[], int & Phase_start){                                             
/* Calculates the average cycle based on the tail and the cycle length. Works 
also for the "one-cycle" special case. */
  //WriteToLog("\t:CORE:\tNow in calc_AvgCycle()",SPE_TempTailEnd);
  /* For now the Phase at t=Tail_start is 0. We might later want to change the 
  Phase_start if the tail changes but the avg_cycle is not recalculated. The 
  Function Phase_Info() allows for a calculation of the correct phase of the avg 
  cycle for any given time.*/
      
    Phase_start = Tail_start; 
    
   
   
    /* We might have partial cycles covers. Therefore we need to count the 
    number of datapoints we have for each bin seperately. Also the arrays need 
    to be initialised. */
     
    int temp_count[Cyclelength];
    for (int lag = 0; lag < Cyclelength; lag++){
          temp_count[lag]=0;               
          AvgCycle_output[lag]=0;  
    }
    
    /* Next we simply sum up over the lag-bins and divide by the lag-count. */
    for (int bin = Tail_start, lag = 0; bin <= Tail_end; bin++,lag++){
      if (lag == Cyclelength){
        lag = 0;             
      }
      temp_count[lag]++;
      AvgCycle_output[lag]+=Data_input[bin]; 
    }
    

    for (int lag = 0; lag < Cyclelength; lag++){
      AvgCycle_output[lag]/=temp_count[lag];
    }

}



void calc_AvgCycle_StdDev (double Data_input[], int Tail_start, int Tail_end, int Cyclelength, double AvgCycle[], int Phase_start, double AvgCycleStdDev_output[]){
/* This algorithm simply calculates the std. dev. of the avg cycle with the 
given data, e.g. for each single bin. Works also for the "one-cycle" special 
case.*/
  //WriteToLog("\t:CORE:\tNow in calc_AvgCycle_StdDev()",SPE_TempTailEnd);
/* In my outputs I printed the std.dev. between the avg cycle and the tail as shaded area. */

    double temp_count[Cyclelength];
    //Initialise array
    for (int lag = 0; lag < Cyclelength; lag++){
    AvgCycleStdDev_output[lag]=0;
    temp_count[lag]=0;
    
    }
    

    for (int bin = Tail_start, lag = Phase_Info(Phase_start, Cyclelength,Tail_start); bin <= Tail_end; bin++,lag++){
      if (lag == Cyclelength){
        lag = 0;             
      }
      temp_count[lag]++;
      AvgCycleStdDev_output[lag]+=pow(AvgCycle[lag]-Data_input[bin],2); 
    }
    
    //normalise
    for (int lag = 0; lag < Cyclelength; lag++){
      AvgCycleStdDev_output[lag]=sqrt(AvgCycleStdDev_output[lag]/temp_count[lag]);
    }

}                                              

double calc_rms(double Input_Data[], int Input_start, int Input_end){
/* Simply calculates the RMS of the given data in the given intervall */  
  //WriteToLog("\t:CORE:\tNow in calc_rms()",SPE_TempTailEnd);  
    double rms = 0;
  for (int bin = Input_start; bin <= Input_end; bin++){
    rms += pow(Input_Data[bin],2);               
  }
   rms /= (Input_end-Input_start+1);
   rms = sqrt(rms);
   
   return rms;
}

double calc_rms_Transform(double AvgCycle_Data[], int Cyclelength, int Phase_start, int Transform_start, int Transform_end){
/* Calculates the RMS of the Transform composed of the Avg Cycle over the range given. Controlls for the Phase. */ 
  //WriteToLog("\t:CORE:\tNow in calc_rms_Transform()",SPE_TempTailEnd); 
  double rms = 0;
  
  for (int bin = Transform_start, lag = Phase_Info (Phase_start, Cyclelength, Transform_start); bin <= Transform_end; bin++, lag++){
    if (lag==Cyclelength){
      lag = 0;
    }
    rms += pow(AvgCycle_Data[lag],2);               
  }
   rms /= (Transform_end-Transform_start+1);
   rms = sqrt(rms);
   
   return rms;

}

void calc_Residual(double Input_Data[], int Data_end, double Input_AvgCyle[], int Cyclelength, int Phase_start, double Output_Residual[]){
/*Simply calculates the residual of the avg cycle and the data.*/

  for (int i = 1, phase = Phase_Info(Phase_start,Cyclelength,1); i<=Data_end; i++, phase++ ){
    if (phase == Cyclelength){
      phase = 0;
    }
    Output_Residual[i]=Input_Data[i]-Input_AvgCyle[phase];  
}

}


void calc_MSER_completeRes(double Input_Data[], int Data_end, double Input_AvgCyle[], int Cyclelength, int Phase_start, double Output_Statistic[]){
/* Calculate the MSER Statistic for given Data */
Output_Statistic[0]=SPE_NotANumber;
double avg_sum=0;
int tail_bins = 0;
double MSER_denom;
double MSER;
double temp_tailMean[Data_end+1];
double temp_resData[Data_end+1];

/* The MSER is the MSE of the constant Regression of the current Tail. Thuss effectively we calculate the tail-mean for all tails and the MSE of according to this tail-avg and save it. The result is weighted by the 1/(Tail_bins^2). */


  for (int i = Data_end, phase = Phase_Info(Phase_start,Cyclelength,Data_end); i>0; i--, phase--){
    if (phase == -1){
      phase = Cyclelength -1; 
    }
    tail_bins++;
    temp_resData[i]=Input_Data[i]-Input_AvgCyle[phase];
    if (i == Data_end){
      temp_tailMean[i]=temp_resData[i];
      avg_sum = temp_resData[i]; 
      Output_Statistic[i] = SPE_NotANumber;
    } else {
    
    avg_sum+=temp_resData[i];
    temp_tailMean[i]=avg_sum/tail_bins;
      
       
    MSER_denom = pow(Data_end-i+1,2);
    MSER = 0;
    for (int bin = i;bin<=Data_end;bin++){
      MSER += pow(temp_resData[bin]-temp_tailMean[i],2);     //12.09.2014 10:13:27 Error found (i/bin switched)
    }
    MSER/=MSER_denom;
    Output_Statistic[i]=MSER;
    }
  }

}
  
void calc_MSER_complete(double Input_Data[], int Data_end, double Output_Statistic[]){
/* Calculate the MSER Statistic for given Data */
Output_Statistic[0]=SPE_NotANumber;
double avg_sum=Input_Data[Data_end];
int tail_bins = 1;
double MSER_denom;
double MSER;
double temp_tailMean[Data_end+1];
temp_tailMean[Data_end]=Input_Data[Data_end];

  for (int i = Data_end-1; i>0; i--){
    tail_bins++;
    avg_sum+=Input_Data[i];
    temp_tailMean[i]=avg_sum/tail_bins;
      
       
    MSER_denom = pow(Data_end-i+1,2);
    MSER = 0;
    for (int bin = i;bin<=Data_end;bin++){
      MSER += pow(Input_Data[bin]-temp_tailMean[bin],2);
    }
    MSER/=MSER_denom;
    Output_Statistic[i]=MSER;
  }
}  

 
double calc_Wilcoxon_signed_rank (double Input_Data[],int Data_start, int Data_end, int *signed_ranks){
  /* Using the Wikipedia set-up. For an odd samplesize last value is neglected. 
  The remaining Data is split in half*/
  
  /* First we split the data in two and calculate the difference between pairs(x2-x1) and 
  the sign of each and store botch in new arrays. */
  int bins_sample = (int) floor(Data_end-Data_start+1); 
  double temp_Data_abs[bins_sample],temp_Data_sign[bins_sample];
  for (int i = 0; i< bins_sample; i++){
    temp_Data_abs[i]=Input_Data[Data_start+i+bins_sample]-Input_Data[Data_start+i];
    if (temp_Data_abs[i]<0) {
      temp_Data_sign[i]=-1;  
    } else if (temp_Data_abs[i]>0){
      temp_Data_sign[i]=1;
    } else {
      temp_Data_sign[i]=0;
    }
    temp_Data_abs[i]=abs(temp_Data_abs[i]);
  }
  
  /* Next we eliminate those elements which have no sign (were the same). 
  We simply add only those pairs (value,sign) to the new arrays which have a sign 
  != 0 and count the number for this is our new array size (the rest we will not 
  make use of). */
  int signed_number = 0;
  double temp_Data_abs2[bins_sample],temp_Data_sign2[bins_sample];
  for (int i = 0; i < bins_sample; i++){
    if (temp_Data_sign[i]!=0){
      temp_Data_abs2[signed_number]=temp_Data_abs[i];
      temp_Data_sign2[signed_number]=temp_Data_sign[i];
      signed_number++;
    } 
  }
  
  /* Now we need to first order the pairs (value,sign) and next rank them. The 
  ordering is done with a switch algorithm, e.g. places of subsequent values are 
  switched as long as necessary. */
  bool no_change_in_order = false;
  double temp_sign, temp_value;
  int bin, compareBin;
  while (!no_change_in_order){
    for (bin = 1, no_change_in_order=true; bin < signed_number; bin++){
      compareBin = bin -1;
      
      if (temp_Data_abs2[bin]<temp_Data_abs2[compareBin]){
        no_change_in_order = false;
        temp_value = temp_Data_abs2[compareBin];
        temp_sign = temp_Data_sign2[compareBin];
        temp_Data_abs2[compareBin]=temp_Data_abs2[bin];
        temp_Data_sign2[compareBin]=temp_Data_sign2[bin];
        temp_Data_abs2[bin]=temp_value;
        temp_Data_sign2[bin]=temp_sign; 
      }   
    
    }
  }
  
  /* Now we need to rank the pairs, taking care of ties. */
  double temp_rank[signed_number];
  
  double rank;
  int with_equal_ranks;
  int rankBin;
  for (bin = 0; bin < signed_number; ){
    rank = bin+1;
    with_equal_ranks = 1;
    
    /*First search which following bins have equal rank*/
    for (compareBin=bin+1;compareBin<signed_number;compareBin++){
      if (temp_Data_abs2[bin] == temp_Data_abs2[compareBin]){
        with_equal_ranks++;
        rank+=compareBin+1;
      } else {
        break;
      }
    }
    /* Calculate average rank */
    rank/=with_equal_ranks;
    /*Update rank statistic. CompareBin is now the first bin not covered.*/
    for (rankBin = bin; rankBin < compareBin;rankBin++){
      temp_rank[rankBin]=rank;
    }
    /* Update the next bin to be checked. */
    bin = compareBin;
  }
  
  /* Calculate Wilcoxon score */
  double Wilcoxon = 0;
  for (bin = 0; bin < signed_number; bin++){
    Wilcoxon += temp_Data_sign2[bin]*temp_rank[bin];
  }
  Wilcoxon = abs (Wilcoxon);
  
  *signed_ranks = signed_number;
  return Wilcoxon;

}



