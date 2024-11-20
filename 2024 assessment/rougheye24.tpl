///////////////////////////////////////////////////////////////////
// BSAI rougheye rockfish model  
//	Has asymptotic fishery selectivity curve
//   
//
// Naming Conventions:
//
//  GENERAL:
//    styr, endyr begining year and ending year of model (catch data available)
//    nages       number of age groups considered
//    nyrs_       number of observations available to specific data set
//
//  DATA SPECIFIC:
//    catch_bio   Observed fishery catch biomass
//    nyrs_srv2   Number of years of slope survey index value (annual)
//    obs_srv2    Observed trawl slope survey index value (annual)
//    nyrs_srv3   Number of years of AI trawl survey index value (annual)
//    obs_srv3   Observed trawl AI trawl survey index value (annual)
//
//    oac_fish    Observed age comp from the fishery
//    oac_srv2 	  Observed age comp from the trawl survey
//
//    pred_catch  Predicted catch from the fishery
//    pred_srv2   Predicted index of the trawl survey
//
//    eac_fish    Expected age comp from fishery data
//    eac_srv2	  Expected age comp from trawl survey.
//
//    sel_fish    selectivity for fishery                
//    sel_srv2    selectivity for the trawl survey
//
//////////////////////////////////////////////////////////////////////////////
DATA_SECTION //------------------------------------------------------------------------------------------------
  init_int styr				// start year for model
  init_int styr_fish			// start year for the fishery
  init_int endyr			// end year for model and fishery
  init_int yrs_r      // number of years to peel for retrospective run
  init_int nages      // number of ages modeled
  init_int nages_dat // number of ages in data, can be less than the modeled ages
  init_ivector ages(1,nages)    // vector of ages (for model)
  init_ivector ages_dat(1,nages_dat)  // vector of ages (for data) 
  init_int nselages			// number of ages with estimated selectivities
  init_int rec_age			// recruitment age
  init_int nlen				// number of length groups
  init_ivector lengths(1,nlen)          // vector of length groups
  init_int nyrs_fish			// number of years in the fishery
  init_ivector yrs_fish(1,nyrs_fish)    // vector of the index years in the fishery
  init_vector catch_bio(styr_fish,endyr)     // observed catch biomass
  init_int nsrv                              // number of surveys (for any type of survey data)
  init_adstring srvnameread;         // string with names of surveys, separated by "%"
  init_ivector nyrs_srv_bio(1,nsrv)                      // vector of number of years in the survey biomass estimates
  init_imatrix yrs_srv_bio(1,nsrv,1,nyrs_srv_bio)            // years of the survey(s) (ragged array)
  init_matrix obs_srv_bio(1,nsrv,1,nyrs_srv_bio)   // observed survey estimates (ragged array)
  init_matrix obs_srv_bio_sd(1,nsrv,1,nyrs_srv_bio)      // observed survey standard deviations (ragged array)
  init_matrix obs_srv_bio_lower(1,nsrv,1,nyrs_srv_bio)   // survey biomass -2SD(ragged array)
  init_matrix obs_srv_bio_upper(1,nsrv,1,nyrs_srv_bio)   // survey biomass +2SD (ragged array)

  init_ivector nyrs_srv_abun(1,nsrv)                      // vector of number of years for survey abundance estimates
  !! cout << nyrs_srv_abun << endl;
  init_imatrix yrs_srv_abun(1,nsrv,1,nyrs_srv_abun)        // years of the abundance estimates for the survey(s) (ragged array)
  init_matrix obs_srv_abun(1,nsrv,1,nyrs_srv_abun)         // observed survey abundance estimates (ragged array)
  init_matrix obs_srv_abun_sd(1,nsrv,1,nyrs_srv_abun)      // observed survey abundance standard deviations (ragged array)
  init_matrix obs_srv_abun_lower(1,nsrv,1,nyrs_srv_abun)   // survey abundance -2SD(ragged array)
  init_matrix obs_srv_abun_upper(1,nsrv,1,nyrs_srv_abun)   // survey abundance +2SD (ragged array)
  
  init_matrix unbiasedages(1,nages,1,nages) // transition age error matrix for unbiased ages
  

  init_matrix translen(1,nlen,1,nages)     // transition matrix from ages to lengths
  init_int nyrs_fish_unbiased_ac			// number of years with  unbiased fishery age comps
  !! cout << " nyrs_fish_unbiased_ac is " << nyrs_fish_unbiased_ac << endl;
  init_ivector yrs_fish_unbiased_ac(1,nyrs_fish_unbiased_ac) 	// vector of index years with unbiased fishery age comps
  !! cout << " yrs_fish_unbiased_ac is " << yrs_fish_unbiased_ac << endl;
  init_matrix oac_fish_unbiased(1,nyrs_fish_unbiased_ac,1,nages_dat) 	// observed unbiased fishery age comps
  init_int nyrs_fish_lc					// number of years with fishery length comps
  init_ivector yrs_fish_lc(1,nyrs_fish_lc) 		// vector of index years with fishery length comps
  init_matrix olc_fish(1,nyrs_fish_lc,1,nlen) 		// observed fishery length comps
  
  init_ivector nyrs_srv_ac(1,nsrv)                           // vector of number of years with age comps, per survey
  init_imatrix yrs_srv_ac(1,nsrv,1,nyrs_srv_ac)              // years with age comp estimates from survey(s) (ragged array)
  !! cout << yrs_srv_ac << endl;
  
  init_3darray obs_srv_ac(1,nsrv,1,nyrs_srv_ac,1,nages_dat)  // observed age comps from the survey(s) (ragged array)
  init_ivector nyrs_srv_lc(1,nsrv)                           // vector of number of years with length comps, per survey
  init_imatrix yrs_srv_lc(1,nsrv,1,nyrs_srv_lc)              // years with length comp estimates from survey(s) (ragged array)
  init_3darray obs_srv_lc(1,nsrv,1,nyrs_srv_lc,1,nlen)       // observed length comps from the survey(s) (ragged array)
  !! cout << obs_srv_lc << endl;
  
  

  init_vector wt_pop(1,nages)   // population weight at age 
  init_vector wt_fsh(1,nages)   // fishery weights at age
  vector wt_pop_bin(1,nages_dat)    // rescale weights to match the data plus group
  vector wt_fsh_bin(1,nages_dat)    // rescale weights to match the data plus group
  //init_vector maturity(1,nages) // Maturity at age
  //vector maturity_bin(1,nages_dat)  // bin maturity to match the data plus group
  //!! cout <<"wt_fsh is " << wt_fsh << endl;
  init_number spawn_mo          // spawning month
  init_int fyear_ac_option	// first year age comp option
  init_int historic_catch	// historic catch level
  init_int sr_type		// option for S-R curve (avg or B-H)
  init_int fixedrec_yrs		// the number of years from the endyr in which we fix the recruitments
  init_number sigr		// the s of the recruitment deviations  
  init_int    sel_option  // option for selectivity functional form (1=logistic, 2=double logistic, 3= bicubic spline)
  init_int    nbins       // the number of bins
  init_int    sel_styr    // the start year for fitting selecitvity
  init_int    sel_fixedyrs    //the number of years from the endyr_r in which we fix the recruitments
  init_ivector binstart(1,nbins)  // the start year for each bin
  init_number sigma_aslp    // sigma for selectivity acsending slope deviations
  init_number sigma_a50     // sigma for selectivity a50 deviations
  init_number sigma_dslp    // sigma for selectivity decsending slope deviations
  init_number sigma_d50     // sigma for selectivity d50 deviations
  init_int    n_yr_nodes      // number of year nodes for bicubic spline selectivity
  init_int    n_age_nodes     // number of age nodes for bicubic spline selectivity
  int         isel_npar       // generalized number of nodes for sel_par
  int         jsel_npar       // generalized number of nodes for sep_par
  init_vector srv_sel_option(1,nsrv)     // option for survey selectivity (vector, by survey)     
  init_vector priormean_q(1,nsrv)     // prior mean of trawl survey q (-9 if not used)
  init_vector priorcv_q(1,nsrv)   // prior cv of trawl survey q (-9 if not used)  
  init_number priormean_M	// prior mean of M
  init_number priorcv_M		// prior cv of M
  //  maturity data
  init_number nages_mat_ogive           // number of ages for maturity ogive
  init_vector matages_ogive(1,nages_mat_ogive) // maturity ages for ogive
  
  init_int  mat_input_switch  // switch for reading in maturity vector
  init_matrix mat_input(1,mat_input_switch,1,nages) // read in matutity curve if not estimated    
  //!! cout << "mat_input is "<< mat_input << endl;

  init_int    nmat_datasets       // number of maturity datasets
  init_ivector nages_mat(1,nmat_datasets)  // number of maturity ages for each data set 
  init_imatrix ages_mat(1,nmat_datasets,1,nages_mat) // maturity ages for each data set
  init_imatrix n_mat(1,nmat_datasets,1,nages_mat)  // fish measured for for each data set, by age
  init_imatrix y_mat(1,nmat_datasets,1,nages_mat)  // mature fish for for each data set, by age  
  init_vector mat_lambda(1,nages_mat_ogive)  // weights for fitting the maturity data 
  !! cout << "mat_lambda is "<< mat_lambda << endl; 




  init_int    fs_option		// option to constrain the fishery sel to close to the survey sel
  				// (1=yes, 0=no)
  init_number priorsdfishslopedev	// prior cv on deviation for fish sel slope parameter 
  init_number priorsdfisha50dev	// prior cv on deviation for fish sel 50% parameter
  !! cout <<" priorsdfisha50dev is  " << priorsdfisha50dev << endl;
  
  



  number spmo_frac		// spawning month as fraction of the year
  int styr_rec			// start year for recruits
  int styr_rec_dev  // start year for which we have recruitment deviations (some may be fixed)
  int lastyr_rec		// last year for which we estimate the recruitment
  int i             // index for years 
  int j				      // index for ages
  int k
  int s             // index for surveys
  int l
  int f
  int bincount       // for looping through the number of bins   
  int styr_fut			// Start year for projections
  int endyr_fut                 // End year for projections
  int num_proj_Fs 		// Number of Fs to evaluate in the future
  
  matrix cv_bio_srv(1,nsrv,1,nyrs_srv_bio)    // vector of CVs for the surveys
  matrix cv_abun_srv(1,nsrv,1,nyrs_srv_abun)    // vector of CVs for the survey abun
  
  int phase_fydev            // phase for first year deviations
  int sr_phase                  // phase for estimation of rzero and steepness
  int mat_phase                // phase for maturity estimation


  matrix tmp(1,nages,1,nages)   // tmp matrix for rescaling the age error matrix
  matrix tmp2(1,nages,1,nlen)   // tmp matrix for rescaling the transition matrix 

  imatrix psrvname(1,nsrv,1,2)  // dimensions of strings for the survey names
  vector ages_dat_mid(1,nages_dat)  // midpoint of the age bins, for Francis
  vector lengths_mid(1,nlen)        // midpoint of the length bins, for Francis
  !! ages_dat_mid = dvector(ages_dat) + 0.5;
  !! lengths_mid = dvector(lengths)+ 0.5; 

 
  
  // things for getting the survey selectivity a10 (as modified by M)
  int firstage;              // the first age that exceeds 10% survey selection
  int excludeage;            // exclude ages at or below this

  int endyr_r   // end year for retrospective run
  int sel_endyr // end year for estimating selectivity

  int comp_yr_count  // for counting the number of years to be used for compositon data in retrospective run
  int surv_yr_count  // for counting the number of years to be used for survey biomass indices
  int nyrs_fish_unbiased_ac_r  // number of years to use for retrospecitve run  
  int nyrs_fish_lc_r  // number of years to use for retrospecitve run
  vector nyrs_srv_bio_r(1,nsrv)  // number of years of the survey biomass indices to use for retrospecitve run
  vector nyrs_srv_abun_r(1,nsrv)  // number of years of the survey abundance indices to use for retrospecitve run

  ivector nyrs_srv_ac_r(1,nsrv)      // number of years of the survey age comps to use for retrospecitve run
  ivector num_srv_ac_resid_r(1,nsrv) // number of residuals of the survey age comps to use for retrospecitve run      
  ivector nyrs_srv_lc_r(1,nsrv)      // number of years of the survey length comps to use for retrospecitve run
  ivector num_srv_lc_resid_r(1,nsrv) // number of residuals of the survey length comps to use for retrospecitve run

 LOCAL_CALCS                               // Get the names of the surveys
  adstring_array CRLF;   // blank to terminate lines
  CRLF+="";
  
  for(f=1;f<=nsrv;f++) {psrvname(f,1)=1; psrvname(f,2)=1;}    // set whole array to equal 1 in case not enough names are read
  f=1;
    for(i=1;i<=strlen(srvnameread);i++)
      if(adstring(srvnameread(i))==adstring("%"))
        {psrvname(f,2)=i-1; f+=1;  psrvname(f,1)=i+1;}
    psrvname(nsrv,2)=strlen(srvnameread);
  for(f=1;f<=nsrv;f++)
  {
    srvname+=srvnameread(psrvname(f,1),psrvname(f,2))+CRLF(1);
  }
  cout<<" surveynames are "<<srvname<<endl;
 END_CALCS

 LOCAL_CALCS
  // calculate endyr_r and sel_year for the surveys
  endyr_r = endyr - yrs_r;
  sel_endyr = endyr_r - sel_fixedyrs;

  // cout the number of years of age comps to be used in the retrospective run
  comp_yr_count = 0;
  for (i=1;i<=nyrs_fish_unbiased_ac;i++)
      if (yrs_fish_unbiased_ac(i) <= endyr_r)  comp_yr_count++;  
  nyrs_fish_unbiased_ac_r = comp_yr_count;

  comp_yr_count = 0;
  for (i=1;i<=nyrs_fish_lc;i++)
      if (yrs_fish_lc(i) <= endyr_r)  comp_yr_count++;  
  nyrs_fish_lc_r = comp_yr_count;

  
  for (s=1;s<=nsrv;s++)
    {
      surv_yr_count = 0;
      for (i=1;i<=nyrs_srv_bio(s);i++)
          if (yrs_srv_bio(s,i) <= endyr_r)  surv_yr_count++;
      nyrs_srv_bio_r(s) = surv_yr_count;
      
      surv_yr_count = 0;
      for (i=1;i<=nyrs_srv_abun(s);i++)
          if (yrs_srv_abun(s,i) <= endyr_r)  surv_yr_count++;
      nyrs_srv_abun_r(s) = surv_yr_count;

      comp_yr_count = 0;
      for (i=1;i<=nyrs_srv_ac(s);i++)
          if (yrs_srv_ac(s,i) <= endyr_r)  comp_yr_count++;
      nyrs_srv_ac_r(s) = comp_yr_count;
      num_srv_ac_resid_r(s) = nyrs_srv_ac_r(s)*nages; 
      
      comp_yr_count = 0;
      for (i=1;i<=nyrs_srv_lc(s);i++)
          if (yrs_srv_lc(s,i) <= endyr_r)  comp_yr_count++;
      nyrs_srv_lc_r(s) = comp_yr_count;
      num_srv_lc_resid_r(s) = nyrs_srv_lc_r(s)*nlen; 

    }

  
 END_CALCS

  matrix rescaled_sel_fish(styr,endyr_r,1,nages_dat)  // rescaled fishery selectivity matrix
  vector rescaled_F(styr,endyr_r)                 // rescaled F values
  vector recent_fish_sel(1,nages_dat)               // the recent fish selectivity (for spr calcs)

  vector yrs_fish_unbiased_ac_r(1,nyrs_fish_unbiased_ac_r)  // years to use for retrospecitve run  
  vector yrs_fish_lc_r(1,nyrs_fish_lc_r)  // years to use for retrospecitve run
  matrix oac_fish_unbiased_r(1,nyrs_fish_unbiased_ac_r,1,nages_dat)  // observed unbiased fishery age comps, retrospective run
  matrix olc_fish_r(1,nyrs_fish_lc_r,1,nlen)     // observed fishery length comps, retrospective run
  matrix yrs_srv_ac_r(1,nsrv,1,nyrs_srv_ac_r)            // years to use for retrospective run, survey age comps
  3darray oac_srv_r(1,nsrv,1,nyrs_srv_ac_r,1,nages_dat)  // observed survey age comps, retrospective run
  matrix yrs_srv_lc_r(1,nsrv,1,nyrs_srv_lc_r)            // years to use for retrospective run, survey length comps
  3darray olc_srv_r(1,nsrv,1,nyrs_srv_lc_r,1,nlen)       // observed survey length comps, retrospective run


 LOCAL_CALCS
  // get the comp data and years for the retrospective run
  yrs_fish_unbiased_ac_r = yrs_fish_unbiased_ac(1,nyrs_fish_unbiased_ac_r);
  yrs_fish_lc_r = yrs_fish_lc(1,nyrs_fish_lc_r);
  
  for (i=1;i<=nsrv;i++){
    yrs_srv_ac_r(i) = yrs_srv_ac(i)(1,nyrs_srv_ac_r(i));
    yrs_srv_lc_r(i) = yrs_srv_lc(i)(1,nyrs_srv_lc_r(i));
  }

  for (i=1;i<=nyrs_fish_unbiased_ac_r;i++)
     oac_fish_unbiased_r(i) = oac_fish_unbiased(i); 
  for (i=1;i<=nyrs_fish_lc_r;i++)
      olc_fish_r(i) =  olc_fish(i);  
  for (i=1;i<=nsrv;i++){
    for (j=1;j<=nyrs_srv_ac_r(i);j++){
        oac_srv_r(i,j) = obs_srv_ac(i,j);
    }
    for (j=1;j<=nyrs_srv_lc_r(i);j++){
        olc_srv_r(i,j) = obs_srv_lc(i,j);
    }
  }


 END_CALCS
    

 LOCAL_CALCS    
  // bin the maturity and weight to match data plus group
  //maturity_bin(1,nages_dat-1) = maturity(1,nages_dat-1);
  //maturity_bin(nages_dat) = mean(maturity(nages_dat,nages));
  wt_pop_bin(1,nages_dat-1) = wt_pop(1,nages_dat-1);
  wt_pop_bin(nages_dat) = mean(wt_pop(nages_dat,nages)); 
  wt_fsh_bin(1,nages_dat-1) = wt_fsh(1,nages_dat-1);
  wt_fsh_bin(nages_dat) = mean(wt_fsh(nages_dat,nages)); 
  jsel_npar = 0;
  isel_npar = 0;
  spmo_frac = (spawn_mo-1)/12.;
  num_proj_Fs = 5;
  styr_fut=endyr_r+1;
  endyr_fut = styr_fut+10;
  // define styr_rec +++++++++++++++++++++++++++++++++++++++

  if (fyear_ac_option == 1) // first year rec are combined with other recruitments  
  {
   styr_rec = styr-nages+1;
   styr_rec_dev = styr-nages_dat+1;    // some cohorts share a recruitment deviation
   phase_fydev = -1;
  }
  else if (fyear_ac_option == 2) // first year recruitments are in equilibrium with historic catch
  { 
   styr_rec = styr;
   styr_rec_dev = styr;
   phase_fydev = -1;
  }
  else if (fyear_ac_option == 3) // first year recruitments are stochastic, but separate from other recs
  { 
   styr_rec = styr;
   styr_rec_dev = styr;
   phase_fydev = 3;
  }  
  lastyr_rec = endyr_r - fixedrec_yrs;   // define the last year for which we estimate recruitment

  if (sr_type==1) sr_phase=-1;
    else sr_phase = 2;

  if (mat_input_switch==1)  mat_phase = -1;
    else mat_phase = 5;


  tmp = trans(unbiasedages);
  for (i=1;i<=nages;i++) tmp(i) = tmp(i)/sum(tmp(i));
  unbiasedages = trans(tmp);
  
  tmp2 = trans(translen);
  for (i=1;i<=nages;i++) tmp2(i) = tmp2(i)/sum(tmp2(i));
  translen = trans(tmp2);
  
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // calculate cv's for the surveys

  
  for (s=1;s<=nsrv;s++)
  {
    cv_bio_srv(s) = elem_div(obs_srv_bio_sd(s),obs_srv_bio(s));
    cv_abun_srv(s) = elem_div(obs_srv_abun_sd(s),obs_srv_abun(s));
  }

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

  // start to read from the control file
  ad_comm::change_datafile_name("bsai_re.ctl");
 END_CALCS
  init_int phase_selcoff
  init_int phase_logist_sel
  init_int phase_f_sel_param
  init_ivector phase_s_sel_aslope(1,nsrv)
  init_ivector phase_s_sel_a50(1,nsrv)
  init_ivector phase_s_sel_mu(1,nsrv)
  init_ivector phase_s_sel_dist(1,nsrv)
  init_ivector phase_s_sel_sig1(1,nsrv)
  init_ivector phase_s_sel_sig2(1,nsrv)
  init_ivector phase_s_sel_q(1,nsrv)
  
  
  init_int phase_proj
  init_int phase_historic_F
  init_int har_flag 
  init_vector srv_bio_flag(1,nsrv)      // flag to fit biomass estimates from the surveys
  init_vector srv_abun_flag(1,nsrv)     // flag to fit abundance estimates from the surveys  
  init_vector srv_age_flag(1,nsrv)      // flag to fit the age comp data from the surveys
  init_vector srv_len_flag(1,nsrv)      // flag to fit the length comp data from the surveys 
  init_int fish_bio_flag                // flag to fit catch biomass  
  init_int fish_age_flag                // flag to fit the proportional catch age comp data from the surveys
  init_int fish_len_flag                // flag to fit the proportional length comp data from the surveys 


  init_vector raw_fish_unbiased_ac_samp(1,nyrs_fish_unbiased_ac)
  init_vector raw_fish_lc_samp(1,nyrs_fish_lc)
  init_matrix raw_srv_ac_samp(1,nsrv,1,nyrs_srv_ac)   // ragged array of survey age sample sizes
  init_matrix raw_srv_lc_samp(1,nsrv,1,nyrs_srv_lc)   // ragged array of survey length sample sizes
  
  init_vector lambda(1,10)
  !! cout << lambda << endl;

 LOCAL_CALCS 
  // start to read from the initial values file 
  ad_comm::change_datafile_name("initvalues.ctl");
 END_CALCS

  // read in the starting parameter values
  init_number M_start
  !! cout << "M_start is "  << M_start << endl;
  init_number mat_beta1_start
  init_number mat_beta2_start
  init_number mean_log_rec_start
  init_number log_rzero_start
  init_number log_rinit_start
  init_number log_avg_fmort_start
  init_number sel_aslope_fish_start
  init_number sel_dslope_fish_start
  init_vector sel_aslope_srv_start(1,nsrv)
  init_number sel_a50_fish_start
  init_number sel_d50_fish_start
  init_vector sel_a50_srv_start(1,nsrv)
  init_number steepness_start
  init_number historic_F_start
  init_vector q_srv_start(1,nsrv)
  init_vector sel_srv_mu_start(1,nsrv)
  init_vector sel_srv_dist_start(1,nsrv)
  init_vector sel_srv_sig1_start(1,nsrv)
  init_vector sel_srv_sig2_start(1,nsrv)
  init_int    jitter_seed   // seed for jittering the starting values
  init_number jitter_stddev // jittering standard deviation
  !! cout << "sel_srv_sig2_start" << endl;
  !! cout << sel_srv_sig2_start << endl;

  vector jitterstdnorm(1,100)         // vector of standard normal errors for jittering
  vector jittererr(1,100)             // vector of jitter errors (expected mean of zero)
  int err_count                       // for jittering
  
  
  // vectors for sample sizes for retrospective runs
  vector fish_unbiased_ac_samp_r(1,nyrs_fish_unbiased_ac_r)
  vector fish_lc_samp_r(1,nyrs_fish_lc_r)
  
  // matrices for sample sizes for retrospective runs
  matrix srv_ac_samp_r(1,nsrv,1,nyrs_srv_ac_r)
  matrix srv_lc_samp_r(1,nsrv,1,nyrs_srv_lc_r)  




 LOCAL_CALCS 
  // start to read from the composition weight file
  // order is fac, flc, sac, slc
  ad_comm::change_datafile_name("compweights1.ctl");
 END_CALCS
  init_vector compweights_fsh(1,2)
  init_vector compweights_sac(1,nsrv)
  init_vector compweights_slc(1,nsrv)
  !! cout <<" compweights_slc are "  << compweights_slc << endl;  
 
 LOCAL_CALCS 
   
    if (nyrs_fish_unbiased_ac_r>0)
        fish_unbiased_ac_samp_r = compweights_fsh(1)*raw_fish_unbiased_ac_samp(1,nyrs_fish_unbiased_ac_r);
    if (nyrs_fish_lc_r>0)
        fish_lc_samp_r = compweights_fsh(2)*raw_fish_lc_samp(1,nyrs_fish_lc_r);
    
    for (i=1;i<=nsrv;i++){
       if(nyrs_srv_ac_r(i)>0)
           srv_ac_samp_r(i) = compweights_sac(i)*raw_srv_ac_samp(i)(1,nyrs_srv_ac_r(i));
       if(nyrs_srv_lc_r(i)>0)
           srv_lc_samp_r(i) = compweights_slc(i)*raw_srv_lc_samp(i)(1,nyrs_srv_lc_r(i));
    }

   
 END_CALCS 
  
  



  int phase_f_sel_ascend;
  int phase_f_sel_descend;
  int phase_f_sel_par;
  int phase_a50_devs;
  int phase_aslope_devs;
  int phase_d50_devs;
  int phase_dslope_devs;
  ivector binindex(sel_styr,sel_endyr);
  vector scal_yr_nodes(1,n_yr_nodes);  // the yr nodes scaled from 0 to 1
  vector scal_age_nodes(1,n_age_nodes);  // the age nodes scaled from 0 to 1

 LOCAL_CALCS  // set phases for fishery and survey selectivity options, and time-varying deviations
  //  phase for the ascending parameters for logistic
  //  phase for the descending parameters for the logistic
  //  phase for the age and year nodes for the bicubic spline
  //  phase for ebs catchabilty and selectivity (turn off if they are not being fit)

  if (fish_age_flag==0 & fish_len_flag==0)  // turn off fishery selectivity parameters if not fitting fishery comp data
     {
      phase_f_sel_param = -1; 
     }
  
  // compute binindex for time-varying selectivity
    if (nbins > 1)
    {
      bincount = 1;
      for (i=sel_styr; i<=sel_endyr; i++)
          {
            if (i == binstart(bincount)) 
            {
             binindex(i) = bincount; 
             if( bincount < nbins ) bincount++;   
            }
            if (i < sel_endyr) binindex(i+1) = binindex(i);
          }
           cout << binindex << endl;
    }
  

  // set phases to -1, and later turn on the ones we need; 
  //   selectivity parameters
  phase_f_sel_ascend = -1;
  phase_f_sel_descend = -1;
  phase_f_sel_par = -1;
  //   time-varying parameters for logistic and double logistic
  phase_a50_devs = -1;
  phase_aslope_devs = -1;
  phase_d50_devs = -1;
  phase_dslope_devs = -1;
  if (sel_option==1)  // logistic fishery selectivity 
  {
    phase_f_sel_ascend = phase_f_sel_param;
    if(nbins>1)
      {
        phase_a50_devs = 4;
        phase_aslope_devs = 4;
      }

  }
  else if (sel_option==2)  // double logistic fishery selectivity 
  {
    phase_f_sel_ascend = phase_f_sel_param;
    phase_f_sel_descend = phase_f_sel_param;
    if(nbins>1)
      {
        phase_a50_devs = 4;
        phase_aslope_devs = 4;
        phase_d50_devs = 4;
        phase_dslope_devs = 4;
      }

  }
  else if (sel_option==3)  // bicubic fishery selectivity 
  {
    phase_f_sel_par = phase_f_sel_param;
    scal_age_nodes.fill_seqadd(0,1./(n_age_nodes-1));
    scal_yr_nodes.fill_seqadd(0,1./(n_yr_nodes-1));
    isel_npar = n_yr_nodes;
    jsel_npar = n_age_nodes; 
  }
  else if (sel_option==4)  // time-varying cubic fishery selectivity 
  {
    phase_f_sel_par = phase_f_sel_param;
    isel_npar = n_age_nodes;
    jsel_npar = nbins; 
  }

  for (i=1;i<=nsrv;i++)   
   {
     if(srv_sel_option(i)==1)   // turn off paramters for double normal if not used
       {phase_s_sel_mu(i)=-1; phase_s_sel_dist(i)=-1; phase_s_sel_sig1(i)=-1; phase_s_sel_sig2(i)=-1; }
     if(srv_sel_option(i)==2)   // turn off logistic survey paramters if not used
       {phase_s_sel_a50(i)=-1; phase_s_sel_aslope(i)=-1;}  
     // turn off all survey catchability and selectivity paramaters if the index (either biomass or abundance index) is not being used
     if(srv_bio_flag(i)==0 && srv_abun_flag(i)==0) 
       {phase_s_sel_aslope(i)=-1; phase_s_sel_a50(i)=-1;phase_s_sel_q(i)=-1;phase_s_sel_mu(i)=-1; phase_s_sel_dist(i)=-1;
          phase_s_sel_sig1(i)=-1; phase_s_sel_sig2(i)=-1; }         
   }

 END_CALCS 

  

INITIALIZATION_SECTION //-------------------------------------------------------------------------------------
  M  M_start
  mat_beta1 mat_beta1_start
  mat_beta2 mat_beta2_start
  mean_log_rec mean_log_rec_start
  log_rzero log_rzero_start
  log_rinit log_rinit_start
  log_avg_fmort log_avg_fmort_start
  sel_aslope_fish sel_aslope_fish_start
  sel_dslope_fish sel_dslope_fish_start
  sel_aslope_srv sel_aslope_srv_start
  sel_a50_fish sel_a50_fish_start
  sel_d50_fish sel_d50_fish_start
  sel_a50_srv sel_a50_srv_start
  steepness steepness_start
  historic_F historic_F_start
  q_srv q_srv_start
  sel_srv_mu sel_srv_mu_start
  sel_srv_dist sel_srv_dist_start
  sel_srv_sig1 sel_srv_sig1_start
  sel_srv_sig2 sel_srv_sig2_start


  

PARAMETER_SECTION //-----------------------------------------------------------------------------------------

 //init_number aa
 // offset parameters
 vector offset(1,6)
 vector sac_offset(1,nsrv)  // offset for survey age comps, by survey
 vector slc_offset(1,nsrv)  // offset for survey length comps, by survey
 // selectivity parameters 
 //  First, the logistic curve parameters for the domestic and foreign fisheries(slope and 50% parameters)
 init_bounded_number sel_aslope_fish(0.1,3.0,phase_f_sel_ascend)  
 init_number sel_dslope_fish(phase_f_sel_descend)

 init_bounded_number sel_a50_fish(0.1,30.0,phase_f_sel_ascend)
 init_number sel_d50_fish(phase_f_sel_descend)
 // fishery bicubic selectivity matrix 
 init_matrix sel_par(1,jsel_npar,1,isel_npar,phase_f_sel_par);
  
 //init_number slopedev(phase_logist_param)  // deviations in the fish sel parameters from the  
 //init_number a50dev(phase_logist_param)    // survey sel parameters

 // now the logistic curve parameters for the surveys
 init_bounded_number_vector sel_aslope_srv(1,nsrv,0.1,3.0,phase_s_sel_aslope)
 init_bounded_number_vector sel_a50_srv(1,nsrv,0.1,30.0,phase_s_sel_a50) 
 
 // parameters for the double normal
 init_bounded_number_vector sel_srv_mu(1,nsrv,0,50,phase_s_sel_mu)
 init_bounded_number_vector sel_srv_dist(1,nsrv,0,50,phase_s_sel_dist)
 init_bounded_number_vector sel_srv_sig1(1,nsrv,0.1,100,phase_s_sel_sig1)
 init_bounded_number_vector sel_srv_sig2(1,nsrv,0.1,100,phase_s_sel_sig2)
 //init_bounded_number_vector sel_dslope_srv(1,nsrv,-3.0,0.1,phase_s_sel_dslope)
 //init_bounded_number_vector sel_d50_srv(1,nsrv,0.1,30.0,phase_s_sel_d50)

 // deviations on fishery selectivity
 init_bounded_dev_vector a50_devs(1,nbins,-10.,10.,phase_a50_devs)
 init_bounded_dev_vector aslope_devs(1,nbins,-10.,10.,phase_aslope_devs)
 init_bounded_dev_vector d50_devs(1,nbins,-10.,10.,phase_d50_devs)
 init_bounded_dev_vector dslope_devs(1,nbins,-10.,10.,phase_dslope_devs)
 
 number a_slptmp
 number a50tmp
 number d_slptmp
 number d50tmp



 //init_bounded_number sel_dslope_srv3(-3.0,0.1,3)
 //init_bounded_number sel_d50_srv3(0.1,30.0,3)
 
 // Or we could estimate the selectivity coefficients directly (not used so far)
 //init_vector log_selcoffs_domfish(1,nselages,phase_selcoff)
 //init_vector log_selcoffs_fish(1,nselages,phase_selcoff)
 //number avgsel_fish;  // (averages are used in the likelihood function (if used))

 // Next, the selectivity values (logged and unlogged) 
 matrix log_sel_fish(sel_styr,sel_endyr,1,nselages)
 matrix sel_fish(styr,endyr_r,1,nages)
 
 matrix log_sel_srv(1,nsrv,1,nselages)     // log of the selectivities for all surveys
 matrix sel_srv(1,nsrv,1,nages)         // unlogged selectivities for all surveys 

 

 // The survival parameters
 init_bounded_number M(0.001,.5,4)
 //init_bounded_number M(0.001,.5,-1)
 number surv
 //init_number log_avg_fmort(-1)
 init_number log_avg_fmort(1)
 init_bounded_dev_vector fmort_dev(styr_fish,endyr_r,-10,10,2)
 //init_bounded_dev_vector fmort_dev(styr_fish,endyr_r,-10,10,-1)
 number avg_fmort_dev
 matrix F(styr_rec,endyr_r,1,nages)
 matrix Z(styr_rec,endyr_r,1,nages)
 matrix S(styr_rec,endyr_r,1,nages)
 matrix mort(styr_rec,endyr_r,1,nages);  // first, the the multiplier for the mean population size
 matrix spawner_S(styr_rec,endyr_r,1,nages);

 // The numbers at age parameters
 init_bounded_dev_vector rec_dev(styr_rec_dev,lastyr_rec,-10,10,2) // recruitment deviations (styr to endyr)
 init_number mean_log_rec(1)			// mean recruitment
 
 //init_bounded_dev_vector rec_dev(styr_rec_dev,lastyr_rec,-10,10,-1) // recruitment deviations (styr to endyr)
 //init_number mean_log_rec(-1)      // mean recruitment

 init_number log_rinit(1)           // initial recruitment to start model
 //init_number log_rinit(-1)           // initial recruitment to start model
 
 matrix natage(styr_rec,endyr_r+1,1,nages)		// numbers at age
 matrix natage_bin(styr_rec,endyr_r+1,1,nages_dat)     // bins the plus group to match the data
 vector natagetmp(1,nages)   // temporary numbers at age (prior to start year)
 init_bounded_dev_vector fydev(2,nages_dat,-10,10,phase_fydev)    // deviations around eq numbers in first year (stage 3)
 //init_bounded_dev_vector fydev(2,nages_dat,-10,10,-1)    // deviations around eq numbers in first year (stage 3)
  

 //vector spawners(styr,endyr) 			// estimated biomass of female spawners
 //vector expbiom(styr,endyr) 			// estimated exploitable biomass
 sdreport_vector totbiom(styr_rec,endyr_r+1)		// total biomass of population
 //vector totbiom9899(styr_rec-rec_age,endyr_r+1)        // total biomass of the 98-99 year classes
 //sdreport_vector totbiomwo9899(styr_rec-rec_age,endyr_r+1)  // total biomass without 98-99 year classes  
//sdreport_number depletion			// change in totbiom over time
 //sdreport_number endbiom			// totbiom in final year
 init_number historic_F(phase_historic_F)	// historic F for computing first year age comps
   
// The parameters for evaluating the objective function
 objective_function_value obj_fun
 vector rec_like(1,3)
 vector srv_bio_like(1,nsrv)
 vector srv_abun_like(1,nsrv)
 number catch_like
 number fpen
 number hf_pen
 vector age_like(1,6)
 vector age_like_sac(1,nsrv)
 vector age_like_slc(1,nsrv)
 matrix fish_effn(1,2,1,80)
 matrix sac_effn(1,nsrv,1,80)
 matrix slc_effn(1,nsrv,1,80)
 number rec_rmse
 vector srv_bio_rmse(1,nsrv)
 vector srv_abun_rmse(1,nsrv)
 vector sel_like(1,9)
 number sprpen
 number prior_M
 vector prior_q(1,nsrv)
 vector survey_bio_sdnr(1,nsrv)
 vector survey_abun_sdnr(1,nsrv)
 vector fish_sdnr(1,2)
 vector sac_sdnr(1,nsrv)
 vector slc_sdnr(1,nsrv)
 vector fish_rmse(1,2)
 vector sac_rmse(1,nsrv)
 vector slc_rmse(1,nsrv)
 vector fish_mean_effn(1,2)
 vector sac_mean_effn(1,nsrv)
 vector slc_mean_effn(1,nsrv)
 vector fish_mean_samp_wts(1,2)
 vector sac_mean_samp_wts(1,nsrv)
 vector slc_mean_samp_wts(1,nsrv)
 number mat_like  
   

  matrix pred_srv_bio(1,nsrv,styr_fish,endyr_r)     // predicted survey biomass
  matrix pred_srv_abun(1,nsrv,styr_fish,endyr_r)     // predicted survey numbers
  init_bounded_number_vector  q_srv(1,nsrv,0.1,10.0,phase_s_sel_q)  // q for each of the surveys
 // The parameters for getting the predicted values for the catch total
  vector pred_catch(styr_fish,endyr_r)  
  number ehc		   			// the estimated historic catch
  matrix catage(styr_fish,endyr_r,1,nages)

// The parameters for getting the predicted age comps
  matrix eac_fish_unbiased_mod(styr_fish,endyr_r,1,nages)
  matrix eac_fish_unbiased_dat(styr_fish,endyr_r,1,nages_dat)
  matrix elc_fish(styr_fish,endyr_r,1,nlen)
  matrix elc_fish_tempages(styr_fish,endyr_r,1,nages)  // the true ages for the year in which we have length comps
  
  3darray eac_srv_mod(1,nsrv,styr,endyr_r,1,nages)       // predicted survey age comps
  3darray eac_srv_dat(1,nsrv,styr,endyr_r,1,nages_dat)   // predicted survey age comps
  3darray elc_srv(1,nsrv,styr,endyr_r,1,nlen)            // predicted survey length comps
  3darray elc_srv_tempages(1,nsrv,styr,endyr_r,1,nages)  // the true survey ages for the year in which we have length comps 



// The parameters for getting the SPR values
   number F40
   number F35
   number F30
   number SB0
   number SBF40
   number SBF35
   number SBF30

// numbers for getting the normalized residuals
 matrix  survey_bio_nr(1,nsrv,1,nyrs_srv_bio)   //  survey biomass normalized residuals
 matrix  survey_abun_nr(1,nsrv,1,nyrs_srv_abun)   //  survey abundance normalized residuals  
 vector  fac_nr(1,nyrs_fish_unbiased_ac_r*nages)   // fishery age comps
 vector  flc_nr(1,nyrs_fish_lc_r*nlen) 		   // fishery length comps
 matrix  sac_nr(1,nsrv,1,num_srv_ac_resid_r)   
 matrix  slc_nr(1,nsrv,1,num_srv_lc_resid_r)

 // number for getting the pearson residuals for the comp data
 3darray sac_pearson(1,nsrv,1,nyrs_srv_ac_r,1,nages_dat)  
 3darray slc_pearson(1,nsrv,1,nyrs_srv_lc_r,1,nlen)           
 matrix  flc_pearson(1,nyrs_fish_lc_r,1,nlen)
 matrix fac_pearson(1,nyrs_fish_unbiased_ac_r,1,nages_dat)
 
 
 vector  fac_mcian_wgt(1,nyrs_fish_unbiased_ac_r)   // McAllister-Ianelli weights (method TA1.1 in Francis 2011)
 vector  flc_mcian_wgt(1,nyrs_fish_lc_r)   
 matrix  sac_mcian_wgt(1,nsrv,1,nyrs_srv_ac_r)
 matrix  slc_mcian_wgt(1,nsrv,1,nyrs_srv_lc_r)



 vector  fac_mcian_wgt_inv(1,nyrs_fish_unbiased_ac_r)   // inverse of McAllister-Ianelli weights for harmonic mean
 vector  flc_mcian_wgt_inv(1,nyrs_fish_lc_r)   
 matrix  sac_mcian_wgt_inv(1,nsrv,1,nyrs_srv_ac_r)
 matrix  slc_mcian_wgt_inv(1,nsrv,1,nyrs_srv_lc_r)

 vector fac_nr_fran(1,nyrs_fish_unbiased_ac_r)       // normalized residuals from the Francis method (method TA1.8 In Francis 2011)
 vector flc_nr_fran(1,nyrs_fish_lc_r)       
 matrix sac_nr_fran(1,nsrv,1,nyrs_srv_ac_r)
 matrix slc_nr_fran(1,nsrv,1,nyrs_srv_lc_r)

 

 vector  fac_resid(1,nyrs_fish_unbiased_ac_r*nages)   // fishery age comps
 vector  flc_resid(1,nyrs_fish_lc_r*nlen)       // fishery length comps
 matrix  sac_resid(1,nsrv,1,num_srv_ac_resid_r)  // survey age comps
 matrix  slc_resid(1,nsrv,1,num_srv_lc_resid_r)  // survey length comps

 vector  fac_pearson_vec(1,nyrs_fish_unbiased_ac_r*nages)   // fishery age comps, pearson resids, as vector
 vector  flc_pearson_vec(1,nyrs_fish_lc_r*nlen)             // fishery length comps, pearson resids, as vector 
 matrix  sac_pearson_vec(1,nsrv,1,num_srv_ac_resid_r)        // survey age comps, pearson resids, as vector
 matrix  slc_pearson_vec(1,nsrv,1,num_srv_lc_resid_r)       // survey length comps, pearson resids, as vector


// Parameters for doing the future projections
 //matrix nage_future(styr_fut,endyr_fut,1,nages)
 //number ftmp
 //number mean_recent_fs
 //matrix Z_future(styr_fut,endyr_fut,1,nages)
 //matrix F_future(styr_fut,endyr_fut,1,nages)
 //matrix S_future(styr_fut,endyr_fut,1,nages)
 //init_vector rec_dev_future(styr_fut,endyr_fut,phase_proj)
 //matrix catage_future(styr_fut,endyr_fut,1,nages)
 //sdreport_matrix catch_future(1,num_proj_Fs-1,styr_fut,endyr_fut)
 //sdreport_matrix biomass_future(1,num_proj_Fs,styr_fut,endyr_fut)
 //sdreport_matrix ssb_future(1,num_proj_Fs,styr_fut,endyr_fut)
 
// Stock recruitment params 
   init_number log_rzero(sr_phase)
   number bzero
   number rzero
   //sdreport_vector sp_biom(styr_rec-rec_age,endyr_r)
   sdreport_vector sp_biom(styr_rec,endyr_r)
   vector exp_biom(styr,endyr_r)
   sdreport_vector est_rec(styr,lastyr_rec-rec_age)    // for SR fit, indexed to year class
   //sdreport_vector est_rec(styr_rec_dev,lastyr_rec)    // for SR fit, indexed to year class
   number alpha   
   number beta
   init_bounded_number steepness(0.20001,1.0,sr_phase)
   vector pred_rec(styr,lastyr_rec-rec_age)   // prediction from SR model, index to year class
   //vector pred_rec(styr_rec_dev,lastyr_rec)   // prediction from SR model, index to year class
   vector est_spb(styr,lastyr_rec-rec_age)
   //vector est_spb(styr_rec_dev,lastyr_rec)
   vector chi(styr,lastyr_rec-rec_age)   // the squared difference between est and pred recruits
   //vector chi(styr_rec_dev,lastyr_rec)   // the squared difference between est and pred recruits
   number sumrecdev // the sum of the lognormal deviations for the recruitments
   vector SRec_spawn(1,20) // data for estimated SR curve
   vector SRec_rec(1,20) //  data for estimated SR curve
   vector xdum2(styr,endyr_r)

// vector for projection data file
   vector Mvec(1,nages)  // natural mortality, repeated for each age

// maturity estimation
   init_bounded_number mat_beta1(-10,2,mat_phase)   // beta0 and beta1 parameters for logistic regression of maturity curve
   init_bounded_number mat_beta2(0,2,mat_phase)
   vector mat_theta(1,nages_mat_ogive)   // theta parameter estimates for logistic regression of maturity curve (uses all  ages for which we have data)
   vector maturity(1,nages) // the maturity for the ages used in the population model
   vector maturity_bin(1,nages_dat) // maturity ogive, binned to match the data plus group


// updated compweights
   vector compweightsnew_ta12_fsh(1,2)     // weight by inverse of variance of normalized resid
   vector compweightsnew_ta11_fsh(1,2)     //  McAllister-Ianelli method
   vector compweightsnew_ta18_fsh(1,2)     // the Francis method

// updated survey age compweights   
   vector compweightsnew_ta12_sac(1,nsrv)     // weight by inverse of variance of normalized resid
   vector compweightsnew_ta11_sac(1,nsrv)     //  McAllister-Ianelli method
   vector compweightsnew_ta18_sac(1,nsrv)     // the Francis method

// updated survey length compweights   
   vector compweightsnew_ta12_slc(1,nsrv)     // weight by inverse of variance of normalized resid
   vector compweightsnew_ta11_slc(1,nsrv)     //  McAllister-Ianelli method
   vector compweightsnew_ta18_slc(1,nsrv)     // the Francis method


       
PRELIMINARY_CALCS_SECTION //-------------------------------------------------------------------------------
  // Compute the offsets for the multinomial distributions

  for ( i=1;i<=nyrs_fish_unbiased_ac_r; i++)  // fishery unbiased age comps
 {
 	oac_fish_unbiased_r(i) = oac_fish_unbiased_r(i)/sum(oac_fish_unbiased_r(i)); // make sure age comps add to 1.0 for each year
 	offset(1)-=fish_unbiased_ac_samp_r(i)*(oac_fish_unbiased_r(i))*log(1.e-13+oac_fish_unbiased_r(i)); //get the negative log like
 }
 for ( i=1;i<=nyrs_fish_lc_r; i++)  // fishery length comps
 {
 	olc_fish_r(i) = olc_fish_r(i)/sum(olc_fish_r(i)); // make sure age comps add to 1.0 for each year
 	offset(2)-=fish_lc_samp_r(i)*(olc_fish_r(i))*log(1.e-13+olc_fish_r(i)); //get the negative log like
 }
 
 for (i=1;i<=nsrv;i++){
  for (j=1;j<=nyrs_srv_ac_r(i);j++){
    oac_srv_r(i,j) =  oac_srv_r(i,j)/sum(oac_srv_r(i,j));  // make sure age comps add to 1.0 for each year
    sac_offset(i) -= srv_ac_samp_r(i,j)*(oac_srv_r(i,j))*log(1.e-13+oac_srv_r(i,j)); // get the negative log like
  }
  for (j=1;j<=nyrs_srv_lc_r(i);j++){
    olc_srv_r(i,j) =  olc_srv_r(i,j)/sum(olc_srv_r(i,j));  // make sure length comps add to 1.0 for each year
    slc_offset(i) -= srv_lc_samp_r(i,j)*(olc_srv_r(i,j))*log(1.e-13+olc_srv_r(i,j)); // get the negative log like            
  }
 }


 
PROCEDURE_SECTION //-----------------------------------------------------------------------------------------
// example of using FUNCTION to structure the procedure section
    get_maturity();
    //cout << " got mat" << endl;
    get_selectivity();
    //cout <<" got sel" << endl;
    get_mortality();
    //cout <<" got mort " << endl;
    get_first_year();
    //cout << " got first year " << endl;
    get_numbers_at_age();
    //cout << " got natage" << endl;
    get_expected_values();
    //cout <<" got expected values " << endl;
    get_sr_inputs();
    //cout << " got sr inputs" << endl;
    get_catch_at_age();
    //cout <<" got catch at age" << endl;
    get_age_comps();
    //cout <<"got age comps"  << endl;
    get_binned();
    //cout << " got binned" << endl;
    
//    if (active(F40))
//      compute_spr_rates();
//    if (current_phase()>6)
//      future_projections();
    evaluate_the_objective_function();
    //cout <<" got the obj funct " << endl;
    comp_metrics();
    get_age10();
    //cout <<" got age 10 " << endl;
    update_compweights();
    //cout <<" got the updatd comp weights" << endl;
  
     if (mceval_phase())
  {
    ofstream evalout("evalout_re.prj", ios::app);
    evalout << obj_fun<<" "<<q_srv(1)<<" "<<M<<" "<<totbiom<<" "<<sp_biom<<" "<<est_rec<<endl;
  } 
  
FUNCTION get_maturity   // get the maturity ogive
 if(mat_input_switch==1)
  { 
   maturity = mat_input(1);
  }
 else
   {  
   mat_theta = elem_div(mfexp(mat_beta1 + mat_beta2*matages_ogive),(1. + mfexp(mat_beta1 + mat_beta2*matages_ogive)));
   maturity = mat_theta(1,nages);
   maturity(nages) = 0.5*(maturity(nages) +  mfexp(mat_beta1 + mat_beta2*100)/(1. + mfexp(mat_beta1 + mat_beta2*100)));
   }


FUNCTION get_selectivity //---------------------------------------------------------------------------------
 // Calculate the selectivity parameters (if they are being used)
 // First the AI logistic selectivity (only if being used)
 dvariable diff;  // for double normal
 dvariable diff2;  // for double normal


  if (current_phase()>=phase_logist_sel)
 {
 
 //  sel_aslope_fish = sel_aslope_srv3*exp(slopedev);  // allow the fish sel parameters to
 //  sel_a50_fish = sel_a50_srv3*exp(a50dev);		// vary a little from the survey sel

  for (s=1;s<=nsrv;s++)
  {
    for (j=1;j<=nselages;j++)  // Get the AI and EBS survey selectivty
      {
        if(srv_sel_option(s)==1)  // logistic survey selectivity
          {
            log_sel_srv(s,j) = -1.*log((1.0+mfexp(-1.0*sel_aslope_srv(s)*(ages(j)-sel_a50_srv(s)))));
          }  
        else if (srv_sel_option(s)==2)   // double normal survey selectivity

          {
            //A = 1./(sqrt(2*PI)*(sel_srv_sig1(s)+sel_srv_sig2(s))/2.0);
            diff = ages(j)-sel_srv_mu(s);  
            diff2 = ages(j)-(sel_srv_mu(s)+sel_srv_dist(s)); 

            if(ages(j)<=sel_srv_mu(s))
              {
               log_sel_srv(s,j)= log(mfexp(-(diff*diff)/(2*sel_srv_sig1(s)*sel_srv_sig1(s))));
              }
            else if(ages(j)>=(sel_srv_mu(s)+sel_srv_dist(s)))
              {
               log_sel_srv(s,j) =log(mfexp(-(diff2*diff2)/(2*sel_srv_sig2(s)*sel_srv_sig2(s))));
              }
            else
              {
               log_sel_srv(s,j) = 0.0; 
              }    
            
            
            
          }       
    }  
  }

     
 if (sel_option==1)
 {
  for (i=sel_styr; i<=sel_endyr; i++)  // Use the selectivity values for the foreign fishery through 1988 
  {
    if (nbins==1) 
    {
     for (j=1;j<=nselages;j++)
       {
        log_sel_fish(i,j) = -1.*log((1.0+mfexp(-1.0*sel_aslope_fish*(ages(j)-sel_a50_fish))));
       }
    }
    else if (nbins>1)
    {
     a_slptmp = sel_aslope_fish*exp(aslope_devs(binindex(i)));
     a50tmp = sel_a50_fish*exp(a50_devs(binindex(i)));
     for (j=1;j<=nselages;j++)
       {
        log_sel_fish(i,j) = -1.*log((1.0+mfexp(-1.0*a_slptmp*(ages(j)-a50tmp))));
       }

    }  
  }
 }
 else if (sel_option==2)
 {
  for (i=sel_styr; i<=sel_endyr; i++)  // Use the selectivity values for the foreign fishery through 1988 
  {
    if (nbins==1) 
    {
     for (j=1;j<=nselages;j++)
       {
        log_sel_fish(i,j) = -1.*log(1.0+mfexp(-1.0*sel_aslope_fish*(ages(j)-sel_a50_fish)))
                 -1.*log(1.0+mfexp(-1.0*sel_dslope_fish*(ages(j)-sel_d50_fish)));
       }
    }
    else if (nbins>1)
    {
     a_slptmp = sel_aslope_fish*exp(aslope_devs(binindex(i)));
     a50tmp = sel_a50_fish*exp(a50_devs(binindex(i)));
     d_slptmp = sel_dslope_fish*exp(dslope_devs(binindex(i)));
     d50tmp = sel_d50_fish*exp(d50_devs(binindex(i)));

     for (j=1;j<=nselages;j++)
       {
         log_sel_fish(i,j) = -1.*log(1.0+mfexp(-1.0*a_slptmp*(ages(j)-a50tmp)))
                   -1.*log(1.0+mfexp(-1.0*d_slptmp*(ages(j)-d50tmp)));
       }
    }  
  }
 }

  else if (sel_option==3)
  {
    bicubic_spline(scal_yr_nodes,scal_age_nodes,sel_par,log_sel_fish);
   
  }

  else if (sel_option==4)
  {
      bincount = 1;
      for (i=sel_styr; i<=sel_endyr; i++)
          {
            if (i == binstart(bincount))
            {
             log_sel_fish(i)=cubic_spline(sel_par(bincount) );
             if( bincount < nbins ) bincount++;   
            }
            if (i < endyr_r) log_sel_fish(i+1) = log_sel_fish(i);
          }
  }
 }

 // exponentionate  log selectivity
 for (i=styr; i<=endyr_r; i++)
  {     
   if (i < sel_styr) 
      { 
      sel_fish(i)(1,nselages) = mfexp(log_sel_fish(sel_styr));
      }
   else if (i > sel_endyr)  
    {
    sel_fish(i)(1,nselages) = mfexp(log_sel_fish(sel_endyr));
    }
   else   { sel_fish(i)(1,nselages) = mfexp(log_sel_fish(i));}

   if (nselages<nages)  
      sel_fish(i)(nselages+1,nages)=sel_fish(i,nselages);
  }

 if (nselages==nages)
 {
   sel_srv = mfexp(log_sel_srv);   
 }
 else
 {
  for (s=1;s<=nsrv;s++)
    {
      sel_srv(s)(1,nselages) = mfexp(log_sel_srv(s)(1,nselages));
      sel_srv(s)(nselages+1,nages) = sel_srv(s)(nselages);  
    }
 }
 
 
 

 
 
 
  

 

FUNCTION get_mortality //----------------------------------------------------------------------------------
 // Calulate the values of F and Z for each age and year

  
   avg_fmort_dev=mean(fmort_dev);
   for (i=styr_rec; i<=styr_fish-1; i++)  // set F to zero for years with no fishery
   {
   	for (j=1;j<=nages;j++)
 	{
          F(i,j) = 0.0;
 	}
   }


   
   for (i=styr_fish; i<=endyr_r; i++)  // use the selectivity values for the foreign fishery through 1988
   {
    for (j=1;j<=nages;j++)
   	{
        F(i,j)=sel_fish(i,j)*mfexp(log_avg_fmort + fmort_dev(i));
  	}
   }

 
  
   Z=F+M;
   S=mfexp(-1.0*Z);
   spawner_S = mfexp(-spmo_frac*Z); 

FUNCTION get_first_year  //---------------------------------------------------------------------------------------------
          // get the first year of the natage matrix. Note that in both cases below, the 'styr_rec' 
          // refers to the first year for which there is an estimated value of recruitment.  For the 
          // first year age comp option 2, this is simply the styr of the model.  
          // For option 1, this is the styr-(nages-1).  The first year of the spawning biomass vector is styr_rec.

 surv = mfexp(-1.0*M);
 if (fyear_ac_option == 1)  // First year stochastic recruitment devs are coupled with regular recruitment deviations
 {
  natage(styr_rec,1) = mfexp(mean_log_rec)*exp((sigr*sigr)/2);
  for (j=2; j<=nages; j++)
    natage(styr_rec,j) = natage(styr_rec,j-1)*surv;
  natage(styr_rec,nages) /= (1.-surv);

  for (j=styr_rec; j<=styr_rec_dev; j++)   // deviations in the cohorts that make up the plus group are shared
    natage(j,1) = mfexp(mean_log_rec+rec_dev(styr_rec_dev));

  for (j=styr_rec_dev+1; j<styr; j++)
    natage(j,1) = mfexp(mean_log_rec+rec_dev(j));  

  for (j=styr_rec; j<styr; j++)   // get sp_biom and natage for years prior to styr
    {
      natage(j+1)(2,nages) = ++elem_prod(natage(j)(1,nages-1),S(j)(1,nages-1));
      natage(j+1,nages) += natage(j,nages)*S(j,nages);
    }  
 }

 else if (fyear_ac_option == 2) {  // Initial age comps are in equilibrium with historic catch
            			
 natage(styr,1) = mfexp(log_rinit)*exp((sigr*sigr)/2);		// first, write the first age, first year as rzero
  					
 for (j=2; j<=nages;j++)			// next, get the first year ages (2,nages)
    natage(styr,j)=natage(styr,j-1)*mfexp(-(historic_F*sel_fish(styr)(j-1)+M));
 natage(styr,nages) /= (1-mfexp(-(historic_F*sel_fish(styr)(nages)+M)));  	// Plus group for first year

    if (historic_catch > 0.) {  // estimate the historical catch
      ehc = 0;
      for (j=1;j<=nages;j++)
          {
          ehc += natage(styr,j)*wt_fsh(j)*(historic_F*sel_fish(styr,j))*
 	      (1.0-mfexp(-(historic_F*sel_fish(styr,j)+M)))/(historic_F*sel_fish(styr,j)+M);
          }
        } 
 }

  else if (fyear_ac_option == 3) {  // Initial age comps are stochastic, but have a different mean than 
  				   //  the other recruitments. fydev is noise around equilibrium decay.					
  for (j=2; j<=nages_dat;j++)			
    natage(styr,j)=mfexp((log_rinit)   -M*double(j-1) + fydev(j));
  if (nages>nages_dat)
    {                                  // for the 'extra' ages needed to account for aging error, 
     for (j=nages_dat+1; j<=nages;j++) // use a single fydev for all cohorts in the plus group                                        
       natage(styr,j)=mfexp((log_rinit)  -M*double(j-1)+fydev(nages_dat));
    }
  natage(styr,nages) = mfexp((log_rinit) - M*double(nages-1) + fydev(nages_dat))/(1-mfexp(-M));  // plus group
 }

FUNCTION get_numbers_at_age //------------------------------------------------------------------------------------------
   // Get numbers for the age of recruitment for all years, and fill out the natage matrix

 // get the recruits  
 for (i=styr;i<=lastyr_rec;i++)  natage(i,1) = mfexp(mean_log_rec+rec_dev(i));

 for (i=lastyr_rec+1;i<=endyr_r+1;i++) natage(i,1) = mfexp(mean_log_rec +sigr*sigr*0.5);   // ages where we fix the recruits
 
 for (i=styr;i<=endyr_r;i++)    // get natage matrix 
  {
   natage(i+1)(2,nages)=++elem_prod(natage(i)(1,nages-1),S(i)(1,nages-1));
   natage(i+1,nages)+=natage(i,nages)*S(i,nages);  // survival of plus group
 }

FUNCTION dvar_matrix get_projection_numbers_at_age() //------------------------------------------------------------------------------------------
  // Get the numbers at age for the projection model, which uses the mean recruitment for year classes which
  // have not exceeded the criteria for the survey and/or fishery selectivity
  RETURN_ARRAYS_INCREMENT();

   // the last year of recruitment to meet the survey a10 condition 
  double lastyr_rec_a10;
  lastyr_rec_a10 = endyr_r - max(fixedrec_yrs,excludeage - (rec_age - 1));

  dvar_matrix natage_mean_tmp(lastyr_rec_a10+1,endyr_r+1,1,nages);    // numbers at age

  // get the first year
  natage_mean_tmp(lastyr_rec_a10+1) = natage(lastyr_rec_a10+1);
  
  // get the recruits  
  for (i=lastyr_rec_a10+1;i<=endyr_r+1;i++)  natage_mean_tmp(i,1) = mfexp(mean_log_rec +sigr*sigr*0.5);   
  
  for (i=lastyr_rec_a10+1;i<=endyr_r;i++)    // get natage matrix 
    {
   natage_mean_tmp(i+1)(2,nages)=++elem_prod(natage_mean_tmp(i)(1,nages-1),S(i)(1,nages-1));
   natage_mean_tmp(i+1,nages)+=natage_mean_tmp(i,nages)*S(i,nages);  // survival of plus group
    }

 RETURN_ARRAYS_DECREMENT();
 return natage_mean_tmp;


FUNCTION get_expected_values //-----------------------------------------------------------------------------------------
   // get reproductive outputs, total biomass, and survey biomass     
   sp_biom.initialize();
   //sp_biom(styr_rec-rec_age,styr_rec-1) = elem_prod(wt_pop,maturity)*elem_prod(natage(styr_rec)/2.,spawner_S(styr_rec));
   //totbiom(styr_rec-rec_age,styr_rec-1) = natage(styr_rec)*wt_pop;
   //totbiom9899 = 0.0;
   //totbiomwo9899 = 0.0;

    for (j=styr_rec; j<=endyr_r; j++)   // get sp_biom and natagetmp for years prior to styr
    {
      sp_biom(j) = elem_prod(wt_pop,maturity)*elem_prod(natage(j)/2.,spawner_S(j));
      totbiom(j) = natage(j)*wt_pop;
      //for (i=1;i<=nages;i++)
      //{
       //  if (j-ages(i) >= 1998 && j-ages(i) <= 1999){
       //    totbiom9899(j) += natage(j,i)*wt_pop(i);}
      
           
      //} 
    }
   totbiom(endyr_r+1) = natage(endyr_r+1)*wt_pop;
   //totbiom9899(endyr_r+1) = natage(endyr_r+1,13)*wt_pop(13) + natage(endyr_r+1,12)*wt_pop(12);
   //totbiomwo9899 = totbiom - totbiom9899;
 

 // compute the predicted values for the surveys, and the exploitable biomass
 for (i=styr;i<=endyr_r;i++){
  mort(i) = elem_div((1.-mfexp(-Z(i))),Z(i));
  }

 
 for (s=1;s<=nsrv;s++){
   for (i=styr_fish;i<=endyr_r;i++) 
     { 
       pred_srv_bio(s,i)= q_srv(s)*elem_prod(natage(i),mort(i))*elem_prod(sel_srv(s),wt_pop);
       pred_srv_abun(s,i)= q_srv(s)*elem_prod(natage(i),mort(i))*sel_srv(s); 
     }
 }
 


 




FUNCTION get_sr_inputs //----------------------------------------------------------------------------------------------
 // get the inputs for the SR curve  
 // first, define rzero and set bzero, if needed
 if (active(log_rzero)) 
 { 
 rzero = mfexp(log_rzero);
 natagetmp = 0.0;
 natagetmp(1) = rzero;
 for (j=2; j<=nages; j++)
  natagetmp(j) = natagetmp(j-1)*surv;
 natagetmp(nages) /= (1.-surv);
 bzero = elem_prod(wt_pop,maturity)*natagetmp*0.5;
 alpha = 0.8*rzero*steepness/(steepness-0.2);
 beta = 0.2*bzero*((1.-steepness)/(steepness-0.2)); 
 }
  

 est_rec(styr,lastyr_rec-rec_age) = column(natage,1)(styr+rec_age,lastyr_rec).shift(styr);  //assign the estimated rec to year class     	// get the estimated recruitment
 //est_rec(styr_rec_dev,lastyr_rec) = column(natage,1)(styr_rec_dev,lastyr_rec);  //assign the estimated rec to year class       // get the estimated recruitment
 dvar_vector Stmp(styr, lastyr_rec-rec_age);    // temporary S
 //dvar_vector Stmp(styr_rec_dev, lastyr_rec);    // temporary S
 Stmp = sp_biom(styr, lastyr_rec-rec_age);  //assign the ssb
 //Stmp = sp_biom(styr_rec_dev-rec_age, lastyr_rec-rec_age).shift(styr_rec_dev);  //assign the ssb
 est_spb = Stmp;				// save the spawning biomass
 pred_rec = SRecruit(Stmp);			// get the predicted recruits

// get the data for the fitted recruitment curve
 dvariable tmpsp=1.1*max(est_spb);
 for (i=1;i<=20;i++)
 {
   SRec_spawn(i)=tmpsp*double(i)/20.;
   SRec_rec(i)=SRecruit(SRec_spawn(i));
 }

FUNCTION get_catch_at_age // by using Baranov's Catch Equation//-----------------------------------------------
 for (i=styr_fish; i<=endyr_r; i++)
 {
    pred_catch(i) = 0.0;
    for (j=1;j<= nages;j++)
    {
 	     catage(i,j) = natage(i,j)*F(i,j)*(1.0-S(i,j))/Z(i,j);
 	     pred_catch(i)+=catage(i,j)*wt_fsh(j);
    }
   }
 
FUNCTION get_age_comps  // need to apply age error matrices  
  for (i=1;i<=nyrs_fish_unbiased_ac_r;i++)
   {
    eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i))=catage(yrs_fish_unbiased_ac_r(i))/sum(catage(yrs_fish_unbiased_ac_r(i)));
    eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i)) = unbiasedages*eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i));
    eac_fish_unbiased_dat(yrs_fish_unbiased_ac_r(i))(1,nages_dat-1) = eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i))(1,nages_dat-1);
    eac_fish_unbiased_dat(yrs_fish_unbiased_ac_r(i))(nages_dat) = sum(eac_fish_unbiased_mod(yrs_fish_unbiased_ac_r(i))(nages_dat,nages));
   }
 
  for (i=1;i<=nyrs_fish_lc_r;i++)
   {
    elc_fish_tempages(yrs_fish_lc_r(i))=catage(yrs_fish_lc_r(i))/sum(catage(yrs_fish_lc_r(i)));
    elc_fish(yrs_fish_lc_r(i))=translen*elc_fish_tempages(yrs_fish_lc_r(i));
   }
   
   for (i=1;i<=nsrv;i++){
    for (j=1;j<=nyrs_srv_ac_r(i);j++){
      eac_srv_mod(i,yrs_srv_ac_r(i,j))=elem_prod(sel_srv(i),elem_prod(  natage(yrs_srv_ac_r(i,j)),  mort(yrs_srv_ac_r(i,j))  ))/
       (sel_srv(i)*elem_prod(  natage(yrs_srv_ac_r(i,j)),  mort(yrs_srv_ac_r(i,j))));
       eac_srv_mod(i,yrs_srv_ac_r(i,j)) = unbiasedages*eac_srv_mod(i,yrs_srv_ac_r(i,j));
       eac_srv_dat(i,yrs_srv_ac_r(i,j))(1,nages_dat-1) = eac_srv_mod(i,yrs_srv_ac_r(i,j))(1,nages_dat-1);
       eac_srv_dat(i,yrs_srv_ac_r(i,j))(nages_dat) = sum(eac_srv_mod(i,yrs_srv_ac_r(i,j))(nages_dat,nages));
    }
   }

   
   for (i=1;i<=nsrv;i++){
    for (j=1;j<=nyrs_srv_lc_r(i);j++){
       elc_srv_tempages(i,yrs_srv_lc_r(i,j))=elem_prod(sel_srv(i),natage(yrs_srv_lc_r(i,j)))/(sel_srv(i)*natage(yrs_srv_lc_r(i,j)));
       elc_srv(i,yrs_srv_lc_r(i,j))=translen*elc_srv_tempages(i,yrs_srv_lc_r(i,j)); 
    }
   }
   

  
 //if (sd_phase())
 //{
 //  depletion = totbiom(endyr)/totbiom(styr);
 //  endbiom=totbiom(endyr);
 //}

FUNCTION get_binned
  // bin the natage matrix to match the plus group for the data
  // also use a weighted average to get the selectivity for the plus group

  for (i=styr_rec;i<=endyr_r+1;i++)
    {
       natage_bin(i)(1,nages_dat-1) = natage(i)(1,nages_dat-1);
       natage_bin(i)(nages_dat) = sum(natage(i)(nages_dat,nages)); 
    }
  maturity_bin(1,nages_dat-1) = maturity(1,nages_dat-1);
  maturity_bin(nages_dat) = mean(maturity(nages_dat,nages));    

FUNCTION dvariable SRecruit(const dvariable& Stmp)
 RETURN_ARRAYS_INCREMENT();
 dvariable RecTmp;
 switch (sr_type)
 {
    case 1:
      RecTmp = mfexp(mean_log_rec +sigr*sigr*0.5);  // average recruitment
      break;
    case 2:
      RecTmp = alpha*Stmp*(1./(beta + Stmp)); // Beverton-Holt form
      break;
 }
 RETURN_ARRAYS_DECREMENT();
 return RecTmp;

FUNCTION dvar_vector SRecruit(const dvar_vector& Stmp)
 RETURN_ARRAYS_INCREMENT();
 dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
 switch (sr_type)
 {
    case 1:
      RecTmp = mfexp(mean_log_rec +sigr*sigr*0.5);  // average recruitment
      break;
    case 2:
      RecTmp = elem_prod(alpha*Stmp,1./(beta + Stmp)); // Beverton-Holt form
      break;
 }
 RETURN_ARRAYS_DECREMENT();
 return RecTmp;

FUNCTION evaluate_the_objective_function //--------------------------------------------------------------------
 if (fyear_ac_option == 2 && historic_catch > 0.)  // computes the ssq for the historical F 
      histFpen();
 mat_likelihood();
 //cout <<" got mat_like"  << endl;
 rec_likelihood();
 //cout <<" got the rec like" << endl;
 srv_likelihood();
 //cout <<" got the srv like" << endl;
 cat_likelihood();
 //cout <<" got the catch like" << endl;
 Fmort_pen();
 //cout <<" got the Fmort pen" << endl;
 age_likelihood();
 //cout <<" got the age like" << endl;
 prior();
 //cout <<" got the prior" << endl;
 sel_likelihood();
 //cout <<" got sel like" << endl;
 
 
 

FUNCTION histFpen  // fit the historical catches if neccessary
 if (active(historic_F))
  {
  hf_pen = 500.*square(ehc - historic_catch);
  obj_fun += hf_pen;
  }
  

FUNCTION rec_likelihood   // fit the recruitment deviations
 rec_like.initialize();
 //chi = log(est_rec+1e-8) - log(pred_rec+1e-8) + (sigr*sigr*0.5);  // correct for bias;
 //rec_rmse = sqrt(norm2(chi)/size_count(chi)+1e-13);
   rec_rmse = sqrt(norm2(log(est_rec+1e-8) - log(pred_rec+1e-8))/size_count(est_rec)+1e-13);
 //rec_rmse = sqrt(norm2(rec_dev)/size_count(rec_dev)+1e-13);
 //sumrecdev = sum(chi); 
 
 //rec_like(1) = norm2(chi)/(2.*sigr*sigr) + size_count(chi)*log(sigr);
 rec_like(1) = norm2(rec_dev)/(2.*sigr*sigr) + size_count(rec_dev)*log(sigr);
 //rec_like(1) = norm2(chi);
 //rec_like(2) = square((rmse - sigr)/(sigr*sigr/size_count(chi)))/2. + log(sigr*sigr/size_count(chi)) +
 //   square(sumrecdev/(sigr/size_count(chi)))/2. + log(sigr/size_count(chi));
 
 if (fyear_ac_option == 3)
  rec_like(2) = norm2(fydev)/(2.*sigr*sigr) + size_count(fydev)*log(sigr);

 // fitting the SR curve
   if (sr_type==2){
     chi = log(est_rec+1e-8) - log(pred_rec+1e-8) + (sigr*sigr*0.5);  // correct for bias;
     rec_like(3) = norm2(chi)/(2.*sigr*sigr) + size_count(chi)*log(sigr); 
   }

 obj_fun +=lambda(1)*sum(rec_like);
 
 
 
 //  rec_like+=1.0*norm2(rec_dev_future); // the deviations for the future recruitments 
  
FUNCTION age_likelihood  // fit the age comps--------------------------------------------------------------------
 dvariable tmp1;
 dvariable tmp2;
 dvariable tmp3;
 dvariable tmp4;
 dvariable tmp5;
 
 age_like=0.;
 age_like_sac=0.;
 age_like_slc=0.;
 fish_effn=0.;
 sac_effn=0.;
 slc_effn=0.;
 int ii;
 int m;

 k=0;
 for (i=1;i<=nyrs_fish_unbiased_ac_r;i++)
 {
   ii=yrs_fish_unbiased_ac_r(i);
   age_like(1)-=fish_unbiased_ac_samp_r(i)*oac_fish_unbiased_r(i)*log(eac_fish_unbiased_dat(ii)+1.e-13);
   fish_effn(1,i) = eac_fish_unbiased_dat(ii)*(1.-eac_fish_unbiased_dat(ii))/(norm2(eac_fish_unbiased_dat(ii)-oac_fish_unbiased_r(i))); 
   for (j=1;j<=nages_dat;j++)
   {
     k = k+1;
     tmp1 = (oac_fish_unbiased_r(i,j)+0.00001) - (eac_fish_unbiased_dat(ii,j) + 0.00001);
     tmp2 = (eac_fish_unbiased_dat(ii,j)+0.00001)*(1.-(eac_fish_unbiased_dat(ii,j)+0.00001)  );
     //fac_nr(k) = tmp1/sqrt(tmp2/(fish_unbiased_ac_samp(i)*(4./6.)));
     fac_nr(k) = tmp1/sqrt(tmp2/(raw_fish_unbiased_ac_samp(i)));
     fac_resid(k) = tmp1;
     fac_pearson(i,j) = tmp1/sqrt(tmp2/fish_unbiased_ac_samp_r(i));
     fac_pearson_vec(k) = fac_pearson(i,j);    
   }
   tmp3 = ages_dat_mid*oac_fish_unbiased_r(i);   // the mean of the observations 
   tmp4 = ages_dat_mid*eac_fish_unbiased_dat(ii);  // the mean of the predictions
   tmp5 = elem_prod(ages_dat_mid,eac_fish_unbiased_dat(ii))*ages_dat_mid - square(tmp4);  // the v term in Francis method
   fac_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_fish_unbiased_ac_samp(i)); 
 }
 age_like(1)-=offset(1);
 
  
 k = 0;
 for (i=1;i<=nyrs_fish_lc_r;i++)
 {
   ii=yrs_fish_lc_r(i);
   age_like(2)-=fish_lc_samp_r(i)*olc_fish_r(i)*log(elc_fish(ii)+1.e-13);
   fish_effn(2,i) = elc_fish(ii)*(1.-elc_fish(ii))/(norm2(elc_fish(ii)-olc_fish_r(i)));
   for (j=1;j<=nlen;j++)
   {
     k = k+1;
     tmp1 = (olc_fish_r(i,j)+0.00001) - (elc_fish(ii,j) + 0.00001);
     tmp2 = (elc_fish(ii,j)+0.00001)*(1.-(elc_fish(ii,j)+0.00001)  );
     //flc_nr(k) = tmp1/sqrt(tmp2/(fish_lc_samp(i)*(8./6.)));
     flc_nr(k) = tmp1/sqrt(tmp2/(raw_fish_lc_samp(i)));
     flc_resid(k) = tmp1;
     flc_pearson(i,j) = tmp1/sqrt(tmp2/fish_lc_samp_r(i));
     flc_pearson_vec(k) = flc_pearson(i,j); 
   }
   tmp3 = lengths_mid*olc_fish_r(i);   // the mean of the observations 
   tmp4 = lengths_mid*elc_fish(ii);  // the mean of the predictions
   tmp5 = elem_prod(lengths_mid,elc_fish(ii))*lengths_mid - square(tmp4);  // the v term in Francis method
   flc_nr_fran(i) = (tmp3-tmp4)/sqrt(tmp5/raw_fish_lc_samp(i));   
 }
 age_like(2)-=offset(2);
 
 for (i=1;i<=nsrv;i++){
  k = 0;
  for (j=1;j<=nyrs_srv_ac_r(i);j++){
    ii = yrs_srv_ac_r(i,j);
    age_like_sac(i) -= srv_ac_samp_r(i,j)*oac_srv_r(i,j)*log(eac_srv_dat(i,ii)+1.e-13);
    sac_effn(i,j) = eac_srv_dat(i,ii)*(1.-eac_srv_dat(i,ii))/(norm2(eac_srv_dat(i,ii)-oac_srv_r(i,j)));
    for (m=1;m<=nages_dat;m++)
    {
    k = k+1;
     tmp1 = (oac_srv_r(i,j,m)+0.00001) - (eac_srv_dat(i,ii,m) + 0.00001);
     tmp2 = (eac_srv_dat(i,ii,m)+0.00001)*(1.-(eac_srv_dat(i,ii,m)+0.00001)  );
     sac_nr(i,k) = tmp1/sqrt(tmp2/(raw_srv_ac_samp(i,j)));
     sac_resid(i,k) = tmp1;
     sac_pearson(i,j,m) = tmp1/sqrt(tmp2/(srv_ac_samp_r(i,j)));
     sac_pearson_vec(i,k) = sac_pearson(i,j,m);  
    }
    tmp3 = ages_dat_mid*oac_srv_r(i,j);
    tmp4 = ages_dat_mid*eac_srv_dat(i,ii);
    tmp5 = elem_prod(ages_dat_mid,eac_srv_dat(i,ii))*ages_dat_mid - square(tmp4);  // the v term in Francis method
    sac_nr_fran(i,j) = (tmp3-tmp4)/sqrt(tmp5/raw_srv_ac_samp(i,j));

  }
  age_like_sac(i)-=sac_offset(i);
  obj_fun += srv_age_flag(i)*age_like_sac(i);

  k = 0;
  for (j=1;j<=nyrs_srv_lc_r(i);j++){
    ii = yrs_srv_lc_r(i,j);
    age_like_slc(i) -= srv_lc_samp_r(i,j)*olc_srv_r(i,j)*log(elc_srv(i,ii)+1.e-13);
    slc_effn(i,j) = elc_srv(i,ii)*(1.-elc_srv(i,ii))/(norm2(elc_srv(i,ii)-olc_srv_r(i,j)));
    for (m=1;m<=nlen;m++)
   {
     k = k+1;
     tmp1 = (olc_srv_r(i,j,m)+0.00001) - (elc_srv(i,ii,m) + 0.00001);
     tmp2 = (elc_srv(i,ii,m)+0.00001)*(1.-(elc_srv(i,ii,m)+0.00001)  );
     slc_nr(i,k) = tmp1/sqrt(tmp2/raw_srv_lc_samp(i,j)); 
     slc_resid(i,k) = tmp1;
     slc_pearson(i,j,m) = tmp1/sqrt(tmp2/srv_lc_samp_r(i,j));
     slc_pearson_vec(i,k) = slc_pearson(i,j,m);
   }
   tmp3 = lengths_mid*olc_srv_r(i,j);
   tmp4 = lengths_mid*elc_srv(i,ii);
   tmp5 = elem_prod(lengths_mid,elc_srv(i,ii))*lengths_mid - square(tmp4);  // the v term in Francis method
   slc_nr_fran(i,j) = (tmp3-tmp4)/sqrt(tmp5/raw_srv_lc_samp(i,j)); 
  }
  age_like_slc(i)-=slc_offset(i);
  obj_fun += srv_len_flag(i)*age_like_slc(i); 

 }
 

 //obj_fun += (fish_age_flag*age_like(1) + fish_len_flag*age_like(2) + srv_age_flag(1)*age_like(3) + srv_len_flag(1)*age_like(4) + srv_age_flag(2)*age_like(5) + srv_len_flag(2)*age_like(6));
 obj_fun += (fish_age_flag*age_like(1) + fish_len_flag*age_like(2));
 //obj_fun += srv_age_flag(1)*age_like_sac(1) + srv_age_flag(2)*age_like_sac(2); 
 
FUNCTION dvariable update_compweight_ta11(const int& nyrs,const dvector& yrs,const dvector& samp, const dmatrix& obs,const dvar_matrix& est)

 RETURN_ARRAYS_INCREMENT();
 dvariable tmp6;
 dvariable tmp7;
 int ii;
 int z;
 dvar_vector wgt(1,nyrs);
 dvar_vector wgt_inv(1,nyrs);
 dvariable newcompweight_ta11;

 for (z=1;z<=nyrs;z++)
 {
   ii=yrs(z);
   tmp6 = (est(ii)+0.00001)*(1.-(est(ii)+0.00001));
   tmp7 = ((obs(z)+0.00001) - (est(ii) + 0.00001))*
             ((obs(z)+0.00001) - (est(ii) + 0.00001));
   wgt(z) = (tmp6/tmp7)/samp(z);
   wgt_inv(z) =  1.0/wgt(z);
   }
      if (har_flag==1) newcompweight_ta11 = 1.0/mean(wgt_inv);
          else   newcompweight_ta11 = mean(wgt);    
        
 RETURN_ARRAYS_DECREMENT();
 return newcompweight_ta11;
 

 
 


FUNCTION prior  // compute the prior parts of the likelood for q_surv3, M, and the fish sel deviations --------------------------------------------------------------------
 prior_M=0.;
 prior_q=0.;

 if (active(M))
      {prior_M = square(log(M) - log(priormean_M) + square(priorcv_M)/2.0)/(2.*square(priorcv_M));}
 for (i=1;i<=nsrv;i++)
   if (active(q_srv(i)))
    {
      if (priormean_q(i)>0)  prior_q(i) = square(log(q_srv(i)) - log(priormean_q(i)) + square(priorcv_q(i))/2.0)/(2.*square(priorcv_q(i)));
    } 
 
 

 //prior_like(1) =  square(log(q_srv(1)) - log(priormeansurvq) + square(priorcvsurvq)/2.0)/(2.*square(priorcvsurvq));
 //prior_like(2) =  square(log(q_srv_ebs) - log(priormeansurvq) + square(priorcvsurvq)/2.0)/(2.*square(priorcvsurvq));
 
 
 //if(fs_option==1)
 //{
 //  prior_like(3) = square(slopedev)/(2.*square(priorsdfishslopedev));
 //  prior_like(4) = square(a50dev)/(2.*square(priorsdfisha50dev));
 //}

 obj_fun += (prior_M + sum(prior_q));

//FUNCTION compute_spr_rates
 // compute SPR rates and add them to the likelihood for females
// SB0=0.;
// SBF40=0.;
// SBF35=0.;
// SBF30=0.;
 // fill in the number of spawners matrix
// for (i=1;i<=4;i++) 
//  Nspr(i,1) = 1.;
  
// for (j=2;j<nages;j++)
 // {
 //  Nspr(1,j)=Nspr(1,j-1)*exp(-M);
 //  Nspr(2,j)=Nspr(2,j-1)*exp(-(M+F40*sel_fish(endyr,j-1)));
 //  Nspr(3,j)=Nspr(3,j-1)*exp(-(M+F35*sel_fish(endyr,j-1)));
 //  Nspr(4,j)=Nspr(4,j-1)*exp(-(M+F30*sel_fish(endyr,j-1)));
 // }
 // Fill in the plus group
 //Nspr(1,nages)=Nspr(1,nages-1)*exp(-M)/(1.-exp(-M));
 //Nspr(2,nages)=Nspr(2,nages-1)*exp(-(M+F40*sel_fish(endyr,nages-1)))/(1.-exp(-(M+F40*sel_fish(endyr,nages))));
 //Nspr(3,nages)=Nspr(3,nages-1)*exp(-(M+F35*sel_fish(endyr,nages-1)))/(1.-exp(-(M+F35*sel_fish(endyr,nages))));
 //Nspr(4,nfages)=Nspr(4,nages-1)*exp(-(M+F30*sel_fish(endyr,nages-1)))/(1.-exp(-(M+F30*sel_fish(endyr,nages))));
 
 // Kill them off until they spawn
 //for (j=1;j<=nages;j++)
 //{
 //SB0 +=Nspr(1,j)*maturity(j)*0.5*wt_pop(j)*exp(-spmo_frac*M);
 //SBF40 += Nspr(2,j)*maturity(j)*0.5*wt_pop(j)*exp(-spmo_frac*(M+F40*sel_fish(endyr,j)));
 //SBF35 += Nspr(3,j)*maturity(j)*0.5*wt_pop(j)*exp(-spmo_frac*(M+F35*sel_fish(endyr,j)));
 //SBF30 += Nspr(4,j)*maturity(j)*0.5*wt_pop(j)*exp(-spmo_frac*(M+F30*sel_fish(endyr,j)));
 //}

 // Use a penalty to find the values of Fxx%
 //sprpen = 10*square(SBF40-0.4*SB0);
  //sprpen +=10*square(SBF35-0.35*SB0);
 //sprpen +=10*square(SBF30-0.30*SB0);
 //obj_fun += sprpen; 

FUNCTION dvariable get_spr(dvariable Ftemp)    // calculation of equilibrium SPR for equilibrium recruitment
  dvariable phi;
  dvar_vector Ntmp(1,nages_dat);
  Ntmp(1)=1.;
  for (j=2;j<=nages_dat;j++)
     Ntmp(j)=Ntmp(j-1)*exp(-(M+Ftemp*recent_fish_sel(j-1)));  // fills in matrix for ages 2 through nages-1           
   Ntmp(nages_dat)=Ntmp(nages_dat-1)*exp(-(M+Ftemp*recent_fish_sel(nages_dat-1)))/(1.-exp(-(M+Ftemp*recent_fish_sel(nages_dat))));
  
   // Kill them off until they spawn
  for (j=1;j<=nages_dat;j++) 
   phi = 0.5*elem_prod(Ntmp,maturity_bin)*elem_prod(wt_pop_bin,exp(-spmo_frac*(M+Ftemp*recent_fish_sel)));

   return(phi);

FUNCTION dvariable get_spr_rates(double spr_percent)
  dvariable df=1.e-8;
  dvariable F1;
  dvariable dd;
  F1 = 0.0;
  if (M<0.2)  
    F1 = 1.2*M*(1-spr_percent);
  else
    F1 = 0.5*M*(1-spr_percent);
  dvariable F2;
  dvariable F3;
  dvariable yld1;
  dvariable yld2;
  dvariable yld3;
  dvariable dyld;
  dvariable dyldp;
  // Newton Raphson stuff to go here
  //for (int ii=1;ii<=12;ii++)
  dd = 1.0;
  while (dd > 1.e-10)   
 {
    F2     = F1 + df;
    F3     = F1 - df;
    yld1   = -1000*square(log(spr_percent/(get_spr(F1)/SB0)));
    yld2   = -1000*square(log(spr_percent/(get_spr(F2)/SB0)));
    yld3   = -1000*square(log(spr_percent/(get_spr(F3)/SB0)));
    dyld   = (yld2 - yld3)/(2*df);                          // First derivative (to find the root of this)
    //dyldp  = (yld3-(2*yld1)+yld2)/(df*df);  // Newton-Raphson approximation for second derivitive
    //F1    -= dyld/dyldp;
    F1 -= yld1/dyld;
    dd = fabs(log(spr_percent/(get_spr(F1)/SB0)));  
  }
  return(F1); 


FUNCTION srv_likelihood  // fit to indices (lognormal) ---------------------------------------------------------
 srv_bio_like=0;
 srv_bio_rmse=0;
 srv_abun_like=0;
 srv_abun_rmse=0;

  int ii;
 
  for (s=1;s<=nsrv;s++)
  {
   for (i=1;i<=nyrs_srv_bio_r(s);i++)  // biomass likelihood, by survey
     {
      ii=yrs_srv_bio(s,i);
      srv_bio_like(s) += square(log(obs_srv_bio(s,i)+1e-13) - log(pred_srv_bio(s,ii)+1e-13))/(2.*cv_bio_srv(s,i)*cv_bio_srv(s,i));
      survey_bio_nr(s,i) = (log(obs_srv_bio(s,i)+1e-13) - log(pred_srv_bio(s,ii)+1e-13))/cv_bio_srv(s,i);
      srv_bio_rmse(s) += square(log(obs_srv_bio(s,i)+1e-13) - log(pred_srv_bio(s,ii)+1e-13));      
     } 
     if(nyrs_srv_bio_r(s)>0)
     {
      srv_bio_rmse(s) = sqrt(srv_bio_rmse(s)/nyrs_srv_bio_r(s));
      survey_bio_sdnr(s) = std_dev(survey_bio_nr(s));
     }
   obj_fun+= srv_bio_flag(s)*srv_bio_like(s);

   for (i=1;i<=nyrs_srv_abun_r(s);i++)  // abundance likelihood, by survey
     {
      ii=yrs_srv_abun(s,i);
      srv_abun_like(s) += square(log(obs_srv_abun(s,i)+1e-13) - log(pred_srv_abun(s,ii)+1e-13))/(2.*cv_abun_srv(s,i)*cv_abun_srv(s,i));
      survey_abun_nr(s,i) = (log(obs_srv_abun(s,i)+1e-13) - log(pred_srv_abun(s,ii)+1e-13))/cv_abun_srv(s,i);
      srv_abun_rmse(s) += square(log(obs_srv_abun(s,i)+1e-13) - log(pred_srv_abun(s,ii)+1e-13));      
     } 
     if(nyrs_srv_abun_r(s)>0)
     {
      srv_abun_rmse(s) = sqrt(srv_abun_rmse(s)/nyrs_srv_abun_r(s));
      survey_abun_sdnr(s) = std_dev(survey_abun_nr(s));
     }
   obj_fun+= srv_abun_flag(s)*srv_abun_like(s); 
  }
  
  //obj_fun+= srv_bio_flag(1)*srv_bio_like(1) + srv_bio_flag(2)*srv_bio_like(2);
  
 
FUNCTION cat_likelihood  // fit the catches -------------------------------------------------------------
 catch_like=norm2(log(catch_bio(styr_fish,endyr_r)+0.0001) - log(pred_catch+0.00001));
 obj_fun += fish_bio_flag*lambda(2)*catch_like;
 

FUNCTION Fmort_pen  // Phases less than 2, penalize high F's ---------------------------------------------
 fpen = 0.0;
 //if (current_phase()<2)
 //  fpen=10.*norm2(mfexp(fmort_dev+log_avg_fmort)-1.0);
 //else
 //  fpen=.01*norm2(mfexp(fmort_dev+log_avg_fmort)-.2);

 //if (active(fmort_dev))
   fpen+= 0.1*norm2(fmort_dev);
 
 //fpen+=100*square(avg_fmort_dev);
 obj_fun+=fpen;
 

FUNCTION sel_likelihood // penalty for selectivity smoothness and dome-shapeness
 sel_like = 0.0;
 dvariable s = 0.0;
 dvar_matrix trans_log_sel_fish = trans(log_sel_fish);

 // for logistic and double logistic curves with time-varying parameters, penalize the param devs
 if (active(a50_devs))
 {
     sel_like(1) = norm2(a50_devs)/(2.*sigma_a50*sigma_a50);   
 }

 if (active(aslope_devs))
 {
     sel_like(2) = norm2(aslope_devs)/(2.*sigma_aslp*sigma_aslp);
 }
 if (active(d50_devs))
 {
     sel_like(3) = norm2(d50_devs)/(2.*sigma_d50*sigma_d50);   
 }

 if (active(dslope_devs))
 {
     sel_like(4) = norm2(dslope_devs)/(2.*sigma_dslp*sigma_dslp);
 }
 

 // for double logistic  and bicubic spline, penalize the dome-shape

  if(sel_option==2 || sel_option==3 || sel_option==4)
   {
     for (i=sel_styr;i<=sel_endyr;i++)
     {
       for(j=1;j<=nselages-1;j++)
        {
         if(log_sel_fish(i,j)>log_sel_fish(i,j+1))
           {
             sel_like(5) += lambda(7)*square(log_sel_fish(i,j)-log_sel_fish(i,j+1));
           }
        }
      }
    }
 
 // for bicubic spline, penalize the smoothness across ages and years, and interannual differences 
 if(sel_option==3 || sel_option==4)
  {                   // the smoothness penalty (across ages)
   for (i=sel_styr;i<=sel_endyr;i++)
     {
       s = mean(log_sel_fish(i));
       sel_like(9) += 10000*s*s;    
       dvar_vector df2 = first_difference(first_difference(log_sel_fish(i)));
       sel_like(6) += lambda(8)/nselages*df2*df2;    
     }
 
   for (j=1;j<=nselages;j++)
     {
       dvar_vector df1 = first_difference(trans_log_sel_fish(j));
       sel_like(7) += lambda(9)/(sel_endyr-sel_styr+1)*df1*df1;

       dvar_vector df2 = first_difference(df1);
       sel_like(8) += lambda(10)/(sel_endyr-sel_styr+1)*df2*df2; 
     }
  } 

 obj_fun += sum(sel_like);
 //cout <<"sel_like is "<< sel_like << endl; 


 //if (active(log_selcoffs_fish))
 //{
 //  sel_like(1)=norm2(first_difference(first_difference(log_sel_fish)));
     //***sel_like(2)=norm2(first_difference(first_difference(log_sel_srv3)));
    
 //  obj_fun+= lambda(7)*square(avgsel_fish);
 //    obj_fun+= lambda(7)*square(avgsel_srv3);
 //}
 
FUNCTION mat_likelihood  // fit the maturity curve
 
 
 if(mat_input_switch!=1)
 {
 mat_like = 0.0;
 int ii;

 
 for (j=1;j<=nmat_datasets;j++)
 {
   for (i=1;i<=nages_mat(j);i++)
   {
    ii = ages_mat(j,i) - rec_age +1;
    mat_like += -0.01*mat_lambda(ii)*(y_mat(j,i)*log(mat_theta(ii)) + (n_mat(j,i) - y_mat(j,i))*log(1.-mat_theta(ii)+1e-15)); // -ln like
   }
 }

 //for (i=1;i<=nages_S;i++)
 //{
 //ii = ages_S(i) - rec_age +1;
 //mat_like += -0.01*mat_lambda(ii)*(S_y(i)*log(mat_theta(ii)) + (S_n(i) - S_y(i))*log(1.-mat_theta(ii)+1e-15)); // -ln like, Shawdata
 //}
 //}
 

 obj_fun += mat_like;
 }

FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
  {
  RETURN_ARRAYS_INCREMENT();
  int nodes=size_count(spline_coffs);
  dvector ia(1,nodes);
  dvector fa(1,nselages);
  ia.fill_seqadd(0,1./(nodes-1));
  fa.fill_seqadd(0,1./(nselages-1));
  vcubic_spline_function ffa(ia,spline_coffs);
  RETURN_ARRAYS_DECREMENT();
  
  return(ffa(fa));
  }


FUNCTION comp_metrics // get metrics like rmse, sndr, effn, and weights sample sizes for the compostion data  
 if (nyrs_fish_unbiased_ac_r>0)
 {
    fish_rmse(1) = sqrt(mean(elem_prod(fac_resid,fac_resid)));
    fish_sdnr(1) = std_dev(fac_pearson_vec);
    fish_mean_effn(1) =  (sum(fish_effn(1)))/nyrs_fish_unbiased_ac_r;
    fish_mean_samp_wts(1) =  (sum(fish_unbiased_ac_samp_r))/nyrs_fish_unbiased_ac_r;
 }
 
 if (nyrs_fish_lc_r>0)
 {
   fish_rmse(2) = sqrt(mean(elem_prod(flc_resid,flc_resid)));
   fish_sdnr(2) = std_dev(flc_pearson_vec);
   fish_mean_effn(2) =  (sum(fish_effn(2)))/nyrs_fish_lc_r;
   fish_mean_samp_wts(2) =  (sum(fish_lc_samp_r))/nyrs_fish_lc_r;
 }
 
 
 for (i=1;i<=nsrv;i++){
  if (nyrs_srv_ac_r(i)>0){
    sac_rmse(i) = sqrt(mean(elem_prod(sac_resid(i),sac_resid(i))));
    sac_sdnr(i) = std_dev(sac_pearson_vec(i));
    sac_mean_effn(i) =  (sum(sac_effn(i)))/nyrs_srv_ac_r(i);
    sac_mean_samp_wts(i) =  (sum(srv_ac_samp_r(i)))/nyrs_srv_ac_r(i);
    }  
 }  
 
 for (i=1;i<=nsrv;i++){
  if (nyrs_srv_lc_r(i)>0){
    slc_rmse(i) = sqrt(mean(elem_prod(slc_resid(i),slc_resid(i))));
    slc_sdnr(i) = std_dev(slc_pearson_vec(i));
    slc_mean_effn(i) =  (sum(slc_effn(i)))/nyrs_srv_lc_r(i);
    slc_mean_samp_wts(i) =  (sum(srv_lc_samp_r(i)))/nyrs_srv_lc_r(i);
    }  
 }



 //FUNCTION future_projections
 // Start calculation for first year of numbers at age matrix in projection
 //nage_future(styr_fut)(2,nages)=++elem_prod(natage(endyr)(1,nages-1),S(endyr)(1,nages-1));
 //nage_future(styr_fut,nages)+=natage(endyr,nages)*S(endyr,nages);
 // Set future total catch biomass to zero
 //catch_future = 0.;
 //biomass_future = 0.;
 //ssb_future = 0.;
 // Compute recent F levels
 //mean_recent_fs = mfexp(sum(fmort_dev(endyr-4,endyr))/5 + log_avg_fmort);

 
// Loop to cycle different fishing mortality values through the projections
 //for (int l=1;l<=num_proj_Fs;l++)
 // {
 //    switch(l)
 //      {
 //        case 1:
 //         ftmp = F40;
 //         break;
 //        case 2:
 //         ftmp = F40/2.0;
 //         break;
 //        case 3:
 //         ftmp = F30;
 //         break;
 //        case 4:
 //         ftmp = mean_recent_fs;
 //         break; 
 // 	 case 5:
 //         ftmp = 0.0;
 //         break;
 	
 //	}
   
  // Calculation of future F's, Z and survival (S)
  //Z_future = M;
  //for (i=styr_fut;i<=endyr_fut;i++)
  //{
  //  F_future(i) = sel_fish*ftmp;
  //  Z_future(i) +=F_future(i);
  //  S_future(i) = exp(-Z_future(i));
  //}
  // Calculation of future recruitment and spawners
  //Mean average recruitment of the time-series is used for projection
  //  dvariable Rectmp=mfexp(mean_log_rec);
  //  for (i=styr_fut;i<endyr_fut;i++)
  //  {
  //   nage_future(i,1) = Rectmp*mfexp(rec_dev_future(i));
    
    // Now graduate for the next year
   //   nage_future(i+1)(2,nages) = ++elem_prod(nage_future(i)(1,nages-1),S_future(i)(1,nages-1));
   //   nage_future(i+1,nages) += nage_future(i,nages)*S_future(i,nages);
   // }
   // nage_future(endyr_fut,1)= Rectmp*mfexp(rec_dev_future(i));
  
// Calculation of catch at predicted future age composition
    //for (i=styr_fut; i<=endyr_fut;i++)
    //{
      //catage_future(i) = 0.;
      //catage_future(i) += elem_prod(nage_future(i),
        //                  elem_prod(F_future(i),
        //                  elem_div ((1.-S_future(i)),Z_future(i))));
      //if (l!=num_proj_Fs) catch_future(l,i) += catage_future(i)*wt_pop;
      //biomass_future(l,i) += nage_future(i)*wt_pop;
      //ssb_future(l,i) += (nage_future(i)/2)*elem_prod(wt_pop,maturity);
    //}
 //}
     
FUNCTION get_age10
 // get the age at which the survey selectivity exceeds 10% (as potentially modified by the natural mortality rate)
  dvector tmp;
  tmp = value(sel_srv(1)) - 0.10;
  
  for (j=2;j<=nages;j++)
    {
     if(tmp(j-1) < 0.0 &  tmp(j) >= 0 )  firstage = ages(j); 
    }

    if (firstage >=25) firstage = 25;

    firstage += round(0.05/value(M));  // modify by the natural morality
    excludeage = firstage -1;          // exclude ages at and below excludeage
    

    
  
FUNCTION update_compweights   // update if the comp is estmated, otherwise carry over previous comp weight
  // *_ta11 -- McAllister-Ianelli weights (method TA11 in Francis 2011)
  // *_ta12 -- weight by inverse of variance of normalized resids (Method TA1.2 in Francis 2011)
  // *_ta18 -- The weights for the Francis method

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) 
     { 
      compweightsnew_ta11_fsh(1) = update_compweight_ta11(nyrs_fish_unbiased_ac_r,yrs_fish_unbiased_ac_r,raw_fish_unbiased_ac_samp,oac_fish_unbiased_r,eac_fish_unbiased_dat);
      compweightsnew_ta12_fsh(1) = 1./var(fac_nr);
        if (nyrs_fish_unbiased_ac_r >1) compweightsnew_ta18_fsh(1) = 1/(var(fac_nr_fran)*((nyrs_fish_unbiased_ac_r - 1.0)/(nyrs_fish_unbiased_ac_r*1.0)));
        else compweightsnew_ta18_fsh(1) = 1/(var(flc_nr_fran)*((nyrs_fish_lc_r - 1.0)/(nyrs_fish_lc_r*1.0)));   // if only one data point, pair with flc       }
     }
   else 
   {
     compweightsnew_ta11_fsh(1) = compweights_fsh(1);
     compweightsnew_ta12_fsh(1) = compweights_fsh(1);
     compweightsnew_ta18_fsh(1) = compweights_fsh(1);
   }  
  
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) 
     { 
      compweightsnew_ta11_fsh(2) = update_compweight_ta11(nyrs_fish_lc_r,yrs_fish_lc_r,raw_fish_lc_samp, olc_fish_r,elc_fish);
      compweightsnew_ta12_fsh(2) = 1./var(flc_nr);
        if (nyrs_fish_lc_r >1) compweightsnew_ta18_fsh(2) = 1/(var(flc_nr_fran)*((nyrs_fish_lc_r - 1.0)/(nyrs_fish_lc_r*1.0)));
        else compweightsnew_ta18_fsh(2) = 1/(var(fac_nr_fran)*((nyrs_fish_unbiased_ac_r - 1.0)/(nyrs_fish_unbiased_ac_r*1.0)));  // if only one data point, pair with fac  
     }
   else 
   {
     compweightsnew_ta11_fsh(2) = compweights_fsh(2);
     compweightsnew_ta12_fsh(2) = compweights_fsh(2);
     compweightsnew_ta18_fsh(2) = compweights_fsh(2);
   }
   
   for (i=1;i<=nsrv;i++){
   if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0) 
     { 
      compweightsnew_ta11_sac(i) =  update_compweight_ta11(nyrs_srv_ac_r(i),yrs_srv_ac_r(i),raw_srv_ac_samp(i),  oac_srv_r(i)  ,eac_srv_dat(i)  );
      compweightsnew_ta12_sac(i) = 1./var(sac_nr(i));
        if (nyrs_srv_ac_r(i) >1) compweightsnew_ta18_sac(i) = 1/(var(sac_nr_fran(i))*((nyrs_srv_ac_r(i) - 1.0)/(nyrs_srv_ac_r(i)*1.0)));
        else compweightsnew_ta18_sac(i) = 1/(var(slc_nr_fran(i))*((nyrs_srv_lc_r(i) - 1.0)/(nyrs_srv_lc_r(i)*1.0)));  // if only one data point, pair with slc  
     }
   else 
   {
     compweightsnew_ta11_sac(i) = compweights_sac(i);
     compweightsnew_ta12_sac(i) = compweights_sac(i);
     compweightsnew_ta18_sac(i) = compweights_sac(i);
   }
  }

  for (i=1;i<=nsrv;i++){
   if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i) >0) 
     { 
      compweightsnew_ta11_slc(i) =  update_compweight_ta11(nyrs_srv_lc_r(i),yrs_srv_lc_r(i),raw_srv_lc_samp(i),  olc_srv_r(i)  ,elc_srv(i)  );
      compweightsnew_ta12_slc(i) = 1./var(slc_nr(i));
        if (nyrs_srv_lc_r(i) >1) compweightsnew_ta18_slc(i) = 1/(var(slc_nr_fran(i))*((nyrs_srv_lc_r(i) - 1.0)/(nyrs_srv_lc_r(i)*1.0)));
        else compweightsnew_ta18_slc(i) = 1/(var(sac_nr_fran(i))*((nyrs_srv_ac_r(i) - 1.0)/(nyrs_srv_ac_r(i)*1.0)));  // if only one data point, pair with slc   
     }
   else 
   {
     compweightsnew_ta11_slc(i) = compweights_slc(i);
     compweightsnew_ta12_slc(i) = compweights_slc(i);
     compweightsnew_ta18_slc(i) = compweights_slc(i);
   }
  }

FUNCTION get_jitter_err
 //--------------------------------------------------------------------------------------------
 // get the errors for jittering
       random_number_generator rng(jitter_seed);  
       jitterstdnorm.fill_randn(rng);      // get the standard normal errors for jittering
       srand(jitter_seed);
       
    for (i=1;i<= 100;i++) 
      {
          jittererr(i) = jitter_stddev*jitterstdnorm(i);
      }


FUNCTION double apply_jitter(const double& start_num, const int& phase,const int& mult_switch,int& err_count)
 RETURN_ARRAYS_INCREMENT();
 double new_start_num;
 if (phase<0)
  {
     new_start_num = start_num;
  } else {
      if (mult_switch==1)
       {
         new_start_num = start_num*exp(jittererr(err_count) - jitter_stddev*jitter_stddev*0.5);
       } else {
         new_start_num = start_num + jittererr(err_count);
       }  
   err_count += 1;
   }    
 RETURN_ARRAYS_DECREMENT();
 return new_start_num;



REPORT_SECTION //-------------------------------------------------------------------------------------------

 get_jitter_err();
 report << "the jitter errors are " << endl;
 report << jittererr << endl;

  int m;
  rescaled_F = value(mfexp(log_avg_fmort + fmort_dev));
  
  for (i=styr;i<=endyr_r;i++)
      {
        rescaled_sel_fish(i) = value(sel_fish(i)(1,nages_dat));
        rescaled_F(i) = rescaled_F(i)*max(rescaled_sel_fish(i));
        rescaled_sel_fish(i) = rescaled_sel_fish(i)/max(rescaled_sel_fish(i)); 
        exp_biom(i) = elem_prod(natage_bin(i),wt_pop_bin)*rescaled_sel_fish(i);
      }

  for (j=1;j<=nages_dat;j++)
       recent_fish_sel(j) = mean(column(rescaled_sel_fish,j)(endyr_r-4,endyr_r));


  // the last year of recruitment to meet the survey a10 condition 
  double lastyr_rec_a10;
  get_age10();
  lastyr_rec_a10 = endyr_r - max(fixedrec_yrs,excludeage - (rec_age - 1));

   // numbers at age with mean recruitments for the survey age10 year classes
   dvar_matrix natage_mean(lastyr_rec_a10+1,endyr_r+1,1,nages);   
   natage_mean = get_projection_numbers_at_age();

 SB0 = get_spr(0.0);
 F40 = get_spr_rates(0.4);
 F35 = get_spr_rates(0.35);
 F30 = get_spr_rates(0.30);
 SBF40 = get_spr(F40);
 SBF35 = get_spr(F35);
 SBF30 = get_spr(F30);


 
 report << "the firstage is "<< firstage <<endl;
 report << "the excludeage is "<< excludeage <<endl;
 report << "the last year of recruitment is "  <<lastyr_rec << endl;
 report << "the last year of recruitment for a10 is "  << lastyr_rec_a10 << endl;
 report << "the first year where we do the new projection is " << lastyr_rec_a10+1 << endl; 
  report << "Total number of fish: years " <<styr<<" to " << endyr_r+1<< endl;
 	report << rowsum(natage) << endl;
 	report << "Numbers of fish: ages "  <<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
 	for (i=styr;i<= endyr_r+1;i++)
      report << i <<" "<<natage_bin(i) << endl;
	report << "Number in first year: ages "<<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
  report << natage_bin(styr) << endl;
  report << "Number in end year: ages "<<ages_dat(1)<<" to " << ages_dat(nages_dat)<< endl;
 	report << natage_bin(endyr_r) << endl;
  report << "Number in year "<<endyr_r<< ", with mean for age10 yc:" << ages(1) <<" to "<< ages(nages) << endl;
  report << natage_mean(endyr_r) << endl;
  report << "SSB and recruitment for SR fitting" << endl;
  report << "Recruitments (age "<<rec_age<<"): years "<<styr+rec_age<<" to " << lastyr_rec << endl;
  report << est_rec(styr,lastyr_rec-rec_age) << endl;
  report << "Spawner biomass: years "<<styr<<" to " << lastyr_rec-rec_age << endl;
  report << est_spb(styr,lastyr_rec-rec_age) << endl;
  report << "time series of recruitment : years "<<styr<<" to " << lastyr_rec << endl;
  report << column(natage,1)(styr,lastyr_rec) << endl;
  report << "time series spawner biomass: years "<<styr<<" to " << endyr_r << endl;
  report << sp_biom(styr,endyr_r) << endl;
 	report << "SR curve SSB: seq(1,20)" << endl;
 	report << SRec_spawn << endl;
 	report << "SR curve recs: seq(1,20)"  << endl;
 	report << SRec_rec << endl;
  

  for (i=1;i<=nsrv;i++){
      report <<srvname(i)<<": survey selectivity" <<ages_dat<< endl;
      report << sel_srv(i)(1,nselages) << endl;
    }
  
  for (i=1;i<=nsrv;i++){
      if(nyrs_srv_bio_r(i)>0)
      {
        report <<srvname(i)<<" observed biomass: years " <<yrs_srv_bio(i)(1,nyrs_srv_bio_r(i)) << endl;
        report << obs_srv_bio(i)(1,nyrs_srv_bio_r(i)) << endl;
        report <<srvname(i)<<" lower CI: years " <<yrs_srv_bio(i)(1,nyrs_srv_bio_r(i)) << endl;
        report << obs_srv_bio_lower(i)(1,nyrs_srv_bio_r(i)) << endl;
        report <<srvname(i)<<" upper CI: years " <<yrs_srv_bio(i)(1,nyrs_srv_bio_r(i)) << endl;
        report << obs_srv_bio_upper(i)(1,nyrs_srv_bio_r(i)) << endl;
        report <<srvname(i)<<" predicted biomass: years "<<styr_fish<<" to " << endyr_r << endl;
        report << pred_srv_bio(i) << endl;
     }
    }  

   for (i=1;i<=nsrv;i++){
      if(nyrs_srv_abun_r(i)>0)
      {
        report <<srvname(i)<<" observed biomass: years " <<yrs_srv_abun(i)(1,nyrs_srv_abun_r(i)) << endl;
        report << obs_srv_abun(i)(1,nyrs_srv_abun_r(i)) << endl;
        report <<srvname(i)<<" lower CI: years " <<yrs_srv_abun(i)(1,nyrs_srv_abun_r(i)) << endl;
        report << obs_srv_abun_lower(i)(1,nyrs_srv_abun_r(i)) << endl;
        report <<srvname(i)<<" upper CI: years " <<yrs_srv_abun(i)(1,nyrs_srv_abun_r(i)) << endl;
        report << obs_srv_abun_upper(i)(1,nyrs_srv_abun_r(i)) << endl;
        report <<srvname(i)<<" predicted abundance: years "<<styr_fish<<" to " << endyr_r << endl;
        report << pred_srv_abun(i) << endl;
     }
    }

  
  report << "Fishing mortality: years "<<styr<<" to " << endyr_r << endl;
 	report << rescaled_F << endl;
 	
 	report << "Total biomass: years "<<styr_rec<<" to " << endyr_r+1 << endl;
 	report << totbiom << endl;
  report << "Exploitable biomass: years "<<styr<<" to " << endyr_r << endl;
  report << exp_biom << endl;
  report << "Observed catch Biomass: years "<<styr_fish<<" to " << endyr_r << endl;
	report << catch_bio(styr_fish,endyr_r) << endl;
 	report << "Predicted catch Biomass: years "<<styr_fish<<" to " << endyr_r << endl;
	report << pred_catch << endl;
 	report << "Estimated historical catch: 'est_hist_catch'"  << endl;
 	report << ehc << endl;

 	report << "Observed Prop(fishery lengths): year, effn, lengths " <<lengths<< endl;
        for (i=1;i<=nyrs_fish_lc_r; i++)
 	{          
 	  if (fish_lc_samp_r(i)>1)
            {
                  report << yrs_fish_lc_r(i)<<" "<<fish_effn(2,i)<<"  "<<olc_fish_r(i)<< endl;
            }
 	}
 	report << "Predicted Prop(fishery lengths): year, effn, lengths " <<lengths<< endl;
  	for (i=1;i<=nyrs_fish_lc_r; i++)
 	{
            report << yrs_fish_lc_r(i)<<" "<<fish_effn(2,i)<<"  "<<elc_fish(yrs_fish_lc_r(i))<< endl;
 	}

  report << "Observed Prop(fishery unbiased): year, effn, ages " <<ages_dat<< endl;
        for (i=1;i<=nyrs_fish_unbiased_ac_r; i++)
 	{
          if (fish_unbiased_ac_samp_r(i)>1)
            {
 	      report << yrs_fish_unbiased_ac_r(i)<<" "<<fish_effn(1,i)<<" "<<oac_fish_unbiased_r(i)<< endl;
 	    }
 	}
 	report << "Predicted Prop(fishery unbiased): year, effn, ages " <<ages_dat<< endl;
  	for (i=1;i<=nyrs_fish_unbiased_ac_r; i++)
 	{
 	  report << yrs_fish_unbiased_ac_r(i)  <<" "<<fish_effn(1,i)<<" "<<eac_fish_unbiased_dat(yrs_fish_unbiased_ac_r(i))<< endl;
 	}

   
   for (i=1;i<=nsrv;i++){
    if(nyrs_srv_ac_r(i)>0 && srv_age_flag(i)>0){
      report <<srvname(i)<<": observed age comps: year, effn, ages " <<ages_dat<< endl;
      for (j=1;j<=nyrs_srv_ac_r(i);j++)
      {
           report << yrs_srv_ac_r(i,j)<<" "<<sac_effn(i,j)<<" "<<oac_srv_r(i,j)<< endl;
      }
      report <<srvname(i)<<": estimated age comps: year, effn, ages " <<ages_dat<< endl;
      for (j=1;j<=nyrs_srv_ac_r(i);j++)
      {
           report << yrs_srv_ac_r(i,j)<<" "<<sac_effn(i,j)<<" "<<eac_srv_dat(i,yrs_srv_ac_r(i,j))<< endl;
      }
    }
    if(nyrs_srv_lc_r(i)>0 && srv_len_flag(i)){
      report <<srvname(i)<<": observed length comps: year, effn, lengths " <<lengths<< endl;
      for (j=1;j<=nyrs_srv_lc_r(i);j++)
      {
           report << yrs_srv_lc_r(i,j)<<" "<<slc_effn(i,j)<<" "<<olc_srv_r(i,j)<< endl;
      }
      report <<srvname(i)<<": estimated age comps: year, effn, lengths " <<ages_dat<< endl;
      for (j=1;j<=nyrs_srv_lc_r(i);j++)
      {
           report << yrs_srv_lc_r(i,j)<<" "<<slc_effn(i,j)<<" "<<elc_srv(i,yrs_srv_lc_r(i,j))<< endl;
      }
    }
  }	

  
  report << " the mean effective N for the";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i) >0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<fish_mean_effn(1)<<" ";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<fish_mean_effn(2)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report <<sac_mean_effn(i)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report <<slc_mean_effn(i)<<" ";
  report << endl;

  report << " the sample weights for the";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i) >0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<fish_mean_samp_wts(1)<<" ";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<fish_mean_samp_wts(2)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report <<sac_mean_samp_wts(i)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report <<slc_mean_samp_wts(i)<<" ";
  report << endl;

  report << " the sdnr for the";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i) >0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<fish_sdnr(1)<<" ";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<fish_sdnr(2)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report <<sac_sdnr(i)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report <<slc_sdnr(i)<<" ";
  report << endl;

  report << " the root mean square error for the";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i) >0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<fish_rmse(1)<<" ";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<fish_rmse(2)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report <<sac_rmse(i)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report <<slc_rmse(i)<<" ";
  report << endl;

  report << "the rmse for: 'rec', ";
  for (i=1;i<=nsrv;i++)  if (srv_bio_flag(i)>0) report << " '"<<srvname(i)<<" biomass',";
  report << endl;

  report << rec_rmse;
  for (i=1;i<=nsrv;i++)  if (srv_bio_flag(i)>0) report << " "<<srv_bio_rmse(i);
  report << endl; 
  
  report << "The likelhood components: 'histFpen', 'selpen(1,9)', 'rec_likelihood','Fmort_pen','mat_like',";
  if (fish_bio_flag>0) report <<" 'cat_likelihood',";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";

  for (i=1;i<=nsrv;i++)  if (srv_bio_flag(i)>0) report << " '"<<srvname(i)<<" biomass',"; 
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  report << hf_pen<<" "<<sel_like<<" "<<sum(rec_like)<<" "<<fpen<<" "<<mat_like;
  if (fish_bio_flag>0) report <<" "<<catch_like;
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" "<<age_like(1);
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" "<<age_like(2);
  for (i=1;i<=nsrv;i++)  if (srv_bio_flag(i)>0) report << " "<<srv_bio_like(i); 
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0) report << " "<<age_like_sac(i);
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report << " "<<age_like_slc(i);
  report << endl;
  
  report << "the prior components of the like: 'prior_M',";
  for (i=1;i<=nsrv;i++)  if (priormean_q(i)>0) report << " '"<<srvname(i)<<" prior_q',";
  report << endl;

  report << prior_M<<" ";
  for (i=1;i<=nsrv;i++)  if (priormean_q(i)>0) report << " "<<prior_q(i);
  report << endl;

  report << " the Var(NR) weights for the";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i) >0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i) >0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<compweightsnew_ta12_fsh(1)<<" ";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<compweightsnew_ta12_fsh(2)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report <<compweightsnew_ta12_sac(i)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report <<compweightsnew_ta12_slc(i)<<" ";
  report << endl;

  report << " the McAllister- Ianelli weights for the";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<compweightsnew_ta11_fsh(1)<<" ";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<compweightsnew_ta11_fsh(2)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report <<compweightsnew_ta11_sac(i)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report <<compweightsnew_ta11_slc(i)<<" ";
  report << endl;

  report << " the Francis weights for the";
  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<" 'fish_ac',";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<" 'fish_lc',";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report << " '"<< srvname(i)<<" age comps',";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report << " '"<< srvname(i)<<" len comps',";
  report << endl;

  if (fish_age_flag>0 && nyrs_fish_unbiased_ac_r >0) report <<compweightsnew_ta18_fsh(1)<<" ";
  if (fish_len_flag>0 && nyrs_fish_lc_r >0) report <<compweightsnew_ta18_fsh(2)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_age_flag(i)>0 && nyrs_srv_ac_r(i)>0 ) report <<compweightsnew_ta18_sac(i)<<" ";
  for (i=1;i<=nsrv;i++)  if (srv_len_flag(i)>0 && nyrs_srv_lc_r(i)>0) report <<compweightsnew_ta18_slc(i)<<" ";
  report << endl;

  report << "estimates of: 'M', ";
  for (i=1;i<=nsrv;i++)  if (srv_bio_flag(i)>0 || srv_abun_flag(i)>0 ) report << " '"<<srvname(i)<<" q',";
  report << endl;

  report << M;
  for (i=1;i<=nsrv;i++)  if (srv_bio_flag(i)>0 || srv_abun_flag(i)>0 ) report << " "<<q_srv(i);
  report << endl;


  
  report << "the fishery selectivity by year and age "<< endl;  
  report << rescaled_sel_fish << endl;
  report << "the recent selectivity (last five years)" << endl;
  report << recent_fish_sel << endl;
  report << "F spr rates " << endl;
  report << F40 <<" "<< F35 << endl;
  report <<"the objective function is "<<endl;
  report <<obj_fun << endl;
  if (sel_option==3){
    report << " the cubic spline parameter matrix is "<< endl;
    report << sel_par << endl;
  }


  report << "the natage is "  << endl;
  report << natage << endl;
  report << " the wt_pop is "  << endl;
  report << wt_pop << endl;
  report << "the maturity is "  << endl;
  report << maturity << endl;
  

  # include "re24-s-report_exp_r.cxx"   // ADMB code to write the S-compatible report
  
  // write the age comp data, sample size, and pearson resids for a flatfile
  ofstream agecomps("agecomps.dat");
  for (i=1;i<=nsrv;i++)
  {
    if(nyrs_srv_ac_r(i)>0)
    {
      if (i==1) agecomps <<" index year n comp_unit obs pred pearson"  << endl;
      for (j=1;j<=nyrs_srv_ac_r(i);j++)
      {
        for (m=1;m<=nages_dat;m++)
        {    
          agecomps << srvname(i) <<" "<<yrs_srv_ac_r(i,j)<<" "<<srv_ac_samp_r(i,j)<<" "<< 
               " "<<ages_dat(m) <<" "<<oac_srv_r(i,j,m)<<" "<<eac_srv_dat(i,yrs_srv_ac_r(i,j),m)<<" "
                  <<sac_pearson(i,j,m) << endl;
        }   
      }
    }
  }
  
  for (j=1;j<=nyrs_fish_unbiased_ac_r; j++)
  {
    if (fish_unbiased_ac_samp_r(j)>1)
    {
      for (m=1;m<=nages_dat;m++)
        {   
          agecomps << "fishery" <<" "<<yrs_fish_unbiased_ac_r(j)<<" "<<fish_unbiased_ac_samp_r(j)<<" "<< 
          " "<<ages_dat(m) <<" "<<oac_fish_unbiased_r(j,m)   <<" "<<eac_fish_unbiased_dat(yrs_fish_unbiased_ac_r(j),m)<<" "
          <<fac_pearson(j,m)<< endl;        
        }
    }
   }

   // write the length comp data, sample size, and pearson resids for a flatfile
  ofstream lengthcomps("lengthcomps.dat");
  for (i=1;i<=nsrv;i++)
  {
    if(nyrs_srv_lc_r(i)>0)
    {
       if (i==1) lengthcomps <<" index year n comp_unit obs pred pearson"  << endl;
      for (j=1;j<=nyrs_srv_lc_r(i);j++)
      {
        for (m=1;m<=nlen;m++)
        {    
          lengthcomps << srvname(i) <<" "<<yrs_srv_lc_r(i,j)<<" "<<srv_lc_samp_r(i,j)<<" "<< 
               " "<<lengths(m) <<" "<<olc_srv_r(i,j,m)<<" "<<elc_srv(i,yrs_srv_lc_r(i,j),m)<<" "
                  <<slc_pearson(i,j,m) << endl;
        }   
      }
    }
  }
  
  
  for (j=1;j<=nyrs_fish_lc_r; j++)
  {
    if (fish_lc_samp_r(j)>1)
    {
      for (m=1;m<=nlen;m++)
        {   
          lengthcomps << "fishery" <<" "<<yrs_fish_lc_r(j)<<" "<<fish_lc_samp_r(j)<<" "<< 
          " "<<lengths(m) <<" "<<olc_fish_r(j,m)   <<" "<<elc_fish(yrs_fish_lc_r(j),m)<<" "
          <<flc_pearson(j,m)<< endl;        
        }
    }
   }
   

   






// write projection model data file

  ofstream projfile("reproj.dat");
  projfile << "Rougheye"  << endl; 
  projfile << 0 <<" # SSL Species???" << endl;
  projfile << 0 <<" # Constant buffer of Dorn" << endl;
  projfile << 1 <<" # Number of fsheries" << endl;
  projfile << 1 <<" # Number of sexes??" << endl;
  projfile << mean(mfexp(log_avg_fmort + fmort_dev(endyr_r-5,endyr_r-1))) << " #  average 5 yr f " << endl;
  projfile << 1 <<" # author f" << endl;
  projfile << 0.4 <<" # ABC SPR" << endl;
  projfile << 0.35 <<" # MSY SPR" << endl;
  projfile << spawn_mo <<" # Spawnmo" << endl;
  projfile << nages_dat <<" # Number of ages" << endl;
  projfile << 1 <<" # Fratio" << endl;
  Mvec = M;
  projfile << " # Natural mortality " << endl; 
  projfile << Mvec(1,nages_dat) << endl;
  projfile << " # Maturity " << endl; 
  projfile << maturity_bin << endl;       
  projfile << " # Wt Spawn " << endl; 
  projfile << wt_pop_bin << endl;	
  projfile << " # Wt Fish " << endl; 
  projfile << wt_fsh_bin << endl;
  projfile << " # selectivity " << endl; 
  projfile << recent_fish_sel << endl;
  projfile << " # natage " << endl; 
  projfile << natage_bin(endyr_r) << endl; 
  projfile << " # Nrec " << endl; 
  projfile << lastyr_rec - (max(1977,styr)+rec_age) +1   << endl; 
  projfile << " # rec " << endl; 
  projfile << column(natage_bin,1)(max(1977,styr) +rec_age,lastyr_rec) << endl;
  projfile << " # ssb " << endl; 
  projfile << sp_biom(max(1977,styr),lastyr_rec-rec_age) << endl;

  
  ofstream projfile2("reproj_age10.dat");
  projfile2 << "Rougheye"  << endl; 
  projfile2 << 0 <<" # SSL Species???" << endl;
  projfile2 << 0 <<" # Constant buffer of Dorn" << endl;
  projfile2 << 1 <<" # Number of fsheries" << endl;
  projfile2 << 1 <<" # Number of sexes??" << endl;
  projfile2 << mean(mfexp(log_avg_fmort + fmort_dev(endyr_r-5,endyr_r-1))) << " #  average 5 yr f " << endl;
  projfile2 << 1 <<" # author f" << endl;
  projfile2 << 0.4 <<" # ABC SPR" << endl;
  projfile2 << 0.35 <<" # MSY SPR" << endl;
  projfile2 << spawn_mo <<" # Spawnmo" << endl;
  projfile2 << nages_dat <<" # Number of ages" << endl;
  projfile2 << 1 <<" # Fratio" << endl;
  Mvec = M;
  projfile2 << " # Natural mortality " << endl; 
  projfile2 << Mvec(1,nages_dat) << endl;
  projfile2 << " # Maturity " << endl; 
  projfile2 << maturity_bin << endl;       
  projfile2 << " # Wt Spawn " << endl; 
  projfile2 << wt_pop_bin << endl; 
  projfile2 << " # Wt Fish " << endl; 
  projfile2 << wt_fsh_bin << endl;
  projfile2 << " # selectivity " << endl; 
  projfile2 << recent_fish_sel << endl;
  projfile2 << " # natage " << endl; 
  projfile2 << natage_bin(endyr_r) << endl;
  projfile2 << " # Nrec " << endl; 
  projfile2 << lastyr_rec_a10 - (max(1977,styr)+rec_age) +1   << endl; 
  projfile2 << " # rec " << endl; 
  
  projfile2 << column(natage_bin,1)(max(1977,styr) +rec_age,lastyr_rec_a10) << endl;

  projfile2 << " # ssb " << endl; 
  projfile2 << sp_biom(max(1977,styr),lastyr_rec_a10-rec_age) << endl;
  

  
  /*ofstream sampsizesnew("sampsizes_new.ctl");
  sampsizesnew << "##### fishery age comp sample sizes" << endl;
  sampsizesnew << "## reweighted number of hauls with read otoliths" << endl;
  for (i=1;i<=nyrs_fish_unbiased_ac;i++)
    sampsizesnew << round(value((1./var(fac_nr))*raw_fish_unbiased_ac_samp)[i]) << " ";
  sampsizesnew << endl;
  sampsizesnew << "##### fishery lengths sample sizes" << endl; 
  sampsizesnew <<  "## reweighted number of hauls w/ lengths" << endl;
  for (i=1;i<=nyrs_fish_lc;i++)
    sampsizesnew << round(value((1./var(flc_nr))*raw_fish_lc_samp)[i]) << " ";
  sampsizesnew <<   endl;
  sampsizesnew << "##### unbiased AI survey age comp sample sizes" << endl;
  sampsizesnew <<  "## reweighted number of hauls with otoliths, 1991 start" << endl;
  for (i=1;i<=nyrs_surv3_unbiased_ac;i++)
    sampsizesnew << round(value((1./var(sac_nr))*raw_surv3_unbiased_ac_samp)[i]) << " ";
  sampsizesnew <<  endl;
  sampsizesnew <<  "###  rewweighted survey lengths sample sizes (hauls with otoliths), 1991 start" << endl; 
  for (i=1;i<=nyrs_surv3_lc;i++)
    sampsizesnew << round(value((1./var(slc_nr))*raw_surv3_lc_samp)[i]) << " ";
  sampsizesnew << endl;
  */
   
  
   
   
  
  ofstream compweightsnew_ta12_file("compweights_new_ta12.ctl");    // new comp weights based on method TA1.2 in Francis 2011
  compweightsnew_ta12_file << compweightsnew_ta12_fsh <<" ";
  compweightsnew_ta12_file << compweightsnew_ta12_sac <<" ";
  compweightsnew_ta12_file << compweightsnew_ta12_slc << endl; 
  
  ofstream compweightsnew_ta18_file("compweights_new_ta18.ctl");    // new comp weights based on Francis method (TA1.8)
  compweightsnew_ta18_file << compweightsnew_ta18_fsh <<" ";
  compweightsnew_ta18_file << compweightsnew_ta18_sac <<" ";
  compweightsnew_ta18_file << compweightsnew_ta18_slc <<endl;  

  ofstream compweightsnew_ta11_file("compweights_new_ta11.ctl");    // new comp weights McAllister-Ianelli method (TA1.1)
  compweightsnew_ta11_file << compweightsnew_ta11_fsh <<" ";
  compweightsnew_ta11_file << compweightsnew_ta11_sac <<" ";
  compweightsnew_ta11_file << compweightsnew_ta11_slc <<endl;

  ofstream newinitvaluesfile("initvalues_new.ctl");
  err_count = 1; 
  //int err_count;
  newinitvaluesfile << "# read in starting parameter values" << endl;
  newinitvaluesfile << "# M" << endl;
  newinitvaluesfile << apply_jitter(M_start,4,1,err_count) << endl;
  newinitvaluesfile << "# mat_beta1 " << endl;
  newinitvaluesfile << mat_beta1_start << endl;
  newinitvaluesfile << "# mat_beta2 " << endl;
  newinitvaluesfile << mat_beta2_start << endl;
  newinitvaluesfile << "# mean_log_rec" << endl;
  newinitvaluesfile << apply_jitter(mean_log_rec_start,1,0,err_count) << endl;     
  newinitvaluesfile << "# log_rzero" << endl;
  newinitvaluesfile << apply_jitter(log_rzero_start,sr_phase,1,err_count) << endl;
  newinitvaluesfile << "# log_rinit" << endl;
  newinitvaluesfile << apply_jitter(log_rinit_start,1,1,err_count) << endl;  
  newinitvaluesfile << "# log_avg_fmort" << endl;
  newinitvaluesfile << apply_jitter(log_avg_fmort_start,1,1,err_count) << endl;   
  newinitvaluesfile << "# sel_aslope_fish" << endl;
  newinitvaluesfile << apply_jitter(sel_aslope_fish_start,phase_f_sel_ascend,1,err_count) << endl;  
  newinitvaluesfile << "# sel_dslope_fish" << endl;
  newinitvaluesfile << apply_jitter(sel_dslope_fish_start,phase_f_sel_descend,1,err_count) << endl;
  newinitvaluesfile << "# sel_aslope_srv" << endl;
  for (i=1;i<=nsrv;i++)
  {
    newinitvaluesfile << apply_jitter(sel_aslope_srv_start(i),phase_s_sel_aslope(i),1,err_count) << " ";
  }
  newinitvaluesfile << endl;
  newinitvaluesfile << "# sel_a50_fish" << endl;
  newinitvaluesfile << apply_jitter(sel_a50_fish_start,phase_f_sel_ascend,1,err_count) << endl;  
  newinitvaluesfile << "# sel_d50_fish" << endl;
  newinitvaluesfile << apply_jitter(sel_d50_fish_start,phase_f_sel_descend,1,err_count) << endl;
  newinitvaluesfile << "# sel_a50_srv" << endl;
  for (i=1;i<=nsrv;i++)
  {
    newinitvaluesfile << apply_jitter(sel_a50_srv_start(i),phase_s_sel_a50(i),1,err_count) << " ";
  }
  newinitvaluesfile << endl;  
  newinitvaluesfile << "# steepness" << endl;
  newinitvaluesfile << apply_jitter(steepness_start,sr_phase,1,err_count) << endl;
  newinitvaluesfile << "# historic_F" << endl;
  newinitvaluesfile << apply_jitter(historic_F_start,phase_historic_F,1,err_count) << endl;
  newinitvaluesfile << "# q_srv" << endl;
  for (i=1;i<=nsrv;i++)
  {
    newinitvaluesfile << apply_jitter(q_srv_start(i),phase_s_sel_q(i),1,err_count) << " ";
  }
  newinitvaluesfile << endl;  
  newinitvaluesfile << "# sel_srv_mu" << endl;
  for (i=1;i<=nsrv;i++)
  {
    newinitvaluesfile << apply_jitter(sel_srv_mu_start(i),phase_s_sel_mu(i),1,err_count) << " ";
  }
  newinitvaluesfile << endl;  
  newinitvaluesfile << "# sel_srv_dist" << endl;
  for (i=1;i<=nsrv;i++)
  {
    newinitvaluesfile << apply_jitter(sel_srv_dist_start(i),phase_s_sel_dist(i),1,err_count) << " ";
  }
  newinitvaluesfile << endl;
  newinitvaluesfile << "# sel_srv_sig1" << endl;
  for (i=1;i<=nsrv;i++)
  {
    newinitvaluesfile << apply_jitter(sel_srv_sig1_start(i),phase_s_sel_sig1(i),1,err_count) << " ";
  }
  newinitvaluesfile << endl;
  newinitvaluesfile << "# sel_srv_sig2" << endl;
  for (i=1;i<=nsrv;i++)
  {
    newinitvaluesfile << apply_jitter(sel_srv_sig2_start(i),phase_s_sel_sig2(i),1,err_count) << " ";
  }
  newinitvaluesfile << endl;
  newinitvaluesfile << "# jitter seed " << endl;
  newinitvaluesfile << jitter_seed +1 << endl;  
  newinitvaluesfile << "# jitter standard deviation " << endl;
  newinitvaluesfile << jitter_stddev << endl;

  ofstream jitteroutfile("jitter_out.dat");
  jitteroutfile << obj_fun <<" "<<M <<" ";
  for (i=1;i<=nsrv;i++)  jitteroutfile <<q_srv(i)<<" ";   
  jitteroutfile <<" "<<mean_log_rec <<" "<<sel_aslope_srv <<" "<<sel_a50_srv;
  jitteroutfile <<" "<<sel_aslope_fish <<" "<<sel_a50_fish <<" "<<totbiom(endyr_r) << endl;


RUNTIME_SECTION //------------------------------------------------------------------------------------------
    convergence_criteria 1.e-4, 1.e-4, 1.e-4, 1.e-7, 1.e-7, 1.e-7
    maximum_function_evaluations 1000, 1000, 1000, 10000, 10000, 10000

TOP_OF_MAIN_SECTION
  arrmblsize = 1000000000;

GLOBALS_SECTION
 # include "admodel.h"          // Include AD class definitions
 //# include "qfclib.h"           // some function, including rounding
 # include "mhp-s-funcs.cpp"    // Include S-compatible output functions (needs preceding)
 adstring_array srvname;  

 void function_minimizer::mcmc_eval(void)
        {
                // |---------------------------------------------------------------------------|
                // | Added DIC calculation.  Martell, Jan 29, 2013                             |
                // |---------------------------------------------------------------------------|
                // | DIC = pd + dbar
                // | pd  = dbar - dtheta  (Effective number of parameters)
                // | dbar   = expectation of the likelihood function (average f)
                // | dtheta = expectation of the parameter sample (average y)

          gradient_structure::set_NO_DERIVATIVES();
          initial_params::current_phase=initial_params::max_number_phases;
          uistream * pifs_psave = NULL;

        #if defined(USE_LAPLACE)
        #endif

        #if defined(USE_LAPLACE)
            initial_params::set_active_random_effects();
            int nvar1=initial_params::nvarcalc();
        #else
          int nvar1=initial_params::nvarcalc(); // get the number of active parameters
        #endif
          int nvar;

          pifs_psave= new
            uistream((char*)(ad_comm::adprogram_name + adstring(".psv")));
          if (!pifs_psave || !(*pifs_psave))
          {
            cerr << "Error opening file "
                    << (char*)(ad_comm::adprogram_name + adstring(".psv"))
               << endl;
            if (pifs_psave)
            {
              delete pifs_psave;
              pifs_psave=NULL;
              return;
            }
          }
          else
          {
            (*pifs_psave) >> nvar;
            if (nvar!=nvar1)
            {
              cout << "Incorrect value for nvar in file "
                   << "should be " << nvar1 << " but read " << nvar << endl;
              if (pifs_psave)
              {
                delete pifs_psave;
                pifs_psave=NULL;
              }
              return;
            }
          }

          int nsamp = 0;
          double sumll = 0;
          independent_variables y(1,nvar);
          independent_variables sumy(1,nvar);

          do
          {
            if (pifs_psave->eof())
            {
              break;
            }
            else
            {
              (*pifs_psave) >> y;
              sumy = sumy + y;
              if (pifs_psave->eof())
              {
                double dbar = sumll/nsamp;
                int ii=1;
                y = sumy/nsamp;
                initial_params::restore_all_values(y,ii);
                initial_params::xinit(y);
                double dtheta = 2.0 * get_monte_carlo_value(nvar,y);
                double pd     = dbar - dtheta;
                double dic    = pd + dbar;
                double dicValue      = dic;
                double dicNoPar      = pd;

                cout<<"Number of posterior samples    = "<<nsamp    <<endl;
                cout<<"Expectation of log-likelihood  = "<<dbar     <<endl;
                cout<<"Expectation of theta           = "<<dtheta   <<endl;
                cout<<"Number of estimated parameters = "<<nvar1    <<endl;
                    cout<<"Effective number of parameters = "<<dicNoPar <<endl;
                    cout<<"DIC                            = "<<dicValue <<endl;
                break;
              }
              int ii=1;
              initial_params::restore_all_values(y,ii);
              initial_params::xinit(y);
              double ll = 2.0 * get_monte_carlo_value(nvar,y);
              sumll    += ll;
              nsamp++;
              // cout<<sumy(1,3)/nsamp<<" "<<get_monte_carlo_value(nvar,y)<<endl;
            }
          }
          while(1);
          if (pifs_psave)
          {
            delete pifs_psave;
            pifs_psave=NULL;
          }
          return;
        }


   
 
 

       
