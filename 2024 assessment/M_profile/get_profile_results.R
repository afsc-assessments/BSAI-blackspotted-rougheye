
####  read into R, make plots, and get model summaries

# get current directory
master="C:/Users/paul.spencer/Work/stock_assess/rougheye/2024/Nov model/m_24_1_reweighted/M_profile"

results_dir <- file.path(master,"results")
setwd(results_dir)

modelnames <- c("rougheye24_profileM_M01","rougheye24_profileM_M02","rougheye24_profileM_M03","rougheye24_profileM_M04","rougheye24_profileM_M05",
                "rougheye24_profileM_M06","rougheye24_profileM_M07","rougheye24_profileM_M08","rougheye24_profileM_M09","rougheye24_profileM_M10",
                "rougheye24_profileM_M11","rougheye24_profileM_M12","rougheye24_profileM_M13")
                
## read in the *.rdat files
for (i in 1:length(modelnames))
{
         eval(parse(text = paste(modelnames[i]," <- dget('",modelnames[i],".rdat')",sep="")))
}

## define the likelihood components for the profile
likeprof_names <- c("obj_fun","fish.ac","fish.lc","AI_survey.biom","AI_survey.sac","AI_survey.slc")
pretty_names <- c("M","obj_fun","fac","flc","AI_srv_biom","AI_srv_ac","AI_srv_lc")



M_vec <- seq(0.01,0.13,0.01)
         
M_profile_results <- as.data.frame(matrix(data=-9,nrow=length(modelnames),ncol=length(likeprof_names)+1))
M_profile_results[,1] <- M_vec 


for (i in 1:length(modelnames))
{
  for (j in 1:length(likeprof_names)){
     M_profile_results[i,j+1] <- as.numeric(eval(parse(text = paste0(modelnames[i],"$datalikecomp","['",likeprof_names[j],"']"))))
  } 
}
names(M_profile_results) <- pretty_names

# rescale results so min likelihood is at 0

for (j in 2:(length(likeprof_names)+1)){
  M_profile_results[,j] <- M_profile_results[,j] - min(M_profile_results[,j])
}

setwd(master)


         

