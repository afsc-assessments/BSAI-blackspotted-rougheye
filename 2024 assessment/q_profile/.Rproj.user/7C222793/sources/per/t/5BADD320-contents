
####  read into R, make plots, and get model summaries

# get current directory
master="C:/Users/paul.spencer/Work/stock_assess/rougheye/2024/Nov model/m_24_1_reweighted/q_profile"

results_dir <- file.path(master,"results")
setwd(results_dir)

modelnames <- c("rougheye24_q05","rougheye24_q06","rougheye24_q07","rougheye24_q08","rougheye24_q09","rougheye24_q10",
                "rougheye24_q11","rougheye24_q12","rougheye24_q13","rougheye24_q14","rougheye24_q15","rougheye24_q16"
                ,"rougheye24_q17","rougheye24_q18")
                
## read in the *.rdat files
for (i in 1:length(modelnames))
{
         eval(parse(text = paste(modelnames[i]," <- dget('",modelnames[i],".rdat')",sep="")))
}

## define the likelihood components for the profile
likeprof_names <- c("obj_fun","fish.ac","fish.lc","AI_survey.biom","AI_survey.sac","AI_survey.slc")
pretty_names <- c("q","obj_fun","fac","flc","AI_srv_biom","AI_srv_ac","AI_srv_lc")



q_vec <- seq(0.5,1.8,0.1)
         
q_profile_results <- as.data.frame(matrix(data=-9,nrow=length(modelnames),ncol=length(likeprof_names)+1))
q_profile_results[,1] <- q_vec 


for (i in 1:length(modelnames))
{
  for (j in 1:length(likeprof_names)){
     q_profile_results[i,j+1] <- as.numeric(eval(parse(text = paste0(modelnames[i],"$datalikecomp","['",likeprof_names[j],"']"))))
  } 
}
names(q_profile_results) <- pretty_names

# rescale results so min likelihood is at 0

for (j in 2:(length(likeprof_names)+1)){
  q_profile_results[,j] <- q_profile_results[,j] - min(q_profile_results[,j])
}

setwd(master)


         

