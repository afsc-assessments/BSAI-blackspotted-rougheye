

postscript(file="q_profile.ps", horizontal=T, family="Times")

par(oma=c(5,5,5,5), mfrow=c(1,1), cex.axis=1.4, cex.lab=1.4, las=1)


#par(las=1)

#	plot the survey biomass
plot(q_profile_results$q,q_profile_results$fac, type="l", lwd=2, xlab="AI Survey catchability", 
     ylab = "Relative Neg. log likelihood", ylim=c(0,30), col="maroon")
lines(q_profile_results$q,q_profile_results$flc, col="red", lwd=2)
lines(q_profile_results$q,q_profile_results$AI_srv_biom, lwd=2, col="green")
lines(q_profile_results$q,q_profile_results$AI_srv_ac, lwd=2, col="blue")
lines(q_profile_results$q,q_profile_results$AI_srv_lc, lwd=2, col="orange")
lines(q_profile_results$q,q_profile_results$obj_fun, lwd=2, col="black", lty=2)
legend(1.3,30,c("Obj. Function", "Fishery ac", "Fishery lc", "AI srv biomass", "AI srv ac", "AI srv lc"),
       text.col=c("black","maroon","red","green","blue","orange"),cex=1.2, bty="n")

dev.off()




