res=read.table("/Users/tianjing/Library/CloudStorage/Box-Box/james_simulation_QMSim_XSim/QMSim/results/res.txt")
names(res)=c("markerID","var_e","var_g","h2")
hist(res$h2)



