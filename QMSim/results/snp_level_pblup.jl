ENV["GKSwstype"]="100" 

# this file is to run JWAS
cd("/group/qtlchenggrp/tianjing/singlestep_nnlmm/qmsim_n114p164/snp_level_pblup")
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
using Plots
Random.seed!(123)
i = parse(Int, ARGS[1]) # 1,2,3
chain_length=10000
@show i

data_folder="jwas_n100p200"
#read pedigree
pedfile    = "/group/qtlchenggrp/tianjing/singlestepdata/QMSim/$data_folder/pedi.csv"
pedigree   = get_pedigree(pedfile,separator=",",header=true)

#read phenotype
genofile = "/group/qtlchenggrp/tianjing/singlestepdata/QMSim/$data_folder/geno.csv"
genotype = CSV.read(genofile,DataFrame)

#PBLUP for each SNP
model_equation = "m$i = intercept + ID";
@show model_equation
model = build_model(model_equation);
set_random(model,"ID",pedigree);
out   = runMCMC(model,genotype,chain_length=chain_length,output_folder="m$i");

mcmc_ebv  = CSV.read("m$i/MCMC_samples_EBV_m$i.txt",DataFrame)
ind_names = names(mcmc_ebv)
mcmc_ebv  = Matrix(mcmc_ebv)
mcmc_mean = zeros(size(mcmc_ebv))
for i in 1:1000
    mcmc_mean[i,:]=mean(mcmc_ebv[1:i,:],dims=1)
end

mcmc_df=DataFrame(mcmc_mean')
insertcols!(mcmc_df,1,:ID => ind_names)
resall=innerjoin(mcmc_df, genotype, on = :ID)

accuracy_datai=ones(1000)*999
for iter in 1:1000
    accuracy_datai[iter]=cor(resall[:,"x$iter"],resall[:,"m$i"])
end

#save results
open("accuracy.m$i.txt", "w") do io
    writedlm(io, accuracy_datai)
end
final_accuracy=round(accuracy_datai[end],digits=3)

#plot results
myfig=plot(collect(chain_length/1000:chain_length/1000:chain_length),accuracy_datai,xlabel="iteration",ylabel="accuracy",title="accu=$final_accuracy")
savefig(myfig,"accuracy.m$i.png")

