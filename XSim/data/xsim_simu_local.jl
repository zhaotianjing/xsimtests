cd("/Users/tianjing/Library/CloudStorage/Box-Box/james_simulation_QMSim_XSim/XSim")

using XSim

println("---------Step0: create founders from cattle genome")
build_genome(species="cattle")
build_phenome(200)


# Step1. create founders, 100 males and 100 females
n_sires = 100
sires   = Founders(n_sires)
n_dams  = 100
dams    = Founders(n_dams)
println("---------END Step1")


println("---------Step2: using base population to generate 5 generations")
######using base population to generate 5 generations (each generation has 200 individuals)
args_mate2 = Dict(:nA         => 10,  #nA individuals will be sampled from Dam
                  :nB_per_A   => 2,   #nB_per_A individuals sampled from Sire will mate with each individual from Dam
                  :replace_A  => false,
                  :replace_B  => false,
                  :n_per_mate => 1,    #n_per_mate progenies will be reproduced from each pair of mating parent
                  :ratio_malefemale=> 1)
#generation1
dams1,sires1 = mate(dams, sires; args_mate2...)

#generation2
dams2,sires2 = mate(dams1, sires1; args_mate2...)

#generation2
dams3,sires3 = mate(dams2, sires2; args_mate2...)

#generation2
dams4,sires4 = mate(dams3, sires3; args_mate2...)

#generation2
dams5,sires5 = mate(dams4, sires4; args_mate2...)

#generation1-5
pedALL5=XSim.get_pedigree(dams1+sires1+dams2+sires2+dams3+sires3+dams4+sires4+dams5+sires5)
genoALL5=XSim.get_genotypes(dams1+sires1+dams2+sires2+dams3+sires3+dams4+sires4+dams5+sires5)
@show size(genoALL5)
#QC
using Statistics
maf=vec(mean(genoALL5,dims=1)/2)
select1=0.01.<maf.<0.99;
fixed=vec(var(genoALL5,dims=1));
select2 = fixed.!=0;
selectAll = select1 .& select2;
@show sum(selectAll)
genoALL5_2=genoALL5[:,selectAll];
@show size(genoALL5_2)

#save data
using DelimitedFiles
writedlm("pedi.txt",pedALL5)

ind_ID=pedALL5[:,1]
geno_df=DataFrame(genoALL5_2[:,1:200],["m$i" for i in 1:200])
insertcols!(geno_df,1,:ID => ind_ID)
CSV.write("geno.csv",geno_df) #save the first 200 SNPs
