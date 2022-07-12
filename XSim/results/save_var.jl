cd("/group/qtlchenggrp/tianjing/xsim_test")
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra

p=200
res=ones(p,4)*999

for i in 1:p
    res[i,1]=i
    if ispath("m$i/residual_variance.txt") && ispath("m$i/genetic_variance.txt") && ispath("m$i/heritability.txt")
        var_e=CSV.read("m$i/residual_variance.txt",DataFrame)[1,2]
        var_g=CSV.read("m$i/genetic_variance.txt", DataFrame)[1,2]
        h2   =CSV.read("m$i/heritability.txt",     DataFrame)[1,2]
        res[i,2]=var_e
        res[i,3]=var_g
        res[i,4]=h2
    end
end

writedlm("res.txt",res)
