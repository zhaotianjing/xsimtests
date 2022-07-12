using DelimitedFiles, DataFrames, CSV,StatsBase

snp_file="/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/QMSim/r_n100p200/p1_mrk_001.txt"

t     = readdlm(snp_file,String,skipstart=1)
indID = parse.(Int64,t[:,1])
nInd=length(indID)
snp = t[:,2]
nSNP=length(snp[1])

snp_matrix=ones(nInd,nSNP)*999
for i in 1:nInd
    snpi=snp[i]
    for j in 1:nSNP
        snp_matrix[i,j]=parse(Int64,snpi[j])
    end
end

snp_matrix
sum(snp_matrix.==999)

#Genotypes (0 = a1,a1; 2 = a2,a2; 3 = a1,a2; 4 = a2,a1; 5 = missing; The first allele is paternal and the second allele is maternal) ...

snp_matrix[snp_matrix .== 3] .= 1.0
snp_matrix[snp_matrix .== 4] .= 1.0
sum(snp_matrix.==5)


#QC
maf=vec(mean(snp_matrix,dims=1)/2)
select1=0.01.<maf.<0.99
fixed=vec(var(snp_matrix,dims=1))
select2 = fixed.!=0
selectAll = select1 .& select2
@show sum(selectAll)
geno=snp_matrix[:,selectAll]
@show size(geno)
p=size(geno,2)
snp_name=["m$i" for i in 1:p]

#calc LD
function calcLD(genotypes)
    corM=cor(genotypes)
    res=0
    for i in 1:(size(corM,1)-1)
        res+=corM[i,i+1]^2
    end
    ld=res/size(corM,1)
    return ld
end

ld=calcLD(geno) #0.3

snp_df=DataFrame(geno,snp_name)
insertcols!(snp_df,1,:ID => indID)
snp_df
CSV.write("/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/QMSim/jwas_n100p200/geno.csv",snp_df)




#y (h2=0.7, 10% QTL)
n=size(geno,1)
p=size(geno,2)

nQTL=Int(floor(p*0.05))
QTL_effect=randn(MersenneTwister(1),nQTL)
QTL_pos=sample(MersenneTwister(2),1:p,nQTL,replace=false, ordered=true)
bv=geno[:,QTL_pos] * QTL_effect
bv=bv/std(bv)*sqrt(0.7)
@show var(bv)

y=bv+randn(MersenneTwister(3),n)*sqrt(0.3)
@show var(y) #~1

@show var(bv)/var(y)

y_df=DataFrame(ID=indID,y=y,bv=bv)
CSV.write("/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/QMSim/jwas_n100p200/y.csv",y_df)



pedfile="/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/QMSim/r_n100p200/pedigree_001.txt"
pedigree= CSV.read(pedfile,DataFrame)
pedigree=pedigree[:,1:3]
CSV.write("/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/QMSim/jwas_n100p200/pedi.csv",pedigree)
