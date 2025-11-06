#################################################################################################################
# Mixed infection estimator                                                                                     #
# An algorithm for detection of mixed infections in bacterial whole genome sequence data.                       #
# The algorithm analyses a set of defined variable sites for evidence of mixed infection with up to 2 strains,  #
# if a mixed infection is present, the relative proportions of the dominant and minor haplotypes are estimated. #
# A database of known haplotypes is used to identify the most likely dominant and minor haplotype.              #
#                                                                                                               #
# david.eyre@ndm.ox.ac.uk                                                                                       #
# 04 March 2013                                                                                                 #
#################################################################################################################

### LIBRARIES ###
libpath = "/nfs/research/goldman/anoufa/Rlibs/" #specify path to R libraries if not using default location
.libPaths(c(libpath, .libPaths()))

# Automatically install 'foreach' if it's missing
if (!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach", repos = "https://cloud.r-project.org", lib=libpath)
}
library("foreach", lib.loc=libpath) #library used for bootstrap loop

# if using Mac or Linux support use of multiple cores for bootstrap
if (Sys.info()['sysname']!="Windows") {
	# please ensure the doMC library (and dependancies) are installed (to support multiple cores for bootstrap calculations)
	if (!requireNamespace("doMC", quietly = TRUE)) {
	  install.packages("doMC", repos = "https://cloud.r-project.org", lib=libpath)
	}
	library("doMC", lib.loc=libpath)
}

### OPTIONS ###

#set working directory to location of dataset
setwd("/")

#file with database of known haplotypes - tab delimited format with 2 columns: id, sequence
#sequences supplied as single string, each must be aligned and of the same length with no missing data
sequences.file = "/nfs/research/goldman/anoufa/src/test_eyre_model/Dataset_S1/haplotype_sequence.txt"

#file with population frequency of each haplotype
#tab delimited with 2 columns: id, frequency 
sequences.freq.file = "/nfs/research/goldman/anoufa/src/test_eyre_model/Dataset_S1/haplotype_frequency.txt"

#path to directory with base counts
path = "/nfs/research/goldman/anoufa/src/test_eyre_model/Dataset_S1/"

#file with list of file names containing the base counts with new line for each file
files = as.vector(read.table("/nfs/research/goldman/anoufa/src/test_eyre_model/Dataset_S1/filenames_list.txt"))

#file to record output
outlog = "/nfs/research/goldman/anoufa/src/test_eyre_model/Dataset_S1/sample_output.txt"

#number of cores to use for bootstrapping on Mac / Linux, ignored on Windows
bs.cores = 4
#number of bootstrap iterations to run for parameter estimates
bs.iter = 1000

#base error probability, probability that any one base call represents an error
p.err = 2e-3 

#deviance statistic threshold, above this perform haplotype matching, below assume not mixed
dev.thres = 19.4

### END OPTIONS ###



### CONSTANTS ###
DIPLOID = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
BASES = c("A", "C", "G", "T")


### FUNCTIONS ###
#functions for defining m, mixture proportion and d, between haplotype diversity by maximum likelihood

#define functions to keep values between 0.5 and 1
get.logit2 = function(p) log((2*p-1)/(1-(2*p-1)))
get.inv.logit2 = function(p) (  ( (exp(p) / (1 + exp(p))) /2 ) + 0.5)

#define functions to keep values between 0 and 1
get.logit = function(p) log(p/(1-p))
get.inv.logit = function(p) exp(p) / (1 + exp(p)) 

#several functions to generate the likelihood of m, the mixture proportion
get.p.b = function(A, m, b, e) {
	# returns p(b|A,m,e) where A = {AA, AC, ..., TT} and b is a base in a single read, and e is the read error probability
	p.b.1 = m * e/3 # p where major ≠ b 
	p.b.2 = (1-m) * e/3 # p where minor ≠ b
	if(b == substring(A,1,1)) p.b.1 = m * (1-e) # p where major = b
	if(b == substring(A,2,2)) p.b.2 = (1-m) * (1-e) # p where minor = b
	p.b = p.b.1 + p.b.2
	return(p.b)
	}

get.p.b.diploid = function(m) {
	# returns a matrix of p(b|A,m,e) for all values of b {A, C, G, T} - rows and values of A - columns
	# e is set by constant above
	mat = matrix(NA, nrow=4, ncol=16)
	for (i in 1:4) mat[i,] = sapply(DIPLOID, get.p.b, m=m, b=BASES[i], e=p.err)
	return(mat)
	}

get.p.site = function(site, m, bc, mat) {
	# returns p(B1 to Bn | A, m, e) where B1 to Bn are all reads at a given site
	# - product of p(b|A,m,e) across all reads - taken from base counts
	# where mat is a 16x4 matrix from get.p.b.diploid
	exp(colSums(log(mat)*bc[site,2:5]))
	}

get.d = function(d) {
	# returns a vector of probabilities for the 16 DIPLOID forms
	aa = (1-d)/4
	ab = d/12
	c(aa, ab, ab, ab, ab, aa, ab, ab, ab, ab, aa, ab, ab, ab, ab, aa)
	}

logliki.m = function(parm, bc) {
	# return the log likihood of m,d - product over all sites [ sum over all values of A [ p(B1-Bn | A, m, e)  * p(A) ] ]
	# p(A) from the vector given by get.d
	m = get.inv.logit2(parm[1])
	d = get.inv.logit(parm[2])
	mat = get.p.b.diploid(m)
	i = 1:nrow(bc)
	sum ( sapply(i, function(i,m,bc,mat,d) log(sum(get.p.site(i,m,bc,mat)*get.d(d))), m=m, bc=bc, mat=mat, d=d) )
	}

get.dev = function(ml, bc) {
	#compare the likelihood obtained to the likelihood under an approximation to the null hypothesis
	h0.liki = logliki.m(c(get.logit2(0.9999),get.logit(.00001)),bc)
	D = -2*(h0.liki-ml)
	return(D)
	}


sample.bc = function(bc) {
	#function for sampling base counts
	bc.samp = bc[sample(seq(1:nrow(bc)), nrow(bc), replace=TRUE),]
	return(bc.samp)
}








### SETUP HAPLOTYPE MATCHING DATABASE ###

#read in haplotypes
sequences = read.table(sequences.file, stringsAsFactors=FALSE)

#read in population frequency of each haplotype
sequences.freq = read.table(sequences.freq.file)
rownames(sequences.freq)=sequences.freq[,1]
#find the least common haplotype
min.freq = min(sequences.freq[,2])

#split sequence strings into matrix of sites
sequences.mat = t(as.matrix(sapply(sequences[,2], function(x) unlist(strsplit(x,"")))))
rownames(sequences.mat)=sequences[,1]

sequences.n = nrow(sequences.mat) #number of haplotypes
sequences.len = ncol(sequences.mat) #number of variable sites
sequences.comb = sequences.n * (sequences.n-1) #number of combinations of different minor / major haplotypes

sequences.pairs = rep(NA, sequences.comb) # generate vector to contain all haplotype pairs

# populate vector of haplotype pair names and generate vector of relative frequencies of different haplotype pairs
sequences.pairs.freq = rep(NA, sequences.comb)
ii=1
for (i in rownames(sequences.mat)) {
	for (j in rownames(sequences.mat)) {
		if (i != j) {
			sequences.pairs[ii] = paste(i,'x',j)	
			#obtain the relative frequency of each pair, if the haplotype frequency is known assume is as common as the least common haplotype
			sequences.pairs.freq[ii] = ifelse(is.na(sequences.freq[i,2]), min.freq, sequences.freq[i,2])  * ifelse(is.na(sequences.freq[j,2]), min.freq, sequences.freq[j,2])
			ii = ii+1
		}
	}
}


# normalise pair frequency data to give prior probability of each pair
sequences.pairs.prior = sequences.pairs.freq/sum(sequences.pairs.freq)

#generate matrix of sequences for all haplotype pairs
sequences.pairs.mat = matrix(NA, ncol=sequences.len, nrow=sequences.comb)
rownames(sequences.pairs.mat) = sequences.pairs
for (i in rownames(sequences.mat)) {
	for (j in rownames(sequences.mat)) {
		if (i != j) {
			sequences.pairs.mat[paste(i,'x',j),] = sapply(paste(sequences.mat[i,],sequences.mat[j,],sep=''), function(x) which(DIPLOID==x))
		}
	}
}





### ANALYSE THE DATA ###



#set up output header
cat("id\tML_m\tML_d\tdeviance\tm0.025\tm0.975\td0.025\td0.975\tp_haplotype_pair\n", file=outlog, append=FALSE)
#id, estimate of m, d, deviance statistic, lower bound m, upper bound, lower bound on d, upper bound, haplotype pairs if mixed

#loop over each file containing base count data
for (id in files$V1) 
{
#set up file path
f = paste(path,id,sep="")
#check file is present
if (is.na(file.info(f)$size)) 
	{
	#print error
	print(paste("Missing file for id",id,"skipping...."))
	}
else
	{
		#print info message to screen
		print(paste("Starting ",id,"....",sep=""))
		
		#read in base count data
		base.count = as.matrix(read.table(f,header=T,sep="\t"))
		
		#if any of the sites have a base count of zero skip the file
		if(min(rowSums(base.count[,2:5]))==0)
			{
			print(paste("Base counts of 0 for at least one site in id",id,"skipping...."))
			next()
			}
		
		#ensure base counts are sorted in site order in genome
		base.count[rank(base.count[,1]),] = base.count[c(1:nrow(base.count)),]
		
		
		#obtain initial esimates of d and m
		#estimate of d from total heterozygous sites
		sites = nrow(base.count)
		d.obs = sum(rowSums(base.count[,2:5]==0)<3)
		if (d.obs==0) d.obs=0.00001
		(d.init = d.obs/sites)
		
		#estimate for starting value of m from mean across heterozygous sites
		m.init = 0.999
		if (sum(rowSums(base.count[,2:5]==0)<3)>0)
			{
			bc.het = base.count[rowSums(base.count[,2:5]==0)<3,2:5] #base counts with at least one heterozygote
			if (sum(rowSums(base.count[,2:5]==0)<3)==1)
				{
				(m.init = max(bc.het)/sum(bc.het))
				}
			else
				{
				bc.het.max = sapply(1:nrow(bc.het), function(i) max(bc.het[i,]))
				bc.het.tot = rowSums(bc.het)
				(m.init = mean(bc.het.max/bc.het.tot))
				}
			}
		
		#avoid exact 50/50 mix
		if (m.init == 0.5)
			{
			m.init = 0.50001
			}
		
		init_parm = c(get.logit2(m.init), get.logit(d.init)) #take initial parameters and express in terms of logit functions above to constrain optimisation values
		
		###numerically optimise the values of m and d and obtain ML
		opt.md = optim(init_parm, logliki.m, control=list(fnscale=-1), bc=base.count)
		opt.md
		#ML value for m
		md.ML.m = get.inv.logit2(opt.md$par[1])
		#ML value for d
		md.ML.d = get.inv.logit(opt.md$par[2])
		
		#get deviance from null hypothesis
		md.dev = get.dev(ml=opt.md$value, bc=base.count)
		
		print(paste("got ML m:", md.ML.m, "and d:", md.ML.d, " (deviance: ",md.dev,")"))
		
		
		#run bootstrap, use multiple cores on Mac / Linux
		if (Sys.info()['sysname']!="Windows") {
			# multiple core support for bootstrap
			registerDoMC(cores=bs.cores)
			bs.md = foreach(i = 1:bs.iter, .combine=rbind) %dopar% {
				bc.samp = sample.bc(base.count)
				opt.md.samp = optim(init_parm, logliki.m, control=list(fnscale=-1), bc=bc.samp)
				bs.m = get.inv.logit2(opt.md.samp$par[1])
				bs.d = get.inv.logit(opt.md.samp$par[2])
				c(bs.m, bs.d)
				}
		}
		else {
			bs.md = foreach(i = 1:bs.iter, .combine=rbind) %do% {
				bc.samp = sample.bc(base.count)
				opt.md.samp = optim(init_parm, logliki.m, control=list(fnscale=-1), bc=bc.samp)
				bs.m = get.inv.logit2(opt.md.samp$par[1])
				bs.d = get.inv.logit(opt.md.samp$par[2])
				c(bs.m, bs.d)
				}		
		}
			
			
		bs.m.ci = quantile(bs.md[,1], probs=c(0.025, 0.975))
		bs.d.ci = quantile(bs.md[,2], probs=c(0.025, 0.975))
		print (paste("m (", bs.m.ci["2.5%"], " - ", bs.m.ci["97.5%"], ")"))
		print (paste("d (", bs.d.ci["2.5%"], " - ", bs.d.ci["97.5%"], ")"))
		
		#if reaches deviance threshold then perform haplotype matching
		sequences.pairs.output = NULL
		if(md.dev>dev.thres) 
		{
			# generate matrix of likelihoods based on md.ML.m - p(b|A,md.ML.m,p.err) - 4 bases (A,C,G,T) by 16 DIPLOID
			mat.ML = get.p.b.diploid(md.ML.m)
			
			# generate p(B1-Bn | A, m, e) over all sites - 150 sites by 16 DIPLOID
			bc.liki = t(sapply(1:nrow(base.count),function(i) get.p.site(i, md.ML.m,base.count,mat.ML)))
			
			# walk through bc.liki using DIPLOID type at each site for each haplotype pair to give overall likelihood for each haplotype pair
			sequences.pairs.liki = sapply(1:sequences.comb, function(x) sum(log(sapply(1:nrow(bc.liki), function(i) bc.liki[i, sequences.pairs.mat[x,i]]))))
			
			# obtain probability of each haplotype pair by multiplying by prior
				sequences.pairs.pre_map = sequences.pairs.liki + log(sequences.pairs.prior)
				#normalise
				mx = -max(sequences.pairs.pre_map)
				sequences.pairs.map = exp(sequences.pairs.pre_map+mx)/sum(exp(sequences.pairs.pre_map+mx))
			
			#return most likely pair
			sequences.pairs[which(sequences.pairs.map==max(sequences.pairs.map))]
			
			#report all pairs upto cumualtive p>0.99
			names(sequences.pairs.map)=sequences.pairs
			
			sequences.pairs.map.sort = sort(sequences.pairs.map,decreasing=T)
			sequences.pairs.map.cumsum = cumsum(sequences.pairs.map.sort)
			
			for (i in 1:length(sequences.pairs.map.sort)) {
				if (i==1) 
					{
					print(sequences.pairs.map.sort[i])
					sequences.pairs.output = paste(sequences.pairs.output, paste(names( sequences.pairs.map.sort[i]), round(sequences.pairs.map.sort[i],5), sep=": "), sep=", ")
					}
				else
				{
					if (sequences.pairs.map.cumsum[i-1]<=0.99)
						{
						print(sequences.pairs.map.sort[i])
						sequences.pairs.output = paste(sequences.pairs.output, paste(names( sequences.pairs.map.sort[i]), round(sequences.pairs.map.sort[i],5), sep=": "), sep=", ")
						}
				}
				
			}
		}
		cat(id, md.ML.m, md.ML.d, md.dev, bs.m.ci["2.5%"], bs.m.ci["97.5%"], bs.d.ci["2.5%"], bs.d.ci["97.5%"], paste(substr(sequences.pairs.output,3,nchar(sequences.pairs.output)),"\n",sep=""), sep="\t", file = outlog, append=TRUE)


	}
	
}