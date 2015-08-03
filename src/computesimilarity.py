from __future__ import division


def get_allpairs_medianic(profile1,profile2,termic,ancestors):
	
	medianic=0
	ic=[]
	for term1 in profile1:
		for term2 in profile2:
			ic.append(getlcsic(term1,term2,ancestors,termic))		
	return np.median(ic)

def check_get_allpairs_medianjaccard(profile1,profile2,ancestorswith,ancestorswithout):
	
	mediansimj=0
	simj=[]
	for term1 in profile1:
		for term2 in profile2:
			simjwith,commonwith,unionwith=checkgetsimj(term1,term2,ancestorswith)
			simjwithout,commonwithout,unionwithout=checkgetsimj(term1,term2,ancestorswithout)
			if simjwithout<simjwith:
				print term1,term2
				print "Common With",commonwith
				print "Common Without",commonwithout
				print "Union With",unionwith
				print "Union Without",unionwithout		
	
def checkgetsimj(term1,term2,ancestors):
	common=set()
	union=set()
	if len(set.union(ancestors[term1],ancestors[term2])) >0:
		simj=len(set.intersection(ancestors[term1],ancestors[term2]))/len(set.union(ancestors[term1],ancestors[term2]))
		common=set.intersection(ancestors[term1],ancestors[term2])
		union=set.union(ancestors[term1],ancestors[term2])
	else:
		simj=0
	return simj,common,union





def get_allpairs_medianjaccard(profile1,profile2,ancestors):
	
	mediansimj=0
	simj=[]
	for term1 in profile1:
		for term2 in profile2:
			simj.append(getsimj(term1,term2,ancestors))		
	return np.median(simj)

def getlcsic(term1,term2,ancestors,termic):
	commonancestors=set.intersection(ancestors[term1],ancestors[term2])
	lcslist=[termic[anc] for anc in commonancestors]
	if len(lcslist)>0:
		return max(lcslist)
	else:
		return 0


def get_bestpairs_medianic(profile1,profile2,termic,ancestors):
	finalsim=0
	bestmatchic=[]
	termmatchic=[]
	for term1 in profile1:
		termmatchic=[]
		for term2 in profile2:
			lcsic=getlcsic(term1,term2,ancestors,termic)
			termmatchic.append(lcsic)
		bestmatchic.append(max(termmatchic))
	sim=np.median(bestmatchic)

	bestmatchic=[]
	termmatchic=[]
	for term1 in profile2:
		termmatchic=[]
		for term2 in profile1:
			lcsic=getlcsic(term1,term2,ancestors,termic)
			termmatchic.append(lcsic)
		bestmatchic.append(max(termmatchic))
	revsim=np.median(bestmatchic)

	finalsim=np.mean([sim,revsim])
	return finalsim


def getsimj(term1,term2,ancestors):
	if len(set.union(ancestors[term1],ancestors[term2])) >0:
		simj=len(set.intersection(ancestors[term1],ancestors[term2]))/len(set.union(ancestors[term1],ancestors[term2]))
	else:
		simj=0
	return simj


def get_bestpairs_medianjaccard(profile1,profile2,ancestors):
	finalsim=0
	bestmatchsimj=[]
	termmatchsimj=[]
	for term1 in profile1:
		termmatchsimj=[]
		for term2 in profile2:
			simj=getsimj(term1,term2,ancestors)
			termmatchsimj.append(simj)
		bestmatchsimj.append(max(termmatchsimj))
	sim=np.median(bestmatchsimj)

	bestmatchsimj=[]
	termmatchsimj=[]
	for term1 in profile2:
		termmatchsimj=[]
		for term2 in profile1:
			simj=getsimj(term1,term2,ancestors)
			termmatchsimj.append(simj)
		bestmatchsimj.append(max(termmatchsimj))
	revsim=np.median(bestmatchsimj)

	finalsim=np.mean([sim,revsim])
	return finalsim


def loadancestors(experimentno,flag):
	if flag=="Real":
		parentshomology=open("../data/"+experimentno+"/SubsumersWithHomology.txt")
	else:
		parentshomology=open("../data/random/"+experimentno+"/SubsumersWithHomology.txt")

	ancestorswith=dict()
	ancestorswithout=dict()

	for line in parentshomology:
		data=line.replace("<http://purl.obolibrary.org/obo/","").replace(">","").split("\t")
		child=data[0].strip()
		parent=data[1].strip()
		if child not in ancestorswith:
			ancestorswith[child]=set()
		if "With" in parent:
			ancestorswith[child].add(parent)

	if flag=="Real":
		parentsnohomology=open("../data/"+experimentno+"/SubsumersWithoutHomology.txt")
	else:
		parentsnohomology=open("../data/random/"+experimentno+"/SubsumersWithoutHomology.txt")
	
	for line in parentsnohomology:
		data=line.replace("<http://purl.obolibrary.org/obo/","").replace(">","").split("\t")
		child=data[0].strip()
		parent=data[1].strip()
		if child not in ancestorswithout:
			ancestorswithout[child]=set()
		if "Without" in parent:
			ancestorswithout[child].add(parent)
	parentsnohomology.close()
	parentshomology.close()

	return ancestorswith,ancestorswithout


def getorthologpairs(geneset,datafile,orthologpairs):
	for line in datafile:
		line=line.replace("NCBIGene_","NCBI_gene:")
		gene1=line.split("\t")[0]
		gene2=line.split("\t")[1].strip()
		if gene1 in geneset and gene2 in geneset:
			temp=(gene1,gene2)
			orthologpairs.add(temp)
	return orthologpairs

def prepareCorpusforIC(ancestorswith,ancestorswithout,flag):
	geneset=set()
	annotationgenecount=dict()
	termic=dict()
	
	geneannotations=dict()
	annotationfile=open("../data/AnnotationCorpus.txt")
	if "flag"=="Real":
		outfile=open("../data/GeneAnnotationDataset.txt",'w')
	else:
		outfile=open("../data/random/GeneAnnotationDataset.txt",'w')
	for line in annotationfile:
		if "phenotype" not in line:
			data=line.split("\t")
			gene=data[0].strip()
			geneset.add(gene)
			annotation=data[1].strip().replace(":","_")
			if gene not in geneannotations:
				geneannotations[gene]=set()
			for anc in ancestorswith[annotation]:
				geneannotations[gene].add(anc)
			for anc in ancestorswithout[annotation]:
				geneannotations[gene].add(anc)


	for gene in geneannotations:
		for annotation in geneannotations[gene]:
			if annotation not in annotationgenecount:
				annotationgenecount[annotation]=0
			annotationgenecount[annotation]+=1
		outfile.write(gene+"\t"+(",".join(geneannotations[gene]))+"\n")

	maxic=	round(-math.log(1/float(len(geneset))),2)
	if flag=="Real":
		geneoutfile=open("../data/annotationgenecounts.txt",'w')
	else:
		geneoutfile=open("../data/random/annotationgenecounts.txt",'w')
	for annotation in annotationgenecount:
		geneoutfile.write(annotation+"\t"+str(annotationgenecount[annotation])+"\n")	
		termic[annotation]=round((-math.log(float(annotationgenecount[annotation])/float(len(geneset))))/maxic,2)
			
	geneoutfile.close()
	outfile.close()
	annotationfile.close()	
	return geneset,termic

def loadzfinprofiles():
	zfinprofiles=dict()
	zfin=open("../data/Annotations/Dr-gene-to-phenotype-BF.txt")
	for line in zfin:
		gene=line.split("\t")[0].strip()
		annotation=line.split("\t")[1].strip().replace(":","_")
		if gene not in zfinprofiles:
			zfinprofiles[gene]=set()	
		zfinprofiles[gene].add(annotation)
	zfin.close()
	zfinprofiles=removeempty(zfinprofiles)
	return(zfinprofiles)


def removeempty(profiles):
	remove=set()
	for gene in profiles:
		if len(profiles[gene])==0:
			remove.add(gene)

	for gene in remove:
		profiles.pop(gene,None)
	return profiles

def loadmgiprofiles():
	mgiprofiles=dict()
	mgi=open("../data/Annotations/Mm-gene-to-phenotype-BF.txt")		
	for line in mgi:
		gene=line.split("\t")[0].strip()
		annotation=line.split("\t")[1].strip().replace(":","_")
		if gene not in mgiprofiles:
			mgiprofiles[gene]=set()
		if "MP_0002169" not in annotation and "MP_0003012" not in annotation:
			mgiprofiles[gene].add(annotation)		
	mgi.close()
	mgiprofilesremoved=removeempty(mgiprofiles)
	return(mgiprofilesremoved)

def loadhumanprofiles():
	humanprofiles=dict()
	human=open("../data/Annotations/Hs-gene-to-phenotype.txt")		
	for line in human:
		gene=line.split("\t")[0].strip()
		annotation=line.split("\t")[1].strip().replace(":","_")
		if gene not in humanprofiles:
			humanprofiles[gene]=set()
		humanprofiles[gene].add(annotation)		
	human.close()
	humanprofiles=removeempty(humanprofiles)
	return(humanprofiles)

def loadorthologpairs(species1,species2,geneset,orthologpairs):
	if "Human" in species1+species2 and "MGI" in species1+species2:
		datafile=open("../data/GeneOrthologs/Human_MGI_orthologs.tsv",'r')
		orthologpairs=getorthologpairs(geneset,datafile,orthologpairs)
		
	
	elif "Human" in species1+species2 and "ZFIN" in species1+species2:
		datafile=open("../data/GeneOrthologs/Human_ZFIN_orthologs.tsv",'r')
		orthologpairs=getorthologpairs(geneset,datafile,orthologpairs)	
	else:
		datafile=open("../data/GeneOrthologs/MGI_ZFIN_orthologs.tsv",'r')
		orthologpairs=getorthologpairs(geneset,datafile,orthologpairs)
	
	return orthologpairs



def addtodict(term,dictionary):
	if term not in dictionary:
		dictionary[term]=dict()
	return dictionary

def summary(ap_medianicwith,ap_medianicwithout,ap_mediansimjwith,ap_mediansimjwithout,bp_medianicwith,bp_medianicwithout,bp_mediansimjwith,bp_mediansimjwithout):

	
	ap_medianic_better=0
	ap_medianic_same=0
	ap_medianic_worse=0
	ap_medianic_improvement=[]

	ap_mediansimj_better=0
	ap_mediansimj_same=0
	ap_mediansimj_worse=0
	ap_mediansimj_improvement=[]

	bp_medianic_better=0
	bp_medianic_same=0
	bp_medianic_worse=0
	bp_medianic_improvement=[]
	
	bp_mediansimj_better=0
	bp_mediansimj_same=0
	bp_mediansimj_worse=0
	bp_mediansimj_improvement=[]
	

	total=0
	

	apmedianicwithlist=apmedianicwithoutlist=bpmedianicwithlist=bpmedianicwithoutlist=bpmediansimjwithlist=bpmediansimjwithoutlist=apmediansimjwithlist=apmediansimjwithoutlist=[]
	for gene1 in ap_medianicwith:
		for gene2 in ap_medianicwith[gene1]:
			total+=1
			if ap_medianicwith[gene1][gene2] > ap_medianicwithout[gene1][gene2]:
				ap_medianic_better+=1
				apmedianicwithlist.append(ap_medianicwith[gene1][gene2])
				apmedianicwithoutlist.append(ap_medianicwithout[gene1][gene2])
				ap_medianic_improvement.append(((ap_medianicwith[gene1][gene2]- ap_medianicwithout[gene1][gene2])/ap_medianicwithout[gene1][gene2])*100)

			elif ap_medianicwith[gene1][gene2] < ap_medianicwithout[gene1][gene2]:
				ap_medianic_worse+=1	
			else:
				ap_medianic_same+=1



			if ap_mediansimjwith[gene1][gene2] > ap_mediansimjwithout[gene1][gene2]:
				ap_mediansimj_better+=1
				ap_mediansimj_improvement.append(((ap_mediansimjwith[gene1][gene2]- ap_mediansimjwithout[gene1][gene2])/ap_mediansimjwithout[gene1][gene2])*100)

			elif ap_mediansimjwith[gene1][gene2] < ap_mediansimjwithout[gene1][gene2]:
				ap_mediansimj_worse+=1
			else:
				ap_mediansimj_same+=1



			if bp_medianicwith[gene1][gene2] > bp_medianicwithout[gene1][gene2]:
				bp_medianic_better+=1
				bp_medianic_improvement.append( ((bp_medianicwith[gene1][gene2]- bp_medianicwithout[gene1][gene2])/bp_medianicwithout[gene1][gene2])*100)
			elif bp_medianicwith[gene1][gene2] < bp_medianicwithout[gene1][gene2]:
				bp_medianic_worse+=1
			else:
				bp_medianic_same+=1


			if bp_mediansimjwith[gene1][gene2] > bp_mediansimjwithout[gene1][gene2]:
				bp_mediansimj_better+=1
				bp_mediansimj_improvement.append(  ((bp_mediansimjwith[gene1][gene2] - bp_mediansimjwithout[gene1][gene2])/bp_mediansimjwithout[gene1][gene2] )*100)
			elif bp_mediansimjwith[gene1][gene2] < bp_mediansimjwithout[gene1][gene2]:
				bp_mediansimj_worse+=1
			else:
				bp_mediansimj_same+=1
			
	print "AP IC, Increased, Unchanged, Decreased ",(ap_medianic_better/total)*100,(ap_medianic_same/total)*100,(ap_medianic_worse/total)*100

	print "AP SimJ, Increased, Unchanged, Decreased ",(ap_mediansimj_better/total)*100,(ap_mediansimj_same/total)*100,(ap_mediansimj_worse/total)*100

	print "BP IC, Increased, Unchanged, Decreased ",(bp_medianic_better/total)*100,(bp_medianic_same/total)*100,(bp_medianic_worse/total)*100

	print "BP SimJ, Increased, Unchanged, Decreased ",(bp_mediansimj_better/total)*100,(bp_mediansimj_same/total)*100,(bp_mediansimj_worse/total)*100			


	print "Percent Improvement AP IC", np.median(ap_medianic_improvement)
	print "Percent Improvement BP IC", np.median(bp_medianic_improvement)
	print "Percent Improvement AP SimJ", np.median(ap_mediansimj_improvement)
	print "Percent Improvement BP SimJ", np.median(bp_mediansimj_improvement)

	
def load_all_othologous_genes():
	orthologousgenes=set()
	datafile=open("../data/GeneOrthologs/Human_MGI_orthologs.tsv",'r')
	for line in datafile:
		gene1,gene2=line.strip().split("\t")
		orthologousgenes.add(gene1.strip())
		orthologousgenes.add(gene2.strip())
	datafile.close()	
	datafile=open("../data/GeneOrthologs/Human_ZFIN_orthologs.tsv",'r')
	for line in datafile:
		gene1,gene2=line.strip().split("\t")
		orthologousgenes.add(gene1.strip())
		orthologousgenes.add(gene2.strip())
	datafile.close()	
	
	datafile=open("../data/GeneOrthologs/MGI_ZFIN_orthologs.tsv",'r')
	for line in datafile:
		gene1,gene2=line.strip().split("\t")
		orthologousgenes.add(gene1.strip())
		orthologousgenes.add(gene2.strip())
	datafile.close()	
	return orthologousgenes


def load_non_ortholog_pairs(species1,species2,orthologpairs,orthologousgenes):
	
	nonorthologpairs=[]
	samplesize=len(orthologpairs)
	if species1 =="MGI":
		file1="../data/Annotations/Mm-gene-to-phenotype-BF.txt"
	elif species1 == "ZFIN":
		file1="../data/Annotations/Dr-gene-to-phenotype-BF.txt"
	else:
		file1="../data/Annotations/Hs-gene-to-phenotype.txt"

	if species2 =="MGI":
		file2="../data/Annotations/Mm-gene-to-phenotype-BF.txt"
	elif species2 == "ZFIN":
		file2="../data/Annotations/Dr-gene-to-phenotype-BF.txt"
	else:
		file2="../data/Annotations/Hs-gene-to-phenotype.txt"

	datafile1=open(file1)
	datafile2=open(file2)
	genelist2=set()
	genelist1=set()
	
	for line in datafile1:
		genelist1.add(line.strip().split("\t")[0])


	for line in datafile2:
		genelist2.add(line.strip().split("\t")[0])

	datafile1.close()
	datafile2.close()
		
	genelist1=set.difference(genelist1,orthologousgenes)
	genelist2=set.difference(genelist2,orthologousgenes)
	genelist1=list(genelist1)
	genelist2=list(genelist2)
	genelist1=random.sample(genelist1, 100)
	genelist2=random.sample(genelist2,100)

	if "MGI" in genelist1[0]:
		if "ZFIN" in genelist2[0]:
			nonorthologpairs=list(itertools.product(genelist1, genelist2))
	if "ZFIN" in genelist1[0]:
		if "MGI" in genelist2[0]:
			nonorthologpairs=list(itertools.product(genelist2, genelist1))


	if "NCBI" in genelist1[0]:
		if "ZFIN" in genelist2[0]:
			genepair=allpairs=list(itertools.product(genelist2, genelist1))
	if "ZFIN" in genelist1[0]:
		if "NCBI" in genelist2[0]:
			nonorthologpairs=list(itertools.product(genelist1, genelist2))


	if "MGI" in genelist1[0]:
		if "NCBI" in genelist2[0]:
			nonorthologpairs=list(itertools.product(genelist1, genelist2))
	if "NCBI" in genelist1[0]:
		if "MGI" in genelist2[0]:
			nonorthologpairs=list(itertools.product(genelist2, genelist1))
	

	
	print samplesize
	randomsample=random.sample(nonorthologpairs, samplesize)
	return randomsample
	

	

def randomize_nonortholog_pairs(orthologpairs):
	genelist1=[]
	genelist2=[]
	nonorthologpairs=set()
	for pair in orthologpairs:
		genelist1.append(pair[0])
		genelist2.append(pair[1])
	

	for pair in orthologpairs:
		temp = list(genelist2)
		if pair[1] in temp:
			temp.remove(pair[1])
		randomgene=random.sample(temp,1)
		newpair=(pair[0],randomgene[0])
		genelist2.remove(randomgene[0])
		nonorthologpairs.add(newpair)	
		if newpair in orthologpairs:
			print "Existing pair"
	return nonorthologpairs


def main():
	
	experimentno="Experiment_"+sys.argv[1]
	species1=sys.argv[2]
	species2=sys.argv[3]
	flag=sys.argv[4]
	print experimentno,species1,species2
	

	ancestorswith,ancestorswithout=loadancestors(experimentno,flag)
	geneset,termic=prepareCorpusforIC(ancestorswith,ancestorswithout,flag)
	
	mgiprofiles=loadmgiprofiles()
	zfinprofiles=loadzfinprofiles()
	humanprofiles=loadhumanprofiles()

	# load ortholog pairs
	orthologpairs=set()
	
	orthologousgenes= load_all_othologous_genes()
	orthologpairs=loadorthologpairs(species1,species2,geneset,orthologpairs)
	
	if flag=="Real":
		outfilename="../results/"+species1+"_"+species2+  "_SimilarityScores.tsv"
		
	else:
		outfilename="../results/random/"+species1+"_"+species2+  "_SimilarityScores.tsv"
		
	outfilename=outfilename.replace("ZFIN","Zebrafish").replace("MGI","Mouse")

	analyze(experimentno,species1,species2, orthologpairs,ancestorswith,ancestorswithout,outfilename,mgiprofiles,zfinprofiles,humanprofiles,termic)
	

	nonorthologpairs=load_non_ortholog_pairs(species1,species2,orthologpairs,orthologousgenes)
	nonorthologfile=open("../data/Nonorthologs/"+"Nonothologs_"+species1+"_"+species2+".tsv",'w')
	
	for pair in nonorthologpairs:
		nonorthologfile.write(pair[0]+"\t"+pair[1]+"\n")
	nonorthologfile.close()

	if flag =='Real':
		outfilename="../results/"+"Nonorthologs_"+species1+"_"+species2+  "_SimilarityScores.tsv"
	else:
		outfilename="../results/random/"+"Nonorthologs_"+species1+"_"+species2+  "_SimilarityScores.tsv"
	outfilename=outfilename.replace("ZFIN","Zebrafish").replace("MGI","Mouse")
	
	analyze(experimentno,species1,species2, nonorthologpairs,ancestorswith,ancestorswithout,outfilename,mgiprofiles,zfinprofiles,humanprofiles,termic)

def analyze(experimentno,species1,species2, orthologpairs,ancestorswith,ancestorswithout,outfilename,mgiprofiles,zfinprofiles,humanprofiles,termic):

	
	outfile=open(outfilename,'w')
	
	outfile.write("Gene1\tGene2\tAll Pairs MedianIC With Homology\tAll Pairs MedianIC Without Homology\t All Pairs MedianSimJ With Homology\t All Pairs MedianSimJ Without Homology \t Best Pairs MedianIC With Homology\t Best Pairs MedianIC Without Homology\t Best Pairs MedianSimJ With Homology\t Best Pairs MedianSimJ Without Homology\n")



	genegenelcs=dict()


	ap_medianicwith=dict()
	ap_medianicwithout=dict()

	ap_mediansimjwith=dict()
	ap_mediansimjwithout=dict()

	bp_medianicwith=dict()
	bp_medianicwithout=dict()


	bp_mediansimjwith=dict()
	bp_mediansimjwithout=dict()



	total=0
	lcswithhomology=dict()
	lcswithouthomology=dict()
	lcssetwithhomology=dict()
	lcssetwithouthomology=dict()
	comparingwithhomology=dict()
	comparingwithouthomology=dict()
	termtermsimilarity=dict()

	for orthologpair in orthologpairs:
		gene1=orthologpair[0]
		gene2=orthologpair[1]

		ap_medianicwith=addtodict(gene1,ap_medianicwith)
		ap_medianicwithout=addtodict(gene1,ap_medianicwithout)
		
		ap_mediansimjwith=addtodict(gene1,ap_mediansimjwith)
		ap_mediansimjwithout=addtodict(gene1,ap_mediansimjwithout)

		bp_medianicwith=addtodict(gene1,bp_medianicwith)
		bp_medianicwithout=addtodict(gene1,bp_medianicwithout)
		
		
		bp_mediansimjwith=addtodict(gene1,bp_mediansimjwith)
		bp_mediansimjwithout=addtodict(gene1,bp_mediansimjwithout)
		profile1=profile2=""
		if gene1 in mgiprofiles:
			profile1=mgiprofiles[gene1]
		if gene1 in zfinprofiles:
			profile1=zfinprofiles[gene1]
		if gene1 in humanprofiles:
			profile1=humanprofiles[gene1]		
		

		if gene2 in mgiprofiles:
			profile2=mgiprofiles[gene2]
		if gene2 in zfinprofiles:
			profile2=zfinprofiles[gene2]
		if gene2 in humanprofiles:
			profile2=humanprofiles[gene2]

		if len(profile1)>0 and len(profile2)>0:	
			total+=1
			ap_medianicwith[gene1][gene2]=get_allpairs_medianic(profile1,profile2,termic,ancestorswith)
			ap_medianicwithout[gene1][gene2]=get_allpairs_medianic(profile1,profile2,termic,ancestorswithout)

			ap_mediansimjwith[gene1][gene2]=get_allpairs_medianjaccard(profile1,profile2,ancestorswith)
			ap_mediansimjwithout[gene1][gene2]=get_allpairs_medianjaccard(profile1,profile2,ancestorswithout)

			

			bp_medianicwith[gene1][gene2]=get_bestpairs_medianic(profile1,profile2,termic,ancestorswith)
			bp_medianicwithout[gene1][gene2]=get_bestpairs_medianic(profile1,profile2,termic,ancestorswithout)


			bp_mediansimjwith[gene1][gene2]=get_bestpairs_medianjaccard(profile1,profile2,ancestorswith)
			bp_mediansimjwithout[gene1][gene2]=get_bestpairs_medianjaccard(profile1,profile2,ancestorswithout)


			outfile.write(gene1+"\t"+gene2+"\t"+  str(ap_medianicwith[gene1][gene2])    +"\t"+  str(ap_medianicwithout[gene1][gene2])+"\t"+ str(ap_mediansimjwith[gene1][gene2])    +"\t"+  str(ap_mediansimjwithout[gene1][gene2])+"\t"+   str(bp_medianicwith[gene1][gene2])   +"\t"+  str(bp_medianicwithout[gene1][gene2]) +"\t"+ str(bp_mediansimjwith[gene1][gene2])  +"\t"+str(bp_mediansimjwithout[gene1][gene2])+ "\n")

	summary(ap_medianicwith,ap_medianicwithout,ap_mediansimjwith,ap_mediansimjwithout,bp_medianicwith,bp_medianicwithout,bp_mediansimjwith,bp_mediansimjwithout)

		
	


	outfile.close()
	cmd='(head -n 1 '+ outfilename+ ' && tail -n +2 ' + outfilename+ '| sort -t$"\t" -k3 -nr) > temp.txt'
	os.system(cmd)
	cmd = 'mv temp.txt '+ outfilename
	os.system(cmd)





if __name__ == "__main__":
	from copy import deepcopy
	import numpy as np 
	import scipy.stats
	import itertools
	import random
	from scipy.stats import rankdata
	import math
	import sys
	import os
	main()