def getsimilarity(profile1,profile2,termic,ancestors):
	
	finalmaxic=0
	finallcs="None"
	lcsset=set()
	comparing=""
	for term1 in profile1:
		for term2 in profile2:
			commonancestors=set.intersection(ancestors[term1],ancestors[term2])
			for anc in commonancestors:	
				if termic[anc]>finalmaxic:
					finalmaxic=termic[anc]
					finallcs=anc
					lcsset = deepcopy(commonancestors)
					comparing=term1+"--"+term2
	return finalmaxic,finallcs,lcsset,comparing


def loadancestors(experimentno):
	parentshomology=open("../data/"+experimentno+"/SubsumersWithHomology.txt")
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

	parentsnohomology=open("../data/"+experimentno+"/SubsumersWithoutHomology.txt")
	
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
			temp=[gene1,gene2]
			orthologpairs.append(temp)
	return orthologpairs

def prepareCorpusforIC(ancestorswith,ancestorswithout):
	geneset=set()
	annotationgenecount=dict()
	termic=dict()
	
	geneannotations=dict()
	annotationfile=open("../data/AnnotationCorpus.txt")
	outfile=open("../data/GeneAnnotationDataset.txt",'w')
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


	geneoutfile=open("../data/annotationgenecounts.txt",'w')
	for annotation in annotationgenecount:
		geneoutfile.write(annotation+"\t"+str(annotationgenecount[annotation])+"\n")	
		termic[annotation]=round(-math.log(float(annotationgenecount[annotation])/float(len(geneset))),2)
			
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
	return(zfinprofiles)

def loadmgiprofiles():
	mgiprofiles=dict()
	mgi=open("../data/Annotations/Mm-gene-to-phenotype-BF.txt")		
	for line in mgi:
		gene=line.split("\t")[0].strip()
		annotation=line.split("\t")[1].strip().replace(":","_")
		if gene not in mgiprofiles:
			mgiprofiles[gene]=set()
		mgiprofiles[gene].add(annotation)		
	mgi.close()
	return(mgiprofiles)

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

def main():
	
	experimentno="Experiment_"+sys.argv[1]
	species1=sys.argv[2]
	species2=sys.argv[3]
	outfile=open("../results/"+experimentno+"/" +species1+"_"+species2+  "_SimilarityScores.tsv",'w')
	outfile.write("Gene1\tGene2\tBetter Similarity With Homology\tSimilarity With Homology\tLCS With Homology\tAnnotation Pair leading to best match\tCommon subsumer set\tSimilarity Without Homology\tLCS Without Homology\tAnnotation Pair leading to best match\tCommon subsumer set\t\n")
	
	genegenelcs=dict()

	ancestorswith,ancestorswithout=loadancestors(experimentno)
	geneset,termic=prepareCorpusforIC(ancestorswith,ancestorswithout)
	
	# load ortholog pairs
	orthologpairs=[]
	orthologpairs=loadorthologpairs(species1,species2,geneset,orthologpairs)
	
	#build profiles of genes
	mgiprofiles=loadmgiprofiles()
	zfinprofiles=loadzfinprofiles()
	humanprofiles=loadhumanprofiles()

	similaritywithhomology=dict()
	similaritywithouthomology=dict()
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
		print gene1,gene2
		if gene1 not in similaritywithhomology:
			similaritywithhomology[gene1]=dict()
		if gene1 not in similaritywithouthomology:
			similaritywithouthomology[gene1]=dict()
	
		if gene1 not in lcswithhomology:
			lcswithhomology[gene1]=dict()
		if gene1 not in lcswithouthomology:
			lcswithouthomology[gene1]=dict()

		if gene1 not in lcssetwithouthomology:
			lcssetwithouthomology[gene1]=dict()
		if gene1 not in comparingwithouthomology:
			comparingwithouthomology[gene1]=dict()
		
		if gene1 not in lcssetwithhomology:
			lcssetwithhomology[gene1]=dict()
		if gene1 not in comparingwithhomology:
			comparingwithhomology[gene1]=dict()

		if gene1 in mgiprofiles:
			profile1=mgiprofiles[gene1]
		elif gene1 in zfinprofiles:
			profile1=zfinprofiles[gene1]
		else:
			profile1=humanprofiles[gene1]		
		
		if gene2 in mgiprofiles:
			profile2=mgiprofiles[gene2]
		elif gene2 in zfinprofiles:
			profile2=zfinprofiles[gene2]
		else:
			profile2=humanprofiles[gene2]

		withhomologyic,withhomologylcs,withhomologylcsset,withhomologycomparing=getsimilarity(profile1,profile2,termic,ancestorswith)

		lcswithhomology[gene1][gene2]=withhomologylcs
		similaritywithhomology[gene1][gene2]= withhomologyic
		lcssetwithhomology[gene1][gene2]=withhomologylcsset
		comparingwithhomology[gene1][gene2]=withhomologycomparing



		withouthomologyic,withouthomologylcs,withhomologylcsset,withhomologycomparing=getsimilarity(profile1,profile2,termic,ancestorswithout)
		lcswithouthomology[gene1][gene2]=withouthomologylcs
		similaritywithouthomology[gene1][gene2]= withouthomologyic	
		lcssetwithouthomology[gene1][gene2]=withhomologylcsset
		comparingwithouthomology[gene1][gene2]=withhomologycomparing	


	same=0
	with_greater=0
	without_greater=0
	for gene1 in similaritywithhomology:
		for gene2 in similaritywithhomology[gene1]:
			lcsstringwith=','.join(lcssetwithhomology[gene1][gene2])
			lcsstringwithout=','.join(lcssetwithouthomology[gene1][gene2])
			better=0
			if similaritywithhomology[gene1][gene2] > similaritywithouthomology[gene1][gene2]:
				better=1
				with_greater+=1
			if similaritywithhomology[gene1][gene2] < similaritywithouthomology[gene1][gene2]:
				better=-1
				without_greater+=1 
			if better==0:
				same+=1	
			outfile.write(gene1+"\t"+gene2+"\t"+ str(better)+"\t"+str(similaritywithhomology[gene1][gene2])+"\t"+lcswithhomology[gene1][gene2]+ "\t" + comparingwithhomology[gene1][gene2]+"\t"+    lcsstringwith+"\t"+str(similaritywithouthomology[gene1][gene2])+ "\t"+lcswithouthomology[gene1][gene2]+ "\t"+comparingwithouthomology[gene1][gene2]+"\t"+lcsstringwithout+"\n")

	print "same, with_greater, without_greater",same,with_greater,without_greater







if __name__ == "__main__":
	from copy import deepcopy 
	import math
	import sys
	main()