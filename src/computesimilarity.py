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





def getsimilarityold(profile1,profile2,termic,ancestors,termtermsimilarity,termtermlcs):
	icscores=[]
	finalmaxic=0
	finallcs="None"
	for term1 in profile1:
		for term2 in profile2:
			present=0
			if term1 in termtermsimilarity: 
				if term2 in termtermsimilarity[term1]:
					icscores.append(termtermsimilarity[term1][term2])
					present=1
					if termtermsimilarity[term1][term2] >finalmaxic:
						finalmaxic=termtermsimilarity[term1][term2]
						finallcs=termtermlcs[term1][term2]

			if present==0:
				maxic=0
				commonancestors=set.intersection(ancestors[term1],ancestors[term2])
				for anc in commonancestors:
					if termic[anc]>maxic:
						maxic=termic[anc]
						if termic[anc]>finalmaxic:
							finalmaxic=termic[anc]
							finallcs=anc
							

				icscores.append(maxic)

				if term1 not in termtermsimilarity:
					termtermsimilarity[term1]=dict()
				if term2 not in termtermsimilarity:
					termtermsimilarity[term2]=dict()					
				termtermsimilarity[term1][term2]=maxic
				termtermsimilarity[term2][term1]=maxic
	
	return max(icscores),termtermsimilarity


def loadancestors():
	parentshomology=open("../data/ParentsWithHomology.txt")
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

	parentsnohomology=open("../data/ParentsWithoutHomology.txt")
	
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


def getorthologpairs(geneset):
	datafile=open("../data/gene-orthologs.txt",'r')
	orthologpairs=[]
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
	zfin=open("../data/Dr-gene-to-phenotype-BF.txt",'r')
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

	mgi=open("../data/Mm-gene-to-phenotype-BF.txt",'r')		
	for line in mgi:
		gene=line.split("\t")[0].strip()
		annotation=line.split("\t")[1].strip().replace(":","_")
		if gene not in mgiprofiles:
			mgiprofiles[gene]=set()
		mgiprofiles[gene].add(annotation)		
	mgi.close()
	return(mgiprofiles)

def main():
	
	
	outfile=open("../results/SimilarityScores.txt",'w')
	outfile.write("Gene1\tGene2\tBetter Similarity With Homology\tSimilarity With Homology\tLCS With Homology\tAnnotation Pair leading to best match\tCommon subsumer set\tSimilarity Without Homology\tLCS Without Homology\tAnnotation Pair leading to best match\tCommon subsumer set\t\n")
	
	genegenelcs=dict()

	ancestorswith,ancestorswithout=loadancestors()
	geneset,termic=prepareCorpusforIC(ancestorswith,ancestorswithout)
	orthologpairs=getorthologpairs(geneset)
	#build profiles of genes
	mgiprofiles=loadmgiprofiles()
	zfinprofiles=loadzfinprofiles()

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


		withhomologyic,withhomologylcs,withhomologylcsset,withhomologycomparing=getsimilarity(mgiprofiles[gene1],zfinprofiles[gene2],termic,ancestorswith)

		lcswithhomology[gene1][gene2]=withhomologylcs
		similaritywithhomology[gene1][gene2]= withhomologyic
		lcssetwithhomology[gene1][gene2]=withhomologylcsset
		comparingwithhomology[gene1][gene2]=withhomologycomparing



		withouthomologyic,withouthomologylcs,withhomologylcsset,withhomologycomparing=getsimilarity(mgiprofiles[gene1],zfinprofiles[gene2],termic,ancestorswithout)
		lcswithouthomology[gene1][gene2]=withouthomologylcs
		similaritywithouthomology[gene1][gene2]= withouthomologyic	
		lcssetwithouthomology[gene1][gene2]=withhomologylcsset
		comparingwithouthomology[gene1][gene2]=withhomologycomparing	



	for gene1 in similaritywithhomology:
		for gene2 in similaritywithhomology[gene1]:
			lcsstringwith=','.join(lcssetwithhomology[gene1][gene2])
			lcsstringwithout=','.join(lcssetwithouthomology[gene1][gene2])
			better=0
			if similaritywithhomology[gene1][gene2] > similaritywithouthomology[gene1][gene2]:
				better=1


			outfile.write(gene1+"\t"+gene2+"\t"+ str(better)+"\t"+str(similaritywithhomology[gene1][gene2])+"\t"+lcswithhomology[gene1][gene2]+ "\t" + comparingwithhomology[gene1][gene2]+"\t"+    lcsstringwith+"\t"+str(similaritywithouthomology[gene1][gene2])+ "\t"+lcswithouthomology[gene1][gene2]+ "\t"+comparingwithouthomology[gene1][gene2]+"\t"+lcsstringwithout+"\n")







if __name__ == "__main__":
	from copy import deepcopy 
	import math
	import sys
	main()