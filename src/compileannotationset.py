def getsimilarity(profile1,profile2,termic,ancestors):
	
	finalmaxic=0
	finallcs="None"
	
	for term1 in profile1:
		for term2 in profile2:
			commonancestors=set.intersection(ancestors[term1],ancestors[term2])
			for anc in commonancestors:
				if termic[anc]>finalmaxic:
					finalmaxic=termic[anc]
					finallcs=anc	

	return finalmaxic,finallcs





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


def loadhomology():
	parentshomology=open("../ParentsWithHomology.txt")
	ancestors=dict()
	for line in parentshomology:
		data=line.replace("<http://purl.obolibrary.org/obo/","").replace(">","").split("\t")
		child=data[0].strip()
		parent=data[1].strip()
		if child not in ancestors:
			ancestors[child]=set()
		if "WithHomology" in parent:
			ancestors[child].add(parent)
	return ancestors
	parentshomology.close()


def loadnohomology():
	parentsnohomology=open("../ParentsWithoutHomology.txt")
	ancestors=dict()
	for line in parentsnohomology:
		data=line.replace("<http://purl.obolibrary.org/obo/","").replace(">","").split("\t")
		child=data[0].strip()
		parent=data[1].strip()
		if child not in ancestors:
			ancestors[child]=set()
		if "Withouthomology" in parent:
			ancestors[child].add(parent)
	parentsnohomology.close()
	return ancestors

def getorthologpairs(geneset):
	datafile=open("./data/gene-orthologs.txt",'r')
	orthologpairs=[]
	for line in datafile:
		#ensembl:ENSG00000005007	PomBase:SPAC16C9.06c
		line=line.replace("NCBIGene_","NCBI_gene:")
		gene1=line.split("\t")[0]
		gene2=line.split("\t")[1].strip()
		if gene1 in geneset and gene2 in geneset:
			temp=[gene1,gene2]
			orthologpairs.append(temp)
	return orthologpairs

def prepareCorpusforIC(ancestors):
	geneset=set()
	annotationgenecount=dict()
	annotationset=set()
	geneannotations=dict()
	annotationfile=open("./data/AnnotationCorpus.txt")
	outfile=open("./data/AnnotationDataset.txt",'w')
	for line in annotationfile:
		data=line.split("\t")
		gene=data[0].strip()
		geneset.add(gene)
		annotation=data[1].strip().replace(":","_")
		annotationset.add(annotation)
		if gene not in geneannotations:
			geneannotations[gene]=set()
		geneannotations[gene].add(annotation)
		for anc in ancestors[annotation]:
			geneannotations[gene].add(anc)
			annotationset.add(anc)
	
	for gene in geneannotations:
		for annotation in geneannotations[gene]:
			if annotation not in annotationgenecount:
				annotationgenecount[annotation]=0
			annotationgenecount[annotation]+=1
		outfile.write(gene+"\t"+(",".join(geneannotations[gene]))+"\n")
	outfile.close()
	annotationfile.close()
	return geneset,annotationgenecount

def loadzfinprofiles():
	zfinprofiles=dict()
	zfin=open("./data/Danio_rerio/Dr-gene-to-phenotype-BF.txt",'r')
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

	mgi=open("./data/Mus_musculus/Mm-gene-to-phenotype-BF.txt",'r')		
	for line in mgi:
		gene=line.split("\t")[0].strip()
		annotation=line.split("\t")[1].strip().replace(":","_")
		if gene not in mgiprofiles:
			mgiprofiles[gene]=set()
		mgiprofiles[gene].add(annotation)		
	mgi.close()
	return(mgiprofiles)

def main():
	
	
	
	ancestors=dict()
	termic=dict()
	genegenelcs=dict()
	ancestors=loadhomology()
	#ancestors=loadnohomology()

	geneset,annotationgenecount=prepareCorpusforIC(ancestors)
	geneoutfile=open("annotationgenecounts.txt",'w')

	
	
	orthologpairs=getorthologpairs(geneset)

	for annotation in annotationgenecount:
		geneoutfile.write(annotation+"\t"+str(annotationgenecount[annotation])+"\n")	
		termic[annotation]=-math.log(float(annotationgenecount[annotation])/float(len(geneset)))
			
	geneoutfile.close()


	
	#build profiles of genes
	mgiprofiles=loadmgiprofiles()
	zfinprofiles=loadzfinprofiles()

	genegenesimilarity=dict()
	termtermsimilarity=dict()

	for orthologpair in orthologpairs:
		gene1=orthologpair[0]
		gene2=orthologpair[1]
	
		if gene1 not in genegenesimilarity:
			genegenesimilarity[gene1]=dict()
		if gene1 not in genegenelcs:
			genegenelcs[gene1]=dict()

		ic,lcs=getsimilarity(mgiprofiles[gene1],zfinprofiles[gene2],termic,ancestors)
		genegenelcs[gene1][gene2]=lcs
		genegenesimilarity[gene1][gene2]=ic
		
	for gene1 in genegenesimilarity:
		for gene2 in genegenesimilarity[gene1]:
			print gene1+"\t"+gene2+"\t"+str(genegenesimilarity[gene1][gene2])+"\t"+genegenelcs[gene1][gene2]







if __name__ == "__main__":
	import math
	import sys
	main()