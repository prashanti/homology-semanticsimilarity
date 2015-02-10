
def monarchhomologs(datafile):
	datafile=open("../data/panther-orthologs.owl",'r')
	for line in datafile:
		if "RO_HOM0000019" in line:
			
			data=line.replace("AnnotationAssertion(","").replace(")","").replace("<http://identifiers.org/","").replace(">","").split()
			outfile.write(data[1]+"\t"+data[2].strip()+"\n")


def getMGI2ZFINorthologs(datafile):
	outfile=open("../data/MGI_ZFIN_orthologs.tsv",'w')
	for line in datafile:
		if "LDO" in line:
			if "MGI" in line and "ZFIN" in line:
				#MOUSE|MGI=MGI=88436|UniProtKB=P29974	DANRE|ZFIN=ZDB-GENE-090312-66|UniProtKB=B8JHY6	LDO	Osteichthyes	PTHR10217
				line=line.replace("MGI=MGI=","MGI:").replace("ZFIN=","ZFIN:")
				data=line.split("\t")
				gene1=data[0].split("|")[1]
				gene2=data[1].split("|")[1]
				outfile.write(gene1+"\t"+gene2+"\n")

def loadNCBIGenes():
	datafile=open("../data/Human_NCBI_GeneIDs.txt")
	NCBIgenes=set()
	for ncbiid in datafile:
		NCBIgenes.add(ncbiid.strip())
	datafile.close()	
	return NCBIgenes


def getHUMAN2MGIorthologs(datafile):
	ENSEMBL2NCBI=loadENSEMBL2NCBI()
	NCBIgenes_annotated=loadNCBIGenes()
	NCBIgenes_orthologs=set()
	outfile=open("../data/Human_MGI_orthologs.tsv",'w')
	for line in datafile:
		if "LDO" in line:
			if "MGI" in line and "HUMAN" in line:
				#HUMAN|Ensembl=ENSG00000165533|UniProtKB=Q8TAM2	MOUSE|MGI:1923510|UniProtKB=Q8VD72	LDO	Euarchontoglires	PTHR23083
				line=line.replace("MGI=MGI=","MGI:")
				data=line.split("\t")
				gene1=data[0].split("|")[1].replace("Ensembl=","")
				gene2=data[1].split("|")[1].replace("Ensembl=","")

				if "ENS" in gene1:
					if gene1 in ENSEMBL2NCBI:
						ncbiid=ENSEMBL2NCBI[gene1]
						NCBIgenes_orthologs.add(ncbiid)
						outfile.write(gene2+"\t"+ncbiid+"\n")
				else:
					if gene2 in ENSEMBL2NCBI:
						ncbiid=ENSEMBL2NCBI[gene2]
						NCBIgenes_orthologs.add(ncbiid)
						outfile.write(gene1+"\t"+ncbiid+"\n")	 
	print len(set.difference(NCBIgenes_annotated,NCBIgenes_orthologs)),set.difference(NCBIgenes_annotated,NCBIgenes_orthologs)



def getHUMAN2ZFINorthologs(datafile):
	ENSEMBL2NCBI=loadENSEMBL2NCBI()
	NCBIgenes_annotated=loadNCBIGenes()
	NCBIgenes_orthologs=set()
	outfile=open("../data/Human_ZFIN_orthologs.tsv",'w')
	for line in datafile:
		if "LDO" in line:
			if "ZFIN" in line and "HUMAN" in line:
				#HUMAN|Ensembl=ENSG00000166454|UniProtKB=O43313	DANRE|ZFIN=ZDB-GENE-030616-220|UniProtKB=F1R9Z7	LDO	Osteichthyes	PTHR10593
				line=line.replace("ZFIN=","ZFIN:")

				#HUMAN|Ensembl=ENSG00000172315|UniProtKB=Q96S44	DANRE|ZFIN:ZDB-GENE-050522-123|UniProtKB=F1QB90	LDO	Osteichthyes	PTHR12209
				data=line.split("\t")
				gene1=data[0].split("|")[1].replace("Ensembl=","")
				gene2=data[1].split("|")[1].replace("Ensembl=","")

				if "ENS" in gene1:
					if gene1 in ENSEMBL2NCBI:
						ncbiid=ENSEMBL2NCBI[gene1]
						NCBIgenes_orthologs.add(ncbiid)
						outfile.write(gene2+"\t"+ncbiid+"\n")
				else:
					if gene2 in ENSEMBL2NCBI:
						ncbiid=ENSEMBL2NCBI[gene2]
						NCBIgenes_orthologs.add(ncbiid)
						outfile.write(gene1+"\t"+ncbiid+"\n")	 
	#print len(set.difference(NCBIgenes_annotated,NCBIgenes_orthologs)),set.difference(NCBIgenes_annotated,NCBIgenes_orthologs)


def loadENSEMBL2NCBI():
	ENSEMBL2NCBI=dict()
	datafile=open("../data/Homo_sapiens-ensembl_gene_id-entrezgene.mart")
	for line in datafile:
		ENSEMBL2NCBI[line.split("\t")[0]]="NCBI_gene:"+line.split("\t")[1].strip()
	datafile.close()	
	return ENSEMBL2NCBI	

def main():
	datafile=open("../data/RefGenomeOrthologs",'r')
	#getHUMAN2MGIorthologs(datafile)
	getHUMAN2ZFINorthologs(datafile)
	datafile.close()




if __name__ == "__main__":
	main()