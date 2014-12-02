
def monarchhomologs(outfile):
	datafile=open("./data/panther-orthologs.owl",'r')
	for line in datafile:
		if "RO_HOM0000019" in line:
			
			data=line.replace("AnnotationAssertion(","").replace(")","").replace("<http://identifiers.org/","").replace(">","").split()
			outfile.write(data[1]+"\t"+data[2].strip()+"\n")

def main():
	datafile=open("./data/RefGenomeOrthologs",'r')
	outfile=open("./data/gene-orthologs.txt",'w')
	for line in datafile:
		if "LDO" in line:
			if "MGI" in line and "ZFIN" in line:
				line=line.replace("MGI=MGI=","MGI:").replace("ZFIN=","ZFIN:")
				data=line.split("\t")
				gene1=data[0].split("|")[1]
				gene2=data[1].split("|")[1]


				#MOUSE|MGI=MGI=88436|UniProtKB=P29974	DANRE|ZFIN=ZDB-GENE-090312-66|UniProtKB=B8JHY6	LDO	Osteichthyes	PTHR10217
				outfile.write(gene1+"\t"+gene2+"\n")
	
	outfile.close()
	datafile.close()




if __name__ == "__main__":
	main()