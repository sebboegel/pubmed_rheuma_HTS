import os,sys

os.system("pysradb search --accession SRP042288 -v 3 --saveto RA.pysradb.tsv")

for line in open("./RA.pysradb.tsv"):
	if "pubmed" in line:
		l=line[0:-1].split("\t")
		column_names=l
		pubmedid=l[column_names.index("DB: pubmed")+1].split(" ")[1]
		print(pubmedid)