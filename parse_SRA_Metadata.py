import sys,os


file_assay = open("sra_assays.tsv", "w")
file_studies = open("sra_studies.tsv", "w")
file_instrument = open("sra_instrument.tsv", "w")
file_layout=open("sra_layout.tsv", "w")
file_source=open("sra_source.tsv", "w")
file_tissue=open("sra_tissue.tsv", "w")
file_org=open("sra_organism.tsv", "w")
file_phenotype=open("sra_phenotype.tsv", "w")
file_srp_pubmed = open("sra_pubmed.tsv", "w")
file_srp_pubmed.write("disease\tsrp\tassay\tpubmedid\n")
file_studies.write("disease\tsrp\tnumber\n")
file_assay.write("disease\tsamples\tassay\n")
file_instrument.write("disease\tsamples\tinstrument\n")
file_layout.write("disease\tsamples\tlayout\n")
file_source.write("disease\tsamples\tsource\n")
file_tissue.write("disease\tsamples\ttissue\n")
file_org.write("disease\tsamples\torganism\n")
file_phenotype.write("disease\tsamples\tphenotype\n")

sra_metadata_files=sys.argv[1]
sra_metadata_file=sra_metadata_files.split(",")
for file in sra_metadata_file:
	i=0
	assay_dict={}
	instrument_dict={}
	layout_dict={}
	source_dict={}
	srp_dict={}
	tissue_dict={}
	phenotype_dict={}
	organism_dict={}
	pubmed_dict={}
	tissue_dict["blood/immune cells"]=0
	tissue_dict["cartilage"]=0
	tissue_dict["joint"]=0
	tissue_dict["knee"]=0
	tissue_dict["synovium"]=0
	phenotype_dict["healthy"]=0
	tissue_dict["spleen"]=0
	tissue_dict["lymph node"]=0
	tissue_dict["retina"]=0
	tissue_dict["bone marrow"]=0
	tissue_dict["lymph node"]=0
	tissue_dict["fibroblasts"]=0
	instrument_dict["illumina HiSeq"]=0
	instrument_dict["genome analyzer"]=0
	tissue_dict["not defined"]=0
	phenotype_dict["not defined"]=0
	tissue_dict["soft tissue"]=0
	tissue_dict["salivary gland"]=0
	organism_dict["mus musculus"]=0
	organism_dict["homo sapiens"]=0
	instrument_dict["454 gs"]=0
	instrument_dict["nextseq 5x0"]=0
	tissue_dict["muscle"]=0
	tissue_dict["stool"]=0
	tissue_dict["muscle"]=0
	tissue_dict["kidney"]=0
	organism_dict["sus scrofa"]=0
	organism_dict["rattus norvegicus"]=0
	tissue_dict["skin"]=0
	organism_dict["not defined"]=0
	tissue_dict["retina"]=0
	phenotype_dict["disease"]=0
	disease=file.split("_")[0]
	pubmed_dict[disease]={}
	for line in open(file):
		l=line[0:-1].lower().split("\t")
		if i==0:
			i=1
			column_names=l
		else:
			#----------------assay----------------------
			assay=l[column_names.index("assay_type")]
			if assay not in assay_dict:
				assay_dict[assay]=1
			else:
				assay_dict[assay]+=1
			#---------------instrument-----------------------
			instrument=l[column_names.index("instrument")].lower()
			if "nextseq 5" in instrument or "hiseq" in instrument or "genome analyzer" in instrument or "454 gs" in instrument:
				if "hiseq" in instrument:
					instrument_dict["illumina HiSeq"]+=1
				if "genome analyzer" in instrument:
					instrument_dict["genome analyzer"]+=1
				if "454 gs" in instrument:
					instrument_dict["454 gs"]+=1
				if "nextseq 5" in instrument:
					instrument_dict["nextseq 5x0"]+=1
			else:
				if instrument not in instrument_dict:
					instrument_dict[instrument]=1
				else:
					instrument_dict[instrument]+=1
			#--------------layout----------------
			layout=l[column_names.index("librarylayout")]
			if layout not in layout_dict:
				layout_dict[layout]=1
			else:
				layout_dict[layout]+=1
			#--------------source----------------------
			source=l[column_names.index("librarysource")]
			if source not in source_dict:
				source_dict[source]=1
			else:
				source_dict[source]+=1
			#--------------srp-----------------------------
			srp=l[column_names.index("sra_study")]
			
			if srp not in srp_dict:
				srp_dict[srp]=1
			else:
				srp_dict[srp]+=1
			#if disease=="SLE":
			#	if not srp in pubmed_dict[disease]:
			#		os.system("pysradb search --accession "+srp+" -v 3 --saveto pysradb.tsv")
				 
			#		for row in open("./pysradb.tsv"):
			#			if "DB: pubmed" in row:
			#				r=row[0:-1].split("\t")
			#				row_names=r
							#print(row)
			#				pubmedid=r[row_names.index("DB: pubmed")+1].split(" ")[1]
			#			else:
			#				pubmedid=""	
			#		pubmed_dict[disease][srp]=assay+"\t"+pubmedid	
				#print(disease,srp,assay,pubmedid)
						
			#--------------organism-----------------------------
			organism=l[column_names.index("organism")].lower()
			if organism=="gut metagenome" or organism=="metagenome" or organism=="metagenomes" or "rat" in organism or "sus scrofa" in organism or "human" in organism or "homo" in organism or "mouse" in organism or "mus" in organism:
				if "sus scrofa" in organism:
					organism_dict["sus scrofa"]+=1
				if "rat" in organism:
					organism_dict["rattus norvegicus"]+=1
				if "human" in organism or "homo" in organism:
					organism_dict["homo sapiens"]+=1
				if "mouse" in organism or "mus" in organism:
					organism_dict["mus musculus"]+=1
				if organism=="gut metagenome" or organism=="metagenome" or organism=="metagenomes":
					organism_dict["not defined"]+=1
			else:
				if organism not in organism_dict:
					organism_dict[organism]=1
				else:
					organism_dict[organism]+=1
			#--------------tissue------------------------
			try:
				if "tissue" in column_names:
					tissue=l[column_names.index("tissue")].lower()
				if "tissue_type" in column_names:
					tissue=l[column_names.index("tissue_type")].lower()
				if "organism_part" in column_names:
					tissue=l[column_names.index("organism_part")].lower()
				if "env_material" in column_names and l[column_names.index("env_material")]!="":
					tissue=l[column_names.index("env_material")].lower()
				if "isolation_source" in column_names and l[column_names.index("isolation_source")]!="":
					tissue=l[column_names.index("isolation_source")].lower()
				if "soft tissue" in tissue:
					print("hallo")			
				if "salivary" in tissue or "muscle" in tissue or "ln" in tissue or "lymph node" in tissue or "skin" in tissue or "retina" in tissue or "spleen" in tissue or "spln" in tissue or "tls" in tissue or "homo sapiens" in tissue or "dendritic" in tissue or "missing" in tissue or "beijing" in tissue or "serum" in tissue or "bal" in tissue or "lung" in tissue or "stool" in tissue or "feces" in tissue or "bone" in tissue or "blood" in tissue or "pbmc" in tissue or "t cell" in tissue or "t-cell" in tissue or "t lymphocyte" in tissue or "synovi" in tissue or "cartilage" in tissue or "catilage" in tissue or "plate" in tissue or "monocyte" in tissue:
					if "spleen" in tissue or "spln" in tissue:
						tissue_dict["spleen"]+=1
					if "ln" in tissue or "lymph node" in tissue:
						tissue_dict["lymph node"]+=1
					if "retina" in tissue:
						tissue_dict["retina"]+=1
					if "skin" in tissue:
						tissue_dict["skin"]+=1
					if "muscle" in tissue:
						tissue_dict["muscle"]+=1
					if "salivary" in tissue:
						tissue_dict["salivary gland"]+=1
					if "synovi" in tissue:
						tissue_dict["synovium"]+=1
					if "bal" in tissue or "lung" in tissue:
						tissue_dict["lung"]+=1
					if "stool" in tissue or "feces" in tissue:
						tissue_dict["stool"]+=1
					
					
					if "cartilage" in tissue or "catilage" in tissue:
						tissue_dict["cartilage"]+=1
					if "plate" in tissue or "bone" in tissue:
						tissue_dict["bone"]+=1
					if "serum" in tissue or "blood" in tissue or "pbmc" in tissue or "dendritic" in tissue or "monocytes" in tissue or "t cell" in tissue or "t-cell" in tissue or "t lymphocyte" in tissue:
						tissue_dict["blood/immune cells"]+=1
					
					if  "homo sapiens" in tissue or "missing" in tissue or "beijing" in tissue or "tls" in tissue:
						tissue_dict["not defined"]+=1
				else:
					if tissue=="":
						tissue_dict["not defined"]+=1
					else:
						if tissue not in tissue_dict:
							tissue_dict[tissue]=1
						else:
							tissue_dict[tissue]+=1
			except:
				tissue_dict["not defined"]+=1
				
			try:
				if "genotype" in column_names and l[column_names.index("genotype")]!="":
					phenotype=l[column_names.index("genotype")].lower()
				if "study_group" in column_names and l[column_names.index("study_group")]!="":
					phenotype=l[column_names.index("study_group")].lower()	
				if "disease" in column_names:
					phenotype=l[column_names.index("disease")].lower()
				if "disease_state" in column_names:
					phenotype=l[column_names.index("disease_state")].lower()
				if "disease_stage" in column_names:
					phenotype=l[column_names.index("disease_stage")].lower()
				if "phenotype" in column_names:	
					phenotype=l[column_names.index("phenotype")].lower()
				if "control" in phenotype or "intact" in phenotype or "healthy" in phenotype or "normal" in phenotype:
					phenotype_dict["healthy"]+=1
				else:
					if phenotype=="":
						phenotype_dict["not defined"]+=1
					else:
						phenotype_dict["disease"]+=1
			except:
				phenotype_dict["not defined"]+=1
			#print(assay,instrument,layout,source,srp,tissue,phenotype)
	for assay in assay_dict:
		out=disease+"\t"+str(assay_dict[assay])+"\t"+assay+"\n"
		file_assay.write(out)
	
	for instrument in instrument_dict:
		out=disease+"\t"+str(instrument_dict[instrument])+"\t"+instrument+"\n"
		file_instrument.write(out)
		
	for layout in layout_dict:
		out=disease+"\t"+str(layout_dict[layout])+"\t"+layout+"\n"
		file_layout.write(out)
		
	for source in source_dict:
		out=disease+"\t"+str(source_dict[source])+"\t"+source+"\n"
		file_source.write(out)
	
	for srp in srp_dict:
		out=disease+"\t"+srp+"\t"+str(srp_dict[srp])+"\n"
		file_studies.write(out)
		
	for tissue in tissue_dict:
		out=disease+"\t"+str(tissue_dict[tissue])+"\t"+tissue+"\n"
		file_tissue.write(out)
	
	for org in organism_dict:
		out=disease+"\t"+str(organism_dict[org])+"\t"+org+"\n"
		file_org.write(out)
		
	for phenotype in phenotype_dict:
		out=disease+"\t"+str(phenotype_dict[phenotype])+"\t"+phenotype+"\n"
		file_phenotype.write(out)
		
	for srp in  pubmed_dict[disease]:
		out=disease+"\t"+srp+"\t"+pubmed_dict[disease][srp]+"\n"
		file_srp_pubmed.write(out)
	
	
	

file_srp_pubmed.close()
file_assay.close()
file_studies.close()
file_instrument.close()
file_layout.close()
file_source.close()
file_tissue.close()
file_org.close()
file_phenotype.close()

