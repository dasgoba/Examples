#!/usr/local/bin/python
# -*- coding: utf-8 -*-

### Python packages required: Biopython

import time, subprocess, os
from Bio import Entrez, SeqIO
from eutils import client
from ftplib import FTP
import gzip
import StringIO
#from StringIO import StringIO
import numpy as np
import time 
import re
#from pybedtools import BedTool




localtime = time.asctime( time.localtime(time.time()) )
print("Starting: "+localtime)


###### Contact Info #####

Entrez.email = "gourab21_sit@jnu.ac.in"
Entrez.tool = "Genome Download"



##### Searching Pubmed and downloading abstracts using Pubmed Query or PMID list #####

def dwnldGenomeFasta(strp):

 search_results = Entrez.read(Entrez.esearch(db="nuccore", term=strp, usehistory="y", retmax=99999))
 gis = search_results["IdList"]
 #print(gis)
 ###### Extracting Refseq Accession from GIs #####
 

 result = Entrez.read(Entrez.epost("nucleotide",id=",".join(map(str,gis))))
 webEnv = result["WebEnv"]
 queryKey = result["QueryKey"] 
 records = list(SeqIO.parse(Entrez.efetch(db="nucleotide",retmode="text", rettype="fasta", webenv=webEnv, query_key=queryKey, retmax=99999), "fasta"))
 return(records)

def dwnldGenomeGB(strp):

 search_results = Entrez.read(Entrez.esearch(db="nuccore", term=strp, usehistory="y", retmax=99999))
 gis = search_results["IdList"]
 
 ###### Extracting Refseq Accession from GIs #####
 

 result = Entrez.read(Entrez.epost("nucleotide",id=",".join(map(str,gis))))
 webEnv = result["WebEnv"]
 queryKey = result["QueryKey"] 
 records = list(SeqIO.parse(Entrez.efetch(db="nuccore",retmode="text", rettype="gbwithparts", webenv=webEnv, query_key=queryKey, retmax=99999), "gb"))
 return(records)


def accessFilter(records): ####### Filtering finished genomes from incomplete ones #####
 
 for index, record in enumerate(records):
  if "NC_" in record.id:
	  print(record.id)
	 # refseqid = record.id.split("|")[3].split(".")[0]
	 # refseqorgn = record.description.split("|")[4].split(",")[0]
	  #print(refseqid+"	"+refseqorgn)
	  print("#########################################################")
	  with open("~/Public/perfect/tmp/temp.gbk", "w") as outgbk_handle:
		  outgbk_handle.write(record.format("gb"))
		  with open("~/Public/perfect/tmp/temp.gff", "w") as outgff_handle:
			GFF.write(record, outgff_handle)
	  with open("~/Public/perfect/tmp/temp.fasta", "w") as outfas_handle:
		  outfas_handle.write(record.format("fasta"))
	  
	  #startup()
	  #microsat(refseqid, refseqorgn)
	  #tempRemove_dir()
	  

def startup():
  print("Making temporary Directories...")
  os.system("~/Public/perfect/Bacteria/bacteria_setup/shell_scripts/startup.sh")
	  
def tempRemove_dir():
  print("Removing temporary Directories...")
  os.system("~/Public/perfect/Bacteria/bacteria_setup/shell_scripts/tmp_remove.sh")
	  

def microsat(refseqid, refseqorgn):
 #subprocess.call("~/Public/perfect/Bacteria/bacteria_setup/shell_scripts/auto_run_troll.sh".format(refseqid))
 os.system("~/Public/perfect/Bacteria/bacteria_setup/shell_scripts/auto_run_troll.sh %s" % (refseqid))
 os.system("~/Public/perfect/Bacteria/bacteria_setup/shell_scripts/troll_retrive.sh %s" % (refseqid))
 os.system("~/Public/perfect/Bacteria/bacteria_setup/shell_scripts/micro_annot.sh %s" % (refseqid))

 print("Microsatellite extraction completed...")


def ftpNCBILogin():
 ftp = FTP('ftp.ncbi.nlm.nih.gov')
 ftp.login() # Username: anonymous password: anonymous@
 return(ftp)


def ftpNCBISequence(ftp):

 sio = StringIO.StringIO()
 def handle_binary(more_data):
    sio.write(more_data)

 resp = ftp.retrbinary("RETR genomes/refseq/plasmid/plasmid.2.genomic.gbff.gz", callback=handle_binary)
 sio.seek(0) # Go back to the start
 zippy = gzip.GzipFile(fileobj=sio)
 uncompressed = zippy.read()
 
 return(uncompressed)


def ftpNCBIFeatureTable(ftp):

 sio = StringIO.StringIO()
 def handle_binary(more_data):
    sio.write(more_data)

 floc = "/genomes/all/GCF_000006945.1_ASM694v1/GCF_000006945.1_ASM694v1_feature_table.txt.gz"
 var1 = 'RETR ' + floc
 resp = ftp.retrbinary(var1 , callback=handle_binary)
 sio.seek(0) # Go back to the start
 zippy = gzip.GzipFile(fileobj=sio)
 uncompressed = zippy.read()
 
 return(uncompressed)


#def genIgen(data):
 

def ftpNCBIGenomes(ftp,strp):
 q = StringIO.StringIO()
 ftp.retrlines('RETR /genomes/GENOME_REPORTS/prokaryotes.txt', lambda line: q.write(line + '\n'))
 q.seek(0)
 baids = q.readlines()
 tmpids = filter(lambda x:re.search(strp, x) , baids) #### Search for species
 tids = filter(lambda x:re.search(r'Complete Genome', x) , tmpids) #### Search for complete genome Tags
 temids = "".join(s for s in tids if "NC_" in s) ##### Search for complete genomic molecule
 refseqids = np.genfromtxt(StringIO.StringIO(temids), dtype=None, delimiter='\t', usecols=(21))
 #print(refseqids)
 aids = "|".join(refseqids)


 p = StringIO.StringIO()
 ftp.retrlines('RETR genomes/refseq/bacteria/Escherichia_coli/assembly_summary.txt', lambda line: p.write(line + '\n'))

 p.seek(0)
 assids = p.readlines()
 #print(assids)


 tids = filter(lambda x:re.findall(aids, x), assids)
 tpids = "".join(tids)


 refsqids = np.genfromtxt(StringIO.StringIO(str(tpids)), dtype=None, delimiter='\t', usecols= (19))
 #print(refsqids)
 np.savetxt('/home/gourab/Public/perfect/tmp/bac_list.csv', refsqids, delimiter='\t', fmt="%s")

# return(refseqids)
 



#strp = "Salmonella enterica subsp. enterica"
strp = "Escherichia coli"

ftp=ftpNCBILogin()
#data=ftpNCBISequence(ftp)
ftpNCBIGenomes(ftp,strp)
#data = ftpNCBIFeatureTable(ftp)
#print(data)






#strp= input("Please enter your query for Genome search: ") ## input your query
#strp = "(((gene) OR protein) AND diabetes mellitus) AND human[MeSH Terms] AND hasabstract[text]"
#strp = '"salmonella enterica"[orgn] AND srcdb_refseq[prop] AND "complete genome"[ti] NOT plasmid[filt]'
strp ='NC_003198'
print("Query is: "+strp)
records = dwnldGenomeFasta(strp)
#records = dwnldGenomeGB(strp)
accessFilter(records)

localtime = time.asctime( time.localtime(time.time()) )
print("Ended: "+localtime)
