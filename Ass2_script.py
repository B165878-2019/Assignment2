#!/usr/bin/python3
import sys
import subprocess
import os
import xml.etree.ElementTree as ET
import re
import pandas as pd
import shutil
##### Entrez utilities
cwd = os.getcwd()
end = (cwd+'/B165878_Programme_test')
if not os.path.exists(end) :
 os.mkdir('B165878_Programme_test')
#
taxonomy = input('\nTaxon name : ')
protname = input('\nProtein name : ')
#
# Check spelling
def spellcheck(taxon, protname):
#Parsing xml file and extracting corrections, if any
 espelltax = subprocess.check_output(['espell', '-db', 'taxonomy', '-query', taxon])
 taxroot = ET.fromstring(espelltax)
 taxquery= taxroot[1].text
 taxcorrection = taxroot[2].text
 espellprot = subprocess.check_output(['espell', '-db', 'protein', '-query', protname])
 protroot = ET.fromstring(espellprot)
 protquery = protroot[1].text
 protcorrection = protroot[2].text
#Long list of loops to catch any corrections found in above xml, and to prompt user if they would like to use a corrected query, then set variables accordingly - they are all the same i will describe one
 if not taxcorrection and not protcorrection :
  spellcheck.tax = taxquery
  spellcheck.prot = protquery  
  return checkhits(spellcheck.tax, spellcheck.prot)
 elif taxcorrection and not protcorrection : 
# e.g. : above "elif" catches taxcorrection = True from xml file 
# catches protcorrection as false so no need to change - set variable 
  spellcheck.prot = protquery
# Lets the user know
  print('\nPossible spelling correction found\n\nQuery : '+taxquery+'\nCorrection : '+taxcorrection)
  while True:
#enter loop, selects for only y/n type answers otherwise prints error
   c = input('\nWould you like to use the corrected query? (yes/no) : ')
   if c.lower() in ['yes','no','y','n'] :
#if input in y/n move in, if "y" to prompted question, use correction if "no" use original query
    if c.lower() in ['yes','y'] :
#set variable and break
     spellcheck.tax = taxcorrection
     break
    if c.lower() in ['no','n'] :
#set variable and break
     spellcheck.tax = taxquery
     break
   else :
#catch incorrect input
    print('\n"'+c+'" not a valid entry, please try again....\n\ntype "yes" to use the spell corrected query or type "no" to use the original query')
 elif protcorrection and not taxcorrection :
  spellcheck.tax = taxquery
  print('\nPossible spelling correction found\n\nQuery : '+protquery+'\nCorrection : '+protcorrection)
  while True:
   c = input('\nWould you like to use the corrected query? (yes/no) : ')
   if c.lower() in ['yes','no','y','n'] :
    if c.lower() in ['yes','y'] :
     spellcheck.prot = protcorrection
     break
    if c.lower() in ['no','n'] :
     spellcheck.prot = protquery
     break
   else :
    print('\n"'+c+'" not a valid entry, please try again....\n\ntype "yes" to use the spell corrected query or type "no" to use the original query')
 elif protcorrection and taxcorrection :
  print('\nPossible spelling correction found\n\nQuery : '+taxquery+'\nCorrection : '+protcorrection+'\n\nQuery : '+protquery+'\nCorrection : '+protcorrection)
  while True:
   c = input('\nWould you like to use the corrected taxon query? (yes/no) : ')
   if c.lower() in ['yes','no','y','n'] :
    if c.lower() in ['yes','y'] :
     spellcheck.tax = taxcorrection
     break
    if c.lower() in ['no','n'] :
     spellcheck.tax = taxquery
     break
   else :
    print('\n"'+c+'" not a valid entry, please try again....\n\ntype "yes" to use the spell corrected query or type "no" to use the original query')
  while True :
   c = input('\nWould you like to use the corrected taxon query? (yes/no) : ')
   if c.lower() in ['yes','no','y','n'] :
    if c.lower() in ['yes','y'] :
     spellcheck.tax = taxcorrection
     break
    if c.lower() in ['no','n'] :
     spellcheck.tax = taxquery
     break
   else :
    print('\n"'+c+'" not a valid entry, please try again....\n\ntype "yes" to use the spell corrected query or type "no" to use the original query')
#Finally, feed variables to checkhits function
 return checkhits(spellcheck.tax, spellcheck.prot)
#Check hits
def checkhits(taxon, protname, manquery=''):
#takes arguments for taxon and protname, optionally takes a manquery (manual-query) argument (assigned in another function)
#Only catches zero values or >10000 values for hitcount
 print('\nNow checking hits for search....')
#use manquery if provided
 if manquery :
  checkhits.query_term = manquery
#otherwise use taxon and protname
 else :
  checkhits.query_term = (f"{taxon}[ORGN] AND {protname} NOT PARTIAL")
#search assigned query term and parse xml for hit count
 esearch = subprocess.check_output(['esearch', '-db', 'protein', '-query', checkhits.query_term])
 root = ET.fromstring(esearch)
 checkhits.hitcount = root[3].text
 print('\nFound '+checkhits.hitcount+' hits')
#check if too many sequences
 if int(checkhits.hitcount) >= 10000 :
  print('\nData sets larger than 10000 sequences are not recommended!')
  while True :
#loop only breaks in yes or no, yes continues, no redirects to search help
   big = input('\nDo you still wish to continue with these query terms? :')
   if big.lower() in ['yes','no','y','n']:
    if big in ['yes', 'y']:
     print('okay tough guy....')
     break
    elif big in ['no','n']:
     return helpsearch()
   else:
    print('\n'+big+' not a valid entry, type "yes" to use oversised sample or "no" to exit the programme and think about what you are actually searching....')
#check if zero hits, likely typo
 elif int(checkhits.hitcount) <= 50 : 
  while True :
#loop broken only with input yes or no, yes takes user to help, no exits
   retry = input('\nerror, sequence number too low, type "yes" to go to search help or "no" to exit the programme : ')
   if retry.lower() in ['yes','no','y','n']:
    if retry.lower() in ['yes', 'y']:
     return helpsearch()
    elif retry.lower() in ['no','n']:
     return sys.exit()
   else :
    print('\n"'+retry+'" not a valid option, please try again...')
 else :
  return fetchFASTA(checkhits.query_term)
#Search help 
def helpsearch(): 
#show what the literal query is
 print('\nThe current search term is findng an unsuitable number of hits within the database :\n\n**"'+checkhits.query_term+'"**')
 while True :
#again, only catch y/n
#"yes" allows manual input of query, also prints NCBI help on searching, then feeds checkhits the optional argument "manquery" to check if hitcount is okay
  man = input('\nTo manually enter a search query type "yes", to try again with new search terms type "no"\n\nType here : ')
  if man.lower() in ['yes','no','y','n']:
   if man.lower() in ['yes', 'y']:
    print('\n\nInformation on NCBI search field descriptions can be found here :\n"https://www.ncbi.nlm.nih.gov/books/NBK49540/"\nProtein nomenclature guidelines can be found here :\nhttps://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/')
    manqueryi = input('\n\nPlease type the EXACT search input here : ')    
    checkhits('none', 'none', manquery=(f"{manqueryi}"))
    break
   elif man.lower() in ['no','n']:
# if answer is no then simply input new terms and start again
    helpsearch.taxonomy = input('\nTaxon name : ')
    helpsearch.protname = input('\nProtein name : ')
    break
  else :
   print('\n"'+man+'" not a valid option, please try again...')
 return spellcheck(helpsearch.taxonomy, helpsearch.protname)
#fetching FASTA
def fetchFASTA(query_term):
 print('\nNumber of hits within acceptable range, moving to fetching stage...')
 while True :
#globally sets a title for all following files 
  global name
  name = input('\nPlease specify a filename, this will be the title to all concurrent files : ')   
#pathname variable includes path for easy sorting of all files to end dir
  global pathname
  pathname = f"{end}/{name}"
#continue only if name is alphanumerical
  if name.isalnum()  :   
   with open(f"{end}/{name}.fa.raw", "w+") as fetchout :
    print('\nFetchng '+checkhits.hitcount+' results, please wait.....')
#fetch fasta files for all hits
    esearch = subprocess.Popen(['esearch', '-db', 'protein', '-query', query_term], stdout=subprocess.PIPE)
    efetch = subprocess.Popen(['efetch', '-db', 'protein', '-format', 	'fasta'],stdin=esearch.stdout, stdout=fetchout)
#wait for efetch to finish otherwise sort can read empty file
    efetch.wait()
    print('\nDone!')
  else :
   print('\nSorry only alphanumeric characters are permitted, please try again')
#call sort function
  return(sort(f"{pathname}.fa.raw"))
#Sorting FASTA
def sort(inf):
 print('\nNow sorting sequences....\n')
#open sequence file for reading
 aves = open(f"{inf}", 'r+')
#create dict
 seqinfo = {"Accession" :[], "Protein_name" :[], "Species" :[]}
 for seq in aves:
#For all lines in sequence file find headers with '>'
  if re.search('>', seq) : 
   #seqinfo["Full_header"].append(seq) #no longer used
#regex term for finding accession, looks backward for '>' and then looks for any repetition of non-whitespace characters, accession should always be preceeded by '>' and followed by space
   accession = re.search('(?<=>)\S*', seq) 
#append to dict
   seqinfo["Accession"].append(accession.group())
#term looks for anything enclosed in square brackets, should be only species name, excludes brackets when appending to dict 
   speciesname = re.search('(?<=\[).*?(?=\])',seq) 
   seqinfo["Species"].append(speciesname.group())
#finds anything preceeded by a whitespace and followed by a "[", works fine as long as headers are formatted uniformly 
   protname = re.search('(?<=\s).*?(?=\[)', seq)
   seqinfo["Protein_name"].append(protname.group())   
#read dict in to dataframe
 df = pd.DataFrame(seqinfo)
#print dataframe because it looks cool
 print(df)
#total sequence count, gathered from total num of accessions extracted 
 count = str(len(seqinfo["Accession"]))
#non redundant list of species names 
 namelist = sorted(set(seqinfo["Species"]))
#number of individual species
 sort.indv = str(len(namelist))
#let user know how many sequences and species are represented 
 print('\nIdentified '+checkhits.hitcount+' sequences in '+sort.indv+' individual species')
 while True : 
#prompt to continue, only catch yes or no
  proceed = input('\nWould you like to continue with this sequence set? : ')
  if proceed.lower() in ['yes','no','y','n'] :
    if proceed.lower() in ['yes','y'] :
#if yes feed variables forward 
     return redun(namelist, df, count, inf)
    if proceed.lower() in ['no','n'] :
#if no, cleanup files unused and reset search
     for f in os.listdir(end) :
      os.remove(f)
     return helpsearch()
  else :
   print('\n"'+proceed+'" not a valid entry, please try again....\n\ntype "yes" to proceed or "no" to input new search terms')
#Sorting sequence redundancy within each species
#fed namelist : non-redundant list of species names, df : dataframe produced earlier, count : num of sequences, inf : fasta file to be processed
def redun(namelist, df, count, inf) :
 print('\nSequences will be sorted for redundancy within each species')
 while True :
#prompt for redundancy threshold
  threshold = input('\nPlease input the percentage similarity threshold for sorting redundancy : ')
#regex term, needs two digits but can be told a 4 digit float
  if re.search("\d\d\.?\d?\d?", threshold) :
#make temporary directory "trash"
   trash = f"{end}/{name}_FASTA"
   if not os.path.exists(trash):
    os.mkdir(trash)
#counter for iteration
   counter = 0
#This loop will iterate through each item in the non-redundant species namelist, and extract all corresponding accessions from the dataframe, writing output to "temp"
#Then will use pullseq to extract the correspondig fasta sequences from the original sequence file, writing output to "temp2"
#FASTA sequences in temp2 are then sorted for redundancy at the input threshold and output in the "trash" directory temporarily, each having a different index denoted by the counter
#Each iteration overwrites the previous "temp" and "temp2" files 
#output should be "trash" dir filled with name_counter.fa.sort files for every individual species
   for item in namelist :
#make temp files, both in write mode to overwrite last entry
    temp = open(f"{trash}/temp.txt", 'w+')
    temp2 = open(f"{trash}/temp2.txt", 'w+')
#path to original sequence file
    sequences = f"{inf}" 
#direct outfile to "trash" and label with counter
    outfile = f"{trash}/{name}_{counter}.fa.sort"
#extract all accessions for given species name 
    individual = df[df['Species'] == item]['Accession']
    b = "\n"
    i = b.join(list(individual.values))
#write accesions to "temp" (overwrites previous contents)
    temp.write(i)
    temp.close()
#set variables for pullseq and skipredundant input 
    targets = f"{trash}/temp.txt"
    targets2 = f"{trash}/temp2.txt"
#run pullseq
    pullseq = subprocess.check_output(['/localdisk/data/BPSM/Assignment2/pullseq', '-i', sequences, '-n', targets])
#decode output and write to temp2
    temp2.write(pullseq.decode('utf-8'))
    #print(pullseq.decode('utf-8')) #TEST
#sort redundant at threshold for pulled sequences
    skipredundant = subprocess.Popen(['skipredundant','-auto','-threshold', threshold, '-sequences', targets2, '-outseq', outfile])
#ammend counter and continue iteration
    counter += 1
#Loop end, wait for process to finish (otherwise next command executes too quickly and last file is lost)
   skipredundant.wait()
#delete previous file if present so duplicates arent made
   if os.path.exists(f"{end}/{name}.fa.sort") :
    os.remove(f"{end}/{name}.fa.sort")
#open final output file for appending 
   with open(f"{end}/{name}.fa.sort", "a") as b :
#loop through all used counter numbers, find corresponding file in "trash", append to name.sort.fa for final output
    for i in range(0,counter) :
     readfile = open(f"{trash}/{name}_{i}.fa.sort", "r")
     shutil.copyfileobj(readfile, b)
#output name.fa.sort now contains all sequences, sorted for redundancy within individual species 
   break
  else :
   print('\n'+threshold+' not a valid input, threshold must be and integer or float value, try again')
 with open(f"{end}/{name}.fa.sort", "r") as redun :
  reduninfo = {"Accession" : [], "Species" : [], "Protein_name" : []}
#loop through non-redundant sequence file and append accessions, species and protein name to dict as before
  for seq in redun:
   if re.search('>', seq) : 
   #seqinfo["Full_header"].append(seq) #no longer used
    accession = re.search('(?<=>)\S*', seq) 
    reduninfo["Accession"].append(accession.group())
    speciesname = re.search('(?<=\[).*?(?=\])',seq) 
    reduninfo["Species"].append(speciesname.group())
    protname = re.search('(?<=\s).*?(?=\[)', seq)
    reduninfo["Protein_name"].append(protname.group())
#get new count as list of accessions found
 red_count = str(len(reduninfo["Accession"]))
#calculate number of sequences dropped
 drops =  int(count) - int(red_count)
#show user
 print('\n'+str(drops)+' sequence(s) dropped after '+threshold+'% species-specific similarity restriction')
 while True : 
#prompt continue with current dataset
  proceed = input('\nContinue with current dataset of '+red_count+' sequences? : ')
  if proceed.lower() in ['yes','no','y','n'] :
#if yes, delete "trash" and call next function (align())
   if proceed.lower() in ['yes','y'] :
    shutil.rmtree(trash)
    return align()
#if no, clear whole programme directory and try new inputs
   elif proceed.lower() in ['no','n'] :
    for f in os.listdir(end) :
      os.remove(f)
    return helpsearch()
   else :
    print('\n"'+proceed+'" not a valid entry, please try again....\n\ntype "yes" to proceed or "no" to input new search terms')
#
def align() :
#quick multiple alignment of whole dataset
#Ask for iterations in later alignment now so that both can run without input interuption
 while True :  
#default iteraton of 2
  align.itera = input('\nInput desired number of iterations for final alignment : ') or "2"
#single digit integer accepted without problem
  if re.search('\d', align.itera) :
   break
#more than one digit may not be a good idea
  elif re.search('\d(\d)+', align.itera) :
   while True :
#even with 250 sequences >10 iterations may take time, and not worth the output, user confirm 
    sure = input('\nIterations of >10 may take a while to complete, are you sure you wish to continue? : ')
    if sure.lower() in ['yes','no','y','n'] :
     if sure.lower() in ['yes','y'] :
      break
     if sure.lower() in ['no','n'] :
      return helpsearch()
    else :
     print('\n"'+itera+'" not a valid entry, input must be an integer')
 print('\nClustering sequences......') 
 #clustalo = subprocess.Popen(['clustalo', '-v', '-i', f"{pathname}.fa.sort", '-o', f"{pathname}.clstl", '--force'])
 #clustalo.wait()
#call BLAST function with non-redundand sequence file and produced alignment file as arguments
 return BLAST(f"{pathname}.fa.sort", f"{pathname}.clstl")
#
def BLAST(seq, alig) :
 blastout = open(f"{pathname}.BLAST", "w")
 print('\nMaking consensus.....')
 blastdb = f"{pathname}db"
 if not os.path.exists(blastdb) :
  os.mkdir(blastdb)
#make consensus from above clustalo alignment
 cons = subprocess.Popen(['cons', '-auto', '-sequence', alig, '-outseq', f"{pathname}.cons"])
 cons.wait()
#make blast database from non-redundant sequence file 
 makedb =subprocess.Popen(['makeblastdb', '-in', seq, '-title', name, '-dbtype', 'prot', '-out', f"{name}.db"], cwd=blastdb)
 makedb.wait()
#Run blastp on db with consensus sequence as query
 blastp =subprocess.Popen(['blastp','-outfmt', '7', '-evalue', '50', '-db', f"{name}.db", '-query', f"{pathname}.cons"], cwd=blastdb, stdout=blastout)
 blastp.wait()
#Read blastp output into dataframe
#comment=# : removing comment lines, column headers copied from blastp output
#sort dataframe by e values and get 250 most similar (smallest evalue)
 df = pd.read_table(blastout.name, sep='\t', comment='#', header=None, names=['query acc.', 'subject acc.', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']).sort_values(axis=0, by='evalue')
 evalue = df[['subject acc.', 'evalue', 'bit score', '% identity']].sort_values(axis=0, by='evalue').head(250)
 print(evalue)
 #bitscore = df[['subject acc.', 'bit score']].sort_values(axis=0, by='bit score', ascending=False).nlargest(250, 'bit score')
 #identity = df[['subject acc.', '% identity']].sort_values(axis=0, by='% identity', ascending=False).nlargest(250, '% identity')
#write dataframe fo file for later 
 with open(f"{pathname}.fa.top.data", "w+") as b :
  df.to_csv(sep='\t', path_or_buf=b)
#write corresponding accessions for 250 lowest e values to "temp" file, use pullseq with temp as reference to pull 250 most similar files from sorted sequences, output to name.fa.top, then remove temp
 with open(f"{pathname}.fa.top", "w+") as top :
  temp = open(f"{end}/temp_acc", "w")
  temp.write('\n'.join(evalue['subject acc.'].values))
  temp.close()
  pullseq = subprocess.check_output(['/localdisk/data/BPSM/Assignment2/pullseq', '-i', f"{pathname}.fa.sort", '-n', f"{end}/temp_acc"])
  top.write(pullseq.decode('utf-8'))
 os.remove(temp.name)
 return mult()
#
def mult() : 
#Running clustalo again to align 250 most similar, this time ask for iterations becuase 250 sequences is not so demanding 
 print('\nClustering sequences......') 
#Run clustalo 
 #clustalo = subprocess.Popen(['clustalo', '-v', f"--iter={align.itera}", '-i', f"{pathname}.fa.top", '-o', f"{pathname}.clstl.top", '--force'])
 #clustalo.wait()
#plotcon for conservation graph
 plotcon = subprocess.Popen(['plotcon', '-auto', '-sequences', f"{pathname}.clstl.top", '-graph', 'cps', '-gtitle', f"conservation_of)_{protname}_in_the_{taxonomy}_taxon"]) 
 plotcon = subprocess.Popen(['plotcon', '-auto', '-sequences', f"{pathname}.clstl.top", '-graph', 'cps', '-gtitle', f"conservation_of)_{protname}_in_the_{taxonomy}_taxon", '-goutfile', f"{pathname}", ]) 
 return PROSITE(f"{pathname}.clstl.top")
#
def PROSITE(seq) :
#search alignment for any motifs
 patmatmotifs = subprocess.Popen(['patmatmotifs', '-auto', '-full', '-sequence', seq, '-outfile', f"{pathname}.PROSITE"])
 #iqtree = subprocess.Popen(['/localdisk/software/iqtree-1.6.6-Linux/bin/iqtree'])
 return 

#RUN
spellcheck(taxonomy, protname)
sys.exit()

