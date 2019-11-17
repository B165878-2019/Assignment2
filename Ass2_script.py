#!/usr/bin/python3
import sys
import subprocess
import os
import xml.etree.ElementTree as ET
import re
import pandas as pd
##### Entrez utilities
cwd = os.getcwd()
end = (cwd+'/B165878_Programme')
if not os.path.exists(end) :
 os.mkdir('B165878_Programme')
#
taxonomy = input('\nTaxon name : ')
protname = input('\nProtein name : ')
#
# Check spelling
#
#check spelling
def spellcheck(taxon, protname):
 espelltax = subprocess.check_output(['espell', '-db', 'taxonomy', '-query', taxon])
 taxroot = ET.fromstring(espelltax)
 taxquery= taxroot[1].text
 taxcorrection = taxroot[2].text
 espellprot = subprocess.check_output(['espell', '-db', 'protein', '-query', protname])
 protroot = ET.fromstring(espellprot)
 protquery = protroot[1].text
 protcorrection = protroot[2].text
#
 if not taxcorrection and not protcorrection :
  spellcheck.tax = taxquery
  spellcheck.prot = protquery  
  return checkhits(spellcheck.tax, spellcheck.prot)
 elif taxcorrection and not protcorrection : 
  spellcheck.prot = protquery
  print('\nPossible spelling correction found\n\nQuery : '+taxquery+'\nCorrection : '+taxcorrection)
  while True:
   c = input('\nWould you like to use the corrected query? (yes/no) : ')
   if c.lower() in ['yes','no','y','n'] :
    if c.lower() in ['yes','y'] :
     spellcheck.tax = taxcorrection
     break
    if c.lower() in ['no','n'] :
     spellcheck.tax = taxquery
     break
   else :
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
 return checkhits(spellcheck.tax, spellcheck.prot)
#Check hits
def checkhits(taxon, protname, manquery=''):
 print('\nNow checking hits for search....')
 if manquery :
  checkhits.query_term = manquery
 else :
  checkhits.query_term = (f"{taxon}[ORGN] AND {protname} NOT PARTIAL")
 esearch = subprocess.check_output(['esearch', '-db', 'protein', '-query', checkhits.query_term])
 root = ET.fromstring(esearch)
 checkhits.hitcount = root[3].text
 print('\nFound '+checkhits.hitcount+' hits')
 if int(checkhits.hitcount) >= 10000 :
  print('\nData sets larger than 10000 sequences are not recommended!')
  while True :
   big = input('\nDo you still wish to continue? :')
   if big.lower() in ['yes','no','y','n']:
    if big in ['yes', 'y']:
     print('okay tough guy....')
     break
    elif big in ['no','n']:
     print('\nOkay, exiting.....') 
     sys.exit() 
   else:
    print('\n'+big+' not a valid entry, type "yes" to use oversised sample or "no" to exit the programme and think about what you are actually searching....')
 elif int(checkhits.hitcount) == 0 : 
  while True :
   retry = input('\nerror, zero hits, type "yes" to go to search help or "no" to exit the programme : ')
   if retry.lower() in ['yes','no','y','n']:
    if retry.lower() in ['yes', 'y']:
     helpsearch()
     break
    elif retry.lower() in ['no','n']:
     sys.exit()
   else :
    print('\n"'+retry+'" not a valid option, please try again...')
 else :
  fetchFASTA(checkhits.query_term)
 return
#Search help 
def helpsearch(): 
 print('\nSorry, the current search term is findng zero matches within the database :\n\n**"'+checkhits.query_term+'"**')
 while True :
  man = input('\nTo manually enter a search query type "yes", to try again with new search terms type "no"\n\nType here : ')
  if man.lower() in ['yes','no','y','n']:
   if man.lower() in ['yes', 'y']:
    print('\n\nInformation on NCBI search field descriptions can be found here :\n"https://www.ncbi.nlm.nih.gov/books/NBK49540/"\nProtein nomenclature guidelines can be found here :\nhttps://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/')
    manqueryi = input('\n\nPlease type the EXACT search input here : ')    
    checkhits('none', 'none', manquery=(f"{manqueryi}"))
    break
   elif man.lower() in ['no','n']:
    helpsearch.taxonomy = input('\nTaxon name : ')
    helpsearch.protname = input('\nProtein name : ')
    break
  else :
   print('\n"'+man+'" not a valid option, please try again...')
 return spellcheck(helpsearch.taxonomy, helpsearch.protname)
#
#fetching FASTA
def fetchFASTA(query_term):
 print('\nNumber of hits within acceptable range, moving to fetching stage...')
 while True :
  global name
  name = input('\nPlease specify a filename, this will be the prefix to all concurrent files : ')   
  if name.isalnum()  :   
   with open(f"{end}/{name}.fa", "w") as fetchout :
    fetchr = open(f"{end}/{name}.fa", "r")
    print('\nFetchng '+checkhits.hitcount+' results.....')
    esearch = subprocess.Popen(['esearch', '-db', 'protein', '-query', query_term],stdout=subprocess.PIPE)
    efetch = subprocess.Popen(['efetch', '-db', 'protein', '-format', 	'fasta'],stdin=esearch.stdout, stdout=fetchout)
    fetchout.close()
    return print('butthole') #sort(f"{end}/{name}.fa")
  else :
   print('\nSorry only alphanumeric characters are permitted, please try again')
#Sorting FASTA
def sort(inf):
 aves = open(f"{inf}", 'r')
 seqinfo = {"Accession" :[], "Protein_name" :[], "Species" :[]}
 for seq in aves:
  if re.search('>', seq) : 
   seqinfo["Full_header"].append(seq)
   accession = re.search('(?<=>)\S*', seq) 
   seqinfo["Accession"].append(accession.group())
   speciesname = re.search('(?<=\[).*?(?=\])',seq) 
   seqinfo["Species"].append(speciesname.group())
   protname = re.search('(?<=\s).*?(?=\[)', seq)
   seqinfo["Protein_name"].append(protname.group())   
 df = pd.DataFrame(seqinfo)
 print(df)
 namelist = sorted(set(seqinfo["Species"]))
 indv = str(len(namelist))
 print('\n'+indv+' individual species found')
 while True :
  threshold = input('\nPlease input the similarity threshold for sorting redundancy : ')
  if re.search("\d\d\.?\d?\d?", threshold) :
   trash = f"{end}/{name}_FASTA"
   if not os.path.exists(trash):
    os.mkdir(trash)
   counter = 0
   for item in namelist :
    temp = open(f"{trash}/temp.txt", 'w+')
    tempr = open(f"{trash}/temp.txt", 'r')
    temp2 = open(f"{trash}/temp2.txt", 'w+')
    temp2r = open(f"{trash}/temp2.txt", 'r') 
    sequences = f"{end}/{inf}" 
    outfile = f"{trash}/{name}_{counter}.fa.sort"
    outfiler = open(outfile, "r")
    individual = df[df['Species'] == item]['Accession']
    b = "\n"
    i = b.join(list(individual.values))
    temp.write(i)
    temp.close()
    targets = f"{trash}/temp.txt"
    targets2 = f"{trash}/temp2.txt"
    pullseq = subprocess.check_output(['/localdisk/data/BPSM/Assignment2/pullseq', '-i', sequences, '-n', targets])
    temp2.write(pullseq.decode('utf-8'))
  #print(temp2r.read())
    skipredundant = subprocess.Popen(['skipredundant','-auto','-threshold', threshold, '-sequences', targets2, '-outseq', outfile])
    outfiler.close()
    counter += 1
   with open(f"{end}/{name}.fa.sort", "a") as outfile :
    for i in range(0,counter) :
     readfile = open(f"{trash}/{name}_{i}.fa.sort", "r")
     shutil.copyfileobj(readfile, outfile)
   break
  else :
   print('\n'+threshold+' not a valid input, threshold must be and integer or float value, try again')
 shutil.rmtree(trash)
 return

 

#def cons(filename) :
 #clustalo = subprocess.Popen(['clustalo', '-i', f"{filename}.fa", '-o', f"{filename}.clstl"], stdout=subprocess.PIPE) 
 #subprocess.Popen(['cons', '-outseq', f"{infile}.cons"], stdin=clustalo.stdout)
#
#def BLAST() :
#makeblastdb -in () -title () -dbtype prot #make db directory because this outputs 3 files
#blastp -db () -query ()
#RUN
spellcheck(taxonomy, protname)


