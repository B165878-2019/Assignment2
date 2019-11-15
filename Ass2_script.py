#!/usr/bin/python3
import sys
import subprocess
import os
import xml.etree.ElementTree as ET
##### Entrez utilities
taxonomy = input('\nTaxon name : ')
protname = input('\nProtein name : ')
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
  checkhits(spellcheck.tax, spellcheck.prot)
 elif taxcorrection and not protcorrection : 
  spellcheck.prot = protquery
  print('\nPossible spelling correction found\n\nQuery : '+taxquery+'\nCorrection : '+taxcorrection)
  while True:
   c = input('\nwould you like to use the corrected query? (yes/no) : ')
   if c.lower() in ['yes','no','y','n'] :
    if c.lower() in ['yes','y'] :
     spellcheck.tax = taxcorrection
     break
    if c.lower() in ['no','n'] :
     spellcheck.tax = taxquery
     break
   else :
    print('\n"'+c+'" not a valid entry, please try again....\n\ntype "yes" to use the spell corrected query or type "no" to use the original query')
  checkhits(spellcheck.tax, spellcheck.prot)
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
  checkhits(spellcheck.tax, spellcheck.prot)
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
  checkhits(spellcheck.tax, spellcheck.prot)
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
     force = True #can be useful later 
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
#Search help 
def helpsearch(): 
 print('\nSorry, the current search term is findng zero matches within the database :\n\n**"'+checkhits.query_term+'"**')
 while True :
  man = input('\nTo manually enter a search query type "yes"\nTo try again with new search terms type "no"\n\nType here : ')
  if man.lower() in ['yes','no','y','n']:
   if man.lower() in ['yes', 'y']:
    print('\n\nInformation on NCBI search field descriptions can be found here :\n"https://www.ncbi.nlm.nih.gov/books/NBK49540/"\nProtein nomenclature guidelines can be found here :\nhttps://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/')
    manqueryi = input('\n\nPlease type the EXACT search input here : ')    
    checkhits('none', 'none', manquery=(f"{manqueryi}"))
    break
   elif man.lower() in ['no','n']:
    helpsearch.taxonomy = input('\nTaxon name : ')
    helpsearch.protname = input('\nProtein name : ')
    try: 
     spellcheck(helpsearch.taxonomy, helpsearch.protname)
    except: 
     print('No spelling errors found')
    else:
     fetchFASTA(checkhits.query_term)
  else :
   print('\n"'+retry+'" not a valid option, please try again...')
#
#fetching FASTA
def fetchFASTA(query_term):
 fetchout = open(f"{spellcheck.tax}.fa", "w+")
 print('\nFetchng '+checkhits.hitcount+' results.....')
 esearch = subprocess.Popen(['esearch', '-db', 'protein', '-query', query_term],stdout=subprocess.PIPE)
 efetch = subprocess.Popen(['efetch', '-db', 'protein', '-format', 	'fasta'],stdin=esearch.stdout, stdout=fetchout)
#Making consessus
def cons(infile) :
 subprocess.Popen(['clustalo', '-i', f"{in}.fa", '-o', f"{in}.clstl"]) 
#
#def BLAST() :
#makeblastdb -in () -title () -dbtype prot #make db directory because this outputs 3 files
#blastp -db () -query ()
#RUN
spellcheck(taxonomy, protname)


