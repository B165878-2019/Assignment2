#!/usr/bin/python3
import sys
import subprocess
import os
import xml.etree.ElementTree as ET
##### Entrez utilities
taxonomy = input('Taxon : ')
protname = input('Protein family : ')

#check spelling
def spellcheck():
 espell = subprocess.check_output(['espell', '-db', 'taxonomy', '-query', taxonomy])
 root = ET.fromstring(espell)
 query = root[1].text
 correction = root[2].text
#
 if not correction :
  term = query
 else : 
  print('\nPossible spelling error found\n\nQuery : '+query+'\nCorrection : '+correction)
  while True:
   c = input('\nwould you like to use the corrected query? (yes/no)')
   if c.lower() in ['yes','no','y','n'] :
    if c.lower() in ['yes','y'] :
     term = correction
     break
    if c.lower() in ['no','n'] :
     term = query
     break
   else :
    print('\n"'+c+'" not a valid entry, please try again....\n\ntype "yes" to use the spell corrected query or type "no" to use the original query')
#Check hits
def checkhits()
 print('\nNow checking hits for search....')
 query_term = (f"{term}[ORGN] AND {protname} NOT PARTIAL")
 checktax = subprocess.check_output(['esearch', '-db', 'protein', '-query', query_term])
 root = ET.fromstring(checktax)
 hitcount = root[3].text
 print('\nFound '+hitcount+' hits')
 if int(hitcount) >= 10000 :
  print('\nData sets larger than 10000 sequences are not recommended!')
  while True :
   big = input('\nDo you still wish to continue? :')
   if big.lower() in ['yes','no','y','n']:
    if big in ['yes', 'y']:
     print('okay tough guy....')
     force = True
     break
    if big in ['no','n']:
     print('\nOkay, exiting.....') 
     sys.exit() #Refinement sub-function instead?? (if i have time)
   else:
    print('\n'+big+' not a valid option, type "yes" to use oversised sample or "no" to exit the programme and think about what you are actually searching....')
#fetching FASTA
with open(f"{term}.fa", "w+") as fetchout:
 query_term = (f"{term}[ORGN] AND {protname}[PROT] NOT PARTIAL")
 esearch = subprocess.Popen(
	['esearch', '-db', 'protein', '-query', query_term],
	stdout=subprocess.PIPE)
 efetch_fasta = subprocess.Popen(
	['efetch', '-db', 'protein', '-format', 'fasta'],
	stdin=esearch.stdout, stdout=fetchout)








