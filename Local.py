import sys, os, subprocess

def generatescript(org):
  "generates a fetch script for all protein seqs of org organism from NCBI"
  DOL='$$'
  AMP='&&&'
  bi=org.split()
  four=bi[0][0]+bi[1][:3]
  four=four.lower()
  template=open('fetchtemp').readlines()
  if not os.path.exists('Proteomes'):
    subprocess.call(['mkdir','Proteomes'])
  genscript=open('Proteomes/fetch'+four,'w')

  for lin in template:
    spldol=lin.split(DOL)
    splamp=lin.split(AMP)
    if len(spldol)==3:
      genscript.write(spldol[0]+org+spldol[2])
      continue
    elif len(splamp)==3:
      genscript.write(splamp[0]+'Proteomes/'+four+splamp[2])
      continue
    else:
      genscript.write(lin)
  genscript.close()
  return four

def fetchfasta(org):
  "runs a fetch script of 4-letter org organism"
  subprocess.call(['perl', 'Proteomes/fetch'+org])
  return 'Proteomes/'+org+'.aa'


