#import Main, multiprocessing
import psutil, subprocess
import sys, shutil

if 'clear' in sys.argv:
  try:
    shutil.rmtree('aligns')
    shutil.rmtree('Data')
    shutil.rmtree('Orthos')
  except:
    pass
  print sys.exit(0)

jobs=open('jobs.txt').readlines()

def launch(query,label,domain):
  "Wrapper method to Main.py execution."
  subprocess.call(['python','Main.py',query,label,domain])

if __name__=='__main__':
  for i in jobs:
    spl=i.split()
    (name,quer)=(spl[0],spl[1])
    launch(quer,name,'all')
#    p = multiprocessing.Process(name=name, target=launch, args=(quer,name,'all'))
#    p.start()
#    p.join()
