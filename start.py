import Main
import multiprocessing, psutil

jobs=open('jobs.txt').readlines()

if __name__=='__main__':
  for i in jobs:
    spl=i.split()
    (name,quer)=(spl[0],spl[1])
    p = multiprocessing.Process(name=name, target=Main, args=(quer,name,'all'))
    p.start()
