import os, Fetchutil, random, psutil, subprocess
bayesmodels=['poisson','jtt','mtrev','mtmam','wag','rtrev','cprev','vt','blosum','dayhoff']

def clusttofasta(clust,out):
    "converts a clustal alignment file to a fasta file out"
    hand=open(clust)

    seqs={}
    hand.readline()
    hand.readline()
    hand.readline()
    while hand:
        line=hand.readline()
        #print line
        if line.split()==[] and len(line)>0:
            continue
        if line=='':
            break
        elif line.split()[0][0] in [':','.','*']:
            #print line
            hand.readline()
            continue
        lin=line.split()
        #print line
        try:
            seqs[lin[0]]+=lin[1]
        except KeyError:
            seqs.update({lin[0]:lin[1]})
    for j in seqs.keys():
        open(out,'a').write('>'+j+'\n'+seqs[j]+'\n\n')


def dicextract(fil):
    "Extracts a dict object from a given file"
    p=open(fil)
    dic={}
    for i in p:
        dicobjects=i.strip()[1:-1]
        ar=dicobjects.split(',')
        #print ar
        for pair in ar:
            pairs=pair.split(': ')
            #print pairs
            dic.update({str(pairs[0]):str(pairs[1])})
    return dic

def rename(fas):
    "Takes in a fasta file fas that has sequences and shortens the names, assigning to a new file"
    infile=open(fas)
    out=fas.split('.')[0]
    outfile=open(out,'w')
    orgs={}
    while infile:
        lin=infile.readline()
        if lin=='':
            break
        elif lin.startswith('>'):
            outfile.write('>')
            linstr=lin[1:].split()
            org=linstr[0][0]+linstr[1][0:3]
            if org in orgs.keys():
                org+=repr(orgs[org]+1)
            else:
                orgs.update({org:1})
            outfile.write(org+'\n')
        else:
            outfile.write(lin)

def modseqs(fil):
    "Modifies the file of aligned sequences to remove gapped positions"
    arr={}
    filin=open(fil)
    while filin:
        title=filin.readline()
        if title.startswith('matrix'):
            while filin:
                title=filin.readline()
                if title.startswith('\n'):
                    continue
                elif title.startswith(';'):
                    break
                spl=title.split()
                try:
                    arr[spl[0]]+=spl[1].strip().split('-')
                except KeyError:
                    arr.update({spl[0]:spl[1].strip().split('-')})
            break
    for k in arr:
        seq=''
        for i in arr[k]:
            seq+=i
        arr[k]=seq
    return arr

def pamlseqfas(seqs,inpath):
    "Converts aligned FASTA file seqs to a ProML input file, inpath, Plus 1 space"
    seqstr=''
    seqsin=open(seqs)
    seqstr+=seqsin.readline()[1:].strip()
    count=1
    length=0
    for i in range(11-len(seqstr)):
        seqstr+=' '
    while seqsin:
        line=seqsin.readline()
        if line.startswith('>'):
            count+=1
            seqstr+='\n'+line[1:].strip()
            for i in range(11-len(line[1:].strip())):
                seqstr+=' '
        else:
            seqstr+=line.strip()
            if count==1:
                length+=len(line.strip())
        if line=='':
            break
    seqsout=open(inpath,'w')
    seqsout.write(repr(count)+'  '+repr(length)+'\n'+
                 seqstr)
    seqsin.close()
    seqsout.close()

def pamlseqnex(seqs,inpath):
    "Converts aligned Nexus file seqs to a ProML input file, inpath"
    seqsin=open(seqs)
    count=0
    length=0
    while seqsin:
        lin=seqsin.readline()
        if lin.startswith('dimensions'):
            count=int(lin.split()[1][5:])
            length=int(lin.split()[2][6:len(lin.split()[2])-1])
            break
    out=open(inpath,'w')
    
    while seqsin:
        lin=seqsin.readline()
        if lin.startswith('matrix'):
            break
    dicto={}
    while seqsin:
        lin=seqsin.readline()
        #print '\''+lin+'\''

        if lin=='\n':
            continue
        elif lin==';\n':
            break
        name=lin.split()[0]
        seq=lin.split()[1].strip()
        #print name, seq
        try:
            dicto[name]+=seq
        except:
            dicto.update({name:seq})
    if not length== len(dicto.values()[0]):
        length=len(dicto.values()[0])
    out.write(repr(count)+'  '+repr(length)+'\n')
    for k in dicto.keys():
        out.write(k)
        for i in range(30-len(k)):
            out.write(' ')
        out.write(dicto[k]+'\n')
def bayesinNex(infile):
    "Modifies a PRANK alignment file to be compatible with MrBayes"
    lines=open(infile).readlines()

    for i in range(len(lines)):
        if 'trees' in lines[i]:
            lines=lines[:i]
            break
        else:
            lines[i]=lines[i].replace('\'','')
    fil=open(infile,'w')
    for j in lines:
        fil.write(j)
    fil.close()
def bayesfile(infile,model,outfile):
    "Writes a Nexus file for use as a MrBayes batch file"
    if not os.path.exists('Bayes'):
      subprocess.call(['mkdir','Bayes'])
    for k in model.keys():
        extra=k.split('+')
        if extra[0].lower()=='jtt':
            extra[0]='jones'
        elif extra[0].lower()=='blosum62':
            extra[0]='blosum'
        if extra[0].lower() in bayesmodels:
            han=open(outfile,'w')
            han.write('#NEXUS\n'+
                      'begin mrbayes;\n'+
                      '\texe '+infile+';\n'+
                      '\tprset aamodelpr=fixed('+extra[0]+');\n')
            if len(extra)>1:
                if 'I' in extra and 'G' in extra:
                    han.write('\tlset rates=Invgamma;\n')
                    han.write('\tprset shapepr=fixed('+model[k][0]+');\n')
                elif 'I' in extra:
                    han.write('\tlset rates=Propinv;\n')
                elif 'G' in extra:
                    han.write('\tlset rates=Gamma;\n')
                    han.write('\tprset shapepr=fixed('+model[k][0]+');\n')

            han.write('\tmcmc ngen=50000 samplefreq=50 file='+outfile+';\n'+
                      '\tsumt burnin=250;\n'+
                      'end;\n\n')
            han.close()

def boot(infile,outfile,norep):
    "takes in an aligned phylip file and creates a file in phylip format with norep bootstrap datasets"
    hand= open(infile)
    lin = hand.readline()
    noseqs=int(lin.split()[0])
    length=int(lin.strip().split()[1])
    seqs={}
    for i in range(noseqs):
        lin=hand.readline().strip().split()
        seqs.update({lin[0]:lin[1]})
   # print 'seqs'+str(seqs)
    hand2=open(outfile,'w')
    for j in range(norep):
        hand2.write(str(noseqs)+'   '+str(length)+'\n')
        newseqs={}
        for i in range(length):
            randInt=random.randint(0,length-1)
            for j in seqs:
                try:
                    newseqs[j]+=seqs[j][randInt]
                except KeyError:
                    newseqs.update({j:seqs[j][randInt]})
        for k in newseqs:
            hand2.write(k+'    '+newseqs[k]+'\n')
        hand2.write('\n\n')
    hand2.close()

def noseqs(fil):
    "returns number of seqs in given fil fasta file"
    try:
        p=open(fil)
    except:
        return 0
    count=0
    while p:
        lin=p.readline()
        if lin=='':
            break
        elif lin.startswith('>'):
            count+=1
    return count

def splicealign(inalign,outalign):
    "Splices out evolutionarily uncertain areas in inalign, realigns the new sequence"
    align=open(inalign)
    dicto={}
    inalign=inalign.split('.')[0]
    noseq=0
    align_size=0
    while align:
        lin=align.readline()
        if lin.startswith('end'):
            break
        if lin.startswith('dimensions'):
            li=lin.strip().split()
            for i in li:
                if i.startswith('ntax'):
                    noseq=int(i.split('=')[1])
                elif i.startswith('nchar'):
                    align_size=int(i.split('=')[1].split(';')[0])
        if lin.startswith('matrix'):
            while align:
                lin=align.readline()
                if lin=='\n':
                    continue
                elif lin==';\n':
                    break
                else:
                    lo=lin.split()
                    try:
                        dicto[lo[0]]+=lo[1]
                    except KeyError:
                        dicto.update({lo[0]:lo[1]})
    ##sequences have been read in; determining the positions to be removed
    pos={}

    for i in dicto.values():
        for c in range(len(i)):
            if i[c]=='-':
                try:
                    pos[c]+=1
                except:
                    pos.update({c:1})
    #print pos
    #processing pos; eliminating key:value pairs whose values are less than half len(i)
    for k in pos.keys():
        if pos[k]<(0.5*noseq):
            pos.pop(k)
    post=list(pos.keys())
    post.sort()
    #initial condition to stop recursion
    if len(post)<=.05*align_size:
        han=open(outalign,'w')
        han.write('#NEXUS\n'+
                  'begin data;\n'+
                  'dimensions ntax='+repr(noseq)+' nchar='+repr(align_size)+';\n'+
                  'format datatype=protein interleave=no gap=-;\n'+
                  'matrix\n\n')

        for k in dicto.keys():
            han.write(k)
            for i in range((Fetchutil.findlonglen(dicto)+6)-len(k)):
                han.write(' ')
            han.write(dicto[k]+'\n')
        han.write(';\nend;\n')
        han.close()
        return
    #adds one residue on either side of the residues indicated in post
    for p in range(len(post)):
        if post[p]>0 and post[p]<align_size-1:
            bef=post[p]-1
            af=post[p]+1
            if bef in post:
                post.append(bef)
            if af in post:
                post.append(af)
    post.sort()
    post.reverse()

    print align_size, len(post)

    #processes dicto; removes residues contained in post
    for i in dicto.keys():
        val = dicto[i]
        new_val=''
        for j in range(len(post)-1):
            new_val+=val[post[j+1]+1:post[j]]
        for k in range(align_size-len(new_val)):
            new_val+='-'
        dicto.update({i:new_val})
    #creates the nex file containing the newly spliced data
    align_size=len(dicto.values()[0])
    han=open(inalign+'.edit','w')
    han.write('#NEXUS\n'+
              'begin data;\n'+
              'dimensions ntax='+repr(len(dicto.keys()))+' nchar='+repr(align_size)+';\n'+
              'format datatype=protein interleave=no gap=-;\n'+
              'matrix\n\n')

    for k in dicto.keys():
        han.write(k)
        for i in range((Fetchutil.findlonglen(dicto)+6)-len(k)):
            han.write(' ')
        han.write(dicto[k])
        #print dicto[k]
        han.write('\n')
    han.write(';\nend;\n')
    han.close()
    print '.edit written. aligning'
    os.system('clustalw -align -infile='+inalign+'.edit -outfile='+inalign+'.ed -output=nexus -quiet')
    print '.ed written splice aligning'
    #splicealign(inalign+'.ed',outalign)

def bestmod(infile):
    "Takes in alignment file runs protTest, and extracts best model(s)"
    if not os.path.exists('Prot'):
      os.mkdir('Prot')
    procs=int(round(psutil.NUM_CPUS/2.0))
    #print procs
    out=infile.split('/')[1].split('.')[0]
    protCmd='java -jar prottest/prottest-3.4.jar -i '+infile+' -o Prot/'+out+'.pro -all-distributions -all -S 1 -threads '+repr(procs)+' -BIC'
    #print protCmd
    os.system(protCmd)
    prot_hand=open('Prot/'+out+'.pro')
    models={}
    ret={}
#reading and processing prottest output
    while prot_hand:
        lin=prot_hand.readline()
        if lin.startswith('Model.'):
            lsplit=lin.split()
            mod= lsplit[2]
            modmod=mod.split('+')
            para=['0','0'] #index 0 is G index 1 is I
            if len(modmod)>1:
                if 'G' in modmod:
                    while 1:
                        lin=prot_hand.readline()
                        #print lin,'g'
                        if lin.split()[0]=='gamma':
                            para[0]= lin.split()[6]
                            break
                if 'I' in modmod:
                    while 1:
                        lin=prot_hand.readline()
                        #print lin.split(),'i'
                        if lin.split()[0]=='proportion':
                            para[1]= lin.split()[5]
                            break
            models.update({mod:para})
            #print models
            continue
        elif lin.startswith('Best model'):
            lsplit=lin.split()
            mod = lsplit[5]
            ret.update({mod:models[mod]})
            #print ret
            if mod.lower().split('+')[0] not in bayesmodels:
                while 1:
                    lin=prot_hand.readline()
                    #print lin,1
                    if lin.endswith('-\n'):
                        while 2:
                            lin=prot_hand.readline().split()
                            #print lin,2
                            if lin[0].split('+')[0].lower() in bayesmodels:
                                ret.update({lin[0]:models[lin[0]]})
                                #print ret
                                break
                        break
        
       # elif lin.startswith('*') and :
        #    print "Saving Tree"
        #    tree=''
         #   while 9:
         #       lin=prot_hand.readline().strip()
          #      tree+=lin
          #      if lin.endswith(';'):
           #         break
           # open('Prot/'+out+'.tre','w').write(tree)
            return ret
