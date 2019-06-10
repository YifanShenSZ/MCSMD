''' User input '''
SMDOrder=2

''' Import library '''
import sys# Standard library
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anm
#sys.path.append('/home-4/yshen57@jhu.edu/Library/Python-Library')# My library path on MARCC
sys.path.append('C:\\Python-Library')# My library path on my Lenovo-E540
#sys.path.append('C:\\Users\\56402\\OneDrive\\Research\\Source\\Python')# My library path on my surface pro
from PythonLibrary import *

''' Auxiliary routine '''
def Get_Array(source):
    with open(source,'r') as f:
        data=f.readlines()
    l=len(data)
    a=numpy.empty(l)
    for i in range(l):
        a[i]=float(data[i].strip())
    return a

def Get_SMD(source,lt):
    with open(source,'r') as f:
        data=f.readlines()
    nSMD=int(len(data)/lt)
    SMD=numpy.empty((nSMD,lt))
    index=0
    for i in range(lt):
        for j in range(nSMD):
            SMD[j,i]=float(data[index].strip())
            index=index+1
    return SMD

''' Do the job '''
t=Get_Array('t.SMD')
SMD=Get_SMD('SMD.SMD',t.shape[0])
with open('SMD.txt','w') as f:
    print('t',end='\t',file=f)
    for i in range(1,SMDOrder+1):
        for j in range(i+1):
            print('x'+str(i-j)+'p'+str(j),end='\t',file=f)
    print(file=f)
    for i in range(t.shape[0]):
        print(t[i],end='\t',file=f)
        for j in range(SMD.shape[0]):
            print(SMD[j,i],end='\t',file=f)
        print(file=f)

Visualization.Plot2D(t,SMD[0,:],title='<q>-t')