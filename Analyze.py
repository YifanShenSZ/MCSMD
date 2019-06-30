''' User input '''
SMDOutputOrder=2

''' Import library '''
import sys# Standard library
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anm
#sys.path.append('/home-4/yshen57@jhu.edu/Library/Python-Library')# My library path on MARCC
#sys.path.append('C:\\Python-Library')# My library path on my Lenovo-E540
sys.path.append('C:\\Users\\56402\\OneDrive\\Research\\Source\\Python')# My library path on my surface pro
from PythonLibrary import *

''' Auxiliary routine '''
def Get_input(source):
    with open(source,'r') as f: data=f.readlines()
    jobtype=str(data[3].strip())
    return jobtype

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

def Get_Wigner(source,lq,lp,lt):
    with open(source,'r') as f: data=f.readlines()
    wigner=numpy.empty((lq,lp,lt))
    index=0
    for i in range(lt):
        for j in range(lp):
            for k in range(lq):
                wigner[k,j,i]=float(data[index].strip())
                index=index+1
    return wigner

''' Do the job '''
jobtype=Get_input('SMD.in')
if(jobtype=='SMD'):
    t=General.Get_Array('t.out')
    SMD=Get_SMD('SMD.out',t.shape[0])
    with open('SMD.txt','w') as f:
        print('t/a.u.','<q>/a.u.','<p>/a.u.','sigmaq/a.u.','correlation','sigmap/a.u.',sep='\t',end='\t',file=f)
        for i in range(3,SMDOutputOrder+1):
            for j in range(i+1):
                print('<q'+str(i-j)+'p'+str(j)+'>',end='\t',file=f)
        print(file=f)
        for i in range(t.shape[0]):
            print(t[i],end='\t',file=f)
            for j in range(int(SMDOutputOrder*(SMDOutputOrder+3)/2)):
                print(SMD[j,i],end='\t',file=f)
            print(file=f)
    Visualization.Plot2D(t,SMD[0,:],title='<q>-t')
elif(jobtype=='Wigner'):
    q=General.Get_Array('Wigner_q.out'); p=General.Get_Array('Wigner_p.out'); t=General.Get_Array('t.out')
    wigner=Get_Wigner('Wigner.out',q.shape[0],p.shape[0],t.shape[0])
    [qtemp,ptemp]=numpy.meshgrid(q,p); Qtemp=qtemp.ravel(); Ptemp=ptemp.ravel()
    Q=numpy.empty((Qtemp.shape[0],t.shape[0])); P=numpy.empty((Ptemp.shape[0],t.shape[0]))
    for i in range(t.shape[0]):
        Q[:,i]=Qtemp; P[:,i]=Ptemp
    WIGNER=numpy.empty((int(q.shape[0]*p.shape[0]),t.shape[0]))
    for j in range(t.shape[0]):
        index=0
        for i in range(p.shape[0]):
            WIGNER[index:index+q.shape[0],j]=wigner[:,i,j]
            index=index+q.shape[0]
    Visualization.Animate3DSurface(t,Q,P,WIGNER,\
        title='Wigner distribution',xlabel='q [a.u.]',ylabel='p [a.u.]',zlabel='density',\
        colormap='seismic',save=True,FileName='Wigner',show=False)
else: print('Unkown job type?')