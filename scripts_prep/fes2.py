import numpy as np
import matplotlib.pyplot as plt
import itertools
import os
import pandas as pd
import argparse, sys
parser = argparse.ArgumentParser()
parser.add_argument("--CV_col", help="CV column in the COLVAR file (starting at 0)",required=True,type=str)
parser.add_argument("--CV_name", help="CVname",required=True,type=str)

args = parser.parse_args()

CVcol =args.CV_col
CVname =args.CV_name
print(int(CVcol))
eq_steps=1

time=np.loadtxt("FULLBIAS")[:, 0]
bias=np.loadtxt("FULLBIAS")[:, 1]
CV=np.loadtxt("COLVAR")[:, int(CVcol)]
KBT = 2.49
weights = np.exp(bias/ KBT)
weights /= weights.sum()

#value, bins=np.histogram(CVname[eq_steps:-1], bins=30, weights=weights[eq_steps:-1], density=True)
value1, bins1=np.histogram(CV[eq_steps:int(len(CV)/5)], bins=30, weights=weights[eq_steps:int(len(CV)/5)], density=True)
value2, bins2=np.histogram(CV[eq_steps:int(2*len(CV)/5)], bins=30, weights=weights[eq_steps:int(2*len(CV)/5)], density=True)
value3, bins3=np.histogram(CV[eq_steps:int(3*len(CV)/5)], bins=30, weights=weights[eq_steps:int(3*len(CV)/5)], density=True)
value4, bins4=np.histogram(CV[eq_steps:int(4*len(CV)/5)], bins=30, weights=weights[eq_steps:int(4*len(CV)/5)], density=True)
value5, bins5=np.histogram(CV[eq_steps:int(5*len(CV)/5)], bins=30, weights=weights[eq_steps:int(5*len(CV)/5)], density=True)

logvalue1=-np.log(value1)
logvalue1=logvalue1-min(logvalue1)
logvalue2=-np.log(value2)
logvalue2=logvalue2-min(logvalue2)
logvalue3=-np.log(value3)
logvalue3=logvalue3-min(logvalue3)
logvalue4=-np.log(value4)
logvalue4=logvalue4-min(logvalue4)

logvalue5=-np.log(value5)
logvalue5=logvalue5-min(logvalue5)
plt.tick_params(labelsize=24)
plt.plot(bins1[:-1],logvalue1,linewidth=3,color='red')
plt.plot(bins2[:-1],logvalue2,linewidth=3,color='green')
plt.plot(bins3[:-1],logvalue3,linewidth=3,color='blue')
plt.plot(bins4[:-1],logvalue4,linewidth=3,color='pink')
plt.plot(bins5[:-1],logvalue5,linewidth=3,color='yellow')
plt.legend( ['t1', 't2','t3','t4','t5'],loc='lower right',ncol=2,prop={'size': 11})
figure = plt.gcf()
figure.set_size_inches(6, 7)
plt.xlabel(str(CVname),fontsize=24 )
plt.ylabel('F (kJ/mol)',fontsize=24 )
plt.savefig('FES'+str(CVname)+'.png',dpi=400,transparent=True, bbox_inches='tight')
