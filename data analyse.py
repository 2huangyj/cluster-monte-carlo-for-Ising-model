import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt 
T = pd.read_csv("result.csv", usecols=[0]) # I want the first col
M = pd.read_csv("result.csv", usecols=[1]) # I want the 3rd col
C = pd.read_csv("result.csv", usecols=[2])
chi = pd.read_csv("result.csv", usecols=[3])
ksai = pd.read_csv("result.csv", usecols=[4])
'''
T =( pd.read_csv("result.csv", usecols=[0])).dropna(thresh = 1) 
M = (pd.read_csv("result.csv", usecols=[1])).dropna(thresh = 1) 
C = pd.read_csv("result.csv", usecols=[2]).dropna(thresh = 1)
chi = pd.read_csv("result.csv", usecols=[3]).dropna(thresh = 1)
ksai = pd.read_csv("result.csv", usecols=[4]).dropna(thresh = 1)'''

T = T.as_matrix().reshape(T.shape[0],)
M = M.as_matrix().reshape(M.shape[0],)
C = C.as_matrix().reshape(C.shape[0],)
chi = chi.as_matrix().reshape(chi.shape[0],)
ksai = ksai.as_matrix().reshape(ksai.shape[0],)


def func(x,A,x0,c_e,const):
	return (A*(np.abs(x-x0))**c_e+const)
'''
def func(x,A,x0,c_e,Const):
	return (A*(x-x0)**c_e+Const)*0.5*(1+np.sign(x0-x))'''

#popt,pcov = curve_fit(func,T[len(T)-12:len(T)],ksai[len(T)-12:len(T)],p0 = [1,2.3,-0.2,-0.55])

#c_e_ksai = popt[2]





l = func(T,0.08,4.25,-0.6,-0.05)
#l = func(T,popt[0],popt[1],popt[2],popt[3])
fig = plt.figure()
plt.scatter(T,C,color = 'r')
plt.plot(T,l,color = 'b')
plt.show()

'''popt,pcov = curve_fit(func,T[len(T)-9:len(T)],C[len(T)-9:len(T)],p0 = [1,2.32,0,0])
c_e_C = popt[2]
l = func(T,popt[0],popt[1],popt[2],popt[3])
fig = plt.figure()
plt.plot(T,l)
plt.show()'''

'''popt,pcov = curve_fit(func,T[len(T)-20:len(T)],chi[len(T)-20:len(T)])
c_e_chi = popt[2]
l = func(T,popt[0],popt[1],popt[2])
fig = plt.figure()
plt.plot(T,l)
plt.show()'''

'''popt,pcov = curve_fit(func,T[len(T)-13:len(T)],ksai[len(T)-13:len(T)])
c_e_ksai = popt[2]
l = func(T,popt[0],popt[1],popt[2])
fig = plt.figure()
plt.plot(T,l)
plt.show()'''

#print c_e_C,c_e_chi,c_e_ksai
#print c_e_ksai
