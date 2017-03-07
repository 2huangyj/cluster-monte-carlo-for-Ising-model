import numpy as np
import matplotlib.pyplot as plt 
import random
from scipy.optimize import curve_fit
import csv

#cluster Monte Carlo
D = 2
l = 100
t_step = 7
J = 1
T_list = [1.7,1.8,2.0,2.2,2.22,2.24,2.26,2.27,2.28,2.29,2.30,2.32,2.34,2.36,2.38,2.40,2.45,2.50,2.55,2.60,2.70,2.80,3.0]
#1.0,1.25,1.5,1.75,2.0,2.1,2.15,2.20,2.24,2.25,2.26,2.27,2.28,2.29,2.30,2.4,2.5,2.75,3.0,3.25,3.75,4.0,4.25,4.5,4.75
#create a stack
M_list = []
C_list = []
chi_list = []
ksai_list = []

class Stack:
	def __init__(self):
		self.items = []

	def isEmpty(self):
		return self.items == []

	def push(self,item):
		self.items.append(item)

	def pop(self):
		return self.items.pop(0)

def func(x,ksai,A,Const):
	return A*np.exp(-x/ksai)+Const

for T in T_list:
	print 'T:', T
	t = 0
	t_max = 980#100000
	beta = 1/T
	p = 1-np.exp(-2*beta*J)
	grid = 2*np.random.randint(2, size=(l,l))-1

	s = Stack()
	sum_M = 0.0
	sum_E = 0.0
	sum_M_2 = 0.0
	sum_E_2 = 0.0
	correlation = np.zeros(shape = (l,))
	
	#flip the spin according to the algorithm
	while t < t_max:
		flipped = np.zeros(shape=(l,l))
		i = np.random.randint(l - 1)
		j = np.random.randint(l - 1)
		s.push((i,j))
		flipped[i,j] = 1
		
		M = 0.0
		E = 0.0

		while not s.isEmpty():
			x,y = s.pop()
			for (x_,y_) in [((x-1)%l,y),(x,(y-1)%l),((x+1)%l,y),(x,(y+1)%l)]:
				if (grid[x,y] == grid[x_,y_]) and (flipped[x_,y_]==0) and (p>random.random()):
					s.push((x_,y_))
					flipped[x_,y_]=1
			grid[x,y] = -grid[x,y] 
			#inner while loop is over

		if (t+1)%t_step == 0:
			print t
			#calculate M
			M = sum(sum(grid))
			M = M/(l**2*1.0)
			sum_M = sum_M + np.abs(M) 
			sum_M_2 = sum_M_2 + M * M
			#calculate E		
			for i in range(l):
				for j in range(l):
					E = E - J*(grid[i,j]*grid[i,(j+1)%l] + grid[i,j]*grid[(i+1)%l,j])
					
			E = E/(l**2*1.0)
			sum_E = sum_E + np.abs(E)
			sum_E_2 = sum_E_2 + E * E

			#calculate the correlation length
			cross  = np.zeros(shape = (l,))
			for i in range(l):
				grid_piece_k = (np.abs(np.fft.fft(grid[i,:])))**2
				cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))
				grid_piece_k = (np.abs(np.fft.fft(grid[:,i])))**2
				cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))

			cross = np.sum(cross,axis = 0)/(1.0*D*l**D) 
			correlation = correlation + cross
		
		t = t + 1
		#outer while loop is over

	ave_M = sum_M/((t_max//t_step)*1.0)
	ave_M_2 = sum_M_2/((t_max//t_step)*1.0)

	ave_E = sum_E/((t_max//t_step)*1.0)
	ave_E_2 = sum_E_2/((t_max//t_step)*1.0)

	correlation = (correlation/((t_max//t_step)*1.0) - ave_M**2)[0:50]
	popt,pcov = curve_fit(func,np.arange(0,50),correlation)
	ksai = popt[0]
	offset = popt[2]

	C = ave_E_2-ave_E**2
	chi = ave_M_2-ave_M**2

	if T == 2.27:
		fig = plt.figure()
		plt.plot(np.arange(len(correlation)),correlation-offset)
		plt.savefig("correlation @ " + str(T) +".png")
		fig.clear()

	M_list.append(ave_M)
	C_list.append(C)
	chi_list.append(chi)
	ksai_list.append(ksai)

with open('result.csv','w') as csvfile:
	fieldnames = ['T','M','C','chi','ksai']
	writer = csv.DictWriter(csvfile,fieldnames = fieldnames)

	writer.writeheader()
	for i in range(len(T_list)):
		writer.writerow({'T':T_list[i],'M':M_list[i],'C':C_list[i],'chi':chi_list[i],'ksai':ksai_list[i]})


