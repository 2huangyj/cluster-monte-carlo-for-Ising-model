import numpy as np
import matplotlib.pyplot as plt 
import random
from scipy.optimize import curve_fit
import csv

#cluster Monte Carlo
D = 4
l = 10
t_step = 7
J = 1
T_list = [2.0,4.0,5.0,5.5,5.7,6.2,6.4,6.5,6.55,6.6,6.63,6.66,6.68,6.70,6.73,6.76,6.80,6.85,7.0,7.3,7.7,8.0,8.4,8.8]
#Tc = 6.68


M_list = []
C_list = []
chi_list = []
ksai_list = []
#create a stack
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
	grid = 2*np.random.randint(2, size=(l,l,l,l))-1

	s = Stack()
	sum_M = 0.0
	sum_E = 0.0
	sum_M_2 = 0.0
	sum_E_2 = 0.0
	correlation = np.zeros(shape = (l,))

	#flip the spin according to the algorithm
	while t < t_max:
		flipped = np.zeros(shape=(l,l,l,l))
		i = np.random.randint(l - 1)
		j = np.random.randint(l - 1)
		k = np.random.randint(l - 1)
		m = np.random.randint(l - 1)

		s.push((i,j,k,m))
		flipped[i,j,k,m] = 1
		
		M = 0.0
		E = 0.0

		while not s.isEmpty():
			x,y,z,u = s.pop()
			for (x_,y_,z_,u_) in [((x-1)%l,y,z,u),(x,(y-1)%l,z,u),(x,y,(z-1)%l,u),((x+1)%l,y,z,u),(x,(y+1)%l,z,u),(x,y,(z+1)%l,u),(x,y,z,(u+1)%l),(x,y,z,(u-1)%l)]:
				if (grid[x,y,z,u] == grid[x_,y_,z_,u_]) and (flipped[x_,y_,z_,u_]==0) and (p>random.random()):
					s.push((x_,y_,z_,u_))
					flipped[x_,y_,z_,u_]=1
			grid[x,y,z,u] = -grid[x,y,z,u] 
			#inner while loop is over


		if (t+1)%t_step == 0:
			print t
			#calculate M
			M = sum(sum(sum(sum(grid))))
			M = M/(l**D*1.0)
			sum_M = sum_M + np.abs(M)
			sum_M_2 = sum_M_2 + M * M

			#calculate E		
			for i in range(l):
				for j in range(l):
					for k in range(l):
						E = E - J*(grid[i,j,k,m]*grid[i,(j+1)%l,k,m] + grid[i,j,k,m]*grid[(i+1)%l,j,k,m] + grid[i,j,k,m] * grid[i,j,(k+1)%l,m]+ grid[i,j,k,m] * grid[i,j,k,(m+1)%l])
						
			E = E/(l**D*1.0)
			sum_E = sum_E + np.abs(E)
			sum_E_2 = sum_E_2 + E * E

			#calculate the correlation length
			cross  = np.zeros(shape = (l,))
			for i in range(l):
				for j in range(l):
					for k in range(l):
						grid_piece_k = (np.abs(np.fft.fft(grid[:,i,j,k])))**2
						cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))
						grid_piece_k = (np.abs(np.fft.fft(grid[i,:,j,k])))**2
						cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))
						grid_piece_k = (np.abs(np.fft.fft(grid[i,j,:,k])))**2
						cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))
						grid_piece_k = (np.abs(np.fft.fft(grid[i,j,k,:])))**2
						cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))

			cross = np.sum(cross,axis = 0)/(1.0*D*l**D) 
			correlation = correlation + cross
		

		#print t
		t = t + 1
		#outer while loop is over

	ave_M = sum_M/((t_max//t_step)*1.0)
	ave_M_2 = sum_M_2/((t_max//t_step)*1.0)

	ave_E = sum_E/((t_max//t_step)*1.0)
	ave_E_2 = sum_E_2/((t_max//t_step)*1.0)

	correlation = (correlation/((t_max//t_step)*1.0))[0:6]
	popt,pcov = curve_fit(func,np.arange(0,6),correlation)
	ksai = popt[0]
	offset = popt[2]

	C = ave_E_2-ave_E**2
	chi = ave_M_2-ave_M**2

	if T == 6.68:
		fig = plt.figure()
		plt.plot(np.arange(len(correlation)),correlation-offset)
		plt.savefig("correlation @ " + str(T) +".png")
		fig.clear()
#T_list
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
