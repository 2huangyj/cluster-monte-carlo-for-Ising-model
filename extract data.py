import numpy as np
import matplotlib.pyplot as plt
data = np.genfromtxt('result.csv', delimiter = ',',skip_header = 1,  names = ['T','M','C','chi','ksai'])
D = 3

fig1 = plt.figure()
plt.scatter(data['T'],data['M'],color = 'r')
plt.plot(data['T'],data['M'],color = 'b')
fig1.suptitle('M-T', fontsize=20)
plt.xlabel('T(/J)', fontsize=18)
plt.ylabel('M/u.c.', fontsize=16)
plt.savefig("M-T in D-" + str(D) +".png")
fig1.clear()

fig2 = plt.figure()
plt.scatter(data['T'],data['C'],color = 'r')
plt.plot(data['T'],data['C'],color = 'b')
fig2.suptitle('C-T', fontsize=20)
plt.xlabel('T(/J)', fontsize=18)
plt.ylabel('C', fontsize=16)
plt.savefig("C-T in D-" + str(D) +".png")
fig2.clear()

fig3 = plt.figure()
plt.scatter(data['T'],data['chi'],color = 'r')
plt.plot(data['T'],data['chi'],color = 'b')
fig3.suptitle('chi-T', fontsize=20)
plt.xlabel('T(/J)', fontsize=18)
plt.ylabel('chi', fontsize=16)
plt.savefig("chi-T in D-" + str(D) +".png")
fig3.clear()

fig4 = plt.figure()
plt.scatter(data['T'],data['ksai'],color = 'r')
plt.plot(data['T'],data['ksai'],color = 'b')
fig4.suptitle('ksai-T', fontsize=20)
plt.xlabel('T(/J)', fontsize=18)
plt.ylabel('ksai', fontsize=16)
plt.savefig("ksai-T in D-" + str(D) +".png")
fig4.clear()