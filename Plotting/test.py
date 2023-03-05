import numpy as np 
from matplotlib import pyplot as plt 
 
a = np.loadtxt('./4.1.1.2.txt')#最普通的loadtxt
#a = np.loadtxt('../EAForEllipse/result0.1.txt')
print(a)

#a_array = np.asarray(a)  

names = ['the first round', 'the second round']
x = range(len(names))


for i in range(10):
    y = []
    y.append(a[0][i])
    y.append(a[1][i])
    plt.plot(x, y, 'bo-')
for i in range(10):
    y = []
    y.append(a[0][i])
    y.append(a[2][i])
    plt.plot(x, y, 'ro-')
plt.ylabel("REU-value") 
plt.xticks(x, names)
plt.show()