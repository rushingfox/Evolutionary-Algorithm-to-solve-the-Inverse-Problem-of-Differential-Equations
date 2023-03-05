import numpy as np 
from matplotlib import pyplot as plt 
 
a = np.loadtxt('../EAForParabola/result0.1.txt')#最普通的loadtxt
#a = np.loadtxt('../EAForEllipse/result0.1.txt')
print(a)

x = np.asarray(a[0])  
y_estimate_q = np.asarray(a[1])
y_real_q = np.asarray(a[2])

plt.title("u_observe and u_estimate") 
plt.xlabel("x axis caption") 
plt.ylabel("y axis caption") 
line1=plt.plot(x,y_estimate_q,color="green",label="y_estimate") 
line2=plt.plot(x,y_real_q,color="red",label="y_real")
plt.legend()
plt.savefig('./figure1.jpg')
plt.show()