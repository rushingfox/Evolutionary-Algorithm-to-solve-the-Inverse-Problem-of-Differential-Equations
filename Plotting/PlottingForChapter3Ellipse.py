import numpy as np 
from matplotlib import pyplot as plt 

a = np.loadtxt('../EAForEllipse/result0.1.txt')

print(a)

x = np.asarray(a[0])  
y_estimate_q = np.asarray(a[1])
y_real_q = np.asarray(a[2])

plt.title("q*(x) and q(x)") 
plt.xlabel("x axis caption") 
plt.ylabel("y axis caption") 
line1=plt.plot(x,y_estimate_q,color="green",label="q*(x)") 
line2=plt.plot(x,y_real_q,color="red",label="q(x)")
plt.legend()
plt.savefig('../EAForEllipse/FigureForEllipseOld.jpg')
plt.show()


a = np.loadtxt('../EAForEllipse/result0.2.txt')

print(a)

x = np.asarray(a[0])  
y_estimate_q = np.asarray(a[1])
y_real_q = np.asarray(a[2])

plt.title("q*(x) and q(x)") 
plt.xlabel("x axis caption") 
plt.ylabel("y axis caption") 
line1=plt.plot(x,y_estimate_q,color="green",label="q*(x)") 
line2=plt.plot(x,y_real_q,color="red",label="q(x)")
plt.legend()
plt.savefig('../EAForEllipse/FigureForEllipseNew.jpg')
plt.show()