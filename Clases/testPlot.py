import matplotlib.pyplot as plt
import numpy as np 


x =	np.arange(1,5,0.5)
y1 = x - 2
y2 = x**2 - 100*np.cos(x)
y3 = x - np.exp(x)

# plot the data
plt.plot(x, y1, color ='blue', label='f1(x)')
plt.plot(x, y2, color ='green', label='f2(x)')
plt.plot(x, y3, color ='red', label='f3(x)')
plt.xlabel('x - axis')
plt.ylabel('y - axis')
# Title of plot 
plt.title('Ejemplo clase')
plt.legend()
# display the plot
plt.show()
