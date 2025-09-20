# -*- coding: utf-8 -*-

# Spyder Editor

# This is a temporary script file.



import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots()
ax.plot([1, 2, 3, 4], [1, 4, 2, 3])
plt.show()


# Simple pyplot figure
fig = plt.figure(figsize=(2, 2), facecolor='lightskyblue',layout='constrained')
fig.suptitle('Figure')
ax = fig.add_subplot()
ax.set_title('Axes', loc='left', fontstyle='oblique', fontsize='medium')
             
  
# This code generates a single figure with multiple Axes corresponding to each plot.
fig, axs = plt.subplots(2, 2, figsize=(4, 3),layout='constrained')


# This code generates a single Figure with multiple Axes corresponding to each plot.
fig, axs = plt.subplot_mosaic([['A', 'right'], ['B', 'right']],
figsize=(4, 3), layout='constrained')
for ax_name, ax in axs.items():
ax.text(0.5, 0.5, ax_name, ha='center', va='center')



# 1D numpy array
np.array([7, 2, 9, 10])

# 2D numpy array
np.array([[5.2, 3.0, 4.5],
[9.1, 0.1, 0.3]])

# 3D numpy array
np.array([[[5.2, 3.0],
[1.9, 4.5],
[9.1, 0.1]],
[[1.3, 4.2],
[2.9, 6.6],
[3.2, 8.2]]])


# number of axes(dimensions) of the array
ndarray.ndim = 

# dimensions of the array, tuple indicating the size in each dimension
ndarray.shape = 

# the total number of elements of the array
ndarray.size = 



# indexing and slicing
a = np.arange(10)**3
a
a[2]
a[2:5]
# equivalent to a[0:6:2] = 1000; from start to position 6, exclusive, set every 2nd element to 1000
a[:6:2] = 1000
a
a[::-1] # reversed a
for i in a:
print(i**(1 / 3.))




# x and y are positions in the array
def f(x, y):
    return 10 * x + y
b = np.fromfunction(f, (5, 4), dtype=int)
b
b[2, 3]
b[0:5, 1] # each row in the second column of b
b[:, 1] # equivalent to the previous example
b[1:3, :] # each column in the second and third row of b




# In-class activity #1
three_d_array = np.array([
    [[1, 2, 3],
     [4, 5, 6]],
    [[7, 8, 9],
     [10, 11, 12]]])




# Using a dicitonary to plot
np.random.seed(19680801) # seed the random number generator.
data = {'a': np.arange(50),
'c': np.random.randint(0, 50, 50),
'd': np.random.randn(50)}
data['b'] = data['a'] + 10 * np.random.randn(50)
data['d'] = np.abs(data['d']) * 100
fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.scatter('a', 'b', c='c', s='d', data=data)
ax.set_xlabel('entry a')
ax.set_ylabel('entry b')




# object oriented programming style of plotting
# all these are basically methods of 'ax'
x = np.linspace(0, 2, 100) # Sample data.
fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.plot(x, x, label='linear') # Plot some data on the Axes.
ax.plot(x, x**2, label='quadratic') # Plot more data on the Axes...
ax.plot(x, x**3, label='cubic') # ... and some more.
ax.set_xlabel('x label') # Add an x-label to the Axes.
ax.set_ylabel('y label') # Add a y-label to the Axes.
ax.set_title("Simple Plot") # Add a title to the Axes.
ax.legend() # Add a legend.
plt.show()



# plotting using pyplot
x = np.linspace(0, 2, 100) # Sample data.
plt.figure(figsize=(5, 2.7), layout='constrained')
plt.plot(x, x, label='linear') # Plot some data on the (implicit) Axes.
plt.plot(x, x**2, label='quadratic') # etc.
plt.plot(x, x**3, label='cubic')
plt.xlabel('x label')
plt.ylabel('y label')
plt.title("Simple Plot")
plt.legend()


# changing color
data1, data2 = np.random.randn(2, 100)
fig, ax = plt.subplots(figsize=(5, 2.7))
ax.scatter(data1, data2, s=50, facecolor='C0', edgecolor='k’)
plt.show()


# marker style
data1, data2, data3, data4 = np.random.randn(4, 100)
fig, ax = plt.subplots(figsize=(5, 2.7))
ax.plot(data1, 'o', label='data1')
ax.plot(data2, 'd', label='data2')
ax.plot(data3, 'v', label='data3')
ax.plot(data4, 's', label='data4')
ax.legend()
plt.show()


# axes scales and tick marks
data1 = np.random.randn(1, 100)[0]
fig, axs = plt.subplots(1, 2, figsize=(5, 2.7), layout='constrained')
xdata = np.arange(len(data1)) # make an ordinal for this
data = 10**data1
axs[0].plot(xdata, data)
axs[1].set_yscale('log')
axs[1].plot(xdata, data)
plt.show()



data1 = np.random.randn(1, 100)[0]
fig, axs = plt.subplots(2, 1, layout='constrained')
axs[0].plot(xdata, data1)
axs[0].set_title('Automatic ticks’)
axs[1].plot(xdata, data1)
axs[1].set_xticks(np.arange(0, 100, 30), ['zero', '30', 'sixty', '90'])
axs[1].set_yticks([-1.5, 0, 1.5]) # note that we don't need to specify labels
axs[1].set_title('Manual ticks’)
plt.show()



# Labeling plots
mu, sigma = 115, 15
x = mu + sigma * np.random.randn(10000)
fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
# the histogram of the data
n, bins, patches = ax.hist(x, 50, density=True, facecolor='C0', alpha=0.75)
ax.set_xlabel('Length [cm]')
ax.set_ylabel('Probability')
ax.set_title('Aardvark lengths\n (not really)')
ax.text(75, .025, r'$\mu=115,\ \sigma=15$')
ax.axis([55, 175, 0, 0.03])
ax.grid(True)
plt.show()



# Adding a legend
data1, data2, data3 = np.random.randn(3, 100)
fig, ax = plt.subplots(figsize=(5, 2.7))
ax.plot(np.arange(len(data1)), data1, label='data1')
ax.plot(np.arange(len(data2)), data2, label='data2')
ax.plot(np.arange(len(data3)), data3, 'd', label='data3')
ax.legend()
plt.show()



# Saving figures
data1, data2, data3 = np.random.randn(3, 100)
fig, ax = plt.subplots(figsize=(5, 2.7))
ax.plot(np.arange(len(data1)), data1, label='data1')
ax.plot(np.arange(len(data2)), data2, label='data2')
ax.plot(np.arange(len(data3)), data3, 'd', label='data3')
ax.legend()
fig.savefig('figure_with_legend.png', dpi=200)