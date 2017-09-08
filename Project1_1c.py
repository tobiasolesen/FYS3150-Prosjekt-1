# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
from numpy import *
from numpy.linalg import solve

N = 1000        #Gridpoints
h = 1/(N+1)    #Step length
a = -1         #Values for Matrix elements
b = 2
c = -1

#Analytic solution
def U(x):
	return 1-(1-exp(-10))*x - exp(-10*x)

#Numerical solution
def f(x):
	f = 100*exp(-10*x) 
	return f 

#Defining all arrays needed for storing values 
#Vektor v for Av = B
#Vektor B og A for Av = B
d_twidd = zeros(N+2)
f_twidd = zeros(N+2)

x = zeros(N+2)
B = zeros(N+2)
A = zeros((N+2,N+2))

#Conditions for x-array and numerical solution Av = B 

x[0] = 0
x[N+1] = 1

v = zeros(N+2)
v[0] = 0  
v[N+1] = 0

#Producing x-array

for j in range(N+1):
	x[j] = j*h 

#Filling A with 2's on diag and -1 on diags below and above diag
#Setting values for B

f_twidd[0] = h**2*f(x[0])  #Init value for row reduction
d_twidd[0] = 2        #Init value for row reduction
d_twidd[N+1] = 2

for i in range(N+2):

	A[i][i] = b
	if i!=0 :
		A[i][i-1] = a
	if i!=(N+2):
		A[i-1][i] = c
	A[N+1][0] = 0
	f_twidd[i] = h**2*f(x[i])

#Diag elements
for k in range(1, N+1):
	d_twidd[k] = (k+1.)/k

#Forward sub
for k in range(2, N+1):
	f_twidd[k] = h**2*f(x[k]) + f_twidd[k-1]/d_twidd[k-1]


v[N] = f_twidd[N]/d_twidd[N]

#Backward sub
for k in reversed(range(2, N+1)):
	v[k-1] = (f_twidd[k-1] + v[k])/d_twidd[k-1]

t = linspace(0,1,100)

plot(t, U(t), label = 'Analytic sol U(x)')
hold('on')
plot(x, v, label = 'Num sol (v(x))')
xlabel('$x$', fontsize = 18)
ylabel('$y$', fontsize = 18)
title('$Possion$ $equation,$ $analytic$ $and$ $numerical$', fontsize = 20)
legend()
show()
