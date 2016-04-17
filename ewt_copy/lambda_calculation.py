import numpy
import func_mp2
import numpy as np
import for_func
f_out1 = open("output1.txt","w")
f_out2 = open("output2.txt","w")
f_out3 = open("output3.txt","w")
f_out4 = open("output4.txt","w")

f_density = open("readfile/CH/fort.50", "r")

f_int = open("readfile/CH/integrals.dat", "r")
f_fock = open("readfile/CH/fock.dat", "r")
size_vir = 14
size_u=4
size_a=20
size_i=1
size_m = 20
gamma_sin = []
gamma_mv = []
eta_sin = []
coeff = []
gamma_sin, gamma_mv, eta_sin, coeff = func_mp2.read_density_file(f_density, size_u, size_a,size_i)
#print 'gamma 1 :\n',gamma_sin, '\neta 1\n',eta_sin, '\ncoeff\n',coeff
#print '\ngamma mv \n', gamma_mv

for i in range(size_a):
    f_out1.write('{}  {}\n'.format(i, coeff[i]))
    if abs(coeff[i]-0.0)<0.00001:
	coeff[i]=0.0
#	print 'changed to 0....'

f_out1.write('single gamma : \n')
for i in range(size_u):
    for j in range(size_u):
	if abs(gamma_sin[i][j]-0.0)>0.00001:
	    f_out1.write('{} {} {}\n'.format(i, j, gamma_sin[i][j]))

f_out1.write('gamma mv : \n')
for i in range(size_u):
    for j in range(size_u):
        for k in range(size_a):
    	    for l in range(size_a):
		if abs(gamma_mv[i][j][k][l]-0.0)>0.00001:
	    	    f_out1.write('{} {} {} {} {}\n'.format(i, j, k, l, gamma_mv[i][j][k][l]))



'''
gamma_sin = np.array(gamma_sin)
gamma_mv = np.array(gamma_mv)
coeff = np.array(coeff)
eta_sin = np.array(eta_sin)
'''


gmv = gamma_mv
g1=gamma_sin

'''
g2 = np.zeros((size_u,size_u,size_u,size_u))
g2mv = np.zeros((size_u,size_u,size_u,size_u,size_a,size_a))
'''
g2 = []
g2mv=[]
g2,g2mv = for_func.make2(gamma_sin, gamma_mv, coeff, size_u, size_a, size_i)
lambda2 = for_func.lambda2(gamma_sin, g2, coeff, size_u, size_a, size_i)

g3 = for_func.make3(gmv, g2, g2mv, coeff, size_u, size_a, size_i)
lambda3 = for_func.lambda3(g1, g3, lambda2, coeff, size_u, size_a, size_i)


gamma2=g2
gamma3=g3
gamma2_mv=g2mv

f_out2.write( '\n gamma2:\n')

for i in range(size_u):
    for j in range(size_u):
	for k in range(size_u):
	    for l in range(size_u):
		#for m in range(size_a):
		    #for n in range(size_a):
		if abs(g2[i][j][k][l]-0.0) >0.00001:
		    f_out2.write('{} {} {} {} {}\n'.format(i, j, k, l, gamma2[i][j][k][l]))

f_out2.write( '\n gamma2_mv:\n')

for i in range(size_u):
    for j in range(size_u):
	for k in range(size_u):
	    for l in range(size_u):
		for m in range(size_a):
		    for n in range(size_a):
		        if abs(gamma2_mv[i][j][k][l][m][n]-0.0) >0.00001:
		            f_out2.write('{} {} {} {} {} {} {}\n'.format(i, j, k, l, m, n, gamma2_mv[i][j][k][l][m][n]))






f_out2.write('\nlambda2 :\n')

for i in range(size_u):
    for j in range(size_u):
	for k in range(size_u):
	    for l in range(size_u):
		if abs(lambda2[i][j][k][l]-0.0)>0.00001:
		    f_out2.write('{} {} {} {} {}\n'.format(i, j, k, l, lambda2[i][j][k][l]))

f_out3.write('\nlambda3:\n')
f_out4.write('\ngamma3:\n')
for i in range(size_u):
    for j in range(size_u):
	for k in range(size_u):
	    for l in range(size_u):
	    	for m in range(size_u):
	    	    for n in range(size_u):

			if abs(lambda3[i][j][k][l][m][n] - 0.0)>0.00001:
			    f_out3.write('{} {} {} {} {} {} {}\n'.format(i, j, k, l, m, n, lambda3[i][j][k][l][m][n]))
			if abs(gamma3[i][j][k][l][m][n] - 0.0)>0.00001:
			    f_out4.write('{} {} {} {} {} {} {}\n'.format(i, j, k, l, m, n, gamma3[i][j][k][l][m][n]))


#calculate energy as 1body term and 2 body term :
energy = 0.0

integrals = func_mp2.readintegrals(f_int, size_vir, size_i, size_u)
fock = func_mp2.readfock(f_fock, size_vir, size_i, size_u)

fock = np.array(fock)
integrals = np.array(integrals)

for i in range(size_i):
    for j in range(size_i):
	energy+= fock[i][i]*2.0
for u in range(size_u):
    for v in range(size_u):
	energy+= fock[u+size_i][v+size_i]*gamma_sin[u][v]
print 'the fock part of energy is :', energy

for u in range(size_u):
    for v in range(size_u):
	for w in range(size_u):
	    for x in range(size_u):

		energy+= integrals[u][v][w+size_i][x+size_i]*gamma2[u][v][w][x]

print 'the CAS energy of the system is : ', energy
#print gamma2

'''
for u in range(size_u):
    for v in range(size_u):
	for w in range(size_u):
	    for x in range(size_u):
#	        for y in range(size_u):
#	    	    for z in range(size_u):
		print u,v,w,x,gamma2[u][v][w][x], gamma2[x][w][v][u]
'''
for u in range(size_u):
    for v in range(size_u):
	for w in range(size_u):
	    for x in range(size_u):
	        for y in range(size_u):
	    	    for z in range(size_u):
			print u,v,w,x,y,z,g2mv[u][v][w][x][y][z], g2mv[u][v][w][x][y][z]

val =for_func.vector_comp(lambda3, gamma3,size_u)
'''

gamma3 = func_mp2.make_gamma3(gamma2_mv, gamma2, gamma_mv, gamma_sin, coeff, size_u, size_a, size_i)
print 'gamma3 calculated'
#print '\ngamma 3:\n', gamma3
lambda3 = func_mp2.make_lambda3(gamma3, gamma_sin, lambda2, coeff, size_u, size_a, size_i)


f_out3.write('\nlambda3:\n')
f_out4.write('\ngamma3:\n')
for i in range(size_u):
    for j in range(size_u):
	for k in range(size_u):
	    for l in range(size_u):
	    	for m in range(size_u):
	    	    for n in range(size_u):

			if abs(lambda3[i][j][k][l][m][n] - 0.0)>0.00001:
			    f_out3.write('{} {} {} {} {} {} {}\n'.format(i, j, k, l, m, n, lambda3[i][j][k][l][m][n]))
			if abs(gamma3[i][j][k][l][m][n] - 0.0)>0.00001:
			    f_out4.write('{} {} {} {} {} {} {}\n'.format(i, j, k, l, m, n, gamma3[i][j][k][l][m][n]))
'''
