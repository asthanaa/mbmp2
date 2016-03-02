import numpy
import func_mp2

f_out1 = open("output1.txt","w")
f_out2 = open("output2.txt","w")
f_out3 = open("output3.txt","w")
f_out4 = open("output4.txt","w")
f_density = open("readfile/fort.50", "r")
size_u=4
size_a=20
size_i=2

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
	print 'changed to 0....'

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

#print gamma_mv
gamma2, gamma2_mv = func_mp2.make_gamma2(gamma_sin, gamma_mv, coeff, size_u, size_a, size_i)
#print '\ngamma 2 :\n', gamma2, '\ngamma mv\n',gamma2_mv, '\n'
f_out2.write( '\n gamma2:\n')
for i in range(size_u):
    for j in range(size_u):
	for k in range(size_u):
	    for l in range(size_u):
		#for m in range(size_a):
		    #for n in range(size_a):
		if abs(gamma2[i][j][k][l]-0.0) >0.00001:
		    f_out2.write('{} {} {} {} {}\n'.format(i, j, k, l, gamma2[i][j][k][l]))



lambda2 = func_mp2.make_lambda2(gamma2, gamma_sin, coeff, size_u, size_a, size_i)
#print '\nlambda 2 :\n', lambda2, '\n'


f_out2.write('\nlambda2 :\n')

for i in range(size_u):
    for j in range(size_u):
	for k in range(size_u):
	    for l in range(size_u):
		if abs(lambda2[i][j][k][l]-0.0)>0.00001:
		    f_out2.write('{} {} {} {} {}\n'.format(i, j, k, l, lambda2[i][j][k][l]))

gamma3 = func_mp2.make_gamma3(gamma2_mv, gamma2, gamma_mv, gamma_sin, coeff, size_u, size_a, size_i)
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

