import ewt
import fix_uv
import func_ewt
import copy
import sys
from collections import deque
import make_c
import func_mp2
import for_func
import numpy as np
fix_temp = fix_uv
func = func_ewt

f = open("tec.txt", "w")

f_density = open("readfile/test1.50", "r")

f_svdev = open("delete1.txt", "w")

f_overlap = open("delete2.txt", "w")
f_svdev.write("start\n")
#...........input for spin free wicks therem
print "\n Many Body Multireference Perturbation Theory\n"

#Read the input orbital indeces Current assumption : 2 i, 2 u, 3 a
#write a function value to read the file and return the value of the density matrix here
size_vir=2
size_u=2
size_a=1
size_i=1
active = [str(i) for i in range(size_u)]
holes = [str(i) for i in range(size_i)]
particles = [str(i) for i in range(size_vir)]

#make operators without dagger or undegged. -1 for the dagger value meaning not assigned
overall_i, overall_u, overall_a = func.make_operators(holes, active, particles)
#make the single excitation density matrixi
print overall_a
gamma_sin = []
gamma_mv = []
eta_sin = []
coeff = []
gamma_sin, gamma_mv, eta_sin, coeff = func_mp2.read_density_file(f_density, size_u, size_a,size_i)

#print gamma_sin, eta_sin, coeff
#test the created operators a little 

'''
test1 = overall_i[0]

test2 = overall_a[0]
i1=copy.deepcopy(test1)
i2=copy.deepcopy(test1)
i3=copy.deepcopy(test1)
i4=copy.deepcopy(test1)
a1=copy.deepcopy(test2)
a2=copy.deepcopy(test2)
a3=copy.deepcopy(test2)
a4=copy.deepcopy(test2)

string1 = func_mp2.make_dagger(i1,i2,a1,a2,0,1) 
string2 = func_mp2.make_dagger(a3,a4,i3,i4,4,2)
print string1, string2, 'here are the two stringi----------------'
print ewt.ewt(string1,string2)

overall_u[0].dag='1'
overall_u[0].string=2
overall_u[0].pos=5
overall_u[1].dag='1'
overall_u[1].string=2
overall_u[1].pos=6
overall_u[2].dag='0'
overall_u[2].string=2
overall_u[2].pos=7

overall_i[1].dag='0'
overall_i[1].name='i'
overall_i[1].string=2
overall_i[1].pos=8
copy_ovi = copy.deepcopy(overall_i[0])
copy_ovi.dag='1'
copy_ovi.string=1

copy_ovi.pos=2
copy_ovi.name='i'
overall_u[3].dag='1'
overall_u[3].string=1
overall_u[3].pos=1
overall_u[4].dag='0'
overall_u[4].string=1
overall_u[4].pos=3
overall_u[5].dag='0'
overall_u[5].string=1
overall_u[5].pos=4



string2 = [overall_u[0], overall_u[1], overall_u[2],overall_i[1]]
string1 = [copy_ovi, overall_u[3], overall_u[4], overall_u[5]]

#string1 = [overall_u[0], overall_u[1], overall_u[2]]
#string2 = [overall_u[3], overall_u[4], overall_u[5]]

for item in string1:
    print 'name', item, 'string',item.string,'dag ', item.dag, 'pos', item.pos, 'type', item.kind
for item in string2:
    print 'name', item, 'string',item.string,'dag ', item.dag, 'pos', item.pos, 'type', item.kind


full_con, const_con = ewt.ewt(string1, string2)
#full_con, const_con = ewt.ewt([copy_ovi], [overall_i[1]])
print 'printing full contrations', full_con, const_con
#full_con_num = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
#print 'the answer is ', full_con_num
#this comment is till the end :



def sum_check(a, b):
    sum1=0.0
    for i in range(len(a)):
        sum1=sum1+a[i][i]
    sum2=0.0
    for i in range(len(b)):
        sum2=sum2+b[i]
    if sum1-sum2<0.00001:
	print 'matched'
        return 1
    else: 
	return 0
'''
#gamma2, gamma2_mv = func_mp2.make_gamma2(gamma_sin, gamma_mv, coeff, size_u, size_a, size_i)
#lambda2 = func_mp2.make_lambda2(gamma2, gamma_sin, coeff, size_u, size_a, size_i)
g2, g2mv = for_func.make2(gamma_sin, gamma_mv, coeff, size_u, size_a, size_i)
lambda2 = for_func.lambda2(gamma_sin,g2, coeff, size_u, size_a, size_i)
#make 1st excitation array of S(u,v)
#Type 1, 2, 3
#Sia=single_exctd_ovrlp(holes, active, particle, 1)
#Siu=func_mp2.single_excited_overlap(overall_i, overall_u, overall_a, 2, gamma_sin, eta_sin, lambda2)
#Siu = np.array(Siu)

#print 'this is the martrix\n', Sua
#print np.linalg.svd(Sui_iu, full_matrices=1, compute_uv=1)
#Sua=func_mp2.single_excited_overlap(overall_i, overall_u, overall_a, 3, gamma_sin, eta_sin, lambda2)
#Sua = np.array(Sua)
#form the 2 body S matrixes

S2_iajb = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 41, gamma_sin, eta_sin, lambda2))

print 'done 1'
S2_iaib = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 42, gamma_sin, eta_sin, lambda2))
print 'done 2'
S2_iaia = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 43, gamma_sin, eta_sin, lambda2))
print 'done 3'

S2_iuju = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 51, gamma_sin, eta_sin, lambda2))
print 'done 4'

S2_iuiu = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 52, gamma_sin, eta_sin, lambda2))
print 'done 5'
S2_aubu = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 61, gamma_sin, eta_sin, lambda2))
print 'done 6'
S2_auau = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 62, gamma_sin, eta_sin, lambda2))
print 'done 7'

S2_iuuu = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 7, gamma_sin, eta_sin, lambda2))
print 'done 8'
S2_uauu = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 8, gamma_sin, eta_sin, lambda2))
print 'done 9'
S2_iauu_iuua = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 9, gamma_sin, eta_sin, lambda2))
print 'done 10'
S2_iuja = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 101, gamma_sin, eta_sin, lambda2))
print 'done 11'
S2_iuia = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 102, gamma_sin, eta_sin, lambda2))
print 'done 12'
S2_iaub = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 111, gamma_sin, eta_sin, lambda2))
print 'done 13'
S2_iaua = np.array(func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 112, gamma_sin, eta_sin, lambda2))
print 'done 14'



#'''
if len(S2_iajb):
    X2_iajb, D2_iajb, Y2_iajb = np.linalg.svd(S2_iajb, full_matrices=1, compute_uv=1)
    if not func_mp2.sum_check(X2_iajb, D2_iajb):
        print 'sum not matched'
    else:
        print 'sum matched'

X2_iaib, D2_iaib, Y2_iaib = np.linalg.svd(S2_iaib, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_iaib, D2_iaib):
    print 'sum not matched'
else:
    print 'sum matched'
X2_iaia, D2_iaia, Y2_iaia = np.linalg.svd(S2_iaia, full_matrices=1, compute_uv=1)

if not func_mp2.sum_check(X2_iaia, D2_iaia):
    print 'sum not matched'
else:
    print 'sum matched'
if len(S2_iuju):
    X2_iuju, D2_iuju, Y2_iuju = np.linalg.svd(S2_iuju, full_matrices=1, compute_uv=1)
    if not func_mp2.sum_check(X2_iuju, D2_iuju):
        print 'sum not matched'
    else:
        print 'sum matched'

X2_iuiu, D2_iuiu, Y2_iuiu = np.linalg.svd(S2_iuiu, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_iuiu, D2_iuiu):
    print 'sum not matched'
else:
    print 'sum matched'
X2_aubu, D2_aubu, Y2_aubu = np.linalg.svd(S2_aubu, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_aubu, D2_aubu):
    print 'sum not matched'
else:
    print 'sum matched'
X2_auau, D2_auau, Y2_auau = np.linalg.svd(S2_auau, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_auau, D2_auau):
    print 'sum not matched'
else:
    print 'sum matched'
#'''
X2_iuuu, D2_iuuu, Y2_iuuu = np.linalg.svd(S2_iuuu, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_iuuu, D2_iuuu):
    print 'sum not matched'
else:
    print 'sum matched'
X2_uauu, D2_uauu, Y2_uauu = np.linalg.svd(S2_uauu, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_uauu, D2_uauu):
    print 'sum not matched'
else:
    print 'sum matched'
X2_iauu_iuua, D2_iauu_iuua, Y2_iauu_iuua= np.linalg.svd(S2_iauu_iuua, full_matrices=1, compute_uv=1)

if not func_mp2.sum_check(X2_iauu_iuua, D2_iauu_iuua):
    print 'sum not matched'
else:
    print 'sum matched'

if len(S2_iuja):
    X2_iuja, D2_iuja, Y2_iuja = np.linalg.svd(S2_iuja, full_matrices=1, compute_uv=1)
    if not func_mp2.sum_check(X2_iuja, D2_iuja):
        print 'sum not matched'
    else:
        print 'sum matched'

X2_iuia, D2_iuia, Y2_iuia = np.linalg.svd(S2_iuia, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_iuia, D2_iuia):
    print 'sum not matched'
else:
    print 'sum matched'
X2_iaub, D2_iaub, Y2_iaub = np.linalg.svd(S2_iaub, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_iaub, D2_iaub):
    print 'sum not matched'
else:
    print 'sum matched'
X2_iaua, D2_iaua, Y2_iaua = np.linalg.svd(S2_iaua, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_iaua, D2_iaua):
    print 'sum not matched'
else:
    print 'sum matched'
#'''





if len(S2_iajb):
    f_svdev.write('41\n')
    f_overlap.write('41\n')

    np.savetxt(f_svdev, D2_iajb, delimiter=' ')   # X is an array
    #for item in S2_iajb:
    #    f_overlap.write(' '.join(map(str, item))) 
    np.savetxt(f_overlap, S2_iajb, delimiter=' ')   # X is an array
    #temp = func_mp2.return_truncatedX(S2_iuju)


f_svdev.write('42\n')
f_overlap.write('42\n')
np.savetxt(f_svdev, D2_iaib, delimiter=' ')   # X is an array

np.savetxt(f_overlap, S2_iaib, delimiter=' ')   # X is an array

f_svdev.write('43\n')
f_overlap.write('43\n')
np.savetxt(f_svdev, D2_iaia, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_iaia, delimiter=' ')   # X is an array



if len(S2_iuju):
    f_svdev.write('51\n')
    f_overlap.write('51\n')
    np.savetxt(f_svdev, D2_iuju, delimiter=' ')   # X is an array
    np.savetxt(f_overlap, S2_iuju, delimiter=' ')   # X is an array




f_svdev.write('52\n')
f_overlap.write('52\n')
np.savetxt(f_svdev, D2_iuiu, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_iuiu, delimiter=' ')   # X is an array

f_svdev.write('61\n')
f_overlap.write('61\n')
np.savetxt(f_svdev, D2_aubu, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_aubu, delimiter=' ')   # X is an array

f_svdev.write('62\n')
f_overlap.write('62\n')
np.savetxt(f_svdev, D2_auau, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_auau, delimiter=' ')   # X is an array

'''
f_overlap.write('00')
f_overlap.close()

print 'started reading'

f_overlap = open("overlap.txt", "r")
temp = func_mp2.readoverlap(f_overlap)
'''
#'''
f_svdev.write('7\n')
f_overlap.write('7\n')
np.savetxt(f_svdev, D2_iuuu, delimiter=' ')   # X is an array


f_svdev.write('8\n')
f_overlap.write('8\n')
np.savetxt(f_svdev, D2_uauu, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_uauu, delimiter=' ')   # X is an array

f_svdev.write('9\n')
f_overlap.write('9\n')
np.savetxt(f_svdev, D2_iauu_iuua, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_iauu_iuua, delimiter=' ')   # X is an array

if len(S2_iuja):
    f_svdev.write('101\n')
    f_overlap.write('101\n')
    np.savetxt(f_svdev, D2_iuja, delimiter=' ')   # X is an array
    np.savetxt(f_overlap, S2_iuja, delimiter=' ')   # X is an array

f_svdev.write('102\n')
f_overlap.write('102\n')
np.savetxt(f_svdev, D2_iuia, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_iuia, delimiter=' ')   # X is an array

f_svdev.write('111\n')
f_overlap.write('111\n')
np.savetxt(f_svdev, D2_iaub, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_iaub, delimiter=' ')   # X is an array

f_svdev.write('112\n')
f_overlap.write('112\n')
np.savetxt(f_svdev, D2_iaua, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_iaua, delimiter=' ')   # X is an array



f_overlap.write('00')
f_overlap.close()

print 'started reading'

f_overlap = open("delete2.txt", "r")
temp = func_mp2.readoverlap(f_overlap)




'''
temp = func_mp2.return_truncatedX(S2_uauu)
print temp


#Schmidt Orthogonalisation :

#Do the singular transformation

#Get the non 0 values
#


'''
