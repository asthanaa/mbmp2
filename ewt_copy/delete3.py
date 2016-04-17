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

f_density = open("readfile/test2.50", "r")
f_int = open("Li2/integrals.dat", "r")
f_fock = open("Li2/fock.dat", "r")

f_svdev = open("delete1.txt", "w")

f_overlap = open("delete2.txt", "w")
f_svdev.write("start\n")
#...........input for spin free wicks therem
print "\n Many Body Multireference Perturbation Theory\n"

#Read the input orbital indeces Current assumption : 2 i, 2 u, 3 a
#write a function value to read the file and return the value of the density matrix here
size_vir=24
size_u=2
size_a=3
size_i=2
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
#print 'gamma\n', gamma_sin
g2, g2mv = for_func.make2(gamma_sin, gamma_mv, coeff, size_u, size_a, size_i)
lambda2 = for_func.lambda2(gamma_sin,g2, coeff, size_u, size_a, size_i)
print 'Lambda\n', np.array(lambda2)

print 'gamma\n', np.array(gamma_sin)
print 'gamma 2\n', np.array(g2)
C = []
#make 1st excitation array of S(u,v)
#Type 1, 2, 3
#Sia=np.array(single_exctd_ovrlp(holes, active, particle, 1))

C.append([])
'''
Siu, Cl=func_mp2.single_excited_overlap(overall_i, overall_u, overall_a, 2, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 2'
Sua, Cl=func_mp2.single_excited_overlap(overall_i, overall_u, overall_a, 3, gamma_sin, eta_sin, lambda2)
print Sua, Cl
C.append(Cl)
print 'done 3'
S2_iajb, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 41, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 4'
S2_iaib, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 42, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 5'
S2_iaia, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 43, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 6'

S2_iuju, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 51, gamma_sin, eta_sin, lambda2)
C.append(Cl)
#print Cl
print 'done 7'
'''
S2_iuiu, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 52, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 8'
'''
S2_aubu, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 61, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print Cl
print 'done 9'

S2_auau, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 62, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 10'
S2_iuuu, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 7, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 11'
S2_uauu, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 8, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 12'
S2_iauu_iuua, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 9, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 13'
S2_iuja, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 101, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 14'
S2_iuia, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 102, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 15'
S2_iaub, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 111, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 16'
S2_iaua, Cl = func_mp2.doubly_excited_overlap(overall_i, overall_u, overall_a, 112, gamma_sin, eta_sin, lambda2)
C.append(Cl)
print 'done 17'
#print C



X_iu, D_iu, Y_iu= np.linalg.svd(Siu, full_matrices=1, compute_uv=1)
X_ua, D_ua, Y_ua= np.linalg.svd(Sua, full_matrices=1, compute_uv=1)
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
    print 'sum not matched', '\nX2\n', X2_iuiu, '\nD2\n', D2_iuiu
else:
    print 'sum matched'
X2_aubu, D2_aubu, Y2_aubu = np.linalg.svd(S2_aubu, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_aubu, D2_aubu):
    print 'sum not matched', '\nX2\n', X2_aubu, '\nD2\n', D2_aubu
else:
    print 'sum matched'
X2_auau, D2_auau, Y2_auau = np.linalg.svd(S2_auau, full_matrices=1, compute_uv=1)
if not func_mp2.sum_check(X2_auau, D2_auau):
    print 'sum not matched', '\nX2\n', X2_auau, '\nD2\n', D2_auau
else:
    print 'sum matched'

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





f_svdev.write('2\n')
f_overlap.write('2\n')
np.savetxt(f_svdev, D_iu, delimiter=' ')   # X is an array
np.savetxt(f_overlap, Siu, delimiter=' ')   # X is an array
f_svdev.write('3\n')
f_overlap.write('3\n')
np.savetxt(f_svdev, D_ua, delimiter=' ')   # X is an array
np.savetxt(f_overlap, Sua, delimiter=' ')   # X is an array

if len(S2_iajb):
    f_svdev.write('41\n')
    f_overlap.write('41\n')
    np.savetxt(f_svdev, D2_iajb, delimiter=' ')   # X is an array
    np.savetxt(f_overlap, S2_iajb, delimiter=' ')   # X is an array
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


#f_overlap.write('00')
#f_overlap.close()

#print 'started reading'

#f_overlap = open("overlap.txt", "r")
#temp = func_mp2.readoverlap(f_overlap)
f_svdev.write('7\n')
f_overlap.write('7\n')

np.savetxt(f_svdev, D2_iuuu, delimiter=' ')   # X is an array
np.savetxt(f_overlap, S2_iuuu, delimiter=' ')   # X is an array

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

#print 'started reading'

#f_overlap = open("delete2.txt", "r")
#temp = func_mp2.readoverlap(f_overlap)




#truncate the mastrics


print C[1], C[10]
mat_S2S11 = func_mp2.orthogonalize(C[1], C[10], Siu, gamma_sin, eta_sin, lambda2)
print 'S11\n', S2_iuuu
S2_iuuu = np.subtract(S2_iuuu, mat_S2S11)
print 'new S11\n', S2_iuuu
print C[2], C[11]
mat_S3S12 = func_mp2.orthogonalize(C[1], C[10], Siu, gamma_sin, eta_sin, lambda2)
print 'S12\n', np.array(S2_uauu)
S2_uauu = np.subtract(S2_uauu, mat_S3S12)
print 'new S12\n', S2_uauu

X=[]
X.append([])
X.append(func_mp2.return_truncatedX(Siu, 2))
X.append(func_mp2.return_truncatedX(Sua, 3))
if len(S2_iajb):
    X1 = func_mp2.return_truncatedX(S2_iajb, 4)
X1 = func_mp2.return_truncatedX(S2_iaib, 5)
X1 = func_mp2.return_truncatedX(S2_iaia, 6)
if len(S2_iuju):
    X1 = func_mp2.return_truncatedX(S2_iuju, 7)
X1 = func_mp2.return_truncatedX(S2_iuiu, 8)
X1 = func_mp2.return_truncatedX(S2_aubu, 9)
X1 = func_mp2.return_truncatedX(S2_auau, 10)
X1 = func_mp2.return_truncatedX(S2_iuuu, 11)
X1 = func_mp2.return_truncatedX(S2_uauu, 12)
X1 = func_mp2.return_truncatedX(S2_iauu_iuua, 13)
if len(S2_iuja):
    X1 = func_mp2.return_truncatedX(S2_iuja, 14)
X1 = func_mp2.return_truncatedX(S2_iuia, 15)
X1 = func_mp2.return_truncatedX(S2_iaub, 16)
X1 = func_mp2.return_truncatedX(S2_iaua, 17)

#integrals = func_mp2.readintegrals(f_int, size_vir, size_i, size_u)
fock = func_mp2.readfock(f_fock, size_vir, size_i, size_u)


#Do the energy calculation

print '----------------- energy contribution ', func_mp2.calculate_F(fock, C[1], Siu, X[1], size_i, size_u)
print '----------------- energy contribution ', func_mp2.calculate_F(fock, C[2], Sua, X[2], size_i, size_u)
#print '----------------- energy contribution ', func_mp2.calculate_F(fock, C[3], Siu, X[3], size_i, size_u)
#make array of fock marix needed

#make 2D martix by multiplying each vector to the destined fock 1D matrix

'''
#use overlap atrix to calculate the total energy contribution 
