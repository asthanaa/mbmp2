import copy
import ewt
import numpy as np
def single_excited_overlap(holes, active, particles, matrix_no):
    #for i_a ia case :
    if matrix_no==2 or matrix_no==3:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	if matrix_no == 2:
	    if particles:

	        orb1 = copy.deepcopy(particles[0])
	        orb2 = copy.deepcopy(particles[0])

	    else :
		print "matrix type 2 cannot be formed : particles or holes not present \n"
		return 0
	elif matrix_no == 3:
	    if holes:

		orb1 = copy.deepcopy(holes[0])
		orb2 = copy.deepcopy(holes[0])
	    else :
		print "matrix type 3 cant be formed : holes or particles not present \n"
		return 0
	#make the daggered matrix of u or the undaggered version: save them in a and b matrix according to the no_matrix

	active_d = copy.deepcopy(active)
	for item in active_d :
    	    item.dag = 1

	matrix = []
	active_n = copy.deepcopy(active)
	for item in active_n :
    	    item.dag = 0
	#run the loop for u:  make a triangle matrix
	for item1 in active_d:
	    matrixtmp = []#important to empty otherwise add extra terms
	    p1 = 1
	    #the two cases of i->u or a->u
	    if matrix_no == 2:
	        item1.dag = '1'
	        orb1.dag = '0'
	        string1 = [item1, orb1]
	    if matrix_no == 3:
	        item1.dag = '0'
	        orb1.dag = '1'
	        string1 = [orb1, item1]
	    for item in string1:
		item.string = 1
		item.pos = p1
		p1+=1
	    for item2 in active_n:
		p2 = p1
		full_con = []
		#the two cases of i->u or a->u
	    	if matrix_no == 2:
	            item2.dag = '0'
	            orb2.dag = '1'
	            string2 = [orb2, item2]
	        if matrix_no == 3:
	            item2.dag = '1'
	            orb2.dag = '0'
	            string2 = [item2, orb2]
		for item in string2:
		    item.string = 2	
		    item.pos = p2
		    p2+=1
		string3 = string1+string2
		for item in string3:
		    print "item position in the operators formed at ewt : ", item.pos
		full_con, const_con = ewt.ewt(string1, string2)
#		full_con_value = evaluate(full_con, const_con)
		matrixtmp.append(copy.deepcopy(full_con))
            matrix.append(copy.deepcopy(matrixtmp))
    print matrix
    print "------------------------------------------------for matrix_no = ", matrix_no
	    #make the triangle matrix : multiply the terms of the same term, add them to the next list by multiplying to the sign and constant
    return matrix	
def read_density_file(f_density, size_u, size_a, size_i):
    mat = [[0.0 for x in range(size_u)] for x in range(size_u)]
    gamma_mv = []
    temp = 0.0
    gamma_mv = [[[[0.0 for x in range(size_a)] for x in range(size_a)] for x in range(size_u)] for x in range(size_u)]
    coeff = []
    flag = 0
    for line in f_density.readlines():
	if flag == 3:
	    try:
		i, j, k, l, den1 = line.split()
		gamma_mv[int(k)-1][int(l)-1][int(i)-1][int(j)-1] = float(den1)
	    except:
	        flag = 0

		print 'end of gamma_mv at ijkl : ', i, j, k, l
	if flag == 2:
	    try:
		coeff.append(float(line))
	    except:
		flag = 0
	if flag == 1:
	    try:
	        i, j, dens = line.split()
	        mat[int(i)-1][int(j)-1] = float(dens)
	    except:
		flag = 0
		print 'end of gamma_sin at ij : ', i, j

	if 'denisty matrix' in line:
	    flag = 1
	elif 'coeff0' in line:
	    flag = 2
	elif 'i    j    k    l    den1' in line:
	    flag = 3
    gamma_sin = copy.deepcopy(mat)
    eta_sin = copy.deepcopy(mat)
    if eta_sin:
	for i in range(len(eta_sin)):
	    for j in range(len(eta_sin)):
		eta_sin[i][j] = 2.0-eta_sin[i][j] 
    return gamma_sin, gamma_mv, eta_sin, coeff


'''
def read_density_file(f_density, size_u, size_a, size_i):
    mat = [[0.0 for x in range(size_u)] for x in range(size_u)]
    gamma_mv = []
    gamma_mv = [[[[1.0 for x in range(size_a)] for x in range(size_a)] for x in range(size_u)] for x in range(size_u)]
    coeff = []
    flag = 0
    for line in f_density.readlines():
	if flag == 2:
	   try:
		mat1 = line.split()
		for item in mat1:
		   coeff.append(float(item))
		flag = 0
	   except:
		print "coeff matrix not found"
	if flag == 1:
	    try:
	        i, j, dens = line.split()
	    except:
		print "not working line split density"
	    mat[int(i)-1][int(j)-1] = float(dens)
	    mat[int(j)-1][int(i)-1] = float(dens)
	if 'u v density(u,v)' in line:
	    flag = 1
	elif 'coeff0' in line:
	    flag = 2
    gamma_sin = copy.deepcopy(mat)
    eta_sin = copy.deepcopy(mat)
    if eta_sin:
	for i in range(len(eta_sin)):
	    for j in range(len(eta_sin)):
		eta_sin[i][j] = 2.0-eta_sin[i][j] 
    return gamma_sin, gamma_mv, eta_sin, coeff
'''


def make_gamma2(gamma_sin, gamma_mv, coeff, size_u, size_a, size_i):
    mat_ele = 0
    matmv_ele = 0
    mat_mat1 = []
    mat_mat2 = []
    mat_mat3 = []
    mat = []
    matmv_mat1 = []
    matmv_mat2 = []
    matmv_mat3 = []
    matmv_mat4 = []
    matmv_mat5 = []
    matmv = []
    for u in range(size_u):
        for v in range(size_u):
            for w in range(size_u):
		for x in range(size_u):
    	    	    for m in range(size_a):
        		for n in range(size_a):
			    for o in range(size_a):
	    	    	        mat_ele+= coeff[m]*coeff[n]*gamma_mv[u][w][m][o]*gamma_mv[v][x][o][n] 
			    if v == w:
			        mat_ele-=coeff[m]*coeff[n]*gamma_mv[u][x][m][n]
		    mat_mat1.append(mat_ele)
	    	    mat_ele = 0
        	mat_mat2.append(mat_mat1)
		mat_mat1 = []
            mat_mat3.append(mat_mat2)
	    mat_mat2 = []
        mat.append(mat_mat3)
	mat_mat3 = []
    for u in range(size_u):
        for v in range(size_u):
            for w in range(size_u):
		for x in range(size_u):
    	    	    for m in range(size_a):
        		for n in range(size_a):
			    for o in range(size_a):
	    	    	        matmv_ele+= gamma_mv[u][w][m][o]*gamma_mv[v][x][o][n] 
			    if v == w:
			        matmv_ele-=gamma_mv[u][x][m][n]
			    matmv_mat1.append(matmv_ele)
			    matmv_ele = 0.0
		        matmv_mat2.append(matmv_mat1)
	    	        matmv_mat1 = []
        	    matmv_mat3.append(matmv_mat2)
		    matmv_mat2 = []
                matmv_mat4.append(matmv_mat3)
	        matmv_mat3 = []
            matmv_mat5.append(matmv_mat4)
	    matmv_mat4 = []
        matmv.append(matmv_mat5)
	matmv_mat5 = []
    return mat, matmv
def make_lambda2(gamma2, gamma_sin, coeff, size_u, size_a, size_i):
    mat = [[[[0.0 for x in range(size_u)] for x in range(size_u)] for x in range(size_u)] for x in range(size_u)]
    for u in range(size_u):    
        for v in range(size_u):    
    	    for w in range(size_u):    
        	for x in range(size_u):   
		     mat[u][v][w][x] = gamma2[u][v][w][x] - gamma_sin[u][w]*gamma_sin[v][x] + 1.0/2.0*gamma_sin[u][x]*gamma_sin[v][w]
    return mat

def make_gamma3(gamma2_mv, gamma2, gamma_mv, gamma_sin, coeff, size_u, size_a, size_i):
    
    mat_ele = 0
    mat_mat1 = []
    mat_mat2 = []

    mat_mat3 = []
    mat_mat4 = []
    mat_mat5 = []
    mat = []
    for u in range(size_u):
        for v in range(size_u):
            for w in range(size_u):
		for x in range(size_u):
    	    	    for y in range(size_u):
        		for z in range(size_u):
			    for m in range(size_a):
			        for n in range(size_a):
			            for o in range(size_a):
	    	    	        	mat_ele+= coeff[m]*coeff[n]*gamma2_mv[u][v][x][y][m][o]*gamma_mv[w][z][o][n] 
				    if w == y:

				        mat_ele+=coeff[m]*coeff[n]*gamma2_mv[u][v][z][x][m][n]
				    if w == x:
				        mat_ele-=coeff[m]*coeff[n]*gamma2_mv[u][v][z][y][m][n]
			    mat_mat1.append(mat_ele)
			    mat_ele = 0.0
		        mat_mat2.append(mat_mat1)
	    	        mat_mat1 = []

        	    mat_mat3.append(mat_mat2)
		    mat_mat2 = []


                mat_mat4.append(mat_mat3)
	        mat_mat3 = []
            mat_mat5.append(mat_mat4)
	    mat_mat4 = []
	
        mat.append(mat_mat5)
	mat_mat5 = []
    return mat
def make_lambda3(gamma3, gamma_sin, lambda2, coeff, size_u, size_a, size_i):
    mat = [[[[[[0.0 for x in range(size_u)] for x in range(size_u)] for x in range(size_u)] for x in range(size_u)] for x in range(size_u)] for x in range(size_u)]
    for u in range(size_u):    
        for v in range(size_u):    
    	    for w in range(size_u):    
        	for x in range(size_u):   
        	    for y in range(size_u):   
        		for z in range(size_u):   
		     	    mat[u][v][w][x][y][z]+= gamma3[u][v][w][x][y][z] - gamma_sin[u][x]*lambda2[v][w][y][z] - gamma_sin[v][y]*lambda2[u][w][x][z] - gamma_sin[w][z]*lambda2[u][v][x][y] 

			    mat[u][v][w][x][y][z]+= 1.0/2.0*gamma_sin[u][y]*lambda2[v][w][x][z] - 1.0/2.0*gamma_sin[u][z]*lambda2[v][w][x][y] + 1.0/2.0*gamma_sin[v][x]*lambda2[u][w][y][z] + 1.0/2.0*gamma_sin[v][z]*lambda2[u][w][x][y] - 1.0/2.0*gamma_sin[w][x]*lambda2[u][v][y][z] + 1.0/2.0*gamma_sin[w][y]*lambda2[u][v][x][z] 

			    mat[u][v][w][x][y][z]-=gamma_sin[u][x]*gamma_sin[v][y]*gamma_sin[w][z] 

			    mat[u][v][w][x][y][z]+= 1.0/2.0*gamma_sin[u][x]*gamma_sin[v][z]*gamma_sin[w][y] + 1.0/2.0*gamma_sin[u][z]*gamma_sin[v][y]*gamma_sin[w][x] + 1.0/2.0*gamma_sin[u][y]*gamma_sin[v][x]*gamma_sin[w][z] 

			    mat[u][v][w][x][y][z]+= -1.0/4.0*gamma_sin[u][y]*gamma_sin[v][z]*gamma_sin[w][x] - 1.0/4.0*gamma_sin[u][z]*gamma_sin[v][x]*gamma_sin[w][y]

    return mat
