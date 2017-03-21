import copy
import ewt
import func_ewt
import numpy as np
from collections import deque
def single_excited_overlap(holes, active, particles, matrix_no, gamma_sin, eta_sin, lambda2):
    #for i_a ia case :
    classmat= []
    if matrix_no==2 or matrix_no==3:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	if matrix_no == 3:
	    if particles:
	        orb1 = copy.deepcopy(particles[0])
	        orb2 = copy.deepcopy(particles[0])
	    else :
		print "matrix type 2 cannot be formed : particles or holes not present \n"
		return 0
	elif matrix_no == 2:
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
	    if matrix_no == 3:
	        item1.dag = '1'
	        orb1.dag = '0'
	        string1 = [item1, orb1]
	        classmat.append([orb1, item1])
	    if matrix_no == 2:
	        item1.dag = '0'
	        orb1.dag = '1'
	        string1 = [orb1, item1]
	        classmat.append([item1, orb1])
	    for item in string1:
		item.string = 1
		item.pos = p1
		p1+=1

	    for item2 in active_n:
		p2 = p1
		full_con = []
		#the two cases of i->u or a->u
	    	if matrix_no == 3:
	            item2.dag = '0'
	            orb2.dag = '1'
	            string2 = [orb2, item2]
	        if matrix_no == 2:
	            item2.dag = '1'
	            orb2.dag = '0'
	            string2 = [item2, orb2]
		for item in string2:
		    item.string = 2	
		    item.pos = p2
		    p2+=1
		full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
		#print full_con, const_con
		full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
		matrixtmp.append(full_val)
            matrix.append(copy.deepcopy(matrixtmp))
	    #make the triangle matrix : multiply the terms of the same term, add them to the next list by multiplying to the sign and constant
    return matrix, classmat
def make_dagger(op1,op2,op3,op4,p,string_no):
    
    op1.spin = 0
    op2.spin = 0
    op3.spin = 0
    op4.spin = 0
    op1.string = string_no
    op2.string = string_no
    op3.string = string_no
    op4.string = string_no
    op1.pos = p+1
    op2.pos = p+2
    op3.pos = p+3
    op4.pos = p+4
    op1.dag = '1'
    op2.dag = '1'
    op3.dag = '0'
    op4.dag = '0'
    return op1,op2,op3,op4

def iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, typ):
    matrix = []

    classmat = []
    a = copy.deepcopy(orb1)
    b = copy.deepcopy(orb2)
    c = copy.deepcopy(orb3)
    d = copy.deepcopy(orb4)
    a1 = copy.deepcopy(orb1)
    b1 = copy.deepcopy(orb2)
    c1 = copy.deepcopy(orb3)
    d1 = copy.deepcopy(orb4)
    if (typ==1):




        for index1 in range(len(a)):
	    for index2 in range(index1, len(b)):
	        for item3 in c:
		    for item4 in d:

			matrixtmp = [] #important
			item4,item3,item2,item1 = make_dagger(item4,item3,b[index2],a[index1],0,1)
			string1 = [item4,item3,item2,item1]

			classmat.append([item1, item2, item3, item4])
		        for index5 in range(len(a1)):
			    for index6 in range(index5, len(b1)):
			        for item7 in c1:
				    for item8 in d1:
			                item5,item6,item7,item8 = make_dagger(a1[index5], b1[index6], item7, item8,4,2)
					string2 = [item5,item6,item7,item8]

					full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
			                full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
                			matrixtmp.append(full_val)
                        matrix.append(copy.deepcopy(matrixtmp))
    elif typ==2:


        for item1 in a:
	    for item2 in b:
        	for index3 in range(len(c)):
	    	    for index4 in range(index3, len(d)):
			matrixtmp = [] #important

			item4,item3,item2,item1 = make_dagger(d[index4], c[index3], item2, item1,0,1)
			string1 = [item4,item3,item2,item1]
			classmat.append([item1, item2, item3, item4])
        		for item5 in a1:
	    		    for item6 in b1:
	    			for index7 in range(len(c1)):
				    for index8 in range(index7, len(d1)):
			                item5,item6,item7,item8 = make_dagger(item5, item6, c1[index7],d1[index8],4,2)
					string2 = [item5,item6,item7,item8]
					full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
			                full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
                			matrixtmp.append(full_val)
                        matrix.append(copy.deepcopy(matrixtmp))
    elif typ == 3:

	matrixtmp = [] #important
	item4,item3,item2,item1 = make_dagger(d[0], c[0], b[0],a[0],0,1)
	string1 = [item4,item3,item2,item1]
	classmat.append([item1, item2, item3, item4])
	item5,item6,item7,item8 = make_dagger(a1[0], b1[0], c1[0], d1[0],4,2)
	string2 = [item5,item6,item7,item8]

	full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
        matrixtmp.append(full_val)
	item4,item3,item2,item1 = make_dagger(d[0], c[0], b[0],a[0],0,1)
	string1 = [item4,item3,item2,item1]
	item5,item6,item7,item8 = make_dagger(b1[0], a1[0], c1[0], d1[0],4,2)
	classmat.append([item1, item2, item3, item4])
	string2 = [item5,item6,item7,item8]
	full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
        matrixtmp.append(full_val)
        matrix.append(copy.deepcopy(matrixtmp))
	matrixtmp = []
	item4,item3,item2,item1 = make_dagger(d[0], c[0], a[0],b[0],0,1)
	string1 = [item4,item3,item2,item1]
	classmat.append([item1, item2, item3, item4])
	item5,item6,item7,item8 = make_dagger(a1[0], b1[0], c1[0], d1[0],4,2)
	string2 = [item5,item6,item7,item8]
	full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
        matrixtmp.append(full_val)

	item4,item3,item2,item1 = make_dagger(d[0], c[0], a[0],b[0],0,1)
	string1 = [item4,item3,item2,item1]
	item5,item6,item7,item8 = make_dagger(b1[0], a1[0], c1[0], d1[0],4,2)
	classmat.append([item1, item2, item3, item4])
	string2 = [item5,item6,item7,item8]
	full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
        matrixtmp.append(full_val)
        matrix.append(copy.deepcopy(matrixtmp))

    elif typ ==4:
	flag1=1
        for i1 in range(2):				
	    if flag1 ==1:
		a2 = copy.deepcopy(orb1)
		b2 = copy.deepcopy(orb2)
		flag1=0
	    else :
		a2=copy.deepcopy(orb2)
		b2=copy.deepcopy(orb1)
            for item1 in a2:
	        for item2 in b2:
        	    for item3 in c:
	    	        for item4 in d:
			    flag2=1
			    matrixtmp = [] #important
			    item4,item3,item2,item1 = make_dagger(item4, item3,item2,item1,0,1)
			    string1 = [item4,item3,item2,item1]
			    classmat.append([item1, item2, item3, item4])
			    for i2 in range(2):
				if flag2 == 1:
				    a1 = copy.deepcopy(orb1)
				    b1 = copy.deepcopy(orb2)
				    flag2=0
			        else :
				    a1=copy.deepcopy(orb2)
				    b1=copy.deepcopy(orb1)
        			for item5 in a1:
	    			    for item6 in b1:
	    				for item7 in c1:
				    	    for item8 in d1:
			                	item5,item6,item7,item8 = make_dagger(item5, item6,item7,item8,4,2)
						string2 = [item5,item6,item7,item8]
						full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
			                	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
                				matrixtmp.append(full_val)
                        matrix.append(copy.deepcopy(matrixtmp))
    return matrix, classmat

def iterate_over_4op1(a, b, c, d, gamma_sin, eta_sin, lambda2):
    matrix = []
    classmat = []
    a1 = copy.deepcopy(a)
    b1 = copy.deepcopy(b)
    c1 = copy.deepcopy(c)
    d1 = copy.deepcopy(d)
    for item1 in a:
	for item2 in b:
	    for item3 in c:
		for item4 in d:
			matrixtmp = [] #important
			item4,item3,item2,item1 = make_dagger(item4,item3,item2,item1,0,1)
			string1 = [item4,item3,item2,item1]
			classmat.append([item1, item2, item3, item4])
    			for item5 in a1:
			    for item6 in b1:
	    			for item7 in c1:
				    for item8 in d1:
			                item5,item6,item7,item8 = make_dagger(item5,item6,item7,item8,4,2)
					string2 = [item5,item6,item7,item8]
					#class_mat.append(copy.deepcopy(string2)
					full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
			                full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
                			matrixtmp.append(full_val)
                        matrix.append(copy.deepcopy(matrixtmp))
    return matrix, classmat

def iterate_over_4op(a, b, c, d, gamma_sin, eta_sin, lambda2, lambda3=None):
    if lambda3 is None :
        lambda3=[]
    matrix = []
    classmat = []
    a1 = copy.deepcopy(a)
    b1 = copy.deepcopy(b)
    c1 = copy.deepcopy(c)
    d1 = copy.deepcopy(d)
    for i1 in range(len(a)):
	for i2 in range(len(b)):
	    for i3 in range(len(c)):
		for i4 in range(len(d)):
			matrixtmp = [] #important
			item4,item3,item2,item1 = make_dagger(d[i4],c[i3],b[i2],a[i1],0,1)
			string1 = [item4,item3,item2,item1]
			classmat.append([item1, item2, item3, item4])
    			for i5 in range(len(a1)):
			    for i6 in range(len(b1)):
	    			for i7 in range(len(c1)):
				    for i8 in range(len(d1)):
			                item5,item6,item7,item8 = make_dagger(a1[i5],b1[i6],c1[i7],d1[i8],4,2)
					string2 = [item5,item6,item7,item8]
					#class_mat.append(copy.deepcopy(string2)
					full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
			                full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2, lambda3)
					#print 'val of box', full_val
                			matrixtmp.append(full_val)
                        matrix.append(copy.deepcopy(matrixtmp))
    #print matrix
    return matrix, classmat

def doubly_excited_overlap(holes, actives, particles, matrix_no, gamma_sin, eta_sin, lambda2, lambda3=None):
    if lambda3 is None :
        lambda3=[]
    #simple i->a i->a S matrix with i->a i->a
    if not particles:
	print "matrix type 4 cannot be formed : particles or holes not present \n"
	return 0
    if not holes:
	print "matrix type 4 cannot be formed : particles or holes not present \n"
	return 0
    if not actives:
	print "matrix type 4 cannot be formed : particles or holes not present \n"
	return 0
    matrix = []
    classmat = []
    if matrix_no==41:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	orb1 = [copy.deepcopy(particles[0])]
	if len(particles)>1:
	    orb2 = [copy.deepcopy(particles[1])]
	else:
	    return matrix, classmat
	orb3 = [copy.deepcopy(holes[0])]
	if len(holes)>1:
	    orb4 = [copy.deepcopy(holes[1])]
	else:
	    return matrix, classmat
	matrix, classmat = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2,3)

    if matrix_no==42:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = [copy.deepcopy(holes[0])]
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)

    if matrix_no==43:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	orb1 = deque([copy.deepcopy(particles[0])])
	if len(particles)>1:
	    orb2 = deque([copy.deepcopy(particles[1])])
	else:
	    return matrix, classmat
	orb3 = deque([copy.deepcopy(holes[0])])
	orb4 = deque([copy.deepcopy(holes[0])])
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 51:
	orb1 = copy.deepcopy(actives)
	orb2 = copy.deepcopy(actives)
	orb3 = deque([copy.deepcopy(holes[0])])
	if len(holes)>1:
	    orb4 = deque([copy.deepcopy(holes[1])])
	else:
	    return matrix, classmat
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 52:
	orb1 = copy.deepcopy(actives)
	orb2 = copy.deepcopy(actives)
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = [copy.deepcopy(holes[0])]
	matrix, classmat = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, 1)
    if matrix_no == 61:
	orb1 = [copy.deepcopy(particles[0])]
	if len(particles)>1:
	    orb2 = [copy.deepcopy(particles[1])]
	else:
	    return matrix, classmat
	orb3 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 62:
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix,classmat = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, 2)
    if matrix_no == 7:
	orb3 = [copy.deepcopy(holes[0])]
	orb1 = copy.deepcopy(actives)
	orb2 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, lambda3)
    if matrix_no == 8:
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = copy.deepcopy(actives)
	orb3 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, lambda3)
    if matrix_no == 9:
	orb1 = [copy.deepcopy(particles[0])]
	orb4 = [copy.deepcopy(holes[0])]
	orb2 = copy.deepcopy(actives)
	orb3 = copy.deepcopy(actives)
	matrix, classmat = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, 4)#special case of supermatrix

    if matrix_no == 101:
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	if len(holes)>1:
	    orb4 = [copy.deepcopy(holes[1])]
	else:
	    return matrix, classmat
	orb1 = copy.deepcopy(actives)
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 102:
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = [copy.deepcopy(holes[0])]
	orb1 = copy.deepcopy(actives)
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)

    if matrix_no == 111:
	orb1 = [copy.deepcopy(particles[0])]
	if len(particles)>1:
	    orb2 = [copy.deepcopy(particles[1])]
	else:
	    return matrix, classmat
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = copy.deepcopy(actives)
	matrix,classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 112:
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = copy.deepcopy(actives)
	matrix, classmat = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    return matrix, classmat
def return_truncatedX(a, no):
    x, d, y = np.linalg.svd(a, full_matrices=1, compute_uv=1)
    y = np.transpose(y)
    for i in range(len(d)):
	if d[i]<=0.0001:
	    y=np.delete(y, np.s_[i:len(d)], axis=1)
	    print 'deleting'
	    break
    f_ind=open("output/Xvec.dat", "a")
    f_ind.write('classno:'+str(no)+'\n')
    np.savetxt(f_ind, y, delimiter=' ')
    return y


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
		gamma_mv[int(l)-1][int(k)-1][int(j)-1][int(i)-1] = float(den1)
	    except:
	        flag = 0

		print 'end of gamma_mv at ijkl : ', i, j, k, l
	if flag == 2:
	    try:
		i, co=line.split()
		coeff.append(float(co))
	    except:
		flag = 0
		print 'end of coeff at i = ', i
	if flag == 1:
	    try:
	        i, j, dens = line.split()
	        mat[int(j)-1][int(i)-1] = float(dens)
	    except:
		flag = 0
		print 'end of gamma_sin at ij : ', i, j

	if 'denisty matrix - ' in line:
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
		if (i==j):
		    eta_sin[i][j] = 2.0-eta_sin[i][j] 
		else:
		    eta_sin[i][j] = -eta_sin[i][j]
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
    ele = 0
    mat_mat1 = np.array([])
    mat_mat2 = np.array([])
    mat_mat3 = np.array([])
    mat = np.array([[[[0.0 for x in range(size_u)] for x in range(size_u)] for x in range(size_u)] for x in range(size_u)])
    matmv = np.array([[[[[[0.0 for x in range(size_a)] for x in range(size_a)] for x in range(size_u)] for x in range(size_u)] for x in range(size_u)] for x in range(size_u)])
    matmv_mat1 = np.array([])
    matmv_mat2 = np.array([])
    matmv_mat3 = np.array([])
    matmv_mat4 = np.array([])
    matmv_mat5 = np.array([])
    for u in range(size_u):
        for v in range(size_u):
            for w in range(size_u):
		for x in range(size_u):
    	    	    for m in range(size_a):
        		for n in range(size_a):
			    for o in range(size_a):

	    	    	        mat_ele+= coeff[m]*coeff[n]*gamma_mv[u, w, m, o]*gamma_mv[v, x, o, n] 
			    if v == w:
			        mat_ele-=coeff[m]*coeff[n]*gamma_mv[u,x,m,n]
			    matmv[u,v,w,x,m,n] = mat_ele
			    ele+=mat_ele
			    mat_ele = 0.0
		    #print u,v,w,x,ele
		    mat[u,v,w,x] = ele
	    	    ele = 0
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
				    check1 = 0
			            for o in range(size_a):

	    	    	        	mat_ele+= coeff[m]*coeff[n]*gamma2_mv[u][v][x][y][m][o]*gamma_mv[w][z][o][n]
	    	    	        	#check1+= coeff[m]*coeff[n]*gamma2_mv[u][v][x][y][m][o]*gamma_mv[w][z][o][n]
				     
				    if w == y:
					

				        #check2=coeff[m]*coeff[n]*gamma2_mv[u][v][z][x][m][n]
				        mat_ele-=coeff[m]*coeff[n]*gamma2_mv[u][v][z][x][m][n]
					if abs( coeff[m]*coeff[n]*gamma2_mv[u][v][z][x][m][n])>0.00001 and u==0 and v==0 and w==0 and x==0 and y==0 and z==0: 
					    print '1', coeff[m]*coeff[n]*gamma2_mv[u][v][z][x][m][n], mat_ele
				    if w == x:

				        #check3=coeff[m]*coeff[n]*gamma2_mv[u][v][z][y][m][n]
				        mat_ele-=coeff[m]*coeff[n]*gamma2_mv[u][v][z][y][m][n]
					if abs(  coeff[m]*coeff[n]*gamma2_mv[u][v][z][y][m][n])>0.00001 and u==0 and v==0 and w==0 and x==0 and y==0 and z==0:
					    print '2', coeff[m]*coeff[n]*gamma2_mv[u][v][z][y][m][n], mat_ele
				    #print 'check all the values :1, 2, 3, val :', check1, check2, check3, mat_ele
			    if u == 0 and v == 0 and w == 0 and x == 0 and y == 0 and z == 0:
				print ' 0 0 0 0 0 0 lambda 3 = ', mat_ele
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

			    mat[u][v][w][x][y][z]+= 1.0/2.0*gamma_sin[u][y]*lambda2[v][w][x][z] + 1.0/2.0*gamma_sin[u][z]*lambda2[v][w][x][y] + 1.0/2.0*gamma_sin[v][x]*lambda2[u][w][y][z] + 1.0/2.0*gamma_sin[v][z]*lambda2[u][w][x][y] + 1.0/2.0*gamma_sin[w][x]*lambda2[u][v][y][z] + 1.0/2.0*gamma_sin[w][y]*lambda2[u][v][x][z] 

			    mat[u][v][w][x][y][z]-=gamma_sin[u][x]*gamma_sin[v][y]*gamma_sin[w][z] 

			    mat[u][v][w][x][y][z]+= 1.0/2.0*gamma_sin[u][x]*gamma_sin[v][z]*gamma_sin[w][y] + 1.0/2.0*gamma_sin[u][z]*gamma_sin[v][y]*gamma_sin[w][x] + 1.0/2.0*gamma_sin[u][y]*gamma_sin[v][x]*gamma_sin[w][z] 

			    mat[u][v][w][x][y][z]+= -1.0/4.0*gamma_sin[u][y]*gamma_sin[v][z]*gamma_sin[w][x] - 1.0/4.0*gamma_sin[u][z]*gamma_sin[v][x]*gamma_sin[w][y]

    return mat


def readintegrals(file1, size_vir, size_i, size_u):
    integrals = [[[[0.0 for x in range(size_u+size_i)] for x in range(size_u+size_i)]for x in range(size_u+size_vir)]for x in range(size_u+size_vir)]
    for line in file1.readlines():
        i, j, a, b, val = line.split() 
	integrals[int(b)-1][int(a)-1][int(j)-1][int(i)-1] = float(val)	
    return integrals
def readfock(f_fock,size_vac, size_i, size_u):
    fock = [[0.0 for x in range(size_u+size_i+size_vac)] for x in range(size_u+size_i+size_vac)]
    for line in f_fock.readlines():
	print line
	a, i, v = line.strip().split() 
	fock[int(i)-1][int(a)-1] = float(v)	
    return fock

def readoverlap(overlap):
    alloverlap = []
    temp_overlap = [] 
    start = 0
    print '---in overlap'
    for line in overlap.readlines():
	mat = []
	mat=line.split()
	for item in mat:
	    item = float(item)
	if len(mat)==1 and start==1:
	    alloverlap.append(temp_overlap)
	    temp_overlap = []
	elif len(mat)>1:
	    start = 1
	    temp_overlap.append(mat)
    return np.array(alloverlap)

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

def prepare(a, b):
    l1 = len(a)/2
    l2 = len(b)/2
    for i in range(l1):
	a[i].dag = '1'
    for i in range(l2):
	b[i].dag = '1'
    for i in range(l1, len(a)):
	a[i].dag = '0'
    for i in range(l2, len(b)):
	b[i].dag = '0'
    for i in range(len(a)):
	a[i].string = 1
	a[i].spin = 0
	a[i].pos = i
    for i in range(len(b)):
	b[i].string = 2
	b[i].spin = 0
	b[i].pos = i+len(a)
def ortho_overlap(c1, c2, xi, gamma_sin, eta_sin, lambda2):
    matrix = []
    matrixtmp = []
    for item2 in c2:
	for item1 in c1:
	    #contract all c1 with c2
	    a= []

	    for i in range(len(item2)-1, -1, -1):
		a.append(item2[i])
	    #a= copy.deepcopy(item2)
	    b= copy.deepcopy(item1)
	    print a, b
	    prepare(a,b)

	    for item in a:
		print 'dag', item.dag, 'pos', item.pos, 'string', item.string, 'kind', item.kind
	    for item in b:
		print 'dag', item.dag, 'pos', item.pos, 'string', item.string, 'kind', item.kind


	    full_con, const_con = ewt.ewt(a, b)
	    full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
	    #store them in a list
	    print full_val
	    matrixtmp.append(full_val)
	#append in another list 
	matrix.append(matrixtmp)
	#empty first list
	matrixtmp = []
    matrixStmp = []
    matrixS = []
    

    for j in range(len(matrix)):
        for i in range(len(xi)):
        #take dot from the xth item in list x1
	    print 'check the dmention:\n', xi[i], matrix[j], '\n size\n', len(xi)
	    dot = np.sum(np.dot(xi[i], matrix[j]))
	    print 'dot\n', dot
	    matrixStmp.append(dot)
	#store in another list
        matrixS.append(matrixStmp)
	print 'matrixS\n', matrixS, 'xi\n', xi, 'matrix\n', matrix
	matrixStmp =[]
    #x+1
    #list store in another
    return matrixS
def orthogonalize(a, b, Sa, gamma_sin, eta_sin, lambda2):
    c= []
    print Sa
    y = return_truncatedX(Sa, 0)

    print 'before transpose\n', y
    y = np.transpose(y)
    print 'after transpose\n', y
    #x, d, y = np.linalg.svd(Sa, full_matrices=0, compute_uv=1)
    c = ortho_overlap(a, b, y, gamma_sin, eta_sin, lambda2)
    print '\nC\n', c
    #c = np.transpose(c) 
    print 'check dim \n', c, y 
    c= np.dot(c, y)
    matrix = []
    matrixtmp=[]
    for item1 in a:
	for item2 in b:
	    c1 = []
	    c2= []
	    for i in range(len(item1)-1, -1, -1):
		c1.append(copy.deepcopy(item1[i]))
	    c2 = copy.deepcopy(item2)
	    prepare(c1, c2)
	    full_con, const_con = ewt.ewt(c1,c2)
	    full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
	    #store them in a list
	    matrixtmp.append(full_val)
	#append in another list 
	matrix.append(matrixtmp)
	matrixtmp = []
    #multiply with [SX]
    submat = np.dot(c, matrix)
    print 'submat\n', submat
    return submat

def calculate_F(fock, Clas, S_ov, X_I, size_i, size_u):
    X_I = np.transpose(X_I)
    iterlist  = []

    fock_use1 = []
    fock_use = []
    energy_add = []
    energy = 0.0
    #build fock matrix in the order of use:
    for exc in Clas:
	for obj in exc:
	    if obj.kind == 'ho':
		print 'hole ind = ', obj.name
		iterlist.append(int(obj.name))
	    elif obj.kind == 'ac':
		iterlist.append(int(obj.name)+size_i)
		print 'active ind = ', obj.name
	    elif obj.kind == 'pa':
		iterlist.append(int(obj.name)+size_i+size_u)
		print 'particle ind = ', obj.name
	    else:
		print 'the object in energy calculation doesnt have a kind'
	fock_use1.append(fock[iterlist[0]][iterlist[1]])
	iterlist = []
    #multiply each row of X_I to produce a different matrix used to calculate energy
    for i in range(len(X_I)):
        fock_use.append(fock_use1)      #This looks wrong. Fock_use1 is empty 
    fock_use = np.transpose(fock_use)
    print '\nFock\n', fock_use, '\nX\n', X_I, '\nS\n', S_ov
    energy_mat = np.dot(fock_use, X_I)
    for i in range(len(energy_mat)):
	energy_add.append(np.sum(np.dot(energy_mat[i], S_ov[i])))
    energy += np.sum(energy_add)
    return energy
def calculate_I(integrals, Clas, S_ov, X_I, size_i, size_u):
    X_I = np.transpose(X_I)
    iterlist  = []
    integ_use1 = []
    integ_use = []
    energy_add = []
    energy = 0.0
    #build fock matrix in the order of use:
    for exc in Clas:
	for obj in exc:
	    if obj.kind == 'ho':
		print 'hole ind = ', obj.name
		iterlist.append(int(obj.name))
	    elif obj.kind == 'ac':
		iterlist.append(int(obj.name)+size_i)
		print 'active ind = ', obj.name
	    elif obj.kind == 'pa':
		iterlist.append(int(obj.name)+size_i+size_u)
		print 'particle ind = ', obj.name
	    else:
		print 'the object in energy calculation doesnt have a kind'
	integ_use1.append(integrals[iterlist[0]][iterlist[1]][iterlist[2]][iterlist[3]])
	iterlist = []
    #multiply each row of X_I to produce a different matrix used to calculate energy
    for i in range(len(X_I)):
        integ_use.append(integ_use1)
    integ_use = np.transpose(integ_use)
    print '\ninteg\n', integ_use, '\nX\n', X_I
    energy_mat = np.dot(integ_use, X_I)
    for i in range(len(energy_mat)):
	energy_add.append(np.sum(np.dot(energy_mat[i], S_ov[i])))
    energy += np.sum(energy_add)
    return energy

