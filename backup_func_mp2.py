import copy
import ewt
import func_ewt
import numpy as np
def single_excited_overlap(holes, active, particles, matrix_no, gamma_sin, eta_sin, lambda2):
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
		full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
		matrixtmp.append(full_val)
            matrix.append(copy.deepcopy(matrixtmp))
	    #make the triangle matrix : multiply the terms of the same term, add them to the next list by multiplying to the sign and constant
    return matrix
def make_dagger(op1,op2,op3,op4,p,string_no):
    
    print op1
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
    			for item5 in a1:
			    for item6 in b1:

        			for index7 in range(len(c1)):
	    			    for index8 in range(index7, len(d1)):
			                item5,item6,item7,item8 = make_dagger(item5,item6,c1[index7],d1[index8],4,2)
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
			item3,item4,item1,item2 = make_dagger(c[index3], d[index4],item1,item2,0,1)
			string1 = [item3,item4,item1,item2]
        		for index5 in range(len(a1)):
	    		    for index6 in range(index5, len(b1)):
	    			for item7 in c1:
				    for item8 in d1:
			                item5,item6,item7,item8 = make_dagger(a1[index5], b1[index6],item7,item8,4,2)
					string2 = [item5,item6,item7,item8]
					full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
			                full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
                			matrixtmp.append(full_val)
                        matrix.append(copy.deepcopy(matrixtmp))
    elif typ == 3:

	matrixtmp = [] #important
	item3,item4,item1,item2 = make_dagger(c[0], d[0], a[0],b[0],0,1)
	string1 = [item3,item4,item1,item2]
	item5,item6,item7,item8 = make_dagger(a1[0], b1[0], c1[0], d1[0],4,2)
	string2 = [item5,item6,item7,item8]

	for item in string1:
	    print 'string', item.string, 'kind ', item.kind, 'pos ', item.pos, 'dag ', item.dag
	for item in string2:
	    print 'string', item.string, 'kind ', item.kind, 'pos ', item.pos, 'dag ', item.dag
	full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
        matrixtmp.append(full_val)
	print matrixtmp
	item3,item4,item1,item2 = make_dagger(c[0], d[0], a[0],b[0],0,1)
	string1 = [item3,item4,item1,item2]
	item5,item6,item7,item8 = make_dagger(b1[0], a1[0], c1[0], d1[0],4,2)
	string2 = [item5,item6,item7,item8]
	full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
        matrixtmp.append(full_val)
        matrix.append(copy.deepcopy(matrixtmp))
	matrixtmp = []
	item3,item4,item1,item2 = make_dagger(c[0], d[0], b[0],a[0],0,1)
	string1 = [item3,item4,item1,item2]
	item5,item6,item7,item8 = make_dagger(a1[0], b1[0], c1[0], d1[0],4,2)
	string2 = [item5,item6,item7,item8]
	full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
	full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
        matrixtmp.append(full_val)

	item3,item4,item1,item2 = make_dagger(c[0], d[0], b[0],a[0],0,1)
	string1 = [item3,item4,item1,item2]
	item5,item6,item7,item8 = make_dagger(b1[0], a1[0], c1[0], d1[0],4,2)
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
			    item3,item4,item1,item2 = make_dagger(item3, item4,item1,item2,0,1)
			    string1 = [item3,item4,item1,item2]
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
    return matrix
def iterate_over_4op(a, b, c, d, gamma_sin, eta_sin, lambda2):
    matrix = []
    a1 = copy.deepcopy(a)
    b1 = copy.deepcopy(b)
    c1 = copy.deepcopy(c)
    d1 = copy.deepcopy(d)
    for item1 in a:
	for item2 in b:
	    for item3 in c:
		for item4 in d:
			matrixtmp = [] #important
			item3,item4,item1,item2 = make_dagger(item3,item4,item1,item2,0,1)
			string1 = [item3,item4,item1,item2]
    			for item5 in a1:
			    for item6 in b1:
	    			for item7 in c1:
				    for item8 in d1:
			                item5,item6,item7,item8 = make_dagger(item5,item6,item7,item8,4,2)
					string2 = [item5,item6,item7,item8]
					full_con, const_con = ewt.ewt(copy.deepcopy(string1), copy.deepcopy(string2))
			                full_val = func_ewt.evaluate(full_con, const_con, gamma_sin, eta_sin, lambda2)
                			matrixtmp.append(full_val)
                        matrix.append(copy.deepcopy(matrixtmp))
    return matrix

def doubly_excited_overlap(holes, actives, particles, matrix_no, gamma_sin, eta_sin, lambda2):
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
    if matrix_no==41:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	orb1 = [copy.deepcopy(particles[0])]
	if len(particles)>1:
	    orb2 = [copy.deepcopy(particles[1])]
	else:
	    return matrix
	orb3 = [copy.deepcopy(holes[0])]
	if len(holes)>1:
	    orb4 = [copy.deepcopy(holes[1])]
	else:
	    return matrix
	matrix = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2,3)

    if matrix_no==42:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = [copy.deepcopy(holes[0])]
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)

    if matrix_no==43:
        #take 1 i or 1 a whatever is required :according to the no_matrix with an if condition
	orb1 = [copy.deepcopy(particles[0])]
	if len(particles)>1:
	    orb2 = [copy.deepcopy(particles[1])]
	else:
	    return matrix
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = [copy.deepcopy(holes[0])]
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 51:
	orb1 = copy.deepcopy(actives)
	orb2 = copy.deepcopy(actives)
	orb3 = [copy.deepcopy(holes[0])]
	if len(holes)>1:
	    orb4 = [copy.deepcopy(holes[1])]
	else:
	    return matrix
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 52:
	orb1 = copy.deepcopy(actives)
	orb2 = copy.deepcopy(actives)
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = [copy.deepcopy(holes[0])]
	matrix = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, 1)
    if matrix_no == 61:
	orb1 = [copy.deepcopy(particles[0])]
	if len(particles)>1:
	    orb2 = [copy.deepcopy(particles[1])]
	else:
	    return matrix
	orb3 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 62:
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, 2)
    if matrix_no == 7:
	orb3 = [copy.deepcopy(holes[0])]
	orb1 = copy.deepcopy(actives)
	orb2 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 8:
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = copy.deepcopy(actives)
	orb3 = copy.deepcopy(actives)
	orb4 = copy.deepcopy(actives)
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 9:
	orb1 = [copy.deepcopy(particles[0])]
	orb4 = [copy.deepcopy(holes[0])]
	orb2 = copy.deepcopy(actives)
	orb3 = copy.deepcopy(actives)
	matrix = iterate_over_4op_special(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2, 4)#special case of supermatrix

    if matrix_no == 101:
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	if len(holes)>1:
	    orb4 = [copy.deepcopy(holes[1])]
	else:
	    return matrix
	orb1 = copy.deepcopy(actives)
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 102:
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = [copy.deepcopy(holes[0])]
	orb1 = copy.deepcopy(actives)
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)

    if matrix_no == 111:
	orb1 = [copy.deepcopy(particles[0])]
	if len(particles)>1:
	    orb2 = [copy.deepcopy(particles[1])]
	else:
	    return matrix
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = copy.deepcopy(actives)
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
    if matrix_no == 112:
	orb1 = [copy.deepcopy(particles[0])]
	orb2 = [copy.deepcopy(particles[0])]
	orb3 = [copy.deepcopy(holes[0])]
	orb4 = copy.deepcopy(actives)
	matrix = iterate_over_4op(orb1,orb2,orb3,orb4, gamma_sin, eta_sin, lambda2)
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
		i, co=line.split()
		coeff.append(float(co))
	    except:
		flag = 0
		print 'end of coeff at i = ', i
	if flag == 1:
	    try:
	        i, j, dens = line.split()
	        mat[int(i)-1][int(j)-1] = float(dens)
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
		    print u,v,w,x,ele
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
    integrals = [[[[0.0 for x in range(size_u+size_vir)] for x in range(size_u+size_vir)]for x in range(size_u+size_i)]for x in range(size_u+size_i)]
    for line in file1.readlines():
        i, j, a, b, val = line.split() 
	integrals[int(i)-1][int(j)-1][int(a)-1][int(b)-1] = float(val)	
    return integrals
def readfock(f_fock,size_vir, size_i, size_u):
    fock = [[0.0 for x in range(size_u)] for x in range(size_u)]
    for line in f_fock.readlines():
	print line
	a, i, v = line.strip().split() 
	fock[int(a)-1][int(i)-1] = float(v)	
    return fock


