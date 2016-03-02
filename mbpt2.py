import ewt
import fix_uv
import func_ewt
import copy
import sys
from collections import deque
import make_c
import func_mp2
fix_temp = fix_uv
func = func_ewt

f = open("tec.txt", "w")
f_density = open("readfile/fort.50", "r")
#...........input for spin free wicks therem
print "\n Many Body Multireference Perturbation Theory\n"

#Read the input orbital indeces Current assumption : 2 i, 2 u, 3 a
#write a function value to read the file and return the value of the density matrix here

size_u=4
size_a=20
size_i=2
active = [str(i) for i in range(size_u)]
holes = [str(i) for i in range(size_i)]
particles = [str(i) for i in range(size_a)]

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


overall_u[0].dag='1'
overall_u[0].string=1
overall_u[0].pos=1
overall_u[1].dag='0'
overall_u[1].string=1
overall_u[1].pos=2
overall_u[2].dag='1'
overall_u[2].string=2
overall_u[2].pos=3
overall_u[3].dag='0'
overall_u[3].string=2
overall_u[3].pos=4


string1 = [overall_u[0], overall_u[1]]
string2 = [overall_u[2], overall_u[3]]
full_con, const_con = ewt.ewt(string1, string2)
print 'printing full contrations', full_con
evaluate(full_con, const_con, 
#this comment is till the end :

#make 1st excitation array of S(u,v)
#Type 1, 2, 3
#Sia_ai=single_exctd_ovrlp(holes, active, particle, 1)
#Sui_iu=func_mp2.single_excited_overlap(overall_i, overall_u, overall_a, 2)
#Sua_au=func_mp2.single_excited_overlap(overall_i, overall_u, overall_a, 3)
'''


#!!!!!!!! I am fixing the commutator to 0 by force here. Remember to remove this line when incorporating the menu !!!!
commutator=0


a = deque([])
i = deque([])
u = deque([])
full = []
full1 = [] #first string
full2 = [] #second string
full_pos = [] #positions of full
p = 0
spin=0

#make 3 lists : a-particle operators, i-hole operators, u-active state operators
for index in range(0, len(string1), 2):
    if (string1[index] == 'o' or string1[index] == 't'):
	x = operator('ac', string1[index+1], p+1, string1[index], 1, -1, 1)
	u.append(x) 
	full1.append(x)
	full.append(x)
	p=p+1
    elif (string1[index] >= 'a') and (string1[index] < 'h') :
	x = operator('pa', string1[index+1], p+1, string1[index], 1, -1, 1)
	a.append(x) 
	full1.append(x)
	full.append(x)
	p=p+1
    elif (string1[index] >= 'i') and (string1[index] < 'u') :
	x = operator('ho', string1[index+1], p+1, string1[index],1, -1, 1)
	i.append(x) 
	full1.append(x)
	full.append(x)
	p=p+1
    elif (string1[index] >='u' and string1[index]<='z'):
	x = operator('ac', string1[index+1], p+1, string1[index], 1, -1, 1)
	u.append(x) 
	full1.append(x)
	full.append(x)
	p=p+1
if string2:
    for index in range(0, len(string2), 2):
        if (string2[index] == 'o' or string2[index] == 't'):
	    x = operator('ac', string2[index+1], p+1, string2[index], 2, -1, 1)
	    u.append(x) 
	    full2.append(x)
	    p=p+1
        elif (string2[index] >= 'a') and (string2[index] < 'h') :
	    x = operator('pa', string2[index+1], p+1, string2[index], 2, -1, 1)
	    a.append(x) 
	    full2.append(x)
	    p=p+1
        elif (string2[index] >= 'i') and (string2[index] < 'u') :
	    x = operator('ho', string2[index+1], p+1, string2[index], 2, -1, 1)
	    i.append(x) 
	    full2.append(x)
	    p=p+1
        elif (string2[index] >='u' and string2[index]<='z'):
	    x = operator('ac', string2[index+1], p+1, string2[index], 2, -1, 1)
	    u.append(x) 
	    full2.append(x)
	    p=p+1

#make list for all possible contractions for any operator
#The commutator here is 1 when menu=3, so 2 set of terms wille be needed (commutator +1)
for i_c in range(commutator+1):
    full = []
    full_pos = []
    store_for_repeat = []
    poss= deque([])
    y = deque([])
    if not i_c:
	full.extend(full1)
	full.extend(full2)
    else :

	for item in full1:
	    item.string=2
	for item in full2:
	    item.string=1
	full.extend(full2)
	full.extend(full1)
    for item in full:
        full_pos.append(item.pos)
    #----------------------------------------Pairing of the operators
    #---------------------------------------Storing the spin of the operators
    is_pair=1
    spin_tracker=0
    if is_pair:
        for index in range(len(full1)/2):
	    full1[index].pair=full1[len(full1)-index-1].pos
	    full1[len(full1)-index-1].pair=full1[index].pos
	    full1[index].spin=spin_tracker
	    full1[len(full1)-index-1].spin=spin_tracker
	    spin_tracker=spin_tracker+1
        for index in range(len(full2)/2):
    	    full2[index].pair=full2[len(full2)-index-1].pos
	    full2[len(full2)-index-1].pair=full2[index].pos
	    full2[index].spin=spin_tracker
	    full2[len(full2)-index-1].spin=spin_tracker
	    spin_tracker=spin_tracker+1

    #-------------------------all the possible contracting operators of each operator is being sored in poss here
    #poss is a matrix 
    if menu=='1':#self normal ordering
        for operator in full:
            y = deque([])
            if operator.kind == 'pa' and operator.dag=='0':
	        for item in a:
	            if operator.pos<item.pos and item.dag=='1':
		        y.append(item)
    
            elif operator.kind == 'ho' and operator.dag=='1':
  	        for item in i:
	            if operator.pos<item.pos and item.dag=='0':
		        y.append(item)
	
            elif operator.kind == 'ac':  #because active states will have eta and gamma
	        for item in u:
	            if operator.pos<item.pos and int(item.dag)!=int(operator.dag):
		        y.append(item)
            poss.append(y) #list of list in dictionary order i.e 1st annhilation -> possible creation then 2nd ...   
    else:
        for operator in full:
            y = deque([])
            if operator.kind == 'pa' and operator.dag=='0':
	        for item in a:
	            if operator.pos<item.pos and item.dag=='1' and operator.string!=item.string:
		        y.append(item)
            elif operator.kind == 'ho' and operator.dag=='1':
  	        for item in i:
	            if operator.pos<item.pos and item.dag=='0' and operator.string!=item.string:
		        y.append(item)
	
            elif operator.kind == 'ac':  #because active states will have eta and gamma
	        for item in u:
	            if operator.pos<item.pos and int(item.dag)!=int(operator.dag) and operator.string!=item.string:
		        y.append(item)
    #if (y): remember that empty strings are also included

            poss.append(y) #list of list in dictionary order i.e 1st annhilation -> possible creation then 2nd ...   
    no = len(full)/2
    contracted = []
    tmp_l = []

    #---------------------The first term of the PDF file is being printed here (not important)

    tmp_l=[]
    tmp_lower=[]
    
    if not i_c:
	if menu == '1' or menu =='2':
	    tmp_l.append("Doing : Normal ordering of String\\\\")
	else :    
	    tmp_l.append("Doing : Commutator expression Generation\\\\")
        tmp_l.append('Here are the operator strings : \[E^{')

	for item in full1:
            if item.dag=='1':
                tmp_upper=item.name
		tmp_l.append(tmp_upper)
        tmp_l.append('}_{')
	for item in full1:
	    if item.dag=='0':
                tmp_lower.append(item.name)

	tmp_lower=tmp_lower[::-1]
	tmp_l=tmp_l+tmp_lower
	tmp_lower=[]
        tmp_l.append('} , ')
	if menu=='2' or menu=='3':
	    tmp_l.append(' E^{')


	    for item in full2:
                if item.dag=='1':
                    tmp_upper=item.name
		    tmp_l.append(tmp_upper)
            tmp_l.append('}_{')
	    for item in full2:
	        if item.dag=='0':
                    tmp_lower.append(item.name)
	    tmp_lower=tmp_lower[::-1]
	    tmp_l=tmp_l+tmp_lower
	    tmp_lower=[]
            tmp_l.append('}')


        tmp_l.append('\]')
        tmp_6 = "Equation : "+'$$'+''.join(tmp_l)+'\\\\'+'$$'+'\n'+"Answer :\n"
        f.write(tmp_6)
    if not i_c and commutator:
        f.write("\nThis is where the first terms start\\\\\n")
    elif commutator : 
        f.write("\nThis is where the second terms start\\\\\n")
    make_c.make_c(len(full), contracted, a, i, u, full, poss, f, store_for_repeat, full_pos, i_c, menu)
print "\n-------------------------------------------------------------------------------------\n     ITS DONE :D Have a look at the tec.txt file !\n     !!CHEERS !!"
'''


