import numpy as np

#Solve system for isolated cell probabilities
single_cell_probabilities = []
for k in range(10):
    single_cell_probabilities.append(np.random.random()/10.)

single_cell_probabilities = np.array(single_cell_probabilities)

Z = 1./10 * 1./(sum(single_cell_probabilities))
#print Z

single_cell_probabilities = np.log(single_cell_probabilities*Z)
#print single_cell_probabilities

single_cell_equations = np.zeros((10,10), dtype=np.float32)-1.
for k in range(10):
    single_cell_equations[k][k]=1

h = np.linalg.solve(single_cell_equations, single_cell_probabilities)

print "Intrinsic cell property values h_i: ", h

#CHECKED already... may wish to recheck

#*****************************************************************
#Solve system for pairwise probabilities
N = 45
pair_cell_probabilities = np.random.random(N)  #Generates random matrix
#pair_cell_probabilities = (pair_cell_probabilities + pair_cell_probabilities.T)/2 #Generates symmetric prob matrix

#Scale LHS by Z
LHS = pair_cell_probabilities*Z

#Take log of LHS:
LHS = np.log(LHS)

#Compute Sum(h_i * sigma_i) for each pairwise values: i, j
sum_h_sigma = np.zeros(45, dtype=np.float32) #Loop over )
counter=0
for i in range(0,10,1):
    for j in range(i+1,10,1):
        temp=0
        for k in range(10):
            if (k==i) or (k==j):
                temp+=h[k] 
            else:
                temp+=-h[k]
        sum_h_sigma[counter] = temp     #**** THIS IS WRONG!  IT SHOULD BE A SINGLE CONSTANT!!!
        counter+=1
        

#Subtract Sum(h_i * sigma_i) from log of LHS:
LHS = LHS - sum_h_sigma

#Multiply LHS x 2
LHS *= 2

#Initialize RHS = Sum(J_ij * sigma_i * sigma_j
RHS = np.zeros((45,45), dtype=np.float32)

#Equation coefficients: J_ij (45,45) matrix
counter2=0
for k in range(0,10,1):
    for p in range(k+1,10,1):
        vals = []
        for i in range(0,10,1):
            for j in range(i+1,10,1):
                temp=1
                if (i!=k) and (i!=p):
                    temp*=-1
                if (j!=k) and (j!=p):
                    temp*=-1
                vals.append(temp)
                counter+=1
        
        RHS[counter2]=vals
        counter2+=1
        
print RHS[0]
print RHS[1]


#Solve system for J values:
J = np.linalg.solve(RHS, LHS)

#********************
#Check results for cells 0 and 1
Check_LHS = np.log(pair_cell_probabilities[0]*Z)
print Check_LHS
        
#RHS h_sigma sum
check_h_sigma = 0
k = 0
p = 1
for i in range(10):
    if (i == k) or (i == p):
        check_h_sigma+=h[i]
    else:
        check_h_sigma+=-h[i]

#RHS J_ij * sigma_i * sigma_j
check_J = 0
counter=0
for i in range(0,10,1):
    for j in range(i+1,10,1):
        temp=J[counter]
        if (i!=k) and (i!=p):
            temp*=-1
        if (j!=k) and (j!=p):
            temp*=-1
            
        check_J+=temp
        counter+=1
#Pring RHS
print check_h_sigma+check_J*0.5


        
        
        
        






