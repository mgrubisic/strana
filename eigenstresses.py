# computes principal stresses and directions in Voigt's notation
# computes stresses in a given direction

import numpy as np
# user input their stress tensor
# stress_input = raw_input('specify sx sy sz txy txz tyz, in that order and format (spaces between values) \n')
# voigt_list = stress_input.split()

# voigt_stresses = []
# for elem in voigt_list:
#     voigt_stresses.append(float(elem))
# print(voigt_list, voigt_stresses)
alpha = 150*np.pi/180
a = 2
b = 1
c = -2
# 45 deg
# voigt_stresses = [a, c, 0, b-0.5*(a+c), 0, 0]
# 60 deg
voigt_stresses = [a, 1/3*(2*(b+c)-a), 0, (b-c)/np.sqrt(3), 0, 0]
# voigt_stresses = [2, -2, 0, 0 ,0 ,0]
direction = [np.cos(60*np.pi/180), np.sin(60*np.pi/180), np.cos(90*np.pi/180)]
versor = direction / np.linalg.norm(direction)
# voigt_stresses = [200, 400, 300, 300, 100, 200]
Stensor = np.zeros([3,3])
for stress in range(len(voigt_stresses)):
    if stress < 3:
        Stensor[stress,stress] = voigt_stresses[stress]
    else:
        if stress == 3:
            Stensor[0, 1] = voigt_stresses[stress]
            Stensor[1, 0] = Stensor [0, 1]
        if stress == 4:
            Stensor[0, 2] = voigt_stresses[stress]
            Stensor[2, 0] = Stensor [0, 2]
        if stress == 5:
            Stensor[1, 2] = voigt_stresses[stress]
            Stensor[2, 1] = Stensor [1, 2]

# print(Stensor)
# take every string sepparated by comma and make it a list
sdet = np.linalg.det(Stensor)
# print(sdet) 
princ_stresses = sorted(np.linalg.eigvals(Stensor), reverse = True)
print(princ_stresses)
# princ_matrix = np.array([[princ_stresses[0],0,0],[0,princ_stresses[1],0],[0,0,princ_stresses[2]]])
princ_stress_matrix = np.diag(princ_stresses)
# make diagon al matrix of eigenvalues in loop with diag.fill or some fast algorithm that
# takes the [princ_stresses] as a list and proceed to do submatrices for final result

print(princ_stresses, princ_stress_matrix)
eigenvectors = np.zeros([3,3])
for eig in range(len(princ_stresses)):

    eigmatrix = Stensor - princ_stresses[eig]*np.eye(len(princ_stresses))
    A = np.linalg.det(eigmatrix[1:,1:]) 
    B = np.linalg.det(eigmatrix[1:,(0,2)]) 
    C = np.linalg.det(eigmatrix[1:,:2])


    K = 1/(A**2 + B**2 + C**2)**(0.5)
    X = A*K
    Y = B*K
    Z = C*K
    eigenvectors[:, eig] = [X,Y,Z]

    # print(eigmatrix[1:,1:])
    # print(eigmatrix[1:,(0,2)])
    # print(eigmatrix[1:,:2])

print(eigenvectors)
# A is det submatrix s11
# B det submatrix s12
# C det submatrix s13

# k is recip sqrt A+B+C

# principal stresses are eigenvalues of matrix


# for a given prinp-stress

# cosa = AK
# cosb = BK
# cosc = CK

stress_vec = np.dot(Stensor, versor)
stress_magnitude = np.linalg.norm(stress_vec)
nstress = np.dot(stress_vec, versor)
tstress = np.sqrt(stress_magnitude**2-nstress**2)

print('n stress is  ' + str(nstress))
print('t stress is  ' + str(tstress))
