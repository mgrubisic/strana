from numpy import array, prod, sin, zeros, sign, dot
from numpy.linalg import norm, eig

def Kronecker(i, j):
        return 1 if i==j else 0

def Levi(*indices):
    return prod([sign(indices[j] - indices[i]) for i in range(len(indices)) for j in range(i+1, len(indices))])

def interior_product(x1, x2):
    return sum([x1[i]*x2[i] for i in range(len(x1))])

def exterior_product(x1, x2):
    dims = range(len(x1))
    return [sum([sum([Levi(i,j,k)*x1[j]*x2[k] for j in dims]) for k in dims]) for i in dims]

def epsilon_delta(*dims):
    epsilons = zeros((dims))
    deltas = zeros((dims))
    for j in range(dims[0]):
        for k in range(dims[1]):
            for m in range(dims[2]):
                for n in range(dims[3]):
                    deltas[j,k,m,n] = Kronecker(j,m)*Kronecker(k,n) - Kronecker(j,n)*Kronecker(k,m)
                    for i in range(dims[0]):
                        epsilons[j,k,m,n] += Levi(i,j,k)*Levi(i,m,n)
    equality = (deltas == epsilons)
    return deltas, epsilons, equality

# def Durelli(tensor):

def stress(tensor, direction):
    direction = array(direction) / norm(direction)
    traction = dot(tensor, direction)
    normal_stress = dot(traction, direction)*direction
    shear_stress = traction - normal_stress
    eigenvalues, eigenvectors = eig(tensor)
    print 'n_j = {}'.format(direction)
    print 't_i = {}, {}'.format(traction, norm(traction))
    print 's_j = {}, {}'.format(normal_stress, norm(normal_stress))
    print 'tau_j = {}, {}'.format(shear_stress, norm(shear_stress))
    print 'T_k = {}'.format(eigenvalues)
    print '{}'.format(eigenvectors)
    return traction, normal_stress, shear_stress, eigenvalues, eigenvectors


ngb = [-2, 0, 4]
neb = [0, -1, 0]
BGE = exterior_product(ngb, neb)
sij = array([[14, 7, -7], [7, 21, 0], [-7, 0, 35]])
sij2 = array([[7,0,-2],[0,5,0],[-2,0,4]])
sij3 = array([[6,5,0],[5, 2*3**.5],[0, 2*3**.5, 0]])
ej1 = stress(sij, [0, 0.5, 3**.5/2])
