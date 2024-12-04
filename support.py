import random as rand
import math

def sign(vector):
    for row in range(len(vector)):
        if(vector[row]>0):
            vector[row]=1
        elif(vector[row]):
            vector[row]=-1
        else:
            vector[row]=0
    return vector

def support(vector):
    support = []
    for i in range(len(vector)):
        if(vector[i]):
            support.append(i)
    return support

def generate_x(n,k, bound = 10):
    x = []
    support = rand.sample(range(n),k)
    for i in range(n):
        if(i in support):
            x.append(rand.uniform(-1*bound,bound))
        else:
            x.append(0)
    return x

def generate_A(n,k,l,a):
    e = math.e
    q = math.ceil((k+l)*(e/a)**2)
    m1 = math.ceil(2/a*(k/l+1)*(math.log(n/(k+1))+e)/math.log(e/a))
    A = [[0 for col in range(n)] for row in range(m1*q)]
    
    # Directly construct A
    for j in range(n):  # Iterate over each column
        for i in range(m1):  # Iterate over each effective row
            k = rand.randint(0, q-1)  # Generate a random value in {0, ..., q-1}
            A[i * q + k][j] = 1  # Set the corresponding row to 1 in column j

    return A

def A_to_B(A):
    B = [set() for _ in range(len(A[0]))]  # Initialize B as a list of sets, one for each column
    for row in range(len(A)):
        for col in range(len(A[row])):
            if A[row][col]:
                B[col].add(row)
    
    return B
    '''
    B = [[] for row in range(len(A))]
    for row in range(len(A)):
        for col in range(len(A[row])):
            if(A[row][col]):
                B[col].append(row+1)
    return B
    '''

def approximate_support_recovery(y,A,k):
    B = A_to_B(A)
    C = []
    y_support = {i for i, val in enumerate(y) if val}
    #print(y_support)
    # Compute intersection
    for i in range(len(A[0])):
        count = len(B[i] & y_support)
        if(count>=len(B[i])/2):
            C.append(i)
    if(len(C)-k>0):
        return C[len(C)-k:]
    return C

def mult(matrix, vector):
    rows = len(matrix)
    cols = len(matrix[0])
    result = [0] * rows
    # Perform the matrix-vector multiplication
    for i in range(rows):
        for j in range(cols):
            result[i] += matrix[i][j] * vector[j]
    
    return result


n = 100
k = 15
#print(generate_x(n,k))
epsilon = .1
l=epsilon*k/2
a=0.5
#d = math.ceil(math.log(n/k)/epsilon)#int(1/(a*math.log(math.e/a))*(k/l+1)*(math.log(n/(k+l)+1)))#
#m = k*d#int((k+l)*(math.e**2/a**3)*(k/l+1)*(math.log(n/(k+l))+1)/math.log(math.e/a))#
A = generate_A(n,k,l,a)
x = generate_x(n,k)
print(approximate_support_recovery(sign(mult(A,x)),A,k))
print(support(x))