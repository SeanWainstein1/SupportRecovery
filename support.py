import random as rand
import math
import galois # type: ignore

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

#Generate a sample k-sparse n-dimension vector
def generate_x(n,k, bound = 10):
    x = []
    support = rand.sample(range(n),k)
    for i in range(n):
        if(i in support):
            x.append(rand.uniform(-1*bound,bound))
        else:
            x.append(0)
    return x

# Probabalistically returns union-free matrix A
# param n: number of columns
# param k, l: determine the family of union free, k is sparsity of x
# param a: alpha, the disjoint factor, used for determining number of rows
def generate_A1(n,k,l,a):
    e = math.e
    q = math.ceil((k+l)*(e/a)**2)
    m1 = math.ceil(2/a*(k/l+1)*(math.log(n/(k+1))+e)/math.log(e/a))
    A = [[0 for col in range(n)] for row in range(m1*q)]
    
    #Set A_ij to a random q-dimensional basis vector 
    for j in range(n):
        for i in range(m1):
            k = rand.randint(0, q-1)
            A[i * q + k][j] = 1
    return A

def generate_A(n,k,m,l,a):
    e = math.e
    q = math.ceil((k+l)*(e/a)**2)
    A = [[0 for col in range(n)] for row in range(m)]
    
    #Set A_ij to a random q-dimensional basis vector 
    for j in range(n):
        for i in range(math.floor(m/q)):
            k = rand.randint(0, q-1)
            A[i * q + k][j] = 1
    return A

#Return the family of sets representation of A, A_ij = 1 iff i in Bj
def A_to_B(A):
    B = [set() for _ in range(len(A[0]))]
    for row in range(len(A)):
        for col in range(len(A[row])):
            if A[row][col]:
                B[col].add(row)
    return B

def bpref(i,l):
    if(i==0):
        return 0
    shift = int(math.log2(i))-l+1
    if(shift>0):
        i>>=int(math.log2(i))-l+1
    return i

def k_independent_hash_family(k, num_functions):
    field_size = 2**(math.ceil(math.log2(k)))

    # Create a list to store hash functions
    hash_functions = []

    for _ in range(num_functions):
        # Generate random coefficients for a polynomial of degree k-1
        coefficients = [rand.randint(0, field_size - 1) for _ in range(k)]
        
        # Define the hash function using the generated coefficients
        def hash_function(x, coeffs=coefficients):
            hash_value = 0
            for i, coeff in enumerate(coeffs):
                hash_value += coeff * (x ** i)
                hash_value %= field_size
            return hash_value
        
        # Append the hash functions to the list
        hash_functions.append(hash_function)

    return hash_functions

def generate_identification_matrix(n, epsilon, k):
    n = 2**(math.ceil(math.log2(n)))
    M=[]
    hash_functions = k_independent_hash_family(k,math.ceil(math.sqrt(k/epsilon)))
    n = 2**(math.ceil(math.log2(n)))
    for l in range(int(math.log2(k)),int(math.log2(n))):
        Ml= [[0 for i in range(n)] for row in range(int(math.sqrt(k/epsilon)))]
        for i in range(n):
            Ml[hash_functions[l](bpref(i,l))%int(math.sqrt(k/epsilon))][i]=1
        M+=Ml
    return M

def B_to_A(B,epsilon,k):
    A=[]
    num_rows = math.ceil(math.sqrt(epsilon*k))
    for z in range(len(B)):
        for i in range(num_rows):
            temp = [0 for _ in range(len(B))]
            num1 = 0
            for j in range(len(B)):
                if(B[z][j]):
                    temp[j]=(j+1)**num1
                    num1+=1
            A.append(temp)
    return A

def approximate_support_recovery(y,A,k):
    B = A_to_B(A)
    C = []
    y_support = {i for i, val in enumerate(y) if val}

    for i in range(len(A[0])):
        if(len(B[i] & y_support)>=len(B[i])/2):
            C.append(i)

    #Remove excess indices (false positives) arbitrarily, here remove from beginning of list in order
    return C[max(len(C)-k,0):]

def mult(matrix, vector):
    rows = len(matrix)
    cols = len(matrix[0])
    result = [0] * rows
    for i in range(rows):
        for j in range(cols):
            result[i] += matrix[i][j] * vector[j]
    return result

def inner_prod(a,b):
    count = 0
    for i in range(len(a)):
        count += a[i]*b[i]
    return count

def superset_recovery(y1,y2,A1,B,epsilon,k):
    B1 = A_to_B(A1)
    C = []
    n = len(A1[0])
    C1 = [i for i in range(n)]
    y1_support = {i for i, val in enumerate(y1) if val}
    for i in range(len(C1)):
        if(len(B1[i] & y1_support)>len(B1[i])/2):
            C.append(i)
    C=C[max(len(C)-k,0):]
    nrows = math.ceil(math.sqrt(k*epsilon))
    for row in range(len(B)):
        #print({i for i, val in enumerate(B[row]) if val})
        if(not ({i for i, val in enumerate(B[row]) if val} & set(C))):
            disjoint = True
            for i in range(nrows*row, nrows*(row+1)):
                if(y2[i]):
                    disjoint = False
                    break
            if(disjoint):
                C1 = [C1[i] for i in range(len(C1)) if not B[row][i]]
    return C

n = 2048
k = 20
#print(generate_x(n,k))
epsilon = 0.1
l=epsilon*k/2
a=0.5
#d = math.ceil(math.log(n/k)/epsilon)#int(1/(a*math.log(math.e/a))*(k/l+1)*(math.log(n/(k+l)+1)))#
#m = k*d#int((k+l)*(math.e**2/a**3)*(k/l+1)*(math.log(n/(k+l))+1)/math.log(math.e/a))#
n = 2**(math.ceil(math.log2(n)))
l = math.ceil(math.sqrt(epsilon*k)/2)
A = generate_A(n,k,5*int(k*math.log(n/k)/math.sqrt(epsilon/k)),l,a)
x = generate_x(n,k)
B = generate_identification_matrix(n,epsilon,k)
A1 = B_to_A(B,epsilon,k)
#print(len({i for i, val in enumerate(B[0]) if val}))'
#print(superset_recovery(sign(mult(A,x)),sign(mult(A1,x)),A,B,epsilon,k))
#print(support(x))

supprecovery = approximate_support_recovery(sign(mult(A,x)),A,k)
suppx = support(x)
count = 0
for i in range(len(suppx)):
    if(suppx[i]-supprecovery[i]):
        count+=1
print("Errors in support recovery:", count)
