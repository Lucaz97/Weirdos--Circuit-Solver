import numpy
import sys


#######################
# # # #FUNCTIONS# # # #
#######################

def checkData(matrix):
    #TODO: chech data consistency (this function should be called passing the full incidence matrix)
    return 1

#######################
# # # #   END   # # # #
#######################

inputFile = "input.txt"
# opens input file
with open(inputFile, "r") as f:
    data = f.readlines()
# read all lines and store them
lines = []
for line in data:
    lines.append(line.split())

l = len(lines)  # instantiate number of edges of the circuit and associated graph
# calculate number of nodes in the circuit and associated graph
nodi = []
for element in lines:
    nodi.append(int(element[2]))
    nodi.append(int(element[3]))
n = max(nodi)
# node zero must be present
if not 0 in nodi:
    print nodi
    sys.exit('Input file is invalid, node zero must be defined!')

MatrixA = numpy.zeros((n + 1, l)) # incidence matrix
MatrixM = numpy.zeros((l, l))     # voltages matrix
MatrixN = numpy.zeros((l, l))     # currents matrix
VectorZ = numpy.zeros((l, 1))     # known terms vector

#feel Matrices and vector
count = 0
for lati in lines:
    MatrixA[int(lati[2])][count] = 1.0
    MatrixA[int(lati[3])][count] = -1.0
    # matrix M andN if element is a resistor
    if lati[0] == 'R':
        MatrixM[count][count] = 1
        MatrixN[count][count] = -1 * int(lati[4])
    # if element is an independent voltage source
    elif lati[0] == 'V':
        MatrixM[count][count] = 1
        VectorZ[count] = int(lati[4])
    #if element is an independent current source
    elif lati[0] == 'I':
        MatrixN[count][count] = 1
        VectorZ[count] = int(lati[4])
    count = count + 1

#debug purpose only
print "matrice A: "
print MatrixA
print "matrice M (tensioni): "
print MatrixM
print "matrice N (correnti): "
print MatrixN
print "vettore Z (termini noti): "
print VectorZ

#check data consistency
if checkData(MatrixA) == 0:
    sys.exit('Invalid input file, program interrupted.')

#reduced incidence matrix
MatrixAr = numpy.delete(MatrixA, 0, 0)
#build three horizontal blocks of the matrix T
Matrix1 = numpy.concatenate((numpy.zeros((n,l+n)), MatrixAr), axis=1)
Matrix2 = numpy.concatenate((MatrixAr.transpose(), numpy.identity(l), numpy.zeros((l,l))), axis = 1)
Matrix3 = numpy.concatenate((numpy.zeros((l,n)), MatrixM, MatrixN), axis = 1)
#Merge the three blocks vertically to build matrix T
MatrixT = numpy.concatenate((Matrix1, Matrix2, Matrix3), axis=0)
#Build vector B
VectorB = numpy.concatenate((numpy.zeros((n+l,1)), VectorZ ), axis=0)

print "matrice T: "
print MatrixT

#find solution of the circuit
solution = numpy.linalg.solve(MatrixT, VectorB)

print "soluzione :"
print solution

if numpy.allclose(numpy.dot(MatrixT, solution), VectorB) == 0:
    sys.exit('Error occurred while solving the circuit. Program interrupted.')

#print results (RESULTS SHOULD BE SHOWN IN A BETTER WAY)
count = 0
print "potenziali ai nodi con nodo di riferimeno nodo 0"
for res in solution:
    if count< n: #stampa potenziali ai nodi
        print "potenziale del nodo ", count+1, ": " , res
    if count>= n and count < n+l: #stampa tensioni
        if lines[count-n][0] == 'R':
            print "tensione ai capi del resistore R" , lines[count-n][1], ": ", res
        elif lines[count-n][0] == 'V':
            print "tensione ai capi del generatore V", lines[count -n][1], ": ", res
        elif lines[count - n][0] == 'I':
            print "tensione ai capi del generatore I", lines[count - n][1], ": ", res
    #stampa correnti
    if count>= n+l:
        if lines[count-(n+l)][0] == 'R':
            print "corrente passante per R" , lines[count - (n+l)][1], ": ", res
        elif lines[count-(n+l)][0] == 'V':
            print "corrente passante per V", lines[count - (n+l)][1], ": ", res
        elif lines[count - (n+l)][0] == 'I':
            print "corrente passante per I", lines[count - (n+l)][1], ": ", res
    count= count+1
