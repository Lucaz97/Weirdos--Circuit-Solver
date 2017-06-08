####################################
# # # #Weirdos Circuit Solver# # # #
#
#
####################################

import numpy
import sys
import matplotlib.pyplot as plt

#######################
# # # #FUNCTIONS# # # #
#######################
def checkNetList(file_name): #check netlist status
    knownSymbols = ["I", "R", "V", "VCCS", "VCVS", "CCVS", "CCCS", "C", "L", "OC", "SC", "AO", "DB"]
    with open(file_name, "r") as f:
        data = f.readlines()
    # read all lines and store them
    lines = []
    for line in data:
        lines.append(line.split())
    count= 0
    for lati in lines:
        #check if all entries are valid
        if lati[0] == "C" or lati[0] == "L":
            dynamic_components.append(count)
            dyn_resV.append([])
            dyn_resC.append([])
        if lati[0] not in knownSymbols:
            return 0
        count+=1

    return 1
    f.close()


def checkData(matrix):  #check if input data is correct by analyzing the corresponding matrices
    # check data consistency (this function should be called passing the full incidence matrix)
    absM = numpy.absolute(matrix)
    (l, n) = matrix.shape
    #check sum of numbers in columns SHOULD BE 0
    rowSum =  numpy.sum(matrix, axis=0)
    for ints in rowSum:
        if ints != 0.0:
            return 0
    #check hanging and isolated nodes
    colSum = numpy.sum(absM , axis= 1)
    for ints in colSum:
        if ints <=1:
            return 0
    return 1
    #TODO:check block matrix

    #check controlled sources


#######################
# # # #   END   # # # #
#######################
has_dynamic_components = 0
dynamic_components = []
dyn_resV = []
dyn_resC = []
h = 0.000001
inputFile = "input.txt"  # open input file
#check netlist
if checkNetList(inputFile) == 0:
    sys.exit("ERROR: the netlist is not valid!")

with open(inputFile, "r") as f:
    data = f.readlines()
# read all lines and store them
f.close()
lines = []
for line in data:
    lines.append(line.split())

# calculate number of edges of the circuit and associated graph
# calculate number of nodes in the circuit and associated graph
nodi = []
l = 0
for element in lines:
    if int(element[2]) not in nodi: #add node if not already in the list
        nodi.append(int(element[2]))
    if int(element[3]) not in nodi:
        nodi.append(int(element[3]))
    if element[0] == "VCCS" or element[0] == "VCVS" or element[0] == "CCCS" or  element[0] == "CCVS" or element[0] ==\
            "C" or element[0] == "L" or element[0] == "AO" or element[0] == "DB":
        l +=2 # add two edges if is a controlled source or  inductor/conductor
    else:
        l +=1 # add one edge for dipoles

n = len(nodi) -1 #minus one because this value is mostly used for matrix dimensions n-1
print "numero di nodi: ", n+1
print "numero di lati: ", l
# node zero must be present
if 0 not in nodi:
    print nodi
    sys.exit('Input file is invalid, node zero must be defined!')

MatrixA = numpy.zeros((n + 1, l)) # incidence matrix
MatrixM = numpy.zeros((l, l))     # voltages matrix
MatrixN = numpy.zeros((l, l))     # currents matrix
VectorZ = numpy.zeros((l, 1))     # known terms vector

#fill Matrices and vector
count = 0
for lati in lines:
    #build matrix A
    MatrixA[int(lati[2])][count] = 1.0
    MatrixA[int(lati[3])][count] = -1.0
    # matrix M andN if element is a resistor
    if lati[0] == 'R':
        MatrixM[count][count] = 1
        MatrixN[count][count] = -1.0 * float(lati[4])
    # if element is an independent voltage source
    elif lati[0] == 'V':
        MatrixM[count][count] = 1
        VectorZ[count] = float(lati[4])
    #if element is an independent current source
    elif lati[0] == 'I':
        MatrixN[count][count] = 1
        VectorZ[count] = float(lati[4])
    #if element is voltage controlled current source:
    elif lati[0] == 'VCCS':
        MatrixN[count][count] = 1
        MatrixM[count][count+1] = -float(lati[6])
        MatrixN[count+1][count+1] = 1
        # remember to add the controlling edge to matrix A!
        MatrixA[int(lati[4])][count + 1] = 1
        MatrixA[int(lati[5])][count + 1] = -1
    #if element is current controlled current source
    elif lati[0] == 'CCCS':
        MatrixN[count][count] = 1
        MatrixN[count][count+1] = -float(lati[6])
        MatrixM[count+1][count+1] = 1
        # remember to add the controlling edge to matrix A!
        MatrixA[int(lati[4])][count + 1] = 1
        MatrixA[int(lati[5])][count + 1] = -1
    #if element is current controlled voltage source
    elif lati[0] == 'CCVS':
        MatrixM[count][count] = 1
        MatrixN[count][count+1] = -float(lati[6])
        MatrixM[count+1][count+1] = 1
        # remember to add the controlling edge to matrix A!
        MatrixA[int(lati[4])][count + 1] = 1
        MatrixA[int(lati[5])][count + 1] = -1
    #if element is voltage controlled voltage source
    elif lati[0] == 'VCVS':
        MatrixM[count][count] = 1
        MatrixM[count][count+1] = -float(lati[6])
        MatrixN[count+1][count+1] =1
        #remember to add the controlling edge to matrix A!
        MatrixA[int(lati[4])][count + 1] = 1
        MatrixA[int(lati[5])][count + 1] = -1
    elif lati[0] == 'SC':   #if short circuit
        MatrixM[count][count] = 1
    elif lati[0] == 'OC':    #if open circuit
        MatrixN[count][count] = 1
    # if element is capacitor
    elif lati[0] == "C":
        has_dynamic_components = 1
        #change with norton equivalent
        #resistor:
        MatrixM[count][count] = 1
        MatrixN[count][count] = float(h/float(lati[4]))
        #current source :
        MatrixN[count+1][count+1] = 1
        VectorZ[count+1] = -1.0*float((float(lati[4])/h)*float(lati[5]))
        #update matrix A
        MatrixA[int(lati[2])][count + 1] = 1
        MatrixA[int(lati[3])][count + 1] = -1
    #if element is inductor
    elif lati[0] == "L":
        has_dynamic_components = 1
        # change with norton equivalent
        # resistor:
        MatrixM[count][count] = 1
        MatrixN[count][count] = float(float(lati[4])/h)
        # current source :
        MatrixN[count+1][count+1] = 1
        VectorZ[count+1] = -1.0*float(lati[5])
        # update matrixA
        MatrixA[int(lati[4])][count + 1] = 1
        MatrixA[int(lati[3])][count + 1] = -1
    #if element is generic double dipole
    if lati[0] == "DB":
        #update matrix A
        MatrixA[int(lati[4])][count + 1] = 1
        MatrixA[int(lati[5])][count + 1] = -1
        if lati[6] == 'G': #if is represented by matrix G
            MatrixN[count][count] = -1
            MatrixN[count + 1][count + 1] = -1
            MatrixM[count][count] = float(lati[7])
            MatrixM[count][count + 1] = float(lati[8])
            MatrixM[count +1][count] = float(lati[9])
            MatrixM[count +1][count +1] = float(lati[10])
    if lati[0] == "VCCS" or lati[0] == "VCVS"or lati[0] == "CCVS"or lati[0] == "CCCS" or lati[0] == "L" or lati[0] ==\
            "L" or lati[0] == "AO" or lati[0] == "DB":
        #skip a line
        #already used for controlling edge of controlled source or second edge of element
        count += 1
    count += 1

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
Matrix2 = numpy.concatenate((numpy.negative(MatrixAr.transpose()), numpy.identity(l), numpy.zeros((l,l))), axis = 1)
Matrix3 = numpy.concatenate((numpy.zeros((l,n)), MatrixM, MatrixN), axis = 1)
#Merge the three blocks vertically to build matrix T
MatrixT = numpy.concatenate((Matrix1, Matrix2, Matrix3), axis=0)
#Build vector B
VectorB = numpy.concatenate((numpy.zeros((n+l,1)), VectorZ ), axis=0)

print "matrice T: "
print MatrixT

#find solution of the circuit
solution = numpy.linalg.solve(MatrixT, VectorB)
#check solution
if numpy.allclose(numpy.dot(MatrixT, solution), VectorB) == 0:
    sys.exit('Error occurred while solving the circuit. Program interrupted.')
#print
print "soluzione :"
print solution
if has_dynamic_components == 0:
    #print results TODO: print results for controlled sources
    count = 0
    lato = 0
    print "potenziali ai nodi con nodo di riferimeno nodo 0"
    while count < (n+l+l):
        if count< n: #stampa potenziali ai nodi
            print "potenziale del nodo ", count+1, ": " , solution[count]
        if count>= n and count < n+l: #stampa tensioni
            if lines[lato-n][0] == 'R':
                print "tensione ai capi del resistore R" , lines[lato-n][1], ": ", solution[count]
            elif lines[lato-n][0] == 'V':
                print "tensione ai capi del generatore indipendente V", lines[lato -n][1], ": ", solution[count]
            elif lines[lato - n][0] == 'I':
                print "tensione ai capi del generatore indipendente I", lines[lato - n][1], ": ", solution[count]
            elif lines[lato-n][0] == "VCCS" or lines[lato-n][0] == "VCVS":
                print "tensione ai capi del generatore pilotato ", lines[lato-n][0], lines[lato -n][1], ": ", solution[count]
                count +=1
            elif lines[lato-n][0] == "CCCS" or lines[lato-n][0] == "CCVS":
                print "tensione ai capi del generatore pilotato ", lines[lato - n][0], lines[lato - n][1], ": ", solution[count]
                count +=1
            elif lines[lato - n][0] == "DB":
                print "tensione ai capi del lato 1 del doppio bipolo ", lines[lato - n][0], lines[lato - n][1], ": ", solution[count]
                print "tensione ai capi del lato 2 del doppio bipolo ", lines[lato - n][0], lines[lato - n][1], ": ", solution[count+1]
                count +=1
        #stampa correnti
        if count>= n+l:
            if lines[lato - (n+l)][0] == 'R':
                print "corrente passante per R" , lines[lato - (n+l)][1], ": ", solution[count]
            elif lines[lato - (n+l)][0] == 'V':
                print "corrente passante per V", lines[lato - (n+l)][1], ": ", solution[count]
            elif lines[lato - (n+l)][0] == 'I':
                print "corrente passante per I", lines[lato - (n+l)][1], ": ", solution[count]
            elif lines[lato - (n+l)][0] == "VCCS" or lines[lato - (n+l)][0] == "VCVS":
                print "corrente passante per ", lines[lato - (n+l)][0], lines[lato - (n+l)][1], ": ", solution[count]
                count+=1
            elif lines[lato - (n+l)][0] == "CCCS" or lines[lato - (n+l)][0] == "CCVS":
                print "corrente passante per ", lines[lato - (n+l)][0], lines[lato - (n+l)][1], ": ", solution[count]
                print "corrente pilotante per ", lines[lato - (n+l)][0], lines[lato - (n+l)][1], ": ", solution[count+1]
                count +=1
            elif lines[lato - (n+l)][0] == "DB":
                print "Corrente passante per il lato 1 del doppio bipolo ", lines[lato - (n+l)][0], lines[lato - (n+l)][1], \
                    ": ", solution[count]
                print "corrente passante per il lato 2 del doppio bipolo ", lines[lato -(n+l)][0], lines[lato - (n+l)][1], \
                    ": ", solution[count+1]
                count +=1
        count= count+1
        lato = lato+1

##end of has not dinamic component
elif has_dynamic_components == 1:
    z = 0
    m = 0
    M=0
    while 1 :
        c = 0
        for comp in dynamic_components:
           if lines[comp][0] == "C": #if ccomp is conductor
               dyn_resV[c].append(solution[n+comp])
               dyn_resC[c].append(solution[n+l+comp] + solution[n+l+comp+1])
               VectorZ[comp+1] = float(solution[n + l + comp] + solution[n + l + comp + 1])
               print" somma correnti: " , solution[n+l+comp] + solution[n+l+comp+1]
               print " corrente resistore: ", solution[n+l+comp]
               print "corrente generatore: ", solution[n+l+comp+1]
               print "valore in vettoreZ: " ,VectorZ[comp+1]
           elif lines[comp][0] == "L":
               dyn_resV[c].append(solution[n + comp])
               dyn_resC[c].append(solution[n + l + comp] + solution[n + l + comp + 1])
               VectorZ[comp+1] = float(solution[n + l + comp] + solution[n + l + comp + 1])
               c += 1
        VectorB = numpy.concatenate((numpy.zeros((n + l, 1)), VectorZ), axis=0)
        solution = numpy.linalg.solve(MatrixT, VectorB)
        z += 1
        print z
        print "C: ", solution[n+l+dynamic_components[0]] + solution[n+l+dynamic_components[0]+1 ]
       # print "C resistore: ", solution[n+l+dynamic_components[0]]
        #print "C generatore: ", solution[n+l+dynamic_components[0]+1 ]
        print "V: ", solution[n+dynamic_components[0]]

        if z == 10000:
            break
    plt.figure()
    plt.plot(dyn_resC[0], label="C")
    M = max([M, max(dyn_resC[0])])
    m = min([m, min(dyn_resC[0])])
    plt.ylim([1.1*m, 1.1 * M])
    plt.show()

else :
    print "ERRORE: VALORE SCONOSCIUTO PER VARIABILE HAS_DINAMIC_COMPONENTS"