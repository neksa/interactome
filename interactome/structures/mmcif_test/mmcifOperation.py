"""
Code from:
http://mmcif.wwpdb.org/docs/sw-examples/python/html/assemblies.html
"""

from pdbx.reader.PdbxContainers import *

def parseOperationExpression(expression) :
    operations = []
    stops = [ "," , "-" , ")" ]

    currentOp = ""
    i = 1
    
    # Iterate over the operation expression
    while i in range(1, len(expression) - 1):
        pos = i

        # Read an operation
        while expression[pos] not in stops and pos < len(expression) - 1 : 
            pos += 1    
        currentOp = expression[i : pos]

        # Handle single operations
        if expression[pos] != "-" :
            operations.append(currentOp)
            i = pos

        # Handle ranges
        if expression[pos] == "-" :
            pos += 1
            i = pos
            
            # Read in the range's end value
            while expression[pos] not in stops :
                pos += 1
            end = int(expression[i : pos])
            
            # Add all the operations in [currentOp, end]
            for val in range((int(currentOp)), end + 1) :
                operations.append(str(val))
            i = pos
        i += 1
    return operations


def prepareOperation(oper_list, op1index, op2index) :
    # Prepare matrices for operations 1 & 2
    op1 = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 1]]
    op2 = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 1]]
    
    # Fill the operation matrices for operations 1 & 2
    for i in range(3) :
        op1[i][3] = float(oper_list.getValue("vector[" + str(i + 1) + "]", op1index))
        
        if (op2index != -1) :
            op2[i][3] = float(oper_list.getValue("vector[" + str(i + 1) + "]", op2index))
        for j in range(3) :
            op1[i][j] = float(oper_list.getValue("matrix[" + str(i + 1) + "][" + str(j + 1) + "]", op1index))
            if (op2index != -1) :
                op2[i][j] = float(oper_list.getValue("matrix[" + str(i + 1) + "][" + str(j + 1) + "]", op2index))
    
    # Handles non-Cartesian product expressions
    if (op2index == -1) :
        return op1

    operation = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 1]]

    # Handles Cartesian product expressions (4x4 matrix multiplication)
    sum = 0.0
    for row in range(4) :
        for col in range(4) :
            sum = 0.0
            for r in range(4) :
                sum += (op1[row][r] * op2[r][col])
            operation[row][col] = sum
    return operation
