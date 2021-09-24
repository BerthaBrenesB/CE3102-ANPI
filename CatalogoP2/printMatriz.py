def matriz2str(matriz):
    cadena = ''
    for i in range(len(matriz)):
        cadena += '['
        for j in range(len(matriz[i])):
            cadena += '{:>4s}'.format(str(matriz[i][j]))
        cadena += ']\n'
    return cadena

A = [[2,-6,12,16],[1,-2,6,6],[-1,3,-3,-7],[0,4,3,-6],[0,4,3,-6]]
s = matriz2str(A)
print(s)
print(len(A))
print(len(A[0]))