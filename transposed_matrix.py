
import numpy as np

matrix=np.random.rand(10000,256**3)#np.loadtxt('noisemap_256_' '.dat',delimiter=' ')

transposed_matrix=np.transpose(matrix)

num=10000 #the number of total step
num_rows=num//10
num_cols=256**3


col_chunk_size=256**3
row_chunk_size=num_rows//10

for i in range(0,num_rows,row_chunk_size):

    chunk =transposed_matrix[i:i+row_chunk_size,:col_chunk_size]
    
    filename=f'tranposed_matrix_part_{i//row_chunk_size}.dat'
    np.savetxt(filename,chunk)

np.savetxt('tranposed_matrix.dat',transposed_matrix)



























