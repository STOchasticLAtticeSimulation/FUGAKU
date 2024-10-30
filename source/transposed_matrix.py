import os
import numpy as np
import time

start_time=time.time()
#directory
#input_file_path='/data/STOLAS/Exercise/source/noisedata/noisemap_0.dat'
#output_dir='/data/STOLAS/Exercise/FUGAKU/practice_matrix'
input_file_path='noisedata/noisemap_256.dat'
output_dir='noisedata/noisetrs_256/'
num_chunk = 16


if not os.path.exists(output_dir):
       os.makedirs(output_dir)

try:
       matrix=np.loadtxt(input_file_path)
       if matrix.size==0:
              raise ValueError("入力ファイルは空です。")
except Exception as e:
       print(f"ファイルの読み込みに失敗しました:{e}")
       raise

rows,cols=matrix.shape
print(f"元の行列のサイズ: {rows} x {cols}")

transposed_matrix=matrix.T
#print(f"転置した行列のサイズ:{transposed_matrix.shape}")
#print("転置行列の一部:", transposed_matrix[:5, :5])


chunk_size=transposed_matrix.shape[0]//num_chunk #転置後の行数で分割
print(f"分割サイズ(行数):{chunk_size}")

for i in range(num_chunk):

    start_row=i*chunk_size
    end_row=start_row+chunk_size if i<(num_chunk-1) else transposed_matrix.shape[0]

    part_matrix =transposed_matrix[start_row:end_row,:]
    
    if part_matrix.size==0:
       print(f"警告: 分割行列 part_matrix_{i} は空です。")

    output_file_path=os.path.join(output_dir,f'part_{i}.dat')
    
    np.savetxt(output_file_path,part_matrix)
   # print(f"分割した行列の部分を {output_file_path} に保存しました。")

print("すべての分割行列が保存されました。")

end_time=time.time()

execution_time = end_time - start_time
print(f"処理にかかった時間: {execution_time:.2f} 秒")


















