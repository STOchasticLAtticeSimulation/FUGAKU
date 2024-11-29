import numpy as np
import os
import math
import time
import sys


import re
# ファイルのパス
parameter_file_path = 'parameters.hpp'

# 読み取りたい変数名totalnoiseNo
NL_name = 'NLnoise'
chunk_name = 'totalnoiseNo'
value = None

# ファイルを開く
with open(parameter_file_path, 'r') as file:
    for line in file:
        count=0
        # 'totalnoiseNo' の定義行を検索
        match = re.search(rf'{NL_name}\s*=\s*pow\((\d+),(\d+)\)', line)
        if match:
            base = int(match.group(1))
            exponent = int(match.group(2))
            NL = base ** exponent
            count += 1

        match = re.search(rf'{chunk_name}\s*=\s*pow\((\d+),(\d+)\)', line)
        if match:
            base = int(match.group(1))
            exponent = int(match.group(2))
            num_chunks = base ** exponent
            count += 1
        
        if count == 2:
            break


# 実行時間の開始時間を記録
start_time = time.time()

# 1. 行列のサイズやデータ型を設定します
# NL = 2**5
print(NL)
print(num_chunks)
dN = 0.01
sigma = 0.1

rows, cols = math.ceil(math.log((NL/2-1)/sigma)/dN), NL**3  # 行数と列数
dtype = 'float64' # データ型

# 出力ファイルのチャンクサイズを指定（例：5000列ごとに分割）
# num_chunks = 2**3
chunk_size = int(cols / num_chunks)

# 2. メモリマップを使用して元データを読み込む
input_file = f'biasdata/biasmap.bin'
output_dir = f'biasdata/biastrs/'

if not(os.path.exists(input_file)):
    print("The noise file couldn't be opened.")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

input_array = np.memmap(input_file, dtype=dtype, mode='r', shape=(rows, cols))

# 各チャンクを転置してファイルに保存
for chunk in range(num_chunks):
    output_file = f'{output_dir}part_{chunk}.bin'

    # 列範囲を計算（最後のチャンクは列数に満たない場合がある）
    start_col = chunk * chunk_size
    end_col = start_col + chunk_size

    if os.path.exists(output_file):
        os.remove(output_file)  # 既存のファイルがあれば削除
    
    # 転置後の形状に合わせてメモリマップファイルを作成
    output_array = np.memmap(output_file, dtype=dtype, mode='w+', shape=(chunk_size, rows))

    # 転置してチャンクごとに保存
    for i in range(rows):
        output_array[:, i] = input_array[i, start_col:end_col]  # 行ごとに転置して列範囲に書き込み

    # データをディスクに確定
    output_array.flush()
    print(f"Transposed chunk {chunk} saved to {output_file}")


# 実行時間の終了時間を記録し、計測
end_time = time.time()
execution_time = end_time - start_time

print(f"Execution time: {execution_time:.2f} seconds")
