import os
import numpy as np
import heapq
import psutil
import time

# 一時ファイルを保存するディレクトリ
TEMP_DIR = 'temp_files_numpy'
os.makedirs(TEMP_DIR, exist_ok=True)

# 使用可能なメモリの割合を指定 (例: 総メモリの80%を使用)
MEMORY_LIMIT_RATIO = 0.8

def get_available_memory():
    """使用可能な物理メモリ量を取得（バイト単位）"""
    mem = psutil.virtual_memory()
    return mem.available * MEMORY_LIMIT_RATIO  # 使用可能なメモリの50%を使用

def save_sorted_chunk(chunk, chunk_index):
    """チャンクをソートして一時ファイルに保存"""
    sorted_chunk = chunk[chunk[:, 0].argsort()]  # 1列目でソート
    temp_file_name = os.path.join(TEMP_DIR, f'sorted_chunk_{chunk_index}.dat')
    np.savetxt(temp_file_name, sorted_chunk)
    return temp_file_name

def calculate_chunk_size(data_sample, available_memory):
    """使用可能なメモリに基づいて適切なチャンクサイズを計算"""
    sample_size = data_sample.nbytes  # データサンプルのメモリサイズを取得
    return int(available_memory // sample_size)  # チャンクに収められる行数を計算

def split_and_sort(file_path):
    """大きなファイルをチャンクごとに読み込みソートして保存"""
    temp_files = []
    chunk_index = 0
    available_memory = get_available_memory()  # 使用可能なメモリ量を取得

    with open(file_path, 'r') as file:
        chunk = []
        sample_line = file.readline().split()
        file.seek(0)  # ファイルポインタを元に戻す

        # サンプル行からNumPy配列に変換してメモリ使用量を計算
        data_sample = np.array([float(x) for x in sample_line])
        chunk_size = calculate_chunk_size(data_sample, available_memory)

        for line in file:
            chunk.append([float(x) for x in line.split()])
            if len(chunk) >= chunk_size:
                chunk_np = np.array(chunk)
                temp_file = save_sorted_chunk(chunk_np, chunk_index)
                temp_files.append(temp_file)
                chunk = []
                chunk_index += 1
        if chunk:  # 最後に残った部分を処理
            chunk_np = np.array(chunk)
            temp_file = save_sorted_chunk(chunk_np, chunk_index)
            temp_files.append(temp_file)

    return temp_files

def merge_sorted_files(output_file, temp_files):
    """ソートされたファイルをマージして最終的なファイルに出力"""
    with open(output_file, 'w') as outfile:
        # 各一時ファイルの先頭行を取得し、ヒープに格納
        open_files = [open(temp_file, 'r') for temp_file in temp_files]
        heap = []
        for file in open_files:
            line = file.readline()
            if line:
                data = [float(x) for x in line.split()]
                heapq.heappush(heap, (data[0], data, file))

        while heap:
            smallest, data, file = heapq.heappop(heap)
            outfile.write(" ".join(map(str, data)) + "\n")
            line = file.readline()
            if line:
                data = [float(x) for x in line.split()]
                heapq.heappush(heap, (data[0], data, file))

        # 全ファイルを閉じる
        for file in open_files:
            file.close()

    # 一時ファイルを削除
    for temp_file in temp_files:
        os.remove(temp_file)


def external_sort(input_file, output_file):
    """外部ソート全体の流れを管理する関数"""
    # 計測開始
    start_time = time.time()

    # 1. ファイルをチャンクごとにソートして保存
    temp_files = split_and_sort(input_file)
    
    # 2. ソート済みのチャンクをマージして出力ファイルに書き出す
    merge_sorted_files(output_file, temp_files)

    # 計測終了
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"time: {elapsed_time:.2f}sec.")

# 実行例
input_file = 'noisedata/noisemap_3.dat'
output_file = 'noisedata/sorted_data_3.dat'
external_sort(input_file, output_file)
