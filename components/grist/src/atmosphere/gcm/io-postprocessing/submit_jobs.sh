#!/bin/bash

# 1. 获取命令行输入的第一个参数作为 -n 的值
N_VAL=$1

# 2. 检查是否输入了参数
if [ -z "$N_VAL" ]; then
    echo "错误: 请输入 -n 参数的值。"
    echo "用法示例: $0 128"
    exit 1
fi

# 设置目标目录 (默认为当前目录)
TARGET_DIR="."

echo "开始遍历目录: $TARGET_DIR"
echo "任务核心数 (-n) 设置为: $N_VAL"

# 遍历目录
for folder in "$TARGET_DIR"/*; do
    # 确保是文件夹
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        echo "----------------------------------------"
        echo "进入文件夹: $folder_name"
        
        # 进入文件夹，如果失败则跳过
        cd "$folder" || continue

        # --- 1. 检查并提交 1D 任务 ---
        if [ -f "scalar_1d_list_info.txt" ]; then
            echo "[检测到] scalar_1d_list_info.txt -> 提交 1d_io 任务..."
            bsub -p -b -o 1dio.log -q q_swgbcm -pr -J 1d_io -n "$N_VAL" -cgsp 64 -host_stack 2048 -share_size 15360 -cache_size 128 ../after_process_1d
        else
            echo "[未检测到] scalar_1d_list_info.txt，跳过 1D 任务。"
        fi

        # --- 2. 检查并提交 2D 任务 ---
        if [ -f "scalar_2d_list_info.txt" ]; then
            echo "[检测到] scalar_2d_list_info.txt -> 提交 2d_io 任务..."
            bsub -p -b -o 2dio.log -q q_swgbcm -pr -J 2d_io -n "$N_VAL" -cgsp 64 -host_stack 2048 -share_size 15360 -cache_size 128 ../after_process_2d
        else
            echo "[未检测到] scalar_2d_list_info.txt，跳过 2D 任务。"
        fi

        # --- 3. 检查并提交 3D 任务 ---
        if [ -f "scalar_3d_list_info.txt" ]; then
            echo "[检测到] scalar_3d_list_info.txt -> 提交 3d_io 任务..."
            bsub -p -b -o 3dio.log -q q_swgbcm -pr -J 3d_io -n "$N_VAL" -cgsp 64 -host_stack 2048 -share_size 15360 -cache_size 128 ../after_process_3d
        else
            echo "[未检测到] scalar_3d_list_info.txt，跳过 3D 任务。"
        fi

        # 返回上一级目录
        cd ..
    fi
done

echo "----------------------------------------"
echo "所有处理结束。"
