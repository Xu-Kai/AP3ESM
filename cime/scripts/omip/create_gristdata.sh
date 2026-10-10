#!/bin/bash

# 定义目录A的路径
DIR_A="/home/export/online1/mdt00/shisuan/swgbcm/xukai/GB25/CESM2.2/GRIST_DATA/g8"

# 获取当前所在目录
DIR_B="$(pwd)"

# 检查目录A是否存在
if [ ! -d "$DIR_A" ]; then
  echo "目录 A 不存在: $DIR_A"
  exit 1
fi

# 遍历目录A下的所有文件
for file in "$DIR_A"/*; do
  # 检查是否为文件（排除目录）
  if [ -f "$file" ]; then
    # 获取文件名
    filename=$(basename "$file")
    # 复制文件到当前目录
    #cp "$file" "$DIR_B/$filename"
    # 创建软链接
    ln -sf "$file" "$DIR_B/$filename"
  fi
done

echo "$DIR_A目录数据准备完毕！"
