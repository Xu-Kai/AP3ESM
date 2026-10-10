#!/bin/bash

# ================= 配置区域 =================
# 设置源目录 (目录1)
SOURCE_DIR="./"

# 设置目标目录 (目录2)
DEST_DIR="../output2"
# ===========================================

# 检查源目录是否存在
if [ ! -d "$SOURCE_DIR" ]; then
    echo "错误: 源目录 $SOURCE_DIR 不存在"
    exit 1
fi

# 创建目标目录（如果不存在）
mkdir -p "$DEST_DIR"

# 查找源目录下所有的 .nc 文件 (包括子目录)
find "$SOURCE_DIR" -type f -name "*.nc" | while read -r file; do
    # 1. 获取文件名 (不带路径) -> GRIST...nc_cldlow.nc
    filename=$(basename "$file")
    
    # 2. 去掉 .nc 后缀 -> GRIST...nc_cldlow
    filename_no_ext="${filename%.nc}"
    
    # 3. 截取最后一个下划线之后的内容 -> cldlow
    # 逻辑说明：${var##*_} 表示从左往右贪婪匹配直到最后一个 _，然后删掉匹配部分，保留剩下的
    var_name="${filename_no_ext##*nc_}"
    
    # 检查提取的变量名是否为空或者和原文件名一样（防止没有下划线的文件被乱移动）
    if [ -z "$var_name" ] || [ "$var_name" == "$filename_no_ext" ]; then
        echo "跳过: 文件 $filename 不符合 '..._变量名.nc' 的命名规则"
        continue
    fi

    # 4. 在目标目录下创建变量名对应的文件夹
    target_folder="$DEST_DIR/$var_name"
    if [ ! -d "$target_folder" ]; then
        echo "创建文件夹: $var_name"
        mkdir -p "$target_folder"
    fi
    
    # 5. 移动文件
    echo "移动: $filename -> $target_folder/"
    mv "$file" "$target_folder/"

done

echo "处理完成。"
