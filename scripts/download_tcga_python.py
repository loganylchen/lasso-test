#!/usr/bin/env python3
"""
使用Python直接从TCGA GDC下载COAD (结直肠癌) RNA-seq数据
不依赖R的TCGAbiolinks/recount3
"""

import requests
import json
import pandas as pd
import gzip
import os
from pathlib import Path

SAVE_DIR = Path.home() / ".openclaw/workspace/lasso-biomarker-tutorial/data/tcga_coad_python"
SAVE_DIR.mkdir(parents=True, exist_ok=True)

def query_tcga_files(project="TCGA-COAD", data_type="Gene Expression Quantification", 
                     workflow="STAR - Counts", max_files=60):
    """查询TCGA数据文件"""
    
    print("1. 查询TCGA-COAD RNA-seq数据...")
    
    # GDC API endpoint
    files_endpoint = "https://api.gdc.cancer.gov/files"
    
    # 构建查询filters
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "data_type", "value": [data_type]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": [workflow]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}}  # 只要公开数据
        ]
    }
    
    # 构建查询参数
    params = {
        "filters": json.dumps(filters),
        "format": "JSON",
        "size": str(max_files),
        "fields": "file_id,file_name,cases.samples.sample_type,cases.case_id"
    }
    
    # 发送请求
    response = requests.get(files_endpoint, params=params)
    
    if response.status_code != 200:
        print(f"❌ 查询失败: {response.status_code}")
        return None
    
    data = response.json()
    hits = data["data"]["hits"]
    
    print(f"   找到 {len(hits)} 个文件")
    
    # 分析样本类型
    sample_types = {}
    for hit in hits:
        samples = hit.get("cases", [{}])[0].get("samples", [])
        if samples:
            sample_type = samples[0].get("sample_type", "Unknown")
            sample_types[sample_type] = sample_types.get(sample_type, 0) + 1
    
    print("   样本类型分布:")
    for st, count in sample_types.items():
        print(f"   - {st}: {count}")
    
    return hits

def select_balanced_samples(hits, n_normal=30, n_tumor=30):
    """选择平衡的Normal和Tumor样本"""
    
    print(f"\n2. 选择样本（{n_normal} Normal + {n_tumor} Tumor）...")
    
    normal_files = []
    tumor_files = []
    
    for hit in hits:
        samples = hit.get("cases", [{}])[0].get("samples", [])
        if not samples:
            continue
        
        sample_type = samples[0].get("sample_type", "")
        
        if "Normal" in sample_type and len(normal_files) < n_normal:
            normal_files.append({
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "case_id": hit["cases"][0]["case_id"],
                "sample_type": "Normal"
            })
        elif "Tumor" in sample_type or "Primary" in sample_type:
            if len(tumor_files) < n_tumor:
                tumor_files.append({
                    "file_id": hit["file_id"],
                    "file_name": hit["file_name"],
                    "case_id": hit["cases"][0]["case_id"],
                    "sample_type": "Tumor"
                })
    
    selected = normal_files + tumor_files
    print(f"   选择了 {len(normal_files)} Normal + {len(tumor_files)} Tumor = {len(selected)} 个样本")
    
    return selected

def download_file(file_id, output_path):
    """下载单个文件"""
    
    data_endpoint = f"https://api.gdc.cancer.gov/data/{file_id}"
    
    response = requests.get(data_endpoint, stream=True)
    
    if response.status_code != 200:
        return False
    
    with open(output_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    
    return True

def download_samples(selected_files):
    """批量下载样本文件"""
    
    print(f"\n3. 下载文件（共 {len(selected_files)} 个）...")
    print("   ⚠️  这可能需要较长时间，每个文件约1-5MB")
    
    download_dir = SAVE_DIR / "raw_files"
    download_dir.mkdir(exist_ok=True)
    
    success_files = []
    
    for i, file_info in enumerate(selected_files, 1):
        file_id = file_info["file_id"]
        file_name = file_info["file_name"]
        output_path = download_dir / file_name
        
        if output_path.exists():
            print(f"   [{i}/{len(selected_files)}] 跳过（已存在）: {file_name}")
            success_files.append({**file_info, "local_path": str(output_path)})
            continue
        
        print(f"   [{i}/{len(selected_files)}] 下载: {file_name}...", end=" ", flush=True)
        
        if download_file(file_id, output_path):
            print("✓")
            success_files.append({**file_info, "local_path": str(output_path)})
        else:
            print("✗")
    
    print(f"\n   成功下载 {len(success_files)} 个文件")
    
    # 保存manifest
    manifest_df = pd.DataFrame(success_files)
    manifest_df.to_csv(SAVE_DIR / "manifest.csv", index=False)
    
    return success_files

def parse_counts_file(file_path):
    """解析STAR counts文件"""
    
    # STAR counts格式: gene_id \t unstranded \t stranded_first \t stranded_second
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None, 
                     names=['gene_id', 'unstranded', 'stranded_first', 'stranded_second'])
    
    # 过滤掉统计行（以N_开头）
    df = df[~df['gene_id'].str.startswith('N_')]
    
    # 使用unstranded counts
    return df[['gene_id', 'unstranded']].set_index('gene_id')

def build_count_matrix(success_files):
    """构建count矩阵"""
    
    print("\n4. 构建count矩阵...")
    
    count_matrices = []
    sample_ids = []
    sample_types = []
    
    for file_info in success_files:
        file_path = file_info["local_path"]
        
        # 解析counts
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                counts = parse_counts_file(f)
        else:
            counts = parse_counts_file(file_path)
        
        # 重命名列
        sample_id = file_info["case_id"]
        counts.columns = [sample_id]
        
        count_matrices.append(counts)
        sample_ids.append(sample_id)
        sample_types.append(file_info["sample_type"])
    
    # 合并所有样本
    count_matrix = pd.concat(count_matrices, axis=1)
    
    print(f"   ✓ Count矩阵: {count_matrix.shape[0]} 基因 × {count_matrix.shape[1]} 样本")
    
    # 保存
    count_matrix.to_csv(SAVE_DIR / "count_matrix.csv")
    
    # 保存样本信息
    sample_info = pd.DataFrame({
        'sample_id': sample_ids,
        'condition': sample_types
    })
    sample_info.to_csv(SAVE_DIR / "sample_info.csv", index=False)
    
    print(f"\n   样本类型:")
    print(sample_info['condition'].value_counts())
    
    return count_matrix, sample_info

def main():
    print("=== 使用Python下载TCGA-COAD数据 ===\n")
    
    # 1. 查询文件
    hits = query_tcga_files(max_files=100)
    if not hits:
        print("❌ 查询失败")
        return
    
    # 2. 选择样本
    selected = select_balanced_samples(hits, n_normal=30, n_tumor=30)
    
    if len(selected) < 20:
        print("❌ 样本数量不足")
        return
    
    # 3. 下载文件
    success_files = download_samples(selected)
    
    if len(success_files) < 10:
        print("❌ 下载失败文件过多")
        return
    
    # 4. 构建count矩阵
    count_matrix, sample_info = build_count_matrix(success_files)
    
    print(f"\n=== 下载完成 ===")
    print(f"数据保存到: {SAVE_DIR}")
    print(f"- count_matrix.csv ({count_matrix.shape[0]} genes × {count_matrix.shape[1]} samples)")
    print(f"- sample_info.csv")
    print(f"- manifest.csv")
    print(f"\n下一步：使用R脚本加载CSV数据进行DESeq2分析")

if __name__ == "__main__":
    main()
