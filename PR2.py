import os
import sys
import argparse
import glob
from collections import defaultdict


'''
密码子第三位碱基含量分析工具

功能描述：计算指定密码子中第三位碱基（A/U/G/C）的比例，用于密码子使用偏好性分析

输入格式：
  基因名1
  AAA 10 AAC 5 AAG 8 AAU 3
  ACA 7 ACC 12 ACG 4 ACU 6
  （密码子-频次对，空格分隔）
  
  基因名2
  GCA 15 GCC 8 GCG 5 GCU 10
  ...

输出文件：
  {输入文件名}_with_labels.txt  : 添加基因标签的中间文件
  {输入文件名}_out.txt          : 结果表格（基因名 A U G C）

目标密码子（32个）：
  GCA GCC GCG GCU  GGA GGC GGG GGU
  CUA CUC CUG CUU  CCA CCC CCG CCU
  CGA CGC CGG CGU  UCA UCC UCG UCU
  ACA ACC ACG ACU  GUA GUC GUG GUU

计算方法：
  对每个基因，提取目标密码子的频次，
  计算第三位碱基为A/U/G/C的比例

使用方法：
  python script.py <输入文件> [输入文件...]

输入方式：
  1. 单文件      : python script.py data.txt
  2. 通配符      : python script.py *.txt
  3. 多文件      : python script.py file1.txt file2.txt

输出文件位置：
  与输入文件相同目录

依赖：
  无需外部依赖（仅使用Python标准库）

输出示例（out.txt）：
  Gene        A       U       G       C
  rps12       0.3245  0.2876  0.1987  0.1892
  rpl16       0.2987  0.3124  0.1876  0.2013
'''


def parse_gene_data(lines):
    """解析基因数据块，返回密码子字典"""
    codon_dict = {}
    for line in lines:
        parts = line.split()
        for i in range(0, len(parts), 2):
            if i + 1 >= len(parts):
                continue
            codon = parts[i]
            try:
                count = int(parts[i+1])
                codon_dict[codon] = codon_dict.get(codon, 0) + count
            except ValueError:
                continue
    return codon_dict

def calculate_third_base_content(target_data):
    """计算目标密码子中第三位碱基A、U、G、C的含量"""
    base_count = {'A': 0, 'U': 0, 'G': 0, 'C': 0}
    total = 0
    
    for codon, count in target_data.items():
        if len(codon) == 3:
            third_base = codon[2]
            if third_base in base_count:
                base_count[third_base] += count
                total += count
    
    if total == 0:
        return {base: 0.0 for base in base_count}
    
    return {base: count / total for base, count in base_count.items()}

# 目标密码子列表
target_codons = [
    'GCA', 'GCC', 'GCG', 'GCU', 'GGA', 'GGC', 'GGG', 'GGU',
    'CUA', 'CUC', 'CUG', 'CUU', 'CCA', 'CCC', 'CCG', 'CCU',
    'CGA', 'CGC', 'CGG', 'CGU', 'UCA', 'UCC', 'UCG', 'UCU',
    'ACA', 'ACC', 'ACG', 'ACU', 'GUA', 'GUC', 'GUG', 'GUU'
]

def process_single_file(input_path):
    """处理单个文件的核心逻辑"""
    try:
        # 获取文件所在目录和基本名称
        file_dir = os.path.dirname(input_path)
        file_basename = os.path.splitext(os.path.basename(input_path))[0]
        
        # 构建中间文件和输出文件路径
        intermediate_path = os.path.join(file_dir, f"{file_basename}_with_labels.txt")
        output_path = os.path.join(file_dir, f"{file_basename}_out.txt")
        
        # 读取文件内容
        with open(input_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        # 步骤1: 反转处理逻辑 - 先识别数据块，再关联标签
        gene_blocks = []
        current_block = []
        pending_labels = []
        
        print(f"\n处理文件: {os.path.basename(input_path)}")
        print("创建带基因标签的中间文件...")
        with open(intermediate_path, 'w', encoding='utf-8') as out_file:
            for line in lines:
                # 检测基因标签 - 单独一行且不包含空格
                if len(line.split()) == 1 and not line.isdigit():
                    pending_labels.append(line)
                    
                    # 如果当前有数据块，将其与之前的标签关联
                    if current_block and pending_labels:
                        gene_label = pending_labels.pop(0)
                        out_file.write(f"#GENE_LABEL: {gene_label}\n")
                        for data_line in current_block:
                            out_file.write(data_line + "\n")
                        out_file.write("\n")
                        gene_blocks.append((gene_label, current_block.copy()))
                        current_block = []
                else:
                    current_block.append(line)
            
            # 处理文件末尾的情况
            if current_block and pending_labels:
                gene_label = pending_labels.pop(0)
                out_file.write(f"#GENE_LABEL: {gene_label}\n")
                for data_line in current_block:
                    out_file.write(data_line + "\n")
                gene_blocks.append((gene_label, current_block.copy()))
            elif current_block and not pending_labels:
                gene_label = f"Unknown_Gene_{len(gene_blocks)+1}"
                print(f"警告: 文件末尾有未标记的数据块，使用默认标签: {gene_label}")
                out_file.write(f"#GENE_LABEL: {gene_label}\n")
                for data_line in current_block:
                    out_file.write(data_line + "\n")
                gene_blocks.append((gene_label, current_block.copy()))
            
            if pending_labels:
                print(f"警告: 有未使用的基因标签: {', '.join(pending_labels)}")
        
        print(f"中间文件已保存到: {intermediate_path}")
        print(f"找到 {len(gene_blocks)} 个基因块")
        
        # 步骤2: 处理每个基因块
        results = {}
        
        print("处理每个基因块...")
        for gene_name, data_lines in gene_blocks:
            print(f"处理基因: {gene_name} ({len(data_lines)}行)")
            codon_dict = parse_gene_data(data_lines)
            target_data = {codon: codon_dict.get(codon, 0) for codon in target_codons}
            third_base_content = calculate_third_base_content(target_data)
            results[gene_name] = third_base_content
        
        # 步骤3: 输出最终结果
        with open(output_path, 'w', encoding='utf-8') as out_file:
            out_file.write("Gene\tA\tU\tG\tC\n")
            for gene, content in results.items():
                out_file.write(f"{gene}\t")
                out_file.write(f"{content['A']:.4f}\t")
                out_file.write(f"{content['U']:.4f}\t")
                out_file.write(f"{content['G']:.4f}\t")
                out_file.write(f"{content['C']:.4f}\n")
        
        print(f"最终结果已保存到: {output_path}")
        return True
        
    except Exception as e:
        print(f"处理文件 {input_path} 时发生错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return False

def expand_file_patterns(patterns):
    """扩展通配符模式为实际文件列表[1,3](@ref)"""
    expanded_files = []
    for pattern in patterns:
        # 处理Windows路径中的特殊字符[1](@ref)
        pattern = pattern.replace('[', '[[]')  # 防止glob将[]视为特殊字符
        
        # 使用glob进行通配符扩展[1,3](@ref)
        matches = glob.glob(pattern, recursive=True)
        
        if not matches:
            print(f"警告: 通配符 '{pattern}' 没有匹配到任何文件")
        else:
            # 过滤掉目录，只保留文件[6](@ref)
            matches = [f for f in matches if os.path.isfile(f)]
            expanded_files.extend(matches)
    
    # 去重处理
    return list(set(expanded_files))

def main():
    """主函数，处理命令行参数和文件遍历[1,3](@ref)"""
    parser = argparse.ArgumentParser(description='基因密码子分析工具（支持通配符）')
    parser.add_argument('files', nargs='+', help='输入文件路径(支持通配符如*.txt)')
    args = parser.parse_args()
    
    # 扩展通配符为实际文件列表[1,3](@ref)
    file_list = expand_file_patterns(args.files)
    
    if not file_list:
        print("错误: 没有找到匹配的文件", file=sys.stderr)
        if sys.platform.startswith('win'):
            input("\n按Enter键退出程序...")
        return
    
    print(f"找到 {len(file_list)} 个匹配的文件:")
    for i, file_path in enumerate(file_list, 1):
        print(f"{i}. {file_path}")
    
    # 处理所有输入文件
    processed_count = 0
    for input_path in file_list:
        if process_single_file(input_path):
            processed_count += 1
    
    print(f"\n=== 执行完成 ===\n成功处理 {processed_count} 个文件")
    
    # 防止窗口关闭（Windows特有）
    if sys.platform.startswith('win'):
        input("\n按Enter键退出程序...")

if __name__ == "__main__":
    main()