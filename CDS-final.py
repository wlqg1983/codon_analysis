import os
import argparse
import glob  # 新增glob模块导入
from Bio import SeqIO
from collections import OrderedDict


'''
FASTA基因序列筛选与去重工具

功能：筛选完整CDS序列（长度>300bp + ATG开头 + 终止密码子），每个基因保留最长序列

用法：
  python script.py <输入> [输入...]

输入方式：
  1. 单文件    : python script.py input.fasta
  2. 通配符    : python script.py *.fasta *.fas
  3. 目录      : python script.py /path/to/dir/

输出：
  {原文件名}_filtered_{唯一基因数}.fasta

筛选条件：
  - 序列长度 > 300 bp
  - 以 ATG 开头（起始密码子）
  - 以 TAA/TAG/TGA 结尾（终止密码子）

去重规则：
  - 按FASTA头部第一个单词作为基因名
  - 同一基因保留最长序列

依赖：
  pip install biopython
'''

def extract_gene_name(header):
    """从FASTA头信息中提取首个单词作为基因名"""
    return header.split()[0]  # 直接取第一个空格前的词（如 "rps12"）

def process_single_fasta(input_file, output_file=None):
    """处理单个FASTA文件：筛选 + 动态去重"""
    # 初始化数据结构
    gene_dict = OrderedDict()  # 保留顺序 {基因名: SeqRecord}
    stats = {
        'total': 0,
        'passed': 0,
        'duplicates': 0,
        'replaced': 0
    }
    
    # 解析并筛选序列
    for record in SeqIO.parse(input_file, "fasta"):
        stats['total'] += 1
        seq = str(record.seq).upper()
        
        # 筛选条件检查
        if len(seq) <= 300: 
            continue
        if not seq.startswith('ATG'):
            continue
        if seq[-3:] not in ['TAA', 'TAG', 'TGA']:
            continue
        stats['passed'] += 1
        
        # 提取基因名
        gene_name = extract_gene_name(record.id)
        
        # 动态去重：保留最长序列
        if gene_name in gene_dict:
            stats['duplicates'] += 1
            existing_seq = str(gene_dict[gene_name].seq)
            if len(seq) > len(existing_seq):
                gene_dict[gene_name] = record
                stats['replaced'] += 1
        else:
            gene_dict[gene_name] = record
    
    # 生成输出文件名
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    if not output_file:
        output_file = f"{base_name}_filtered_{len(gene_dict)}.fasta"
    
    # 写入结果
    with open(output_file, 'w') as f:
        SeqIO.write(gene_dict.values(), f, 'fasta')
    
    # 返回统计信息
    return stats, len(gene_dict), output_file

def batch_process(input_dir):
    """批量处理目录下所有FASTA文件"""
    fasta_files = [f for f in os.listdir(input_dir) 
                  if f.endswith(('.fasta', '.fas', '.fa'))]
    print(f"发现 {len(fasta_files)} 个FASTA文件待处理...")
    results = []
    
    for filename in fasta_files:
        input_path = os.path.join(input_dir, filename)
        print(f"\n处理: {filename}")
        stats, unique_count, output_path = process_single_fasta(input_path)
        
        # 打印单文件报告
        print(f"├─ 原始序列: {stats['total']}")
        print(f"├─ 合格序列: {stats['passed']} (长度>300+ATG+有效终止)")
        print(f"├─ 重复基因: {stats['duplicates']}组 | 更新: {stats['replaced']}条")
        print(f"└─ 保存: {os.path.basename(output_path)} (唯一基因: {unique_count})")
        results.append((filename, unique_count))
    
    # 打印批量总结
    print("\n" + "="*50)
    print("批量处理完成！汇总统计:")
    for file, count in results:
        print(f"· {file}: 保留 {count} 条唯一序列")
    print("="*50)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FASTA基因序列筛选与去重工具')
    # 修改参数定义：使用nargs='*'接收多个参数，包括通配符扩展结果[6](@ref)
    parser.add_argument('input_paths', nargs='*', type=str, help='输入文件/目录路径，支持通配符如*.fas')
    args = parser.parse_args()
    
    # 如果没有提供任何参数
    if not args.input_paths:
        print("请提供输入文件或目录路径，支持通配符如: python script.py *.fas")
        print("或处理单个文件: python script.py input.fasta")  
        print("或处理目录: python script.py /path/to/directory")
        exit(1)
    
    # 使用glob模块处理通配符，展开所有文件模式[1,4](@ref)
    expanded_files = []
    for pattern in args.input_paths:
        # 对每个模式使用glob进行文件匹配，支持Unix shell风格的通配符[2,5](@ref)
        matched_files = glob.glob(pattern)
        if matched_files:
            expanded_files.extend(matched_files)
        else:
            # 如果没有匹配到文件，保留原模式（可能是目录或不存在的路径）
            expanded_files.append(pattern)
    
    # 去重并处理文件列表
    unique_paths = list(OrderedDict.fromkeys(expanded_files))
    
    # 处理逻辑分为三种情况：多文件、单文件、单目录
    if len(unique_paths) > 1:
        print(f"发现 {len(unique_paths)} 个路径待处理...")
        results = []
        valid_count = 0
        
        for input_path in unique_paths:
            if os.path.isfile(input_path):
                valid_count += 1
                stats, count, out_file = process_single_fasta(input_path)
                print(f"处理: {os.path.basename(input_path)} → 输出基因: {count}条")
                results.append((os.path.basename(input_path), count))
            elif os.path.isdir(input_path):
                print(f"跳过目录（需单独处理）: {input_path}")
            else:
                print(f"警告: {input_path} 不是有效文件，跳过")
        
        # 打印批量处理总结
        if valid_count > 0:
            print("\n" + "="*50)
            print("通配符批量处理完成！汇总统计:")
            for file, count in results:
                print(f"· {file}: 保留 {count} 条唯一序列")
            print("="*50)
        else:
            print("未找到任何有效文件进行处理")
    
    # 单个路径的情况
    else:
        input_path = unique_paths[0]
        
        if os.path.isfile(input_path):
            # 单文件模式
            stats, count, out_file = process_single_fasta(input_path)
            print("\n" + "="*50)
            print(f"处理报告: {os.path.basename(input_path)}")
            print("="*50)
            print(f"原始序列: {stats['total']} | 合格序列: {stats['passed']}")
            print(f"重复基因: {stats['duplicates']}组 | 更新: {stats['replaced']}条")
            print(f"输出文件: {out_file} (唯一基因: {count})")
            print("="*50)
        elif os.path.isdir(input_path):
            # 目录批量模式
            batch_process(input_path)
        else:
            print(f"错误：路径 '{input_path}' 不存在！")
            print("请检查路径是否正确，或使用通配符如: python script.py *.fas")