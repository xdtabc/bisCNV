# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 11:32:52 2024

@author: zsh
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from brokenaxes import brokenaxes

data = pd.read_csv('D:/My_WorkSpace/CNV_program/徐丹同-毕业生材料/1.毕业论文研究资料/1.2 毕业论文的数据、图表电子文档/1.2.1 数据/CNVnumber-data/human_simuDATA_wgbs_wgs_cnvnum.csv')
depth = []
for item in data['group']:
    parts = item.split("_")
    middle_part = parts[1]
    depth.append(middle_part)   
data['group'] = depth

data = data.rename(columns={'group': 'depth'})

data = data.groupby(['depth','CNVER', 'mapper', 'cnvtype'], as_index=False).agg({'number':'mean'})

data = data[~((data['mapper'] == 'bwa') & (data['depth'] != '15x'))]
data['CNVER'] = data['CNVER'].replace('cn_mops', 'cn.mops') 

data['z'] = data['CNVER'] + '-' + data['mapper'] + '-' + data['depth'] + '-' + data['cnvtype']
data['strategy'] = data['mapper'] + '-' + data['CNVER']
data['strategy1'] = data['mapper'] + '-' + data['CNVER'] + '-' + data['cnvtype']
data.sort_values(by='z', inplace=True)

data['check_bwa'] = data['strategy1'].str.startswith('bwa-')
data = data.sort_values(by=['check_bwa', 'z'], ascending=[True, True])

# data = data[data['CNVER'] != 'Pindel']
# data = data[data['CNVER'] != 'GASV']

data['number'] = np.log2(data['number'] + 1)
# data['number'] = np.where(data['number'] == -np.inf, 0, data['number'])
# data['number'] = np.where(data['number'] <= 0, 0, data['number'])


colors = {'10x': 'DodgerBlue', '15x': '#EE3B3B', '20x': 'DarkOrange', '30x': '#8A2BE2'}
cnver_colors = {'BreakDancer': '#9400D3', 'CNVkit': '#FFA500', 'CNVnator': '#EE3B3B', 'DELLY':'#8B2323', 'GASV':'#CD8162', 'Pindel':'#006400', 'cn.mops': '#27408B'}

data = data[data['mapper'] != 'bwa']


del_data = data[data['cnvtype'] == 'DEL']
dup_data = data[data['cnvtype'] == 'DUP']

## ======================================================================================================================================================================
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))

# 设置柱状图的基本参数
bar_width = 0.2  # 柱宽
depths = del_data['depth'].unique()
num_depths = len(depths)
total_width = bar_width * num_depths
# offsets = np.linspace(-total_width / 2, total_width / 2, num_depths)
# offsets = np.linspace(-bar_width/2, bar_width/2, num=num_depths)

# 计算组内偏移量，其中没有组内空隙
offsets = bar_width * np.arange(num_depths) - bar_width * (num_depths - 1) / 2

s_values = del_data['strategy'].unique()
s_indexes = {s: index for index, s in enumerate(s_values)}
# 绘制并列柱状图
for i, depth in enumerate(depths):
    # 获取特定cnvtype的数据
    depth_data = del_data[del_data['depth'] == depth]
    
    # 根据`s`值的索引和偏移量来确定x位置
    xs = np.array([s_indexes[row['strategy']] for _, row in depth_data.iterrows()]) + offsets[i]
    
    # # 为每个s值和对应的cnvtype计算位置
    # xs = np.arange(len(cnvtype_data['s'])) + offsets[i]

    # 绘制柱状图
    ax.bar(xs, depth_data['number'], width=bar_width, label=depth,
           color=colors[depth], alpha=0.8, linewidth=0.05)

# 设置x轴刻度位置到组的中心
ax.set_xticks(np.arange(len(s_values)))

ax.set_xticklabels(del_data['strategy'].unique(), rotation=90, fontsize=10, fontweight= 'bold')

# 更改x轴标签颜色
for label in ax.get_xticklabels():
    text = label.get_text()
    if text.startswith('bwa-'):
        label.set_color('black')
    else:
    # 获取CNVER值以决定使用哪种颜色，这需要从标签的文本中解析CNVER
        cnver = label.get_text().split('-')[1]  # 假设格式始终是 'mapper-CNVER'
        label.set_color(cnver_colors.get(cnver, 'black'))  # 使用映射的颜色，如果没有匹配则默认为黑色

# 设置图表的其他属性
# ax.set_ylim(0, 40)
# ax.set_yticks(y_ticks)  # 假定 y_ticks 已经定义

ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
ax.set_ylabel("CNVs number (log2)", fontsize=12, fontweight='bold')
fig.suptitle("Human_simuDATA (DEL)", fontsize=16, fontweight='bold', color='black')

# 图例
patch1 = mpatches.Patch(color='DodgerBlue', label='10x')
patch2 = mpatches.Patch(color='#EE3B3B', label='15x')
patch3 = mpatches.Patch(color='DarkOrange', label='20x')
patch4 = mpatches.Patch(color='#8A2BE2', label='30x')


fig.legend(handles=[patch1, patch2, patch3, patch4], ncol=1, frameon=True, borderaxespad=2, fontsize=10)

plt.tight_layout()
plt.savefig('D:/My_WorkSpace/fig/CNVnum/人模拟数据-DEL.tiff', format='tiff', dpi=300)
plt.show()


## ===================================================================================================================================================================

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))

# 设置柱状图的基本参数
bar_width = 0.2  # 柱宽
depths = dup_data['depth'].unique()
num_depths = len(depths)
total_width = bar_width * num_depths
# offsets = np.linspace(-total_width / 2, total_width / 2, num_depths)
# offsets = np.linspace(-bar_width/2, bar_width/2, num=num_depths)

# 计算组内偏移量，其中没有组内空隙
offsets = bar_width * np.arange(num_depths) - bar_width * (num_depths - 1) / 2

s_values = dup_data['strategy'].unique()
s_indexes = {s: index for index, s in enumerate(s_values)}
# 绘制并列柱状图
for i, depth in enumerate(depths):
    # 获取特定cnvtype的数据
    depth_data = dup_data[dup_data['depth'] == depth]
    
    # 根据`s`值的索引和偏移量来确定x位置
    xs = np.array([s_indexes[row['strategy']] for _, row in depth_data.iterrows()]) + offsets[i]
    
    # # 为每个s值和对应的cnvtype计算位置
    # xs = np.arange(len(cnvtype_data['s'])) + offsets[i]

    # 绘制柱状图
    ax.bar(xs, depth_data['number'], width=bar_width, label=depth,
           color=colors[depth], alpha=0.8, linewidth=0.05)

# 设置x轴刻度位置到组的中心
ax.set_xticks(np.arange(len(s_values)))

ax.set_xticklabels(dup_data['strategy'].unique(), rotation=90, fontsize=10, fontweight= 'bold')

# 更改x轴标签颜色
for label in ax.get_xticklabels():
    text = label.get_text()
    if text.startswith('bwa-'):
        label.set_color('black')
    else:
    # 获取CNVER值以决定使用哪种颜色，这需要从标签的文本中解析CNVER
        cnver = label.get_text().split('-')[1]  # 假设格式始终是 'mapper-CNVER'
        label.set_color(cnver_colors.get(cnver, 'black'))  # 使用映射的颜色，如果没有匹配则默认为黑色

# 设置图表的其他属性
# ax.set_ylim(0, 40)
# ax.set_yticks(y_ticks)  # 假定 y_ticks 已经定义

ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
ax.set_ylabel("CNVs number (log2)", fontsize=12, fontweight='bold')
fig.suptitle("Human_simuDATA (DUP)", fontsize=16, fontweight='bold', color='black')

# 图例
patch1 = mpatches.Patch(color='DodgerBlue', label='10x')
patch2 = mpatches.Patch(color='#EE3B3B', label='15x')
patch3 = mpatches.Patch(color='DarkOrange', label='20x')
patch4 = mpatches.Patch(color='#8A2BE2', label='30x')


fig.legend(handles=[patch1, patch2, patch3, patch4], ncol=1, frameon=True, borderaxespad=2, fontsize=10)

plt.tight_layout()
plt.savefig('D:/My_WorkSpace/fig/CNVnum/人模拟数据-DUP.tiff', format='tiff', dpi=300)
plt.show()

































