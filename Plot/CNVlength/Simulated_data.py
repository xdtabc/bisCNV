# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 10:36:12 2024

@author: zsh
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import re
from matplotlib.ticker import FuncFormatter
import matplotlib.lines as lines

data = pd.read_csv('D:/My_WorkSpace/CNV_program/徐丹同-毕业生材料/1.毕业论文研究资料/1.2 毕业论文的数据、图表电子文档/1.2.1 数据/CNVlength-data/human_simuDATA_wgbs_wgs_cnvlen.percent.csv')
data = data.drop(columns='CNVnum_perc')
data = data[data['CNVER'] != 'CNVnator']
data = data[data['CNVER'] != 'CNVkit']
data = data[data['CNVER'] != 'Pindel']

CNVnator_data = pd.read_csv('D:/My_workspace/CNV_program/CNVlen/human/simuDATA_CNVnator.csv')

CNVkit_data = pd.read_csv('D:/My_workspace/CNV_program/CNVlen/human/simuDATA_CNVkit.csv')

Pindel_data = pd.read_csv('D:/My_workspace/CNV_program/CNVlen/human/simuDATA_Pindel.csv')

data = pd.concat([data, CNVnator_data, CNVkit_data, Pindel_data], ignore_index=True)
# data['CNVnum'] = np.log2(data['CNVnum']).replace(float('-inf'),0)



data_grouped = data.groupby(['CNVER', 'mapper', 'CNVtype', 'CNVlen'], as_index=False).agg({'CNVnum':'mean'})
# data_grouped = pd.concat([data_grouped, new_rows_df], ignore_index=True)

data_grouped['z'] = data_grouped['CNVER'] + '-' + data_grouped['CNVtype']
data_grouped['strategy'] = data_grouped['CNVER'] + '-' + data_grouped['mapper'] + '-' + data_grouped['CNVtype'] + '-' + data_grouped['CNVlen']
data_grouped.sort_values(by='strategy', inplace=True)

data_grouped['s'] = data_grouped['CNVER'] + '-' + data_grouped['mapper']+ '-' + data_grouped['CNVtype']
total_number = data_grouped.groupby('s')['CNVnum'].sum()

data_grouped['CNV_percent'] = data_grouped.apply(lambda x: x['CNVnum'] / total_number[x['s']], axis=1)
data_grouped['CNV_percent'] = data_grouped['CNV_percent'].fillna(0)

y_max = data_grouped['CNVnum'].max()
y_ticks = list(range(0, 45, 5))

mappers = ['bwa', 'bismarkbt2', 'bsbolt', 'bsmap', 'bwameth', 'walt']

cnver_colors = {'BreakDancer': '#9400D3', 'CNVkit': '#FFA500', 'CNVnator': '#EE3B3B', 'DELLY':'#8B2323', 'GASV':'#CD8162', 'Pindel':'#006400', 'cn.mops': '#27408B'}

del_data = data_grouped[data_grouped['CNVtype'] == 'DEL']
dup_data = data_grouped[data_grouped['CNVtype'] == 'DUP']

# data_grouped.to_csv('D:\\My_WorkSpace\\CNV_program\\CNVlen\\human_len_simu.csv', index=False)

fig, axes = plt.subplots(nrows=2, ncols=6, figsize=(18, 5))
for i, mapper in enumerate(mappers):
    # 获取特定mapper的数据
    mapper_data = del_data[del_data['mapper'] == mapper]
    row = 0
    col = i
    # row = i // 2
    # col = i % 2

    # 在子图中创建条形图
    ax = axes[row, col]
    
    for cnver in mapper_data['CNVER'].unique():
        cnver_data = mapper_data[mapper_data['CNVER'] == cnver]
        ax.plot(cnver_data['CNVlen'], cnver_data['CNV_percent'], label=cnver, marker='o', color=cnver_colors[cnver], 
                linestyle=':', 
                alpha=0.9,
                markersize=7)
        
    ax.set_xticks(cnver_data['CNVlen']) 
    ax.set_xticklabels(cnver_data['CNVlen'], rotation=0,  fontsize=8, fontweight='bold')
    
    if i == 0:
        ax.set_ylabel("Percentage", fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('')

        
    ax.set_title(f"{mapper}", fontsize=12, fontweight='bold')
    ax.set_xlabel("", fontsize=14, fontweight='bold')
    # 标签格式化
    # ax.get_xaxis().set_major_formatter(FuncFormatter(lambda x, _: f'  {x}' if x == min(mapper_data['z']) else f'{x}'))
    
    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
    ax.set_ylim(0, 1.05)



for i, mapper in enumerate(mappers):
    # 获取特定mapper的数据
    mapper_data = dup_data[dup_data['mapper'] == mapper]
    row = 1
    col = i
    # row = i // 2
    # col = i % 2

    # 在子图中创建条形图
    ax = axes[row, col]
    
    for cnver in mapper_data['CNVER'].unique():
        cnver_data = mapper_data[mapper_data['CNVER'] == cnver]
        ax.plot(cnver_data['CNVlen'], cnver_data['CNV_percent'], label=cnver, marker='v', color=cnver_colors[cnver], 
                linestyle=':', 
                alpha=0.9,
                markersize=7)
        
    ax.set_xticks(cnver_data['CNVlen']) 
    ax.set_xticklabels(cnver_data['CNVlen'], rotation=0,  fontsize=8, fontweight='bold')
    
    if i == 0:
        ax.set_ylabel("Percentage", fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('')

        
    ax.set_title(f"{mapper}", fontsize=12, fontweight='bold')
    ax.set_xlabel("", fontsize=14, fontweight='bold')
    # 标签格式化
    # ax.get_xaxis().set_major_formatter(FuncFormatter(lambda x, _: f'  {x}' if x == min(mapper_data['z']) else f'{x}'))
    
    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
    ax.set_ylim(0, 1.05)


    
   
line1 = lines.Line2D([0],[0], label='BreakDancer', marker='o', linestyle='None', markerfacecolor='#9400D3', markeredgecolor='white', markersize=12)
line2 = lines.Line2D([0],[0], label='CNVkit', marker='o', linestyle='None', markerfacecolor='#FFA500', markeredgecolor='white', markersize=12)
line3 = lines.Line2D([0],[0], label='CNVnator', marker='o', linestyle='None', markerfacecolor='#EE3B3B', markeredgecolor='white', markersize=12)
line4 = lines.Line2D([0],[0], label='DELLY', marker='o', linestyle='None', markerfacecolor='#8B2323', markeredgecolor='white', markersize=12)
line5 = lines.Line2D([0],[0], label='GASV', marker='o', linestyle='None', markerfacecolor='#CD8162', markeredgecolor='white', markersize=12)
line6 = lines.Line2D([0],[0], label='Pindel', marker='o', linestyle='None', markerfacecolor='#006400', markeredgecolor='white', markersize=12)
line7 = lines.Line2D([0],[0], label='cn.mops', marker='o', linestyle='None', markerfacecolor='#27408B', markeredgecolor='white', markersize=12)


line8 = lines.Line2D([0],[0], label='DUP', marker='v', linestyle='None', markerfacecolor='white', markeredgecolor='black', markersize=14)
line9 = lines.Line2D([0],[0], label='DEL', marker='o', linestyle='None', markerfacecolor='white', markeredgecolor='black', markersize=14)


handles = [line1, line2, line3, line4, line5, line6 ,line7]
handles2 = [line8, line9]

fig.legend(handles=handles, ncol=4, frameon=False, borderaxespad=0, fontsize=10)
fig.legend(handles=handles2, ncol=2, frameon=True, borderaxespad=0, loc='upper left', fontsize=10)

fig.suptitle("Human_simuDATA", fontsize=14, fontweight='bold', color='black')
plt.tight_layout()
plt.show()
plt.savefig('D:/My_WorkSpace/fig/CNV长度分布/模拟数据/折线图/人类模拟数据-1.tiff', format='tiff', dpi=300)



















