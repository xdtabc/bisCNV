# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:30:59 2024

@author: zsh
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from brokenaxes import brokenaxes

###############################################################################################################################################
data = pd.read_csv('D:/My_WorkSpace/CNV_program/徐丹同-毕业生材料/1.毕业论文研究资料/1.2 毕业论文的数据、图表电子文档/1.2.1 数据/CNVnumber-data/human_simuDATA_wgbs_wgs_cnvnum.csv')
depth = []
for item in data['group']:
    parts = item.split("_")
    middle_part = parts[1]
    depth.append(middle_part)   
data['group'] = depth

data = data.rename(columns={'group': 'depth'})

data['CNVER'] = data['CNVER'].replace({'cn_mops':'cn.mops'}) 


data = data.groupby(['depth','CNVER', 'mapper', 'cnvtype'], as_index=False).agg({'number':'mean'})

data = data[~((data['mapper'] == 'bwa') & (data['depth'] != '15x'))]
data['CNVER'] = data['CNVER'].replace('cn_mops', 'cn.mops') 

data['z'] = data['CNVER'] + '-' + data['mapper'] + '-' + data['depth'] + '-' + data['cnvtype']
data['strategy'] = data['mapper'] + '-' + data['CNVER']
data['strategy1'] = data['mapper'] + '-' + data['CNVER'] + '-' + data['cnvtype']
data.sort_values(by='z', inplace=True)

data['check_bwa'] = data['strategy1'].str.startswith('bwa-')
data = data.sort_values(by=['check_bwa', 'z'], ascending=[True, True])

data['species'] = 'Human'

simu_data = data[data['mapper'] == 'bwa']

simu_data['number'] = np.log2(simu_data['number'] + 1).fillna(0)
###############################################################################################################################################
data = pd.read_csv('D:/My_WorkSpace/CNV_program/徐丹同-毕业生材料/1.毕业论文研究资料/1.2 毕业论文的数据、图表电子文档/1.2.1 数据/CNVlength-data/human_realDATA_wgbs_wgs_cnvlen.percent.csv')
data = data.drop(columns=['depth', 'CNVnum_perc'])
data = data[data['CNVER'] != 'CNVnator']
data = data[data['CNVER'] != 'CNVkit']
data = data[data['CNVER'] != 'Pindel']

CNVnator_data = pd.read_csv('D:/My_workspace/CNV_program/CNVlen/human/realDATA_CNVnator.csv')
CNVnator_data = CNVnator_data.rename(columns={'group': 'sample'})

CNVkit_data = pd.read_csv('D:/My_workspace/CNV_program/CNVlen/human/realDATA_CNVkit.csv')
CNVkit_data = CNVkit_data.rename(columns={'group': 'sample'})

Pindel_data = pd.read_csv('D:/My_workspace/CNV_program/CNVlen/human/realDATA_Pindel.csv')
Pindel_data = Pindel_data.rename(columns={'group': 'sample'})
Pindel_data['CNVtype'] = Pindel_data['CNVtype'].str.replace('DUP:TANDEM', 'DUP') 
Pindel_data = Pindel_data[(Pindel_data['CNVtype'] != 'INV') & (Pindel_data['CNVtype'] != 'INS') & (Pindel_data['CNVtype'] != 'RPL')]

data = pd.concat([data, CNVnator_data, CNVkit_data, Pindel_data], ignore_index=True)

new_rows = [{'CNVER':'GASV', 'mapper':'bsbolt', 'CNVtype':'DEL', 'CNVnum':0},
           {'CNVER':'GASV', 'mapper':'bsbolt', 'CNVtype':'DUP', 'CNVnum':0},

           {'CNVER':'GASV', 'mapper':'bwameth', 'CNVtype':'DEL', 'CNVnum':0},
           {'CNVER':'GASV', 'mapper':'bwameth', 'CNVtype':'DUP', 'CNVnum':0},

           {'CNVER':'GASV', 'mapper':'walt', 'CNVtype':'DEL', 'CNVnum':0},
           {'CNVER':'GASV', 'mapper':'walt', 'CNVtype':'DUP', 'CNVnum':0},
           ]

new_rows_df = pd.DataFrame(new_rows)

data_grouped = data.groupby(['CNVER', 'mapper', 'CNVtype', 'CNVlen'], as_index=False).agg({'CNVnum':'mean'})
data_grouped = data.groupby(['CNVER', 'mapper', 'CNVtype'], as_index=False).agg({'CNVnum':'sum'})

data_grouped = pd.concat([data_grouped, new_rows_df], ignore_index=True)

data_grouped['strategy'] = data_grouped['CNVER'] + '-' + data_grouped['mapper'] + "-" + data_grouped['CNVtype']
data_grouped['strategy1'] = data_grouped['CNVER'] + '-' + data_grouped['mapper'] 
data_grouped.sort_values(by='strategy', inplace=True)

data_grouped['CNVnum'] = np.log2(data_grouped['CNVnum'] + 1).fillna(0)
# data_grouped['CNVnum'] = np.where(data_grouped['CNVnum'] == -np.inf, 0, data_grouped['CNVnum'])
# data_grouped['CNVnum'] = np.where(data_grouped['CNVnum'] <= 0, 0, data_grouped['CNVnum'])

# colors = {'DEL': '#CD6889', 'DUP': '#4169E1'}
colors = {'DEL': 'RoyalBlue', 'DUP': 'Tomato'}
cnver_colors = {'BreakDancer': '#9400D3', 'CNVkit': '#FFA500', 'CNVnator': '#EE3B3B', 'DELLY':'#8B2323', 'GASV':'#CD8162', 'Pindel':'#006400', 'cn.mops': '#27408B'}


data_grouped['z'] = data_grouped['mapper'] + '-' + data_grouped['CNVER'] + '-' + data_grouped['CNVtype']
data_grouped['s'] = data_grouped['mapper'] + '-' + data_grouped['CNVER']
data_grouped['check_bwa'] = data_grouped['z'].str.startswith('bwa-')
data_grouped = data_grouped.sort_values(by=['check_bwa', 'strategy'], ascending=[True, True])

data_grouped = data_grouped[data_grouped['mapper'] == 'bwa']
data_grouped['species'] = 'Human'

human_data = data_grouped

real_data = human_data
###############################################################################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
bar_width = 0.2  # 柱宽
depths = simu_data['cnvtype'].unique()
strategies = simu_data['strategy'].unique()

num_strategies = len(strategies)
num_depths = len(depths)

total_group_width = num_depths * bar_width  # 所有柱子占用的总宽度
# total_width = bar_width * num_depths
# 计算组与组之间的间隙，这里可根据需要做调整
inter_group_gap = bar_width * 1.2  # 组间间隙设置为柱宽的一半
# 基于策略计算每个组的x坐标起始点
group_starts = np.arange(num_strategies) * (total_group_width + inter_group_gap)
# 最后一组数据的起始位置
last_group_start = group_starts[-1]
# 根据最后一组的起点和组内柱状图的总宽度，计算x轴的结束点
x_axis_end = last_group_start + total_group_width

# 设置x轴的界限，这里给界限留出一点额外的空间，以确保所有的柱状图及标签都能清晰显示
# 调整此处的margin值, 进一步微调空白区域的大小
margin = bar_width
ax.set_xlim(-margin, x_axis_end + margin)

s_values = simu_data['strategy'].unique()
s_indexes = {s: index for index, s in enumerate(s_values)}

# 绘制并列柱状图
for j, depth in enumerate(depths):
    # 获取特定cnvtype的数据
    depth_data = simu_data[simu_data['cnvtype'] == depth]
    
    # 计算当前深度内所有策略的x位置
    xs = group_starts + j * bar_width
    
    # # 为每个s值和对应的cnvtype计算位置
    # xs = np.arange(len(cnvtype_data['s'])) + offsets[i]

    # 绘制当前深度的柱状图
    ax.bar(xs, depth_data['number'], width=bar_width, label=depth,
           color=colors[depth], alpha=0.8, linewidth=0.05)
    
    

# 设置x轴刻度位置到组的中心
ax.set_xticks(group_starts + total_group_width / 2 - bar_width / 2)

ax.set_xticklabels(simu_data['strategy'].unique(), rotation=90, fontsize=14, fontweight= 'bold')

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
y_ticks = [0, 5, 10, 15, 20]
ax.set_ylim(0, 22)
ax.set_yticks(y_ticks)  # 假定 y_ticks 已经定义

ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
ax.set_ylabel("CNVs number (log2)", fontsize=16, fontweight='bold')
# ax.set_ylim(0, 25)
# ax.set_yticks([0, 5, 10 ,15 ,20])  # 假定 y_ticks 已经定义

fig.suptitle("Simulated WGS data", fontsize=20, fontweight='bold', color='black')

# 图例
# patch1 = mpatches.Patch(color='DodgerBlue', label='Human')
# patch2 = mpatches.Patch(color='#EE3B3B', label='Pig')
# patch3 = mpatches.Patch(color='DarkOrange', label='Cattle')
# patch4 = mpatches.Patch(color='#8A2BE2', label='30x')


# fig.legend(handles=[patch1, patch2, patch3], ncol=1, frameon=True, borderaxespad=2, fontsize=14)

plt.tight_layout()
plt.savefig('D:/My_WorkSpace/fig/CNVnum/WGS模拟-human.tiff', format='tiff', dpi=300)
plt.show()

###############################################################################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
bar_width = 0.2  # 柱宽
depths = real_data['CNVtype'].unique()
strategies = real_data['s'].unique()

num_strategies = len(strategies)
num_depths = len(depths)

total_group_width = num_depths * bar_width  # 所有柱子占用的总宽度
# total_width = bar_width * num_depths
# 计算组与组之间的间隙，这里可根据需要做调整
inter_group_gap = bar_width * 1.2  # 组间间隙设置为柱宽的一半
# 基于策略计算每个组的x坐标起始点
group_starts = np.arange(num_strategies) * (total_group_width + inter_group_gap)
# 最后一组数据的起始位置
last_group_start = group_starts[-1]
# 根据最后一组的起点和组内柱状图的总宽度，计算x轴的结束点
x_axis_end = last_group_start + total_group_width

# 设置x轴的界限，这里给界限留出一点额外的空间，以确保所有的柱状图及标签都能清晰显示
# 调整此处的margin值, 进一步微调空白区域的大小
margin = bar_width
ax.set_xlim(-margin, x_axis_end + margin)

s_values = real_data['s'].unique()
s_indexes = {s: index for index, s in enumerate(s_values)}

# 绘制并列柱状图
for j, depth in enumerate(depths):
    # 获取特定cnvtype的数据
    depth_data = real_data[real_data['CNVtype'] == depth]
    
    # 计算当前深度内所有策略的x位置
    xs = group_starts + j * bar_width
    
    # # 为每个s值和对应的cnvtype计算位置
    # xs = np.arange(len(cnvtype_data['s'])) + offsets[i]

    # 绘制当前深度的柱状图
    ax.bar(xs, depth_data['CNVnum'], width=bar_width, label=depth,
           color=colors[depth], alpha=0.8, linewidth=0.05)
    
    

# 设置x轴刻度位置到组的中心
ax.set_xticks(group_starts + total_group_width / 2 - bar_width / 2)

ax.set_xticklabels(real_data['s'].unique(), rotation=90, fontsize=14, fontweight= 'bold')

# 更改x轴标签颜色
for label in ax.get_xticklabels():
    text = label.get_text()
    if text.startswith('bwa-'):
        label.set_color('black')
    else:
    # 获取CNVER值以决定使用哪种颜色，这需要从标签的文本中解析CNVER
        cnver = label.get_text().split('-')[1]  # 假设格式始终是 'mapper-CNVER'
        label.set_color(cnver_colors.get(cnver, 'black'))  # 使用映射的颜色，如果没有匹配则默认为黑色

y_ticks = [0, 5, 10, 15, 20]
ax.set_ylim(0, 22)
ax.set_yticks(y_ticks)  # 假定 y_ticks 已经定义

ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
ax.set_ylabel("CNVs number (log2)", fontsize=16, fontweight='bold')
# ax.set_ylim(0, 25)
# ax.set_yticks([0, 5, 10 ,15 ,20])  # 假定 y_ticks 已经定义

fig.suptitle("Real WGS data", fontsize=20, fontweight='bold', color='black')

# 图例
# patch1 = mpatches.Patch(color='DodgerBlue', label='Human')
# patch2 = mpatches.Patch(color='#EE3B3B', label='Pig')
# patch3 = mpatches.Patch(color='DarkOrange', label='Cattle')
# patch4 = mpatches.Patch(color='#8A2BE2', label='30x')


# fig.legend(handles=[patch1, patch2, patch3], ncol=1, frameon=True, borderaxespad=2, fontsize=14)

plt.tight_layout()
plt.savefig('D:/My_WorkSpace/fig/CNVnum/WGS真实-human.tiff', format='tiff', dpi=300)
plt.show()
























