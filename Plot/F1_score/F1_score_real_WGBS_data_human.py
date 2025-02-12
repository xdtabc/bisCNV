# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 11:23:43 2024

@author: zsh
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

human_data = pd.read_csv('D:/My_WorkSpace/CNV_program/徐丹同-毕业生材料/1.毕业论文研究资料/1.2 毕业论文的数据、图表电子文档/1.2.1 数据/CNV_F1score-data/human.realDATA.alltype.f1_average.csv')
human_data = human_data[human_data['mapper'] == 'bwa']
human_data['F1_score'] = human_data['F1_score'].str.strip('%').astype('float') / 100
human_data['species'] = 'Human'

data = human_data
data['z'] = data['CNVER'] + '-' + data['mapper'] + '-' + data['cnvtype']
data['strategy'] = data['mapper'] + '-' + data['CNVER']
data['strategy1'] = data['mapper'] + '-' + data['CNVER'] + '-' + data['cnvtype']
data.sort_values(by='z', inplace=True)

real_data = data
## ===================================================================================================================================================
colors = {'del': 'RoyalBlue', 'dup': 'Tomato'}
cnver_colors = {'BreakDancer': '#9400D3', 'CNVkit': '#FFA500', 'CNVnator': '#EE3B3B', 'DELLY':'#8B2323', 'GASV':'#CD8162', 'Pindel':'#006400', 'cn.mops': '#27408B'}
###############################################################################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
bar_width = 0.2  # 柱宽
depths = real_data['cnvtype'].unique()
strategies = real_data['strategy'].unique()

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

s_values = real_data['strategy'].unique()
s_indexes = {s: index for index, s in enumerate(s_values)}

# 绘制并列柱状图
for j, depth in enumerate(depths):
    # 获取特定cnvtype的数据
    depth_data = real_data[real_data['cnvtype'] == depth]
    
    # 计算当前深度内所有策略的x位置
    xs = group_starts + j * bar_width
    
    # # 为每个s值和对应的cnvtype计算位置
    # xs = np.arange(len(cnvtype_data['s'])) + offsets[i]

    # 绘制当前深度的柱状图
    ax.bar(xs, depth_data['F1_score'], width=bar_width, label=depth,
           color=colors[depth], alpha=0.8, linewidth=0.05)
    
    

# 设置x轴刻度位置到组的中心
ax.set_xticks(group_starts + total_group_width / 2 - bar_width / 2)

ax.set_xticklabels(real_data['strategy'].unique(), rotation=90, fontsize=14, fontweight= 'bold')

# 更改x轴标签颜色
for label in ax.get_xticklabels():
    text = label.get_text()
    if text.startswith('bwa-'):
        label.set_color('black')
    else:
    # 获取CNVER值以决定使用哪种颜色，这需要从标签的文本中解析CNVER
        cnver = label.get_text().split('-')[1]  # 假设格式始终是 'mapper-CNVER'
        label.set_color(cnver_colors.get(cnver, 'black'))  # 使用映射的颜色，如果没有匹配则默认为黑色

y_ticks = [0, 0.2, 0.4, 0.6, 0.8]
ax.set_ylim(0, 1)
ax.set_yticks(y_ticks)  # 假定 y_ticks 已经定义

ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.5)
ax.set_ylabel("F1 score", fontsize=16, fontweight='bold')
# ax.set_ylim(0, 25)
# ax.set_yticks([0, 5, 10 ,15 ,20])  # 假定 y_ticks 已经定义

# fig.suptitle("Real WGS data", fontsize=20, fontweight='bold', color='black')

# 图例
# patch1 = mpatches.Patch(color='DodgerBlue', label='Human')
# patch2 = mpatches.Patch(color='#EE3B3B', label='Pig')
# patch3 = mpatches.Patch(color='DarkOrange', label='Cattle')
# patch4 = mpatches.Patch(color='#8A2BE2', label='30x')


# fig.legend(handles=[patch1, patch2, patch3], ncol=1, frameon=True, borderaxespad=2, fontsize=14)

plt.tight_layout()
plt.savefig('D:/My_WorkSpace/fig/F1分数-修正/WGS真实-human.tiff', format='tiff', dpi=300)
plt.show()
