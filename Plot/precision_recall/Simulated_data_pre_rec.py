# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 15:22:50 2024

@author: zsh
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.lines as lines
###
## ==========================================================================================================================================
data = pd.read_csv("D:/My_WorkSpace/CNV_program/徐丹同-毕业生材料/1.毕业论文研究资料/1.2 毕业论文的数据、图表电子文档/1.2.1 数据/CNVpre_rec-data/human.simuDATA.alltype.pre-rec.csv")
data['depth'] = data['sampleNAME'].apply(lambda j: j.split('_')[1])
#data = data.drop(columns=['sampleNAME', 'positiveNUM', 'CNVNUM', 'dataTYPE'])

## 计算平均值
data = data.groupby(['mapper', 'CNVER', 'cnvtype', 'rateTYPE'], as_index=False).agg({'positiveRATE':'mean'})
data['positiveRATE'] = data['positiveRATE'].fillna(0)

filter_data = data

filter_data_dup = filter_data[filter_data['cnvtype'] == 'dup']
filter_data_del = filter_data[filter_data['cnvtype'] == 'del']

mappers = ['bwa', 'bismarkbt2', 'bsbolt', 'bsmap', 'bwameth', 'walt']

human_data = data
## =============================================================================================================================================

def F1_values(grid_x, grid_y):
    return 2 * (grid_x * grid_y) / (grid_x + grid_y + 1e-6)


grid_x, grid_y = np.meshgrid(np.linspace(0, 1.05, 100), 
                              np.linspace(0, 1.05, 100))

##======================================================================================================================================================
fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(24, 4))
for i, mapper in enumerate(mappers):
    mapper_data = filter_data_dup[filter_data_dup['mapper'] == mapper]
    row = 0
    col = i
    ax = axes[col]
    colors = {'CNVnator': '#AA0000', 'Pindel': '#CC6600', 'CNVkit': '#88AA00', 'cn_mops': '#0088A8', 'BreakDancer': '#990099', 'DELLY': '#330066', 'GASV':'#2F4F4F'}
        
    for CNVER in mapper_data['CNVER'].unique():
        cnver_data = mapper_data[mapper_data['CNVER'] == CNVER]
        ax.scatter(cnver_data[cnver_data['rateTYPE'] == 'recall']['positiveRATE'].values, 
                   cnver_data[cnver_data['rateTYPE'] == 'precision']['positiveRATE'].values, 
                   color=colors[CNVER], 
                   marker='v', 
                   s=150)
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, 1.05)
    
for i, mapper in enumerate(mappers):
    mapper_data_1 = filter_data_del[filter_data_del['mapper'] == mapper]
    row = 0
    col = i
    ax = axes[col]
    colors1 = {'CNVnator': '#AA0000', 'Pindel': '#CC6600', 'CNVkit': '#88AA00', 'cn_mops': '#0088A8', 'BreakDancer': '#990099', 'DELLY': '#330066', 'GASV':'#2F4F4F'}
        
    for CNVER1 in mapper_data_1['CNVER'].unique():
        cnver_data1 = mapper_data_1[mapper_data_1['CNVER'] == CNVER1]
        ax.scatter(cnver_data1[cnver_data1['rateTYPE'] == 'recall']['positiveRATE'].values, 
                   cnver_data1[cnver_data1['rateTYPE'] == 'precision']['positiveRATE'].values, 
                   color=colors1[CNVER1], 
                   marker='o', 
                   s=150)
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, 1.05)
             
    contour = ax.contour(grid_x, grid_y, F1_values(grid_x,grid_y), 10, linestyles='--', colors='black')
        
    plt.clabel(contour, inline=True, fontsize=8, colors='black')
        
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.linspace(0, 1, 6))     
    ax.set_title(f"{mapper}", fontsize=16, fontweight='bold')
    ax.set_xlabel("Recall", fontsize=16, fontweight='bold')
    if i == 0:
        ax.set_ylabel("Precision", fontsize=16, fontweight='bold')
        
        
line1 = lines.Line2D([0],[0], label='DUP', marker='v', linestyle='None', markerfacecolor='white', markeredgecolor='black', markersize=14)
line2 = lines.Line2D([0],[0], label='DEL', marker='o', linestyle='None', markerfacecolor='white', markeredgecolor='black', markersize=14)
patch1 = mpatches.Patch(color='#AA0000', label='CNVnator')
patch2 = mpatches.Patch(color='#CC6600', label='Pindel')
patch3 = mpatches.Patch(color='#88AA00', label='CNVkit')
patch4 = mpatches.Patch(color='#0088A8', label='cn_mops')
patch5 = mpatches.Patch(color='#990099', label='BreakDancer')
patch6 = mpatches.Patch(color='#330066', label='DELLY')
patch7 = mpatches.Patch(color='#2F4F4F', label='GASV')
            
handles1 = [patch1, patch2, patch3, patch4, patch5, patch6, patch7]
handles2 = [line1, line2]
            
fig.legend(handles=handles1, ncol=7, frameon=False, borderaxespad=1.5, fontsize=12)
fig.legend(handles=handles2, ncol=2, frameon=True, borderaxespad=1, loc='upper left', fontsize=13)
    
fig.suptitle("Human_simuDATA", fontsize=16, fontweight='bold', color='black')

plt.subplots_adjust(hspace=0.5, wspace=0.5)

plt.tight_layout()

plt.savefig('D:/My_WorkSpace/fig/准确率-召回率曲线-修正/模拟数据/Human_simuDATA_F1_平均-1.tiff', format='tiff', dpi=300)

plt.show()

































