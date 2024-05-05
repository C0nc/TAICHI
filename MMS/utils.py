import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
 
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu


def plot_violin_plot(train_adata, ct_obs, batch_obs):
    label = ct_obs
    df = train_adata.obs

    # 计算每个new_labels和batch_obs组合中所有cell_type的总计数
    total_counts = df.groupby(['new_labels', batch_obs]).size().reset_index(name='total_counts')

    # 计算每个cell_type的计数
    counts = df.groupby(['new_labels', batch_obs, label]).size().reset_index(name='counts')

    # 将总计数合并回每个cell_type的计数DataFrame
    counts = pd.merge(counts, total_counts, on=['new_labels', batch_obs])

    # 计算百分比
    counts['percentage'] = (counts['counts'] / counts['total_counts']) * 100

    # 绘制violinplot
    plt.figure(figsize=(12, 6))
    sns.violinplot(x='new_labels', y='percentage', hue=label, data=counts, split=True)
    plt.title('Cell Type Percentage by Condition and Slice ID')

    # 对每种cell_type组合进行Mann-Whitney U检验，并输出p-value（可选）
    cell_types = counts[label].unique()
    for cell_type in cell_types:
        group_a = counts[(counts[label] == cell_type) & (counts['new_labels'] == 0)]['percentage']
        group_b = counts[(counts[label] == cell_type) & (counts['new_labels'] == 1)]['percentage']
        u_stat, p_value = mannwhitneyu(group_a, group_b, alternative='two-sided', nan_policy='omit')
        print(f'{cell_type} P-value: {p_value:.5f}')

    plt.tight_layout()
    plt.show()



import numpy as np
import tensorly as tl
from tensorly.decomposition import non_negative_tucker, non_negative_tucker_hals
import time
from tensorly.metrics.regression import RMSE
import matplotlib.pyplot as plt
import seaborn as sns


def tensor_decompose(train_adata):

    df1 = train_adata[train_adata.obs['new_labels'] == 0].obs 

    grouped = df1.groupby(['cellAnnotation', 'Label']).size().unstack(fill_value=0)

    # Calculate percentages
    percentages1 = grouped.div(grouped.sum(axis=1), axis=0) * 100

    df2 = train_adata[train_adata.obs['new_labels'] == 1].obs 

    grouped = df2.groupby(['cellAnnotation', 'Label']).size().unstack(fill_value=0)

    # Calculate percentages
    percentages2 = grouped.div(grouped.sum(axis=1), axis=0) * 100


    values = np.stack([percentages1.values.T,
    np.concatenate([np.zeros(percentages2.shape[1]).reshape(-1,1), percentages2.values.T], axis=1)], axis=0)


    tensor = tl.tensor(values, dtype='float')


    core, tucker_factors = non_negative_tucker(tensor, rank=[2, 4, 2], tol=1e-12, n_iter_max=100)

    plt.figure(figsize=(10, 8))

    row_labels = list(percentages2.columns)
    col_labels = ['Factor_{}'.format(i) for i in range(4)]

    #sns.heatmap(tucker_factors[1], annot=False, fmt=".2f", xticklabels=col_labels, yticklabels=row_labels, cmap='coolwarm')


    plt.figure(figsize=(10, 8))

    row_labels = list(percentages1.index.values)
    col_labels = ['Factor_{}'.format(i) for i in range(2)]

    #sns.heatmap(tucker_factors[2], annot=False, fmt=".2f", xticklabels=col_labels, yticklabels=row_labels, cmap='coolwarm')


    plt.figure(figsize=(10, 8))

    row_labels = ['0', '1']
    col_labels = ['Factor_{}'.format(i) for i in range(2)]

    #sns.heatmap(tucker_factors[0], annot=False, fmt=".2f", xticklabels=col_labels, yticklabels=row_labels, cmap='coolwarm')


    plt.figure(figsize=(10, 8))

    row_labels = ['CT_Factor_{}'.format(i) for i in range(4)]
    col_labels = ['CN_Factor_{}'.format(i) for i in range(2)]

    #sns.heatmap(core[0], annot=False, fmt=".2f", xticklabels=col_labels, yticklabels=row_labels, cmap='coolwarm')


    plt.figure(figsize=(10, 8))

    row_labels = ['CT_Factor_{}'.format(i) for i in range(4)]
    col_labels = ['CN_Factor_{}'.format(i) for i in range(2)]

    #sns.heatmap(core[1], annot=False, fmt=".2f", xticklabels=col_labels, yticklabels=row_labels, cmap='coolwarm')
    return core, tucker_factors


import pandas as pd
import matplotlib.pyplot as plt


def plot_perc(df, condition, label):

    # Group by 'condition' and 'cell_type', count occurrences, unstack for easier percentage calculation
    grouped = df.groupby([condition, label]).size().unstack(fill_value=0)

    # Calculate percentages
    percentages = grouped.div(grouped.sum(axis=1), axis=0) * 100

    # Set figure size (e.g., 10 inches wide by 6 inches tall)

    # Plot
    percentages.plot(kind='bar', stacked=True, figsize=(20, 10))
    plt.ylabel('Percentage')
    plt.title('Percentage of Cell Types by Condition')
    plt.legend(title=label)
    plt.xticks(rotation=0)  # Keeps the condition labels horizontal for readability
    plt.show()
    plt.close()
