import matplotlib.pyplot as plt

from palettable.cartocolors.sequential import agSunset_7

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 20  # Adjust title size
plt.rcParams['axes.labelsize'] = 15

def plot_categorical_scatter(ax, x, y, groups, colors_markers, figsize=(5, 3), label=None):
    """
    Plots a scatter graph where color and marker type are determined by categorical groups.

    :param x: List or array of x coordinates
    :param y: List or array of y coordinates
    :param groups: List or array of group labels (categorical data)
    :param colors_markers: Dictionary mapping groups to (color, marker) tuples
    """
    
    unique_groups = set(groups)

    for group in unique_groups:
        # Select data for the current group
        ix = [i for i, g in enumerate(groups) if g == group]
        color, marker,s = colors_markers[group]
        sc = ax.scatter([x[i] for i in ix], [y[i] for i in ix], label=group, color=color, marker=marker, alpha=0.7, s=s)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis(False)
    #plt.legend().remove()
    return sc


def draw_cons_plot(adata_plot, slice_id ='slice_id', key='Taichi'):

    all_handles_labels = []

    for s in adata_plot.obs['slice_id'].unique():
        adata = adata[adata_plot.obs['slice_id'] == s]
        num_methods = 1
        # Create a subplot grid
        fig, axs = plt.subplots(1, num_methods, figsize=(6 * num_methods, 5))
        axs = [axs]  # Make axs iterable if there's only one subplot

        x = adata.obsm['spatial'][:,0]
        y = adata.obsm['spatial'][:,1]
        groups = adata.obs[f'Key_cons'].values
    
        # Define colors and markers for each group
        colors_markers = {
            'TP': ('#FEDFB8', 'o', 10),  
            'TN': ('#B2E08A', 'o', 10), 
            'FP': ('#E31A1C', '^', 100),
            'FN':('#6A3D9A', '*', 100) 
        }
        sc = plot_categorical_scatter(ax, x, y, groups, colors_markers=colors_markers, label=key)

        g = [ax]

        for i, ax in enumerate(g):
            #ax.get_legend().remove()
            handles, labels = ax.get_legend_handles_labels()
            for handle, label in zip(handles, labels):
                if label not in [l for _, l in all_handles_labels]:
                    all_handles_labels.append((handle, label))