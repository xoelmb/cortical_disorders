import os
import scanpy as sc
import numpy as np
import pandas as pd

from tqdm.auto import tqdm

from celloracle.applications import Pseudotime_calculator


def preprocess_adata(
    adata,
    min_total_counts_per_cell=1,
    key_n_counts='total_counts',
    n_top_genes=3000,
    verbose=True,
    copy=True,
    genes_to_keep=None
):
    
    v = verbose
    if v: pbar = tqdm(total=8)
    if copy: adata = adata.copy()
    
    if v: print('[1] Saving raw counts')
    # Keep raw
    adata.raw = adata.copy()
    adata.layers["raw_count"] = adata.raw.X.copy()
    if v: pbar.update(1)
    
    if v: print('[2] Computing QC')
    # Compute automatic metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    if v: pbar.update(1)
    
    if v: print(f'[3] Filtering genes: total_counts > {min_total_counts_per_cell}')
    # Only consider genes with more than 1 count
    # sc.pp.filter_genes(adata, min_counts=1)  # 20 min
    adata = adata[:,adata.var['total_counts'] > min_total_counts_per_cell]
    if v: pbar.update(1)
    
    if v: print(f'[4] Normalizing per cell using {key_n_counts}')
    # Normalize
    sc.pp.normalize_per_cell(adata, 
                             key_n_counts=key_n_counts)
    if v: pbar.update(1)
    
    if v: print(f'[5] Computing highly variable genes and filtering to top {n_top_genes}')
    # Select top 2000 highly-variable genes
    print(adata.shape) # rem
    filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                  flavor='cell_ranger',
                                                  n_top_genes=n_top_genes,
                                                  log=False)
    sc.pl.filter_genes_dispersion(filter_result, log=True)
    filter_mask = pd.DataFrame(filter_result, index=adata.var_names)['gene_subset']

    # If genes to keep:
    if genes_to_keep is None:
        pass
    else: 
        if v: print(f'Keeping provided genes if available...')
        filter_mask = filter_mask | filter_mask.index.isin(genes_to_keep) 
    
    adata = adata[:,filter_mask]
    if v: pbar.update(1)
    print(adata.shape) # rem

    if v: print(f'[6] Normalizing per cell again using {key_n_counts}')
    # Renormalize after filtering
    sc.pp.normalize_per_cell(adata, key_n_counts=key_n_counts)
    if v: pbar.update(1)
    
    if v: print(f'[7] Log1p transformation')
    # Log transformation and scaling
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    if v: pbar.update(1)
    
    if v: print(f'[8] Scaling')
    sc.pp.scale(adata) # max_value 10000 not in original cell oracle
    adata.layers["scaled"] = adata.X.copy()
    if v: pbar.update(1)
    
    return(adata)


def run_phate(adata,
              run_on='X',
              copy=True,
              **operator_kwargs):
    
    if copy:
        adata = adata.copy()
        
    import phate
    phate_operator = phate.PHATE(n_jobs=-2)
    phate_operator.set_params(**{k: v for k, v in operator_kwargs.items() if not v is None})
    
    if run_on == 'X':
        print('running on X layer')
        adata.obsm['X_phate'] = phate_operator.fit_transform(adata.X)
    elif run_on in adata.obsm.keys():
        print(f'Running on {run_on} obsm')
        adata.obsm['X_phate'] = phate_operator.fit_transform(adata.obsm[run_on])
        
    elif run_on in adata.layers.keys():
        print(f'Running on {run_on} layer')
        adata.obsm['X_phate'] = phate_operator.fit_transform(adata.layers[run_on])
        
    return(adata)


def embed_adata(
    adata,
    nn_1_n_neighbors=5, 
    nn_1_use_rep='X_pca',
    nn_1_n_pcs=30,
    dm_n_comps=20,
    nn_2_n_neighbors=50,
    phate_run_on='X',
    phate_knn=None,
    phate_n_pca=None,
    phate_decay=None,
    phate_t=None,
    group_by='cell.type',
    skip_umap=True,
    skip_paga=True,
    skip_phate=False,
    random_state=1997,
    verbose=True,
    copy=True
):

    if copy: adata = adata.copy()
    
    v = verbose
    if v: pbar = tqdm(total=8)
    
    if v: print('[1] Computing PCA')
    # PCA
    sc.tl.pca(adata, svd_solver='arpack', random_state=random_state)
    sc.pl.pca(adata, color=group_by, legend_loc='on data')
    if v: pbar.update(1)
    
    if v: print('[2] Computing neighbors')
    # Diffusion map
    sc.pp.neighbors(adata, 
                    n_neighbors=nn_1_n_neighbors,
                    use_rep=nn_1_use_rep,
                    n_pcs=nn_1_n_pcs, 
                    random_state=random_state)
    if v: pbar.update(1)
    
    if v: print('[3] Computing diffusion map')
    sc.tl.diffmap(adata, 
                  n_comps=dm_n_comps)
    if v: pbar.update(1)
    
    if v: print('[4] Computing neighbors in diffusion map')
    # Calculate neighbors again based on diffusionmap 
    sc.pp.neighbors(adata, 
                    n_neighbors=nn_2_n_neighbors,
                    use_rep='X_diffmap',
                    random_state=random_state)
    if v: pbar.update(1)
    
    if not skip_umap:
        if v: print('[5] Constructing UMAP')
        # PAGA graph construction z
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=group_by, legend_loc='on data')
    if v: pbar.update(1)
    
    if not skip_paga:
        if v: print('[6] Constructing PAGA graph (1/2)')
        # PAGA graph construction z
        sc.tl.paga(adata, groups=group_by)
        sc.pl.paga(adata, color=group_by)
    if v: pbar.update(1)

    if not skip_paga:
        if v: print('[7] Constructing PAGA graph (2/2)')
        sc.tl.draw_graph(adata, init_pos='paga', random_state=random_state)
        sc.pl.draw_graph(adata, color=group_by, legend_loc='on data')
    if v: pbar.update(1)
        
    if not skip_phate:
        if v: print('[8] Constructing PHATE graph')
        ## PHATE
        print(phate_run_on)
        adata = run_phate(adata, run_on=phate_run_on, knn=phate_knn, t=phate_t, decay=phate_decay, n_pca=phate_n_pca, verbose=verbose, copy=False)
        sc.external.pl.phate(adata, color=group_by, legend_loc='on data')
    if v: pbar.update(1)
    
    return(adata)

# Pseudotime
## Root cell selection

### Manual root cell selection
def manual_root_cell_selection(
    adata,
    embedding_key,
    group_by,
    palette_dict=None,
    save_dir=None):
    
    import plotly.express as px
    import pandas as pd

    def plot(adata, embedding_key, cluster_column_name, palette=None):
        
        embedding = adata.obsm[embedding_key]
        df = pd.DataFrame(embedding, columns=["x", "y"])
        df["cluster"] = adata.obs[cluster_column_name].values
        df["label"] = adata.obs.index.values
        
        fig = px.scatter(df, x="x", y="y", hover_name=df["label"],
                         color="cluster", color_discrete_map=palette)
        fig.show()
        return(fig)

    fig = plot(adata=adata,
         embedding_key=embedding_key,
         cluster_column_name=group_by,
         palette=palette_dict)
    
    if save_dir:
        fig.write_html(os.path.join(save_dir, 'RootCellSelection.html'))

        

def prepare_lineages(adata, root_cells=None, lineages=None, 
                     group_by='cell.type', lineage_name='Lin_1',
                     embedding_key='X_phate', save_dir=None):
    
    if lineages is None:
        lineages = {lineage_name: adata.obs[group_by].unique()}

    if root_cells is None:
        root_cells = {k: '' for k in lineages.keys()}
        
    if isinstance(root_cells, str):
        root_cells = {k: root_cells for k in lineages.keys()}
        
    for lin, lin_groups in lineages.items():
        
        if lin not in root_cells.keys():
            lin_cell = ''
        else: 
            lin_cell = root_cells[lin]
        
        plotted = False 
        
        while not lin_cell in adata.obs_names[adata.obs[group_by].isin(lin_groups)] or lin_cell == '':
            if not plotted: 
                manual_root_cell_selection(adata, 
                                           embedding_key=embedding_key,
                                           palette_dict=dict(
                                               zip(adata.obs[group_by].cat.categories, 
                                                   adata.uns[group_by+'_colors'])),
                                           group_by= group_by,
                                           save_dir=save_dir)
                plotted=True
            lin_cell = input(f'Please select valid cell for {lin}:\n{lin_groups}')
        
        root_cells[lin] = lin_cell
    return(lineages, root_cells)

def compute_pseudotime(
    adata, root_cells=None, lineages=None, 
    group_by='cell.type', lineage_name='Lin_1',
    embedding_key='X_phate', save_dir=None):

    lineages, root_cells = prepare_lineages(
        adata, root_cells=root_cells, lineages=lineages, 
        group_by=group_by, lineage_name=lineage_name,
        embedding_key=embedding_key, save_dir=save_dir)
    
    # Instantiate pseudotime object using anndata object.
    pt = Pseudotime_calculator(adata=adata,
                               obsm_key=embedding_key, # Dimensional reduction data name
                               cluster_column_name=group_by # Clustering data name
                              )

    # Input lineage information into pseudotime object
    pt.set_lineage(lineage_dictionary=lineages)
    # Estimated root cell name for each lineage
    pt.set_root_cells(root_cells=root_cells)

    # Check root cell and lineage
    pt.plot_root_cells()

    # Calculate pseudotime
    pt.get_pseudotime_per_each_lineage()

    # Check results
    pt.plot_pseudotime(cmap="rainbow")

    if sum(pt.adata.obs['Pseudotime'] > 0.9) < 10:
        print("There are fewer than 10 cells with high values of pseudotime. They are outliers most likely.")

    sc.pl.violin(pt.adata, 'Pseudotime', groupby=group_by)
    
    # Add calculated pseudotime data to the oracle object
    adata.obs = pt.adata.obs
    
    return(adata, pt)

def prepare_adata(
    
    adata,
    group_by='cell.type',
    save_dir=None,
    prefix=None,
    
    min_total_counts_per_cell=1,
    key_n_counts='total_counts',
    n_top_genes=3000,
    genes_to_keep=None,

    nn_1_n_neighbors=5, 
    nn_1_use_rep='X_scvi',
    nn_1_n_pcs=30,
    dm_n_comps=20,
    nn_2_n_neighbors=50,
    phate_run_on='X_scvi',
    phate_knn=None,
    phate_n_pca=None,
    phate_decay=None,
    phate_t=None,
    skip_umap=True,
    skip_paga=True,
    skip_phate=False,
    random_state=1997,    

    root_cells='hft_w16_p7_r2_GGTGATTCATGACTGT',
    lineages=None, 
    lineage_name='Lin_1',
    embedding_key='X_phate',
    verbose=True):
    
    
    
    adata = adata.copy()
    
    if verbose: print('### PIPE [1] PREPROCESSING DATA')
    adata = preprocess_adata(
        adata,
        min_total_counts_per_cell=min_total_counts_per_cell,
        key_n_counts=key_n_counts,
        n_top_genes=n_top_genes,
        verbose=verbose,
        copy=False,
        genes_to_keep=genes_to_keep
    )

    if verbose: print('### PIPE [2] EMBEDDING DATA')
    adata = embed_adata(
        adata,
        group_by=group_by,
        
        nn_1_n_neighbors=nn_1_n_neighbors,
        nn_1_use_rep=nn_1_use_rep,
        nn_1_n_pcs=nn_1_n_pcs,
        dm_n_comps=dm_n_comps,
        nn_2_n_neighbors=nn_2_n_neighbors,
        
        phate_run_on=phate_run_on,
        phate_knn=phate_knn,
        phate_n_pca=phate_n_pca,
        phate_decay=phate_decay,
        phate_t=phate_t,
        
        skip_paga=skip_paga,
        skip_phate=skip_phate,
        skip_umap=skip_umap,
        
        random_state=random_state,
        verbose=verbose,
        copy=False
    )

    if verbose: print('### PIPE [3] COMPUTING PSEUDOTIME')
    adata, pt = compute_pseudotime(
        adata,
        root_cells=root_cells, 
        lineages=lineages, 
        group_by=group_by,
        lineage_name=lineage_name,
        embedding_key=embedding_key,
        save_dir=save_dir
    )
    

    if save_dir:

        os.makedirs(save_dir, exist_ok=True)
        if verbose: print('### PIPE [4] SAVING DATA')
        
        fname = 'RNA.processed.h5ad'
        if prefix: fname = prefix+'.'+fname

        fname = os.path.join(save_dir, fname)
        adata.write_h5ad(fname)

    return(adata)