### Load libs

import shutil
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import celloracle as co
import os
from tqdm.auto import tqdm



## Load input data
def load_data(adata_fname, base_GRN_fname, verbose=True):
    '''
    Loads data and base GRN
    
    Parameters
    ----------
    adata_fname : str
        Path to adata file
    base_GRN_fname : str
        Path to base GRN file
    
    Returns
    -------
    adata : AnnData
        Annotated data
    base_GRN : pd.DataFrame
        Base GRN
    '''
    if verbose: print('Loading adata')
    adata = sc.read_h5ad(adata_fname)
    if verbose: print('Loading GRN')
    base_GRN = pd.read_parquet(base_GRN_fname)
    return (adata, base_GRN)

### Prepare Oracle
## Oracle object
def prepare_oracle(adata,
                   base_GRN,
                   cluster_column_name='cell.type',
                   embedding_name='X_draw_graph_fa', 
                   max_n_pcs= 50,
                   seed=1997,
                   verbose=True):
    '''
    Prepare oracle object.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    base_GRN: pandas.DataFrame
        Dataframe with the links between genes and cell lines.
    cluster_column_name: str, optional
        Name of the column in the dataframe that contains the cell type.
    embedding_name: str, optional
        Name of the embedding in the AnnData object to use in the network fit.
    max_n_pcs: int, optional
        Maximum number of PCs to use in the network fit. Default is 50.
    seed: int, optional
        Seed for reproducibility. Unused

    Returns
    -------
    oracle: celloracle.Oracle
        Oracle object.

    '''
    
    
    v = verbose
    
    if v: pbar = tqdm(total=4)
    
    
    
    # In this notebook, we use the unscaled mRNA count for the input of Oracle object.
    if v: print(f'''[1] Instatiating Oracle with raw RNA-Seq data
    Expression range: {pd.DataFrame(adata.layers["raw_count"].sum(1)).describe()}''')
    
    adata.X = adata.layers["raw_count"].copy()

    # Instantiate Oracle object.
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata,
                                       cluster_column_name=cluster_column_name,
                                       embedding_name=embedding_name)
    if v: pbar.update(1)

    
    
    # You can load TF info dataframe with the following code.
    if v: print(f'[2] Adding GRN data to Oracle')
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    if v: pbar.update(1)

    
    # Perform PCA
    if v: print(f'[3] Computing and selecting PCs and k for neighbors')
    oracle.perform_PCA()
    # Select important PCs
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    n_comps = min(n_comps, max_n_pcs)
    if v: pbar.update(1)
    
    
    n_cell = oracle.adata.shape[0]
    k = int(0.025*n_cell)
    if v: print(f'''[4] Computing KNNs:
    n cells: {n_cell}
    k: {k}
    oracle.knn_imputation(n_pca_dims=n_comps, k={k}, balanced=True, b_sight={k*8}, b_maxl={k*4}, n_jobs=-1)''')
    
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                          b_maxl=k*4, n_jobs=-1)
    if v: pbar.update(1)
    
    return(oracle)    

### Carto scores


def run_r(r_cmd, env_path='~/venvs/cicero'):
    from IPython import get_ipython
    
    r = get_ipython().run_cell_magic(
        'bash', 
        None,
        f'''\
eval "$(micromamba shell hook --shell=bash)"
micromamba activate {env_path}
which R
''' + r_cmd)

    return(r)

carto_r_script = """\
library(igraph)
library(rnetcarto)

get_scores_cartography <- function(g){
  adj <- get.adjacency(g,sparse=FALSE)
  score <- netcarto(adj)[[1]]
  rownames(score) <- score$name
  return(score)
}

carto_score <- function(d){
    g <- graph.data.frame(d[1:2], directed = T)
    E(g)$weight <- d[[3]]
    return(get_scores_cartography(g))
}

d <- data.table::fread('links.carto.tsv', sep='\t', data.table=FALSE)
# print(head(d))
rownames(d) <- d[,1]
d[,1] <- NULL
# print(head(d))
carto.res <- carto_score(d)
tryCatch({
    data.table::fwrite(carto.res, './results.carto.tsv', sep='\t')
}, error = function(cond){
    print(as.character(cond))
    print(getwd())
    
    write.csv(carto.res, paste0(getwd(), '/results.carto.tsv'), sep='\t')
})
"""

def cluster_carto(arg):
    
    wdir = str(np.random.choice(list(range(0,100))))
    while os.path.exists(wdir):
        wdir = str(np.random.choice(list(range(0,100))))
        
    os.makedirs(wdir, exist_ok=True)
    os.chdir(wdir)
    
    
    cluster, df = arg
    df.to_csv('links.carto.tsv', sep='\t', index=True, header=True)
    with open('carto.script.R', 'wt') as f:
        f.writelines(carto_r_script)
        
    print(os.getcwd())
    run_r(r_cmd='Rscript --vanilla carto.script.R')
    res = pd.read_csv('results.carto.tsv', sep='\t', index_col=0)
    order = ["Ultra peripheral", "Peripheral", "Connector","Kinless","Provincical Hub","Connector Hub", "Kinless Hub"]
    res['Role name'] = res['role'].replace(dict([(i+1, o) for i,o in enumerate(order)]))
    res['cluster'] = cluster
    os.chdir('..')
    shutil.rmtree(wdir)
    
    return(res)

def get_links_carto(links):
    from tqdm.contrib.concurrent import process_map  # or thread_map
    res = pd.concat(process_map(cluster_carto,
                                links.filtered_links.items(), 
                                max_workers=20), 0)

    if 'merged_score' in dir(links):
        mscore = links.merged_score
        mscore['name'] = mscore.index
        
        return(pd.merge(mscore, res, on=['cluster', 'name'], how='outer').set_index('name'))
    else:
        return(res)

### Network scores

def get_all_link_scores(links, tempdir='/scratch/xoel/tmp/'):

    old_dir = os.getcwd()
    wdir = os.path.join(tempdir, str(np.random.choice(list(range(0,100)))))
    # Check preexisstintg temp directory
    while os.path.exists(wdir):
        wdir = os.path.join(tempdir, str(np.random.choice(list(range(0,100)))))
    
    # Create temp directory
    os.makedirs(wdir, exist_ok=True)
    os.chdir(wdir)

    # Calculate network scores. It takes several minutes.
    links.get_network_score()
    links.merged_score = get_links_carto(links)
    
    
    os.chdir(old_dir)
    
    return(links)

### Pipe


def network_fit_pipe(
    
    adata_fname, base_GRN_fname,
    
    cluster_column_name='cell.type', 
    
    embedding_name='X_draw_graph_fa',
    max_n_pcs=50,
    
    links_alpha=10,
    links_pval=0.001,
    links_top_n_coef=2000,
    
    save_dir='./',
    prefix=None,
    tempdir='/scratch/xoel/tmp/',
    
    exclude=None,
    
    verbose=True,
    test_mode=False,
    seed=1997 ):
    
    
    
    v = verbose
    if v: pbar = tqdm(total=6 if not save_dir else 7)
    
    # Load data
    if v: print(f'[1] Reading data')
    adata, base_GRN = load_data(adata_fname, base_GRN_fname)
    if v: pbar.update(1)
    
    
    # Exclude groups of cells (if failed previously)
    if not exclude is None:
        if v: print(f'- Excluding {exclude}')
        adata = adata[~adata.obs[cluster_column_name].isin(exclude)].copy()
        
        
    
    # Prepare oracle
    if v: print(f'[2] Prepairing oracle object')
    oracle = prepare_oracle(adata,
                            base_GRN, 
                            embedding_name=embedding_name,
                            cluster_column_name=cluster_column_name, 
                            max_n_pcs=max_n_pcs,
                            seed=seed)
    if v: print(oracle)
    if v: pbar.update(1)


    # Compute links
    if v: print(f'[3] Getting links')
    links = oracle.get_links(
        cluster_name_for_GRN_unit=cluster_column_name,
        alpha=links_alpha,
        verbose_level=2*int(v), 
        test_mode=test_mode, 
        n_jobs=-1)
    if v: print(links)
    if v: pbar.update(1)

    
    # Check result or recall
    if v: print(f'[4] Checking links')
    
    no_fits = np.array(links.cluster)[[d['p'].isna().all() for _, d in links.links_dict.items()]].tolist()
    if len(no_fits) > 0:
        if v: print(f'WARNING: Repeating analysis without {no_fits}')
        return( network_fit_pipe( adata_fname = adata_fname,
                                  base_GRN_fname = base_GRN_fname,
                                  cluster_column_name = cluster_column_name, 
                                  
                                  max_n_pcs = max_n_pcs,
                                  links_alpha = links_alpha,
                                  links_pval = links_pval,
                                  links_top_n_coef = links_top_n_coef,
                                  
                                  verbose = verbose,
                                  seed = seed,
                                  tempdir = tempdir,
 
                                  exclude = no_fits ))
    else:
        if v: print('All fits are valid.')

    if v: pbar.update(1)

    
    # Post-process links
    if v: print(f'[5] Postprocessing links')
    ## Filtering 
    links.filter_links(p=links_pval, weight="coef_abs", threshold_number=links_top_n_coef)    
    if v: pbar.update(1)

    
    if v: print(f'[6] Scoring genes')
    links = get_all_link_scores(links, tempdir=tempdir)
    if v: pbar.update(1)

    
    # Save results
    if save_dir:
        if v: print(f'[7] Saving results')
        cluster_dir = os.path.join(save_dir, 'cluster_GRN')
        os.makedirs(cluster_dir, exist_ok=True)

        ## Save raw links table per cluster
        if v: print('- Saving raw links per cluster')
        for cluster, links_df in tqdm(links.links_dict.items()):

            fname = f"{cluster_dir}/" + (prefix+'.' if prefix else '')
            fname += f"{cluster}.raw_GRN.alpha={links_alpha}.csv"
            links_df.to_csv(fname)


        ## Save filtered links table per cluster
        if v: print('- Saving filtered links per cluster')
        for cluster, links_df in tqdm(list(links.filtered_links.items())):

            fname = f"{cluster_dir}/" + (prefix+'.' if prefix else '')
            fname += f"{cluster}.filtered_GRN.pval={links_pval}.top={links_top_n_coef}.csv"
            links_df.to_csv(fname)

        ## Save links object
        if v: print('- Saving processed links')
        links_fname = (prefix+'.' if prefix else '')+ 'Links.celloracle.links'
        links.to_hdf5(file_path=os.path.join(save_dir, links_fname))

        ## Save oracle object
        if v: print('- Saving oracle')
        oracle_fname = (prefix+'.' if prefix else '')+ 'Oracle.celloracle.oracle'
        oracle.to_hdf5(file_path=os.path.join(save_dir, oracle_fname))
        
        
        if v: pbar.update(1)

    return(oracle, links)


### Plots

def plot_cartography_term(links, goi, save=None, plt_show=True):
    """
    Plot the summary of gene network cartography like a heatmap.
    Please read the original paper of gene network cartography for detail.
    https://www.nature.com/articles/nature03288

    Args:
        links (Links object): See network_analisis.Links class for detail.
        gois (list of srt): List of Gene name to highlight.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    plt.close('all') 
    order = ["Ultra peripheral", "Peripheral", "Connector","Kinless","Provincical Hub","Connector Hub", "Kinless Hub"]

    # print(goi)
    
    tt = pd.get_dummies(links.merged_score[["cluster", "role"]],columns=["role"])
    tt = tt.loc[goi].set_index("cluster")
    tt.columns = [i.replace("role_", "") for i in tt.columns]
    tt = tt.reindex(index=links.palette.index.values,columns=[str(i+1) for i in range(len(order))]).fillna(0)
    
    tt.columns = order
    plt.rcParams["figure.figsize"] = [5, 15]
    sns.heatmap(data=tt, cmap="Blues", cbar=False)
    if not save is None:
        os.makedirs(save, exist_ok=True)
        path = os.path.join(save,
                           f"cartography_role_in_{links.name}_{links.threshold_number}_{goi}.pdf")
        plt.savefig(path, transparent=True, bbox_inches='tight')

    if plt_show:
        plt.show()
        

def plot_network_entropy_distributions(links, update_network_entropy=False, save=None, show=False):
    """
    Plot the distribution of network entropy.
    See the CellOracle paper for the detail.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        values (list of str): The list of netwrok score type. If it is None, all network score (listed above) will be used.
        update_network_entropy (bool): Whether to recalculate network entropy.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    if links.entropy is None:
        links.get_network_entropy()

    if update_network_entropy:
        links.get_network_entropy()

    # fig = plt.figure()

    ax = sns.boxplot(data=links.entropy, x="cluster", y="entropy_norm",
                palette=links.palette.palette.values,
                order=links.palette.index.values, fliersize=0.0)

    ax.tick_params(axis="x", rotation=90)
    ax.set_ylim([links.entropy.groupby('cluster').quantile(q=0.01)['entropy_norm'].min(), 1.0])

    if not save is None:
        os.makedirs(save, exist_ok=True)
        path = os.path.join(save, f"network_entropy_in_{links.name}_{links.threshold_number}.png")
        ax.set_ylabel("normalized\nentropy")
        plt.savefig(path, transparent=True)
    if show:
        plt.show()
    return


def network_plots(
    links,
    
    n_top_scores=30,
    
    save_dir='./network_plots/',
    verbose=True
):
    
    maxes = links.merged_score.groupby('cluster').quantile(q=0.99).max(0)
    
    os.makedirs(save_dir, exist_ok=True)
    v = verbose
    if v: pbar= tqdm(total=5)
    
    # Degree distribution
    plt.rcParams["figure.figsize"] = [9, 4.5]
    links.plot_degree_distributions(plot_model=True, 
                                    save=f"{save_dir}degree_distribution/")
    plt.close('all') 
    if v: pbar.update(1)
    
    # Top scored genes per cluster
    ## Visualize top n-th genes with high scores.
    plt.rcParams["figure.figsize"] = [6, 6]
    for cluster in tqdm(list(links.links_dict.keys())):
        # print(cluster)
        links.plot_scores_as_rank(cluster=cluster, n_gene=n_top_scores, save=f"{save_dir}ranked_score/{cluster}/")

        plt.close('all') 
    if v: pbar.update(1)

    # Boxplots per cluster
    ## Degree centrality
    plt.rcParams["figure.figsize"] = [6, 6]
    plt.subplots_adjust(left=0.15, bottom=0.3)
    plt.ylim([0,maxes['degree_centrality_all']])
    links.plot_score_discributions(values=["degree_centrality_all"], 
                                   method="boxplot", 
                                   save=f"{save_dir}")
    plt.close('all') 
    if v: pbar.update(1)

    ## Eigenvector centrality
    plt.subplots_adjust(left=0.15, bottom=0.3)
    plt.ylim([0, maxes['eigenvector_centrality']])
    links.plot_score_discributions(values=["eigenvector_centrality"], method="boxplot", save=f"{save_dir}")
    plt.close('all') 
    if v: pbar.update(1)

    ## Network entropy
    plt.subplots_adjust(left=0.15, bottom=0.3)
    plot_network_entropy_distributions(links=links, save=f"{save_dir}")
    plt.close('all') 
    if v: pbar.update(1)


    return()

def network_gene_plots(
    links,
    
    genes_of_interest='all',
    save_dir='./network_gene_plots/',
    verbose=True
):

    v = verbose
    if v: pbar = tqdm(total=3)
    
    os.makedirs(save_dir, exist_ok=True)
    if save_dir[-1] != '/': save_dir+='/'
    
    # Plot genes
    gois = list(set([g for _, d in links.filtered_links.items() for g in d[['source', 'target']].values.flatten()]))
    if v: print(f'Genes in filtered network: {len(gois)}')
    
    if isinstance(genes_of_interest, list):
        gois = [g for g in gois if g in genes_of_interest]
        if v: print(f'''
        Requested genes: {len(genes_of_interest)}
        Available: {len(gois)}''')
        
    elif genes_of_interest != 'all':
        raise()
    
    plt.rcParams["figure.figsize"] = [11, 9]
    # Visualize network score dynamics
    for goi in tqdm(sorted(gois)):
        links.plot_score_per_cluster(goi=goi, save=f"{save_dir}network_score_per_gene/")
        plt.close('all') 
    if v: pbar.update(1)


    # Cartography
    ## Per cluster
    ### Plot cartography as a scatter plot
    plt.rcParams["figure.figsize"] = [16, 9]
    links.plot_cartography_scatter_per_cluster(scatter=True,
                                               kde=False,
                                               gois=gois, # Highlight genes of interest
                                               auto_gene_annot=False,
                                               args_dot={"n_levels": 105},
                                               args_line={"c":"gray"}, 
                                               save=f"{save_dir}cartography_scatter/")
    plt.close('all') 
    if v: pbar.update(1)


    ## Per gene
    ### Plot the summary of cartography analysis
    plt.rcParams["figure.figsize"] = [4, 15]
    for goi in tqdm(gois):
        plot_cartography_term(links=links, goi=[goi], 
                              save=f"{save_dir}cartography/",
                              plt_show=False)
        plt.close('all') 
    if v: pbar.update(1)

    return()

def pipe_plots(links,
               palette=None,
    
               n_top_scores=25,
               genes_of_interest='all',
   
               do_network_plots=True,
               do_gene_plots=True,
   
               save_dir='./',
   
               verbose=True, 
               plot_format='pdf'):

    
    if not save_dir: return
   
    if plot_format:
        plt.rcParams['savefig.format'] = plot_format
        sc.set_figure_params(format=plot_format)
        co.network_analysis.gene_analysis.settings['save_figure_as'] = plot_format

    if palette is not None:
        if not links.palette.index.isin(palette.index).all():
            print('Invalid palette provided. Missing: ', links.palette.index[~links.palette.index.isin(palette.index)])

        else:
            links.palette = palette[palette.index.isin(links.palette.index)].rename({palette.columns[0]:'palette'},axis=1)

    if do_network_plots:
        network_plots(links, n_top_scores=n_top_scores,
                      save_dir=os.path.join(save_dir, 'network_plots/'),
                      verbose=verbose)
        
    if do_gene_plots:
        network_gene_plots(links, genes_of_interest=genes_of_interest,
                           save_dir=os.path.join(save_dir, 'network_plots/'),
                           verbose=verbose)
    return()