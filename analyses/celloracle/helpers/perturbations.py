import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.backends.backend_pdf import PdfPages

from tqdm.auto import tqdm

import scanpy as sc

import celloracle as co
from celloracle.applications import Oracle_development_module, Oracle_systematic_analysis_helper
from celloracle.applications import Gradient_calculator

import pandas as pd
import numpy as np

### Prepare data

### Update palette
def update_palette(oracle, links, palette, group_by):
    
    # Update oracle
    if not all(oracle.adata.obs[group_by].unique().isin(palette.index)):
        raise('Not all groups are in palette')
        
    else:
        oracle.update_cluster_colors(palette[palette.index.isin(oracle.adata.obs[group_by].cat.categories)].iloc[:,0].to_dict())
        
        oracle.adata.obs[group_by] = oracle.adata.obs[group_by].astype('str')
        oracle.adata.obs[group_by] = pd.Categorical(oracle.adata.obs[group_by], 
                                                categories=palette.loc[palette.index.isin(oracle.adata.obs[group_by].unique()),:].index, ordered=True)
        oracle.adata.uns[f'{group_by}_colors'] = palette[palette.index.isin(oracle.adata.obs['cell.type'].cat.categories)].iloc[:,0].tolist()
        
    # Update links
    if not links.palette.index.isin(palette.index).all():
        raise('Links: Invalid palette provided. Missing: ', links.palette.index[~links.palette.index.isin(palette.index)])

    else:
        links.palette = palette[palette.index.isin(links.palette.index)].rename({palette.columns[0]:'palette'},axis=1)

    return(oracle, links)

### Check pseudotime
def plot_pseudotime(oracle):
    

    # Visualize pseudotime
    fig, ax = plt.subplots(figsize=[6,6])

    sc.pl.embedding(adata=oracle.adata,
                    edges=False, # connecting cells edges network neighbors
                    arrows=False, # requires `'velocity_X_phate'` from scvelo or `'Delta_X_phate'` from velocyto
                    basis=oracle.embedding_name, 
                    ax=ax, 
                    size=50,
                    # frameon=True,
                    cmap="RdBu_r",
                    color=["Pseudotime"], 
                    add_outline=True,
                    outline_width=(0.01,0.07),
                    show=True
                   )
    return(fig)
    



### Create grid
def compute_p_mass(gradient,

                   p_mass_smooth = 0.8,
                   p_mass_n_grid = 40,
                   p_mass_n_neighbors = 200,
                   p_mass_filter_min = None      
                   
                  ):
    ans = ''
    
    while not p_mass_filter_min:
        print(f'Min mass suggestions for `smooth`: {p_mass_smooth} `n_grid`: {p_mass_n_grid} `n_neighbors`: {p_mass_n_neighbors}')
        gradient.calculate_p_mass(smooth=p_mass_smooth, n_grid=p_mass_n_grid, n_neighbors=p_mass_n_neighbors)
        gradient.suggest_mass_thresholds(n_suggestion=12)
        plt.show()
    
        ans = input('Select min mass or leave empty to repeat suggestions.')
        if not ans:
            p_mass_smooth = float(input(f'Select new `smooth`: {p_mass_smooth}') or p_mass_smooth)
            p_mass_n_grid = int(input(f'Select new `n_grid`: {p_mass_n_grid}') or p_mass_n_grid)
            p_mass_n_neighbors = int(input(f'Select new `n_neighbors`: {p_mass_n_neighbors}') or p_mass_n_neighbors )
        else:    p_mass_filter_min = float(ans)
    
    gradient.calculate_p_mass(smooth=p_mass_smooth, n_grid=p_mass_n_grid, n_neighbors=p_mass_n_neighbors)
    gradient.calculate_mass_filter(min_mass=p_mass_filter_min, plot=True)
    
    return(gradient)




### Transfer pseudotime values to the grid points
def transfer_pseudotime(gradient,
                        method='both',
                        n_knn=50,
                        n_poly=3,
                        scale_dev = 40,
                        s=5
                        ):
    
    if method in ['knn', 'both']:
        print('KNN pseudotime-transfer model')
        gradient.transfer_data_into_grid(args={"method": "knn",
                                               "n_knn": n_knn}, plot=True)
    
        plt.show()
        
    if method in ['polynomial', 'both']: # Please use this method if knn method does not work.
        print('Polynomial pseudotime-transfer model')
        gradient.transfer_data_into_grid(args={"method": "polynomial",
                                               "n_poly": n_poly}, plot=True)
        plt.show()
    
    if method != 'both':
        # Calculate graddient
        gradient.calculate_gradient()
        
        # Show results
        gradient.visualize_results(scale=scale_dev, s=s)

        # print('Developmental flow')
        fig, ax = plt.subplots(figsize=[6, 6])
        gradient.plot_dev_flow_on_grid(scale=scale_dev, ax=ax) 
        plt.show()
        
        return(gradient)
    
    
    else: 
        method = None
        while not method in ['knn', 'polynomial']:
            method = input('Select method [ knn | polynomial ]')
        return(transfer_pseudotime(gradient,
                                   method=method,
                                   n_knn=n_knn,
                                   n_poly=n_poly))

def prepare_data(
    oracle_fname, 
    links_fname,
    
    palette = None,
    group_by = None,
    
    grn_fit_ridge_alpha = 10,
    
    p_mass_smooth = 0.8,
    p_mass_n_grid = 40,
    p_mass_n_neighbors = 200,
    p_mass_filter_min = None,   

    method='both',
    n_knn=50,
    n_poly=3,
    scale_dev = 40,
    s=5,
    
    save_dir=None,
    prefix=None,
    verbose=True):
    
    v = verbose
    if v: pbar=tqdm(total=6+int(not save_dir is None))
    
    if v: print('[1] Loading data')
    oracle = co.load_hdf5(oracle_fname)
    links = co.load_hdf5(links_fname)
    
    group_by = group_by or oracle.cluster_column_name
    
    
    if palette:
        if v: print('Updating palette')
        oracle, links = update_palette(oracle=oracle, links=links, palette=palette, group_by=group_by)
        
    if v: pbar.update(1)
        
        
        
    if v: print('[2] Re-filtering links')    
    links.filter_links()
    if v: pbar.update(1)
    if v: print('[3] Getting TF dict')    
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    if v: pbar.update(1)
    if v: print('[4] Fitting for simulation')    
    oracle.fit_GRN_for_simulation(alpha=grn_fit_ridge_alpha, use_cluster_specific_TFdict=True)
    if v: pbar.update(1)

    
    plot_pseudotime(oracle)
    
    if v: print('[5] Computing grid p_mass')    
    # Instantiate Gradient calculator object
    gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
    
    gradient = compute_p_mass(
        gradient,
        p_mass_smooth=p_mass_smooth,
        p_mass_n_grid=p_mass_n_grid,
        p_mass_n_neighbors=p_mass_n_neighbors,
        p_mass_filter_min=p_mass_filter_min
    )
    if v: pbar.update(1)
    
    if v: print('[6] Transferring pseudotime to grid')    
    gradient = gradient = transfer_pseudotime(
        gradient, 
        method=method, 
        n_knn=n_knn, 
        n_poly=n_poly,
        scale_dev=scale_dev,
        s=s)
    if v: pbar.update(1)
    
    
    if save_dir:
        if v: print('[7] Saving data')    
       
        fname = os.path.join(save_dir, prefix+'.' if prefix else '')
            
        if v: print('\tLinks')
        links.to_hdf5(fname+'Links.Perturbation.celloracle.links')
        if v: print('\tOracle')
        oracle.to_hdf5(fname+'Oracle.Perturbation.celloracle.oracle')
        if v: print('\tGradient')
        gradient.to_hdf5(fname+'Gradient.Perturbation.celloracle.gradient')   
        if v: pbar.update(1)
    
    
    return(oracle, links, gradient)

### Computations

def compute_perturbation(oracle, gene, exp_value, 
                         
                         shift_n_propagation=3,
                         transprob_n_neighbors=None,
                         embedding_shift_sigma=0.05,
                         
                         p_mass_smooth=0.8,
                         p_mass_n_grid=40,
                         p_mass_n_neighbors=200,
                         p_mass_filter_min=None,
                         
                         markov_n_steps=500,
                         markov_n_duplications=5,
                         
                         n_cores=16,
                         verbose=True
                         ):
    

    v = verbose
    if v: pbar = tqdm(total=5, colour='lightgreen')
    
    if v: print('Calculating expression shift')
    oracle.simulate_shift(perturb_condition={gene: exp_value},
                          GRN_unit='cluster',
                          n_propagation=shift_n_propagation,
                          ignore_warning=True)
    if v: pbar.update(1)

    if v: print('Estimating transition probabilities')
    oracle.estimate_transition_prob(n_neighbors=transprob_n_neighbors,  ## CUSTOMIZABLE
                                    knn_random=True, 
                                    sampled_fraction=1,
                                    n_jobs=n_cores,
                                    threads=n_cores)
    if v: pbar.update(1)

    if v: print('Calculating embedding shift')
    oracle.calculate_embedding_shift(sigma_corr=embedding_shift_sigma)  ## CUSTOMIZABLE
    if v: pbar.update(1)

    
    if v: print('Calculating p_mass')
    # oracle.calculate_p_mass(smooth=p_mass_smooth, n_grid=p_mass_grid, n_neighbors=p_mass_n_neighbors)
    # oracle.calculate_mass_filter(min_mass=p_mass_min, plot=True)
    oracle = compute_p_mass(gradient=oracle, 
                            p_mass_smooth=p_mass_smooth,
                            p_mass_n_grid=p_mass_n_grid,
                            p_mass_n_neighbors=p_mass_n_neighbors,
                            p_mass_filter_min=p_mass_filter_min)
    if v: pbar.update(1)

    if v: print('Running Markov chain simulation')
    oracle.run_markov_chain_simulation(n_steps=markov_n_steps, n_duplication=markov_n_duplications)
    if v: pbar.update(1)
    
    return(oracle)

### Perturbation plots

### Sankey

def plot_markov_celltype_simulation(oracle, 
                                    order=None, 
                                    cluster_use='cell.type'):
    
    order = order if not order is None else oracle.adata.obs[cluster_use].cat.categories.tolist()
    
    plt.figure(figsize=[5,6])
    plt.subplots_adjust(left=0.3, right=0.7)
    oracle.plot_mc_results_as_sankey(cluster_use=cluster_use,
                                     order = order)

    return plt.gcf()



### Comparison

def compare_perturbation(
    oracle, gradient,
    s=5, 
    scale_for_simulation=0.5,
    s_grid=50,
    scale_for_pseudotime=80, 
    vm=0.02):
    
    plots = {}
    
    ## INIT
    
    # print('# Make Oracle_development_module to compare two vector field')
    dev = Oracle_development_module()

    # print('# Load development flow')
    dev.load_differentiation_reference_data(gradient_object=gradient)

    # print('# Load simulation result')
    dev.load_perturb_simulation_data(oracle_object=oracle)


    ## COMPUTE
    # print('# Calculate inner produc scores')
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)


    ## VIZ
    # print('# Let\'s visualize the results ')
    dev.visualize_development_module_layout_0(s=s, 
                                              scale_for_simulation=scale_for_simulation,
                                              s_grid=s_grid,
                                              scale_for_pseudotime=scale_for_pseudotime, 
                                              vm=vm)
    
    # plt.savefig(f'{plots_folder}Perturbation summary.pdf')
    plots['Perturbation summary'] = plt.gcf()
    
    # print('# Show perturbation scores')
    # fig, ax = plt.subplots(figsize=[6, 6])
    # dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax)

    # print('# Show perturbation scores with perturbation simulation vector field')
    fig, ax = plt.subplots(figsize=[6, 6])
    dev.plot_inner_product_on_grid(vm=vm, s=s_grid, ax=ax)
    dev.plot_simulation_flow_on_grid(scale=scale_for_simulation, show_background=False, ax=ax)
    # plt.savefig(f'{plots_folder}Perturbation score and vector map.pdf')
    plots['Perturbation score and vector map'] = plt.gcf()

    return plots



### Expression map

import seaborn as sns

def expression_plots(oracle, gene, group_by):
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(20,12))
    sc.pl.embedding(oracle.adata, basis = oracle.embedding_name, 
                    legend_loc= 'on data', color = group_by, 
                    ax = ax[0][0], title=group_by,
                    show=False, save=False)
    sc.pl.embedding(oracle.adata, basis = oracle.embedding_name, 
                    color = gene, layer="imputed_count", cmap="viridis",
                    ax = ax[0][1], title=f'{gene} (Imputed count)',
                    show=False, save=False)
    sc.pl.embedding(oracle.adata, basis = oracle.embedding_name, 
                    color = gene, use_raw=True, vmax='p90', cmap="viridis", 
                    ax = ax[0][2], title=f'{gene} (Raw count)',
                    show=False, save=False)

    ax[1][0].axis('off')

    raw_dist = sc.get.obs_df(oracle.adata, keys=[gene, group_by], use_raw=True)
    raw_dist[gene]+=1

    sns.kdeplot(sc.get.obs_df(oracle.adata, keys=[gene, group_by], layer="imputed_count"),
                x=gene, hue=group_by, 
                legend=False,
                fill=True, ax=ax[1,1])
    sns.kdeplot(raw_dist,
                x=gene, hue=group_by, 
                # kind="kde", 
                fill=True, log_scale=[True, False], ax=ax[1,2])
    return(fig)



### Quivers

def quiver_plots(oracle, gene, exp_value, scale_for_simulation=0.5, quiver_scale=20, quiver_scale_random=10):
    # print('## CELL LEVEL PLOT')
    fig, ax = plt.subplots(2, 2,  figsize=[13, 13])

    # print('# Show quiver plot')
    oracle.plot_quiver(scale=quiver_scale, ax=ax[0][0])
    ax[0][0].set_title(f"Perturbation simulation results: {gene} {exp_value}")

    # print('# Show quiver plot that was calculated with randomized GRN.')
    oracle.plot_quiver_random(scale=quiver_scale_random, ax=ax[0][1])
    ax[0][1].set_title(f"Perturbation simulation with randomized GRNs")

    # print('# Show quiver plot')
    oracle.plot_simulation_flow_on_grid(scale=scale_for_simulation, ax=ax[1][0])
    # ax[1][0].set_title(f"Perturbation simulation results: {gene} {exp_value}")

    # print('# Show quiver plot that was calculated with randomized GRN.')
    oracle.plot_simulation_flow_random_on_grid(scale=scale_for_simulation, ax=ax[1][1])
    # ax[1][1].set_title(f"Perturbation simulation with randomized GRNs")

    return(fig)



### Vector map on cell types

def plot_simulation_on_cell_types(oracle, scale_for_simulation=0.5, s=10):
    
    # print('## SHOW ON GRID COLORED PLOT')
    # Plot vector field with cell cluster 
    fig, ax = plt.subplots(figsize=[8, 8])

    oracle.plot_cluster_whole(ax=ax, s=10)
    oracle.plot_simulation_flow_on_grid(scale=scale_for_simulation, ax=ax, show_background=False)
    return(fig)

def pertubation_plots(oracle, gradient, group_by, 
                      gene, exp_value, 
                      save_dir='./perturbation_plots',

                      s=5, 
                      scale_for_simulation=0.5,
                      s_grid=50,
                      scale_for_pseudotime=80, 
                      vm=0.02):

        
        
        
    plots = {}
    pbar=tqdm(total=5, colour='darkgreen')
    
    plots['Expression'] = expression_plots(oracle=oracle, gene=gene, group_by=group_by)
    pbar.update(1)
    plots['Perturbation per cell'] = quiver_plots(
        oracle=oracle, gene=gene, exp_value=exp_value, 
        scale_for_simulation=scale_for_simulation,
        quiver_scale=s_grid//2, quiver_scale_random=s_grid//4
    )
    pbar.update(1)
    
    plots['Vector map on cell types'] = plot_simulation_on_cell_types(
        oracle=oracle,
        scale_for_simulation=scale_for_simulation,
        s=s
    )
    pbar.update(1)
    
    plots['Cell transitions'] = plot_markov_celltype_simulation(oracle=oracle, cluster_use=group_by)
    pbar.update(1)

    plots.update(compare_perturbation(
        oracle=oracle, gradient=gradient,
        scale_for_simulation=scale_for_simulation, 
        scale_for_pseudotime=scale_for_pseudotime,
        s_grid=s_grid, s=s, vm=vm
    ))
    pbar.update(1)
    
    
    if save_dir:
        
        save_dir = os.path.join(save_dir, f'{gene}_{exp_value}')
        os.makedirs(save_dir, exist_ok=True)
        for k, p in plots.items():
            p.savefig(os.path.join(save_dir, f'{k}.pdf'))

    return(plots)
    

### Group perturbation scores

def get_cell_groups(oracle, group_dict=None,
                    group_by='final.label',
                    add_individuals=True,
                    add_wholedata=True,
                    add_antigroups=True,
                    min_cells=0):
    
    groups = {k: v for k, v in group_dict.items()} if group_dict else {}
    
    if add_individuals:
        groups_individual = {k: [k] for k in oracle.adata.obs[group_by].cat.categories}
        groups.update(groups_individual)
    
    if add_wholedata:
        groups['Whole data'] = oracle.adata.obs[group_by].cat.categories.tolist()
        
    index_dictionary = {k: np.where(oracle.adata.obs[group_by].isin(v))[0] for k, v in groups.items()}
    
    if add_antigroups:
        anti_dictionary = {f'not_{k}': np.where(~oracle.adata.obs[group_by].isin(v))[0] for k, v in groups.items()}
        index_dictionary.update(anti_dictionary)
        
    # print('Cells per group')
    # print({k: len(v) for k,v in index_dictionary.items()})

    
    if isinstance(min_cells, int):
        print('Excluding:', [k for k, v in index_dictionary.items() if len(v) < min_cells])

        index_dictionary = {k: v for k, v in index_dictionary.items() if len(v) > min_cells}
        
    return(index_dictionary)

def compute_group_ps(lineage_name, cell_idx, gradient, oracle, file_path, gene, verbose=True):
    
    # if verbose: print(lineage_name)

    dev = Oracle_development_module()
    # Load development flow
    dev.load_differentiation_reference_data(gradient_object=gradient)
    # if verbose: print('\tloaded ref', end='\t')
    # Load simulation result
    dev.load_perturb_simulation_data(oracle_object=oracle, 
                                     cell_idx_use=cell_idx, 
                                     name=lineage_name)
    # if verbose: print('loaded pert', end='\t')
    # Calculate inner product
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)
    # if verbose: print('Setting dir', end='\t')
    # Save results in a hdf5 file.
    dev.set_hdf_path(path=file_path) 
    # if verbose: print('Dumping')
    dev.dump_hdf5(gene=gene, misc=lineage_name)

    return()

def get_all_ps_p_values(file_path, negative=True):

    # Load data with Oracle_systematic_analysis_helper.
    helper = Oracle_systematic_analysis_helper(hdf5_file_path=file_path)

    pss = []
    for group in helper.hdf5_info['misc_list']:
        if negative:
            ps = helper.calculate_negative_ps_p_value(misc=group)
        if not negative:
            ps = helper.calculate_positive_ps_p_value(misc=group)
        ps['group'] = group
        ps['log1p'] = np.log1p(ps['ps_sum'])
        pss.append(ps)
        
    pss = pd.concat(pss, axis=0, ignore_index=True)
    return(pss)

### Get scales

def get_scales(file_path):

    # Load data with Oracle_systematic_analysis_helper.
    helper = Oracle_systematic_analysis_helper(hdf5_file_path=file_path)


    scales = helper.estimate_scale_for_visualization()
    return(scales)

### PS plots

def plot_ps_comparison(ps_sums):

    ps_table = ps_sums.pivot(index='gene', columns='group', values='log1p')

    fig, ax = plt.subplots(nrows = ps_table.shape[1], ncols = ps_table.shape[1],
                           figsize=[6*ps_table.shape[1]]*2, constrained_layout=True)
    
    for i, i_c in enumerate(ps_table.columns):

        for j, j_c in enumerate(ps_table.columns):

            # if i == j:
            #     ax[i][j].axis('off')
            # else:
            comp_ps = ps_table.iloc[:,[j,i]].reset_index()
            ax[i][j].scatter(x = comp_ps.iloc[:,1], y = comp_ps.iloc[:,2])
            
            # ax[i][j].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            # ax[i][j].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            
            for idx, (gene, x, y) in comp_ps.iterrows():
                ax[i][j].text(s=gene, x=x, y=y)

            if j == 0:
                ax[i][j].set_ylabel(i_c, fontsize=6*ps_table.shape[1])
            if i == ps_table.shape[1]-1:
                ax[i][j].set_xlabel(j_c, fontsize=6*ps_table.shape[1])

    return(fig,ax)

def run_systemic_simulation(
    
    oracle, gradient,
    
    group_by = None,
    genes = None,

    groups = None, 

    exp_fct = 0,
    prefix='KO',
    save_dir='./perturbations',
    transprob_n_neighbors=None,
    p_mass_smooth=0.8,
    p_mass_n_grid=40,
    p_mass_n_neighbors=200,
    p_mass_filter_min=0.001,
    
    s_grid=50,

    markov_n_steps=500,
    markov_n_duplications=5,

    overwrite=False,
    n_cores=16,
    verbose=True):
    



    ## Prepare arguments and files

    v= verbose
    # Output directories
    save_dir = save_dir or './perturbations'
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(os.path.join(save_dir,'perturbation_plots/'), exist_ok=True)
    os.makedirs(os.path.join(save_dir,'perturbation_transitions/'), exist_ok=True)


    # Store results
    file_path = f"Systematic_simulation_results.celloracle.hdf5" # Please use .hdf5 for extension.
    if prefix: file_path = prefix+'.'+file_path
    if save_dir: file_path = os.path.join(save_dir, file_path)


    # Groups and genes
    group_by = group_by or oracle.cluster_column_name
    genes = sorted([g for g in genes if g in oracle.active_regulatory_genes] if genes else oracle.active_regulatory_genes)

    plot_groups = {'Cell types': oracle.adata.obs[group_by].cat.categories.tolist() + ['Whole data']}
    if groups:
        plot_groups['Lineages & Groups'] = list(groups.keys())+['Whole data']

    # Cells per group
    cell_id_groups = get_cell_groups(oracle, groups, group_by=group_by, add_wholedata=True)

    # Adjust number of neighbors
    min_neigh = min([min([len(v) for v in cell_id_groups.values()])-1, 200])
    transprob_n_neighbors = min([transprob_n_neighbors, min_neigh]) if transprob_n_neighbors else min_neigh
    p_mass_n_neighbors = min([p_mass_n_neighbors, min_neigh]) if p_mass_n_neighbors else min_neigh
    
    if v: print(f'Selected number of neighbors: {transprob_n_neighbors}')
    

    if os.path.exists(file_path) and not overwrite :
        # Get table of genes and groups present in the perturbation file
        done_table = pd.DataFrame(
            Oracle_systematic_analysis_helper(hdf5_file_path=file_path).hdf5_info['gene_misc_lists'],
            columns=['Gene', 'Group']
        ).assign(Done=True).pivot(index='Gene', columns='Group', values= 'Done')
        done_genes = done_table.index[done_table.all(1)]
        genes = [g for g in genes if g not in done_genes]

    ## Run simulations
    genebar = tqdm(genes)

    for gene in genebar:
        
                
        # Expression value for simulation
        max_exp = max(oracle.adata[:,gene].X)[0]
        exp_value = max_exp * exp_fct

        genebar.set_description(f"Simulating gene: {gene} with expression {round(exp_value, 2)}")

        # Run perturbation
        pert_oracle = compute_perturbation(oracle, 
                                           gene, exp_value,
                                           transprob_n_neighbors=transprob_n_neighbors,
                                           p_mass_smooth=p_mass_smooth,
                                           p_mass_n_grid=p_mass_n_grid,
                                           p_mass_n_neighbors=p_mass_n_neighbors,
                                           p_mass_filter_min=p_mass_filter_min,
                                           markov_n_steps=markov_n_steps,
                                           markov_n_duplications=markov_n_duplications,

                                           n_cores=n_cores)
        # Run MCM 
        ct_transitions = pert_oracle.summarize_mc_results_by_cluster(group_by)
        # Save transitions
        ct_transitions.to_csv(os.path.join(save_dir, 'perturbation_transitions', 
                                           f'{prefix+"." if prefix else ""}{gene}_{exp_value}.csv'),
                              index=True)


        if v: print('Getting PS for all groups')
        # Do simulation for all conditions.
        for lineage_name, cell_idx in tqdm(cell_id_groups.items(), colour='green'):
            P = compute_group_ps(
                lineage_name=lineage_name, cell_idx=cell_idx, 
                gradient=gradient, oracle=oracle, 
                file_path=file_path, gene=gene, verbose=True)
            print(lineage_name, P)


        # Do plots
        if v: print('Plotting')
        scales = get_scales(file_path)
        # print(scales.loc['Whole data',:])
        plots = pertubation_plots(pert_oracle,
                                  gradient=gradient,
                                  group_by=group_by,
                                  gene = gene,
                                  exp_value = exp_value,
                                  scale_for_simulation=scales.loc['Whole data', 'simulation'],
                                  s_grid=s_grid,
                                  scale_for_pseudotime=scales.loc['Whole data', 'pseudotime'],
                                  vm=scales.loc['Whole data', 'vm'],
                                  save_dir = os.path.join(save_dir, 'perturbation_plots'))


    if v: print('Gathering PS scores and plotting')

    # Compute negative score sums
    ps_sums = get_all_ps_p_values(file_path=file_path, negative=exp_fct<1)
    ps_sums.to_csv(os.path.join(save_dir,  
                                f'{prefix+"." if prefix else ""}PerturbationScores.csv'))

    # Get plot of negative scores
    with PdfPages(os.path.join(
        save_dir, 'perturbation_plots',
        f'{prefix+"." if prefix else ""}PerturbationScoreComparison.pdf')
                 ) as pdf:

        for k, v in plot_groups.items():
            sub_pss = ps_sums[ps_sums['group'].isin(v)]
            sub_pss['group'] = pd.Categorical(sub_pss['group'], categories=v)
            fig, ax = plot_ps_comparison(sub_pss)
            fig.suptitle(k.upper(), fontsize=12*len(v))

            pdf.savefig(fig)
    plt.close()

    return(file_path, ps_sums)



from tqdm.contrib.concurrent import thread_map as process_map
import shutil



# V1 of async, not workable
# def run_systemic_simulation_async(
    
#     oracle, gradient,
    
#     group_by = None,
#     genes = None,

#     groups = None, 

#     exp_fct = 0,
#     prefix='KO',
#     save_dir='./perturbations',
#     transprob_n_neighbors=None,
#     p_mass_smooth=0.8,
#     p_mass_n_grid=40,
#     p_mass_n_neighbors=200,
#     p_mass_filter_min=0.001,
    
#     s_grid=50,

#     markov_n_steps=500,
#     markov_n_duplications=5,

#     overwrite=False,
#     n_cores=16,
#     n_parallel_genes=5,
#     verbose=True,
#     bck_dir='/scratch/bck'
# ):
    



#     ## Prepare arguments and files

#     v= verbose
#     # Output directories
#     save_dir = save_dir or './perturbations'
#     os.makedirs(save_dir, exist_ok=True)
#     os.makedirs(bck_dir, exist_ok=True)
#     os.makedirs(os.path.join(save_dir,'perturbation_plots/'), exist_ok=True)
#     os.makedirs(os.path.join(save_dir,'perturbation_transitions/'), exist_ok=True)


#     # Store results
#     file_path = f"Systematic_simulation_results.celloracle.hdf5" # Please use .hdf5 for extension.
#     if prefix: file_path = prefix+'.'+file_path
#     if save_dir: file_path = os.path.join(save_dir, file_path)


#     # Groups and genes
#     group_by = group_by or oracle.cluster_column_name
#     genes = sorted([g for g in genes if g in oracle.active_regulatory_genes] if genes else oracle.active_regulatory_genes)

#     plot_groups = {'Cell types': oracle.adata.obs[group_by].cat.categories.tolist() + ['Whole data']}
#     if groups:
#         plot_groups['Lineages & Groups'] = list(groups.keys())+['Whole data']

#     # Cells per group
#     cell_id_groups = get_cell_groups(oracle, groups, group_by=group_by, add_wholedata=True)

#     # Adjust number of neighbors
#     min_neigh = min([min([len(v) for v in cell_id_groups.values()]), 200])
#     transprob_n_neighbors = min([transprob_n_neighbors, min_neigh]) if transprob_n_neighbors else min_neigh
#     if v: print(f'Selected number of neighbors: {transprob_n_neighbors}')
    
#     if os.path.exists(file_path) and not overwrite :
#         # Get table of genes and groups present in the perturbation file
#         done_table = pd.DataFrame(
#             Oracle_systematic_analysis_helper(hdf5_file_path=file_path).hdf5_info['gene_misc_lists'],
#             columns=['Gene', 'Group']
#         ).assign(Done=True).pivot(index='Gene', columns='Group', values= 'Done')
#         done_genes = done_table.index[done_table.all(1)]
#         genes = [g for g in genes if g not in done_genes]

#     ## Run simulations
#     def parallel_wrapper(gene):
#         return((
#             gene,
#             perturb_and_transitions(
#                 oracle=oracle.copy(), gene=gene, exp_fct=exp_fct,
#                 transprob_n_neighbors=transprob_n_neighbors,
#                 p_mass_smooth=p_mass_smooth,
#                 p_mass_n_grid=p_mass_n_grid,
#                 p_mass_n_neighbors=p_mass_n_neighbors,
#                 p_mass_filter_min=p_mass_filter_min,
#                 markov_n_steps=markov_n_steps,
#                 markov_n_duplications=markov_n_duplications,
#                 group_by=group_by,
#                 prefix=prefix,
#                 save_dir=save_dir,
#                 n_cores=n_cores,
#                 verbose=False)
#         ))
    
#     if n_parallel_genes == 1:
#         pert_oracles = [parallel_wrapper(g) for g in tqdm(genes)]
#     else:
#         pert_oracles = process_map(parallel_wrapper, genes, max_workers=n_parallel_genes)


#     if v: print('Getting PS for all groups')
#     file_path_bck = os.path.join(bck_dir, os.path.basename(file_path)+'.bck.hdf5')
#     # Compute and dump negative scores
#     for gene, pert_ocle in tqdm(pert_oracles):
        
#         # this should not be parallelized
#         # Do simulation for all conditions.
#         for lineage_name, cell_idx in tqdm(cell_id_groups.items(), colour='green'):

#             try:

#                 compute_group_ps(
#                     lineage_name=lineage_name, cell_idx=cell_idx, 
#                     gradient=gradient, oracle=pert_ocle, 
#                     file_path=file_path_bck, gene=gene, verbose=False)
                
#                 shutil.copyfile(file_path_bck, file_path)

#             except:
#                 print(f'Failed, retrieving from {file_path}')
#                 shutil.copyfile(file_path, file_path_bck)
                
#                 compute_group_ps(
#                     lineage_name=lineage_name, cell_idx=cell_idx, 
#                     gradient=gradient, oracle=pert_ocle, 
#                     file_path=file_path_bck, gene=gene, verbose=False)
                
#                 shutil.copyfile(file_path_bck, file_path)

#     if v: print('Gathering PS scores and plotting')

#     # Compute negative score sums
#     ps_sums = get_all_ps_p_values(file_path=file_path, negative=exp_fct<1)
#     ps_sums.to_csv(os.path.join(save_dir,  
#                                 f'{prefix+"." if prefix else ""}PerturbationScores.csv'))
    
#     # Get plot of negative scores
#     with PdfPages(os.path.join(
#         save_dir, 'perturbation_plots',
#         f'{prefix+"." if prefix else ""}PerturbationScoreComparison.pdf')
#                  ) as pdf:

#         for k, v in plot_groups.items():
#             sub_pss = ps_sums[ps_sums['group'].isin(v)]
#             sub_pss['group'] = pd.Categorical(sub_pss['group'], categories=v)
#             fig, ax = plot_ps_comparison(sub_pss)
#             fig.suptitle(k.upper(), fontsize=12*len(v))

#             pdf.savefig(fig)
#     plt.close()

#     # Do plots
#     if v: print('Plotting')
#     scales = get_scales(file_path)
#     def parallel_plot_wrapper(pert_res):
        
#         gene, pert_ocle = pert_res
        
#         # Expression value for simulation
#         max_exp = max(oracle.adata[:,gene].X)[0]
#         exp_value = max_exp * exp_fct

#         return(pertubation_plots(
#                pert_ocle,
#                gradient=gradient,
#                group_by=group_by,
#                gene = gene,
#                exp_value = exp_value,
#                scale_for_simulation=scales.loc['Whole data', 'simulation'],
#                s_grid=s_grid,
#                scale_for_pseudotime=scales.loc['Whole data', 'pseudotime'],
#                vm=scales.loc['Whole data', 'vm'],
#                save_dir = os.path.join(save_dir, 'perturbation_plots')
#         ))
                                
    
#     pert_plots = process_map(parallel_plot_wrapper, pert_oracles, max_workers=n_cores)

#     return(file_path, ps_sums)

# def perturb_and_transitions(
#     oracle, gene, exp_fct,
#     transprob_n_neighbors,
#     p_mass_smooth,
#     p_mass_n_grid,
#     p_mass_n_neighbors,
#     p_mass_filter_min,
#     markov_n_steps,
#     markov_n_duplications,
#     group_by,
#     prefix,
#     save_dir,
#     n_cores, 
#     verbose):         
    
    
#     # Expression value for simulation
#     max_exp = max(oracle.adata[:,gene].X)[0]
#     exp_value = max_exp * exp_fct


#     # Run perturbation
#     pert_oracle = compute_perturbation(oracle, 
#                                        gene, exp_value,
#                                        transprob_n_neighbors=transprob_n_neighbors,
#                                        p_mass_smooth=p_mass_smooth,
#                                        p_mass_n_grid=p_mass_n_grid,
#                                        p_mass_n_neighbors=p_mass_n_neighbors,
#                                        p_mass_filter_min=p_mass_filter_min,
#                                        markov_n_steps=markov_n_steps,
#                                        markov_n_duplications=markov_n_duplications,

#                                        n_cores=n_cores,
#                                        verbose=verbose)
#     # Run MCM 
#     ct_transitions = pert_oracle.summarize_mc_results_by_cluster(group_by)
#     # Save transitions
#     ct_transitions.to_csv(os.path.join(save_dir, 'perturbation_transitions', 
#                                        f'{prefix+"." if prefix else ""}{gene}_{exp_value}.csv'),
#                           index=True)

#     return(pert_oracle)



# V2 

from celloracle.applications import Oracle_systematic_analysis_helper

from scipy.stats import wilcoxon
import numpy as np
from celloracle.applications import Oracle_development_module
import pandas as pd

from tqdm.contrib.concurrent import thread_map as process_map
import shutil



def run_systemic_simulation_async(
    
    oracle, gradient,
    
    group_by = None,
    genes = None,

    groups = None, 

    exp_fct = 0,
    prefix='KO',
    save_dir='./perturbations',
    transprob_n_neighbors=None,
    p_mass_smooth=0.8,
    p_mass_n_grid=40,
    p_mass_n_neighbors=200,
    p_mass_filter_min=0.001,
    
    s_grid=50,

    markov_n_steps=500,
    markov_n_duplications=5,

    overwrite=False,
    n_cores=16,
    n_parallel_genes=5,
    verbose=True
):
    



    ## Prepare arguments and files

    v= verbose
    # Output directories
    save_dir = save_dir or './perturbations'
    os.makedirs(save_dir, exist_ok=True)
    # os.makedirs(bck_dir, exist_ok=True)
    os.makedirs(os.path.join(save_dir,'perturbation_plots/'), exist_ok=True)
    os.makedirs(os.path.join(save_dir,'perturbation_transitions/'), exist_ok=True)
    os.makedirs(os.path.join(save_dir,'perturbation_hdf5/'), exist_ok=True)


    # Groups and genes
    group_by = group_by or oracle.cluster_column_name
    
    if genes:
        na_genes = [g for g in genes if g not in oracle.active_regulatory_genes]
        if na_genes:
            print(f'Not available genes: {sorted(na_genes)}')
            
        genes = sorted([g for g in genes if g in oracle.active_regulatory_genes])
    else:
        genes = oracle.active_regulatory_genes

    plot_groups = {'Cell types': oracle.adata.obs[group_by].cat.categories.tolist() + ['Whole data']}
    if groups:
        plot_groups['Lineages & Groups'] = list(groups.keys())+['Whole data']

    # Cells per group
    cell_id_groups = get_cell_groups(oracle, groups, group_by=group_by, add_wholedata=True)
    cell_id_groups = {k: v for k,v in cell_id_groups.items() if len(v)>0}

    # Adjust number of neighbors
    min_neigh = min([min([len(v) for v in cell_id_groups.values()]), 200])
    transprob_n_neighbors = min([transprob_n_neighbors, min_neigh]) if transprob_n_neighbors else min_neigh
    if v: print(f'Selected number of neighbors: {transprob_n_neighbors}')
    
    # if os.path.exists(file_path) and not overwrite :
    #     # Get table of genes and groups present in the perturbation file
    #     done_table = pd.DataFrame(
    #         Oracle_systematic_analysis_helper(hdf5_file_path=file_path).hdf5_info['gene_misc_lists'],
    #         columns=['Gene', 'Group']
    #     ).assign(Done=True).pivot(index='Gene', columns='Group', values= 'Done')
    #     done_genes = done_table.index[done_table.all(1)]
    #     genes = [g for g in genes if g not in done_genes]

    ## Run simulations
    def parallel_wrapper(gene):

        # 1. Perturb gene -> save transitions
        # Expression value for simulation
        max_exp = max(oracle.adata[:,gene].X)[0]
        exp_value = max_exp * exp_fct

        # Store results
        file_path = f"perturbation_hdf5/{prefix+'.' if prefix else ''}{gene}.celloracle.hdf5" # Please use .hdf5 for extension.
        print(file_path)
        if save_dir: file_path = os.path.join(save_dir, file_path)

        
        # Run perturbation
        pert_oracle = compute_perturbation(oracle.copy(), 
                                           gene, exp_value,
                                           transprob_n_neighbors=transprob_n_neighbors,
                                           p_mass_smooth=p_mass_smooth,
                                           p_mass_n_grid=p_mass_n_grid,
                                           p_mass_n_neighbors=p_mass_n_neighbors,
                                           p_mass_filter_min=p_mass_filter_min,
                                           markov_n_steps=markov_n_steps,
                                           markov_n_duplications=markov_n_duplications,

                                           n_cores=n_cores,
                                           verbose=verbose)
        
        try:
            # Run MCM 
            ct_transitions = pert_oracle.summarize_mc_results_by_cluster(group_by)
            # Save transitions
            ct_transitions.to_csv(os.path.join(save_dir, 'perturbation_transitions', 
                                               f'{prefix+"." if prefix else ""}{gene}_{exp_value}.csv'),
                                  index=True)
        except:
            print(f'Could not compute transitions for gene {gene}')
            ct_transitions = None
        
        
        # 2. compute ps scores -> save scores
        try:       
            pss = compute_ps(pert_oracle, gradient=gradient, gene=gene, cell_groups=cell_id_groups, file_path=file_path)
            pss.to_csv(os.path.join(save_dir, 'perturbation_hdf5', 
                                    f'{prefix+"." if prefix else ""}NegPSsum.{gene}_{exp_value}.csv'),
                       index=True)
        except:
            print(f'Could not compute PSs for gene {gene}')
            pss = None
                    
        # 3. get scales
        try:
            scales = get_scales(file_path)
            # scales = perturbations.get_scales(file_path)
        
            pertubation_plots(
                           pert_oracle,
                           gradient=gradient,
                           group_by=group_by,
                           gene = gene,
                           exp_value = exp_value,
                           scale_for_simulation=scales.loc['Whole data', 'simulation'],
                           s_grid=s_grid,
                           scale_for_pseudotime=scales.loc['Whole data', 'pseudotime'],
                           vm=scales.loc['Whole data', 'vm'],
                           save_dir = os.path.join(save_dir, 'perturbation_plots')
                    )        
        except:
            print(f'Could not compute scales for gene {gene}')
            scales = None
            try:
                pertubation_plots(
                   pert_oracle,
                   gradient=gradient,
                   group_by=group_by,
                   gene = gene,
                   exp_value = exp_value,
                   # scale_for_simulation=scales.loc['Whole data', 'simulation'],
                   s_grid=s_grid,
                   # scale_for_pseudotime=scales.loc['Whole data', 'pseudotime'],
                   # vm=scales.loc['Whole data', 'vm'],
                   save_dir = os.path.join(save_dir, 'perturbation_plots')
                )   
            except:
                print(f'Could not plot gene {gene}')
        return(pss)      
        
    
    if n_parallel_genes == 1:
        pss = [parallel_wrapper(g) for g in tqdm(genes)]
    else:
        pss = process_map(parallel_wrapper, genes, max_workers=n_parallel_genes)

        
    # 5. plot all scores        
    if v: print('Getting PS for all groups')
    ps_sums = pd.concat(pss, axis=0, ignore_index=True)
    ps_sums.to_csv(os.path.join(save_dir,  
                                f'{prefix+"." if prefix else ""}PerturbationScores.csv'))
    
    ps_sums['log1p'] = np.log1p(ps_sums['score'])
    
    # Get plot of negative scores
    with PdfPages(os.path.join(
        save_dir, 'perturbation_plots',
        f'{prefix+"." if prefix else ""}PerturbationScoreComparison.pdf')
                 ) as pdf:

        for k, v in plot_groups.items():
            sub_pss = ps_sums[ps_sums['group'].isin(v)]
            sub_pss['group'] = pd.Categorical(sub_pss['group'], categories=v)
            fig, ax = plot_ps_comparison(sub_pss)
            fig.suptitle(k.upper(), fontsize=12*len(v))

            pdf.savefig(fig)
    plt.close()

    return(ps_sums)


def perturb_and_transitions(
    oracle, gene, exp_fct,
    transprob_n_neighbors,
    p_mass_smooth,
    p_mass_n_grid,
    p_mass_n_neighbors,
    p_mass_filter_min,
    markov_n_steps,
    markov_n_duplications,
    group_by,
    prefix,
    save_dir,
    n_cores, 
    verbose):         
    
    
    # Expression value for simulation
    max_exp = max(oracle.adata[:,gene].X)[0]
    exp_value = max_exp * exp_fct


    # Run perturbation
    pert_oracle = compute_perturbation(oracle, 
                                       gene, exp_value,
                                       transprob_n_neighbors=transprob_n_neighbors,
                                       p_mass_smooth=p_mass_smooth,
                                       p_mass_n_grid=p_mass_n_grid,
                                       p_mass_n_neighbors=p_mass_n_neighbors,
                                       p_mass_filter_min=p_mass_filter_min,
                                       markov_n_steps=markov_n_steps,
                                       markov_n_duplications=markov_n_duplications,

                                       n_cores=n_cores,
                                       verbose=verbose)
    # Run MCM 
    ct_transitions = pert_oracle.summarize_mc_results_by_cluster(group_by)
    # Save transitions
    ct_transitions.to_csv(os.path.join(save_dir, 'perturbation_transitions', 
                                       f'{prefix+"." if prefix else ""}{gene}_{exp_value}.csv'),
                          index=True)
    
    

    return(pert_oracle)


def compute_ps(oracle, gradient, gene, cell_groups, file_path=None): 
    
    print(f'File path: {file_path}')
    pss = {}
    
    for lineage_name, cell_idx in cell_groups.items():
    

        try:
            dev = Oracle_development_module()
            # Load development flow
            dev.load_differentiation_reference_data(gradient_object=gradient)
            # if verbose: print('\tloaded ref', end='\t')
            # Load simulation result
            dev.load_perturb_simulation_data(oracle_object=oracle, 
                                             cell_idx_use=cell_idx, 
                                             name=lineage_name)
            # if verbose: print('loaded pert', end='\t')
            # Calculate inner product
            dev.calculate_inner_product()
            dev.calculate_digitized_ip(n_bins=10)


            df = dev.inner_product_df
            x = df["score"]
            # Clipping positive value to focus on negative IP
            x = np.clip(x, -np.inf, 0)

            pss[lineage_name] = -x.sum()
            
            if file_path:
                dev.set_hdf_path(path=file_path) 
                print(f'Dumping to {file_path}')
                dev.dump_hdf5(gene=gene, misc=lineage_name)


        except:
            pss[lineage_name] = np.nan
    
    
    return pd.DataFrame(pss, index=['score']).T.assign(gene=gene).reset_index().rename({'index':'group'}, axis=1)
