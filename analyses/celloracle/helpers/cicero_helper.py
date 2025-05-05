def create_cicero_cmd( 
    atac_path,
    peaks_bed_path,
    cell_metadata_path, 
    plot_by=None, 
    save_dir=None, 
    prefix=None, 
    verbose=True,
    R_quote="\"",
    str_quote="'",
    cicero_src='/users/genomics/xoel/canonades/bioinforgalician/src/python/celloracle/cicero.R'):
    

    if not plot_by: plot_by='NULL'
    else: plot_by = str_quote+plot_by+str_quote
    
    if not save_dir: save_dir='NULL'
    else: save_dir = str_quote+save_dir+str_quote
    
    if not prefix: prefix='NULL'
    else: prefix = str_quote+prefix+str_quote
    
    verbose = 'TRUE' if verbose else 'FALSE'
    
    r_cmd = f'''\
source({str_quote}{cicero_src}{str_quote},echo={verbose})

pipe_cicero_conns(atac_path={str_quote}{atac_path}{str_quote},peaks_bed_path={str_quote}{peaks_bed_path}{str_quote},cell_metadata_path={str_quote}{cell_metadata_path}{str_quote},plot_by={plot_by},save_dir={save_dir},prefix={prefix},verbose={verbose})'''
    
    return(r_cmd)