import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

def check_bedtools(
    bedtools_path='/soft/system/software/BEDTools/2.30.0-GCC-10.2.0/bin'
):
    
    import shutil
    shutil.which('bedtools')

    bed_path = shutil.which('bedtools')

    if bed_path:
        print(f'Detected BEDTools in: {bed_path}')
        
    elif bedtools_path:
        print(f'Adding to PATH: {bedtools_path}')
        os.environ.update({
            'PATH': bedtools_path+':'+os.environ.get('PATH')
           })
    else:
        raise
        
    import shutil
    shutil.which('bedtools')
    
    return(bed_path)

check_bedtools()

from celloracle import motif_analysis as ma

#######################################################################
#### global functions
#######################################################################

# Define quality check fuction
def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'
        
    Returns:
        tuple: chromosome name, start position, end position
    """
    
    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end


def check_peak_format(peaks_df, ref_genome):
    """
    Check peak format. 
     (1) Check chromosome name. 
     (2) Check peak size (length) and remove sort DNA sequences (<5bp)
    
    """
    from genomepy import Genome

    df = peaks_df.copy()
    
    n_peaks_before = df.shape[0]
    
    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed))
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(np.int)
    df_decomposed["end"] = df_decomposed["end"].astype(np.int)
    
    # Load genome data
    genome_data = Genome(ref_genome)
    all_chr_list = list(genome_data.keys())
    
    
    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])
    
    
    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]
    
    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])
    
    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]
    
    #
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)
    
    return df

#######################################################################

#######################################################################
#### functions to incorporate into class MotifAnalysis
#######################################################################
def check_peak_sep(peaks):
    '''
    Check if peak separator is '_' in all peaks.

    Args:
        peaks (list): list of peak strings

    Returns:
        bool: True if peak separator is '_' in all peaks.
    '''
    return(all([p.split('_').__len__() >= 3 for p in peaks]))




def load_data(peaks_csv, connections_csv):
    '''
    Load data from peaks and connections csv files.
    Checks if peak separator is '_' in both files.

    Args:
        peaks_csv (str): path to peaks csv file
        connections_csv (str): path to connections csv file

    Returns:
        tuple: peaks, connections
    '''
    peaks = pd.read_csv(peaks_csv, index_col=0).x.values
    connections = pd.read_csv(connections_csv, index_col=0)
    
    if not check_peak_sep(peaks):
        seps = set((p for peak in peaks for p in peak if not p.isalnum()))
        for sep in seps:
            subs_peaks = np.array([p.replace(sep, '_') for p in peaks])
            if check_peak_sep(subs_peaks):
                peaks = subs_peaks
                connections.iloc[:,0] = connections.iloc[:,0].str.replace(sep, '_')
                connections.iloc[:,1] = connections.iloc[:,1].str.replace(sep, '_')
                break
        
    return(peaks, connections)




def check_genome(ref_genome='hg38', provider='UCSC', verbose=True):
    '''
    Check if the genome is available.

    Args:
        ref_genome (str): genome name
        provider (str): genome provider
        verbose (bool): print out the status

    Returns:
        bool: True if genome is available.
    '''

    if not ma.is_genome_installed(ref_genome=ref_genome):
        if verbose:
            print(f"Genome {ref_genome} is not installed. Installing.")
        import genomepy
        genomepy.install_genome(ref_genome, provider)
    elif verbose:
        print(ref_genome, "is installed.")
    return(True)




def tss_annotation(peaks, ref_genome, cicero_connections, extra_peaks=None):
    '''
    Annotate cicero connections with gene names based on TSS annotation.

    Args:
        peaks (list): list of peak strings
        ref_genome (str): genome name
        cicero_connections (pandas.DataFrame): cicero connections
        extra_peaks (dict): per gene lists of extra peaks to annotate genes to
                            {'GENE': ['chr_start_end']}
    Returns:
        pandas.DataFrame: cicero connections with gene names
    '''
    ##!! Please make sure to specify the correct reference genome here
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=ref_genome) 
    
    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                                   cicero_connections=cicero_connections)
    return(integrated)




def filter_coaccessibility(integrated, coacc_thres, verbose=True, plot=True):
    '''
    Filter connections by coaccessibility.
    
    Args:
        integrated (pandas.DataFrame): integrated cicero connections
        coacc_thres (float): coaccessibility threshold
        verbose (bool): print out the status
        plot (bool): plot the histogram of coaccessibility values and threshold

    Returns:
        pandas.DataFrame: filtered cicero connections
    '''

    if verbose:
        print('Over threshold:', 
              round((integrated.coaccess >= coacc_thres).mean(), 4)*100, 
              '%')
    if plot:
        integrated.coaccess.hist(bins=40)
        plt.axvline(x=coacc_thres, color='red')

    return(integrated[integrated.coaccess >= coacc_thres])
    



def peaks_from_coaccs(connections):
    '''
    Get peak table from coaccessibility-filtered peaks.

    Args:
        connections (pandas.DataFrame): coaccessibility-filtered cicero connections

    Returns:
        pandas.DataFrame: peak table
    '''

    return(connections[["peak_id", "gene_short_name"]].reset_index(drop=True))




def motif_scan(peak_df, ref_genome='hg38', frp=0.02, motifs=None, verbose=True):
    '''
    Scan motifs in peak regions.

    Args:
        peak_df (pandas.DataFrame): peak table
        ref_genome (str): genome name
        frp (float): motif scan frp
        motifs (str): motif database
        verbose (bool): print out the status

    Returns:
        celloracle.motif_analysis.TFInfo: TFInfo object
        '''

    tfi = ma.TFinfo(peak_data_frame=peak_df, 
                    ref_genome=ref_genome) 

    # Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
    tfi.scan(fpr=frp, 
             motifs=motifs,  # If you enter None, default motifs will be loaded.
             verbose=verbose)
    return(tfi)




def motif_filtering(tfi, motif_score_thres, verbose=True, plot=True):
    '''
    Filter motifs by score.
    
    Args:
        tfi: (celloracle.motif_analysis.TFInfo) object
        motif_score_thres (float): motif score threshold
        verbose (bool): print out the status
        plot (bool): plot the histogram of motif scores and threshold

    Returns:
        celloracle.motif_analysis.TFInfo: TFInfo object
    '''

    if verbose:
        print('Over threshold:', 
              round((tfi.scanned_df.score > motif_score_thres).mean(), 4)*100, 
              '%')
    if plot:
        tfi.scanned_df.score.hist(bins=40)
        plt.axvline(x=motif_score_thres, color='red')
        plt.show()
    # Reset filtering 
    tfi.reset_filtering()

    # Do filtering
    tfi.filter_motifs_by_score(threshold=motif_score_thres)

    # Format post-filtering results.
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=verbose)
    
    return(tfi)




#######################################################################

#######################################################################
# class MotifAnatqdm.notebook######################################################################
def run_motif_analysis(
    peaks_csv, connections_csv, extra_peaks=None,
    
    ref_genome='hg38', provider='UCSC',
    
    frp=0.02, coacc_thres=0.8, 
    motifs=None, motif_score_thres = 8,
    
    save_dir=None, prefix=None,
    
    verbose=True):
    
    '''
    Motif analysis runner.

    Args:
        peaks_csv (str): peaks csv file
        connections_csv (str): connections csv file
        extra_peaks (dict): per gene lists of extra peaks to annotate genes to
                            {'GENE': ['chr_start_end', ...]}
        ref_genome (str): genome name
        provider (str): genome provider
        frp (float): motif scan frp
        coacc_thres (float): coaccessibility threshold
        motifs (str): motif database
        motif_score_thres (float): motif score threshold
        verbose (bool): print out the status
        save_dir (str): dir for saving results
        prefix (str): basename for saving results

    Returns:
        Dictionary of peaks filtered by coaccessibility (peaks) and score (base_grn) and celloracle.tfi object
    '''

    results = {}
    v=verbose

    if v: print('[1] Loading peaks and connections')
    peaks, connections = load_data(peaks_csv, connections_csv)

    if v: print('[2] Loading genome')
    check_genome(ref_genome, provider)

    if v: print('[3] Integrating peaks, connections and genome')
    integrated = tss_annotation(peaks=peaks,
                                ref_genome=ref_genome, 
                                cicero_connections=connections)

    if v: print(f'[4] Filtering peaks by coaccessibility ({coacc_thres})')
    filtered_integrated = filter_coaccessibility(integrated, coacc_thres, verbose=verbose)

    if v: print('[5] Extracting peaks')
    peak_df = peaks_from_coaccs(filtered_integrated)

    if extra_peaks:

        extra_peaks_df = pd.DataFrame(
            [{'peak_id':p, 'gene_short_name': g} for g, v in extra_peaks.items() for p in v]
        )[peak_df.columns]
        full_peak_df = pd.concat([peak_df, extra_peaks_df], ignore_index=True)

        if v:
            print(f'[EXTRA] Adding {extra_peaks_df.shape[0]} extra peak to TSS annotations')

        if full_peak_df.duplicated().any():
            if v:
                print(f'[EXTRA] Duplicated peaks (only kept once):\n{full_peak_df[full_peak_df.duplicated()].to_string()}')
            peak_df = full_peak_df[~full_peak_df.duplicated()]

        else: 
            peak_df = full_peak_df

    results[f'coaccessible_peaks'] = peak_df.copy()

    if v: print(f'[6] Performing motif scan (frp: {frp})')
    tfi = motif_scan(peak_df, ref_genome=ref_genome, frp=frp, motifs=motifs, verbose=verbose)

    if v: print(f'[7] Filtering motif score (motif_score_thres)')
    tfi = motif_filtering(tfi, motif_score_thres, verbose=verbose)

    results[f'base_grn'] = tfi.to_dataframe()
    results[f'celloracle_tfinfo'] = tfi

    if save_dir:

        os.makedirs(save_dir, exist_ok=True)
        
        if v: print('[8] Saving results')
        for k, v in tqdm(results.items()):

            fname = k
            if prefix: fname = prefix+'.'+fname

            fname = os.path.join(save_dir, fname)

            if k == 'coaccessible_peaks':
                v.to_csv(f'{fname}.coaccessibility={coacc_thres}.csv')
            elif k == 'base_grn':
                v.to_parquet(f'{fname}.coaccessibility={coacc_thres}.score={motif_score_thres}.base_grn.parquet')
            elif k == 'celloracle_tfinfo':
                v.to_hdf5(f'{fname}.coaccessibility={coacc_thres}.score={motif_score_thres}.celloracle.tfinfo')





    return(results)

