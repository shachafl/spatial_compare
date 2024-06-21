import pandas as pd
from scanpy import read_h5ad
from scipy.spatial import cKDTree
import numpy as np
import matplotlib.pyplot as plt 
from pathlib import Path
import plotly.graph_objs as go


def transcripts_per_cell(cell_x_gene):
    return cell_x_gene.sum(axis=1).to_numpy()

    
def create_seg_comp_df(barcode, seg_name, base_path, min_transcripts):
    """
    Gathers related segmentation results to create descriptive dataframe 
    INPUTS
        barcode: unique identifier, for locating specific segmentation path. String.
        seg name: descriptor of segmentation, i.e. algo name. String.
        base path: path to segmentation results
        min transcripts: minimum number of transcripts needed to define a cell too low quality to be considered for mapping
    OUTPUTS
        seg_df: dataframe with segmentation cell xy locations, low quality cell identifier, and source segmentation ID. index is cell IDs + segmentation identifier string
    """

    seg_path = Path(base_path).joinpath(barcode)
    cell_check = [x for x in seg_path.glob('*.csv') if 'cellpose' in x.stem]
    if cell_check: 
        cxg = pd.read_table(str(seg_path)+'/cellpose-cell-by-gene.csv', index_col=0, sep=',')
        metadata = pd.read_table(str(seg_path)+'/cellpose_metadata.csv', index_col=0, sep=',')
    else:
        cxg = pd.read_table(str(seg_path)+'/cell_by_gene.csv', index_col=0, sep=',')
        meta = read_h5ad(str(seg_path)+'/metadata.h5ad')
        metadata = meta.obs
    
    #assemble seg df with filt cells col, xy cols, and seg name
    high_quality_cells = cxg.index[np.where((transcripts_per_cell(cxg) >= min_transcripts))[0]].astype(str).values.tolist()
    seg_df = metadata[['center_x', 'center_y']].copy()
    seg_df.index = seg_df.index.astype(str)
    seg_df.loc[:, 'source'] = seg_name
    seg_df.loc[:, 'low_quality_cells'] = np.where(seg_df.index.isin(high_quality_cells), False, True)
    seg_df.index= seg_name+'_'+ seg_df.index
    return seg_df


def get_segmentation_data(barcode, seg_name_a, seg_name_b, seg_a_path, seg_b_path, save, reuse_saved=True, min_transcripts=40):
    """
    Loads or calls function to collect segmentation data
    Inputs:
        Barcode: unique section identifier, string
        seg_name_a, seg_name_b: names of segmentations to be compared, string
        seg_a_path, seg_b_path: path to segmentation results, string
        save: whether to save the intermediary dataframe (useful if running functions individually), bool
        reuse_saved: load in the previously saved intermediary dataframe (useful if not making changes to segmentation metadata), bool, default True
        min_transcripts: minimum number of transcripts needed to define a cell too low quality to be considered for mapping
    RETURNS
        seg_comp_df: dataframe describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, with unique index
    """
    
    filepath = seg_b_path+barcode+'_seg_comp_df_'+seg_name_a+'_and_'+seg_name_b+'.csv'
    if Path(filepath).exists() and reuse_saved: 
        seg_comp_df = pd.read_csv(filepath, index_col=0)
    else:
        seg_a_df = create_seg_comp_df(barcode, seg_name_a, seg_a_path, min_transcripts)
        seg_b_df = create_seg_comp_df(barcode, seg_name_b, seg_b_path, min_transcripts)
        seg_comp_df = pd.concat([seg_a_df, seg_b_df], axis=0)
    if save:
        seg_comp_df.to_csv(filepath)
    return seg_comp_df


def get_mutual_matches(dfa, dfb, nn_dist): 
    """
    Compare x and y locations using nearest neighbor function of each segmentation result to identify same cells across segmentations
    INPUTS
        dfa, dfb: subset dataframe for each segmentation
        nn_dist: distance cutoff for identifying likely same cell between segmentations. float.
    OUTPUT:
        mutual_match: Dictionary of cell IDs to matched cell ids. Keys are stacked source a and source b cell ids, with values being match from other source (so source a cell id keys have source b cell id values, and vice versa)
    """
    
    #get single nearest neighbor for each
    tree_a = cKDTree(dfa[["center_x", "center_y"]].copy())
    dists_a, inds_a = tree_a.query(dfb[['center_x', 'center_y']].copy(), 1) #gives me index locs for dfa with neighbor dfb
    tree_b = cKDTree(dfb[["center_x", "center_y"]].copy())
    dists_b, inds_b = tree_b.query(dfa[['center_x', 'center_y']].copy(), 1) #vice versa

    #get mutually matching pairs and save (index of seg_comp_df is str)
    match_to_dfb = pd.DataFrame(data =dfa.iloc[inds_a].index.values.tolist(), 
                            index=pd.Index(dfb.index, name='match_dfb_index'), 
                            columns = ['match_dfa_index'])
    match_to_dfb['same_cell'] = np.where(dists_a <= nn_dist, True, False)
    match_to_dfa = pd.DataFrame(data = dfa.index,
                                index=pd.Index(dfb.iloc[inds_b].index.values.tolist(), name='match_dfb_index'),
                                columns = ['match_dfa_index'])
    match_to_dfa['same_cell'] = np.where(dists_b <= nn_dist, True, False)
    mutual_matches = pd.merge(match_to_dfa, match_to_dfb, how='outer')
    mutual_matches = mutual_matches.set_index('match_dfa_index').join(match_to_dfa.reset_index().set_index('match_dfa_index'), how='left', rsuffix='_match')
    mutual_match_dict = mutual_matches[mutual_matches['same_cell']==True].drop(['same_cell', 'same_cell_match'], axis=1).to_dict()
    inv_mutual_match_dict = {v:k for k, v in mutual_match_dict['match_dfb_index'].items()}
    # added union operwtor to dictionaries in 2020 for 3.9+ (pep 584), discussed https://stackoverflow.com/questions/38987/how-do-i-merge-two-dictionaries-in-a-single-expression-in-python
    mutual_matches_stacked = mutual_match_dict['match_dfb_index'] | inv_mutual_match_dict
    return mutual_matches_stacked


def collect_mutual_match_and_doublets(barcode, seg_name_a ='VPT', seg_name_b ='SIS', save = True, nn_dist = 2.5, reuse_saved= True,
                                   seg_a_path = '//allen/programs/celltypes/workgroups/hct/emilyg/manuscript_23_completed_files/resegd_nas05_files/',
                                   seg_b_path='/allen/programs/celltypes/workgroups/hct/emilyg/reseg_project/new_seg/'):
    """
    Runs all relevant functions in required order, generating dataframe summarizing comparisons between two segmentations for a single section.
    INPUTS
        Barcode: unique identifier for section segmentation comparison was run on. String.
        seg_name_a, seg_name_b: names of segmentations to be compared, string
        save: whether to save the intermediate and final results. Bool, default true
        nn_dist: distance cutoff for identifying likely same cell between segmentations. float, default 2.5 (um)
        seg_a_path, seg_b_path: path to segmentation results, string
    OUTPUTS:
        seg_comp_df: dataframe with unique index describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, mutual matches, and putative doublets. 
    """
    
    #grab base comparison df
    seg_comp_df = get_segmentation_data(barcode, seg_name_a, seg_name_b, seg_a_path, seg_b_path, save=save, reuse_saved=reuse_saved) 
    #mutual matches with both segmentations filtered
    dfa = seg_comp_df[(seg_comp_df['source']==seg_name_a)&(seg_comp_df['low_quality_cells']==False)]
    dfb = seg_comp_df[(seg_comp_df['source']==seg_name_b)&(seg_comp_df['low_quality_cells']==False)]
    col_name_both_filt = 'match_'+seg_name_a+'_filt_'+seg_name_b+'_filt'
    seg_comp_df[col_name_both_filt] = seg_comp_df.index.map(get_mutual_matches(dfa, dfb, nn_dist)) #filt both
    
    #mutual matches with seg a filtered only
    dfa = seg_comp_df[(seg_comp_df['source']==seg_name_a)&(seg_comp_df['low_quality_cells']==False)]
    dfb = seg_comp_df[seg_comp_df['source']==seg_name_b]
    col_name_afilt = 'match_'+seg_name_a+'_filt_'+seg_name_b+'_unfilt'
    seg_comp_df[col_name_afilt] = seg_comp_df.index.map(get_mutual_matches(dfa, dfb, nn_dist)) #a filt only

    #mutual matches with seg b filtered only
    dfa = seg_comp_df[seg_comp_df['source']==seg_name_a]
    dfb = seg_comp_df[(seg_comp_df['source']==seg_name_b)&(seg_comp_df['low_quality_cells']==False)]
    col_name_bfilt = 'match_'+seg_name_a+'_unfilt_'+seg_name_b+'_filt'
    seg_comp_df[col_name_bfilt] = seg_comp_df.index.map(get_mutual_matches(dfa, dfb, nn_dist)) #b filt only

    #find doublets
    seg_comp_df['nota'] = seg_comp_df['match_VPT_filt_SIS_filt'].equals(seg_comp_df['match_VPT_filt_SIS_unfilt'])
    seg_comp_df['notb']= seg_comp_df['match_VPT_filt_SIS_filt'].equals(seg_comp_df['match_VPT_unfilt_SIS_filt'])
    seg_comp_df['putative doublets'] = ~(seg_comp_df['nota'] | seg_comp_df['notb'])
    seg_comp_df.drop(['nota', 'notb'], axis=1, inplace=True)

    #save results
    if save:
        seg_comp_df.to_csv(seg_b_path+barcode+'_seg_comp_df_'+seg_name_a+'_and_'+seg_name_b+'_populated.csv', index= True)
        print('Saved to: '+seg_b_path+barcode+'_seg_comp_df_'+seg_name_a+'_and_'+seg_name_b+'_populated.csv')
    return seg_comp_df


def create_node_df_sankey(seg_comp_df, barcode, save=True, savepath = '/allen/programs/celltypes/workgroups/hct/emilyg/reseg_project/new_seg/'):
    """
    Dataframe describing each node to be included in sankey diagram. Additionally contains data helpful for generating required links dataframe. 
    INPUTS
        seg_comp_df: dataframe with unique index describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, mutual matches, and putative doublets. 
        barcode: unique identifier for section, used for creating the save name for the dataframe, string.
        save: whether to save the final results. Bool, default true
        savepath: path to which results should be saved, string
    OUTPUTS
        nodes_df: dataframe containing node label, color, and value for sankey diagram creation, with additional level and source keys for ease of link dataframe creation
    """
    
    color_dict = {seg_comp_df.source.unique().tolist()[0]: 'red', seg_comp_df.source.unique().tolist()[1]:'blue'}
    nodes_df = pd.DataFrame(columns = ['Label','Color', 'Level', 'Source', 'Value']) 
    
    for source, g in seg_comp_df.groupby('source'):
        new_rows = []
        low_q_and_match = g.groupby('low_quality_cells').agg('count') # gives me low q t/f counts, and matched. 
        total_cells = low_q_and_match.iloc[:,0].sum()
        nodes_df.loc[len(nodes_df)] = {'Label': 'total <br>'+str(total_cells), 'Color': color_dict[source], 
                                      'Level': 0, 'Source':source, 'Value': total_cells}
        #low and normal quality cells
        nodes_df.loc[len(nodes_df)] = {'Label': 'low quality cells <br>'+str(low_q_and_match.iloc[1, 1]),'Color': color_dict[source], 
                                      'Level': 1, 'Source':source, 'Value': low_q_and_match.iloc[1, 1]}
        nodes_df.loc[len(nodes_df)] = {'Label': 'normal quality cells <br>'+str(low_q_and_match.iloc[0, 1]),'Color': color_dict[source], 
                                      'Level': 1, 'Source':source, 'Value': low_q_and_match.iloc[0, 1]}
        #matched and unmatched cells
        unmatched_cells = low_q_and_match.iloc[0, 1] - low_q_and_match.iloc[0, 3] #normal minus matched = unmatched
        nodes_df.loc[len(nodes_df)] = {'Label': 'matched cells <br>'+str(low_q_and_match.iloc[0, 3]),'Color': color_dict[source], 
                                      'Level': 2, 'Source':source, 'Value': low_q_and_match.iloc[0, 3]}
        nodes_df.loc[len(nodes_df)] = {'Label': 'unmatched cells <br>'+str(unmatched_cells),'Color': color_dict[source], 
                                      'Level': 2, 'Source':source, 'Value': unmatched_cells}
        
        rem_col_names = [x for x in g.columns[5:] if source+'_unfilt' not in x] #get remaining columns 
        rem_col_counts = g[g['low_quality_cells']==False].loc[:, rem_col_names].isna().value_counts().reset_index() #want only normal quality cells
        name = ' '.join(rem_col_names[0].split('_'))+ '<br>'
        nodes_df.loc[len(nodes_df)] = {'Label': name+str(rem_col_counts.iloc[1,2]),'Color': color_dict[source], 
                                              'Level': 3, 'Source':source, 'Value': rem_col_counts.iloc[1,2]}
        rem_unmatched_cells = unmatched_cells-rem_col_counts.iloc[1,2]
        if rem_col_counts.loc[:, rem_col_names[-1]].unique() == False: #if no remaining cells
            name = 'unknown unmatched '+source+' cells <br>'
        else:
            name = ' '.join(rem_col_names[-1].split('_'))
        nodes_df.loc[len(nodes_df)] = {'Label': name+str(rem_unmatched_cells),'Color': color_dict[source], 
                                'Level': 3, 'Source':source, 'Value': rem_unmatched_cells}
    if save:
        nodes_df.to_csv(savepath+'/sankey_nodes_df'+barcode+'.csv')
    return nodes_df


def create_link_df_sankey(nodes_df, barcode, save=True, savepath='/allen/programs/celltypes/workgroups/hct/emilyg/reseg_project/new_seg/'):
    """
    Generates links dataframe from nodes dataframe, connecting relevant nodes to one another and ensuring distinct segmentations are separated. 
    INPUTS
        nodes_df: dataframe containing information on each node for sankey diagram.
        barcode: unique identifier for section, used for creating the save name for the dataframe, string.
        save: whether to save the final results. Bool, default true
        savepath: path to which results should be saved, string
    OUTPUTS
        links_df: dataframe containing source, target, value, and colors for sankey diagram creation
    """
    
    link_color_dict = {'red': 'rgb(205, 209, 228)', 'blue': 'rgb(205, 209, 228)'} 
    links_df = pd.DataFrame(columns=['Source', 'Target', 'Value', 'Link Color'])
    for i, row in nodes_df.iterrows():
        if 'total' in row.Label:
            continue #skip 0 since all comes from zero
        else: 
            target = i
            value = row['Value']
            source = nodes_df[(nodes_df['Level']==row['Level']-1) & (nodes_df['Source']==row['Source'])].index.values.tolist()
            if len(source)> 1 and row['Level'] != 3:
                source = [i for i, x in nodes_df.iloc[source, :].iterrows() if 'un' not in x['Label'] if 'low' not in x['Label']][0]
            elif len(source)> 1 and row['Level'] == 3:
                source = [i for i, x in nodes_df.iloc[source, :].iterrows() if 'un' in x['Label']][0]
            else:
                source = source[0]
            prev_row = nodes_df.iloc[source, :]
            link_color = link_color_dict[prev_row['Color']] 
            links_df.loc[len(links_df)] = {'Source':source, 'Target': target, 'Value': value, 'Link Color': link_color}
    
    #special case: matched cells should have same name, new source, new color
    matched_cell_rows = nodes_df[(nodes_df['Label'].str.contains('matched cells'))&(~nodes_df['Label'].str.contains('un'))]
    matched_cell_val = matched_cell_rows['Value'].values.tolist()[0]
    nodes_df.loc[len(nodes_df)]= {'Label': 'matched cells <br> '+ str(matched_cell_val), 'Color': 'violet', 
                                  'Level': matched_cell_rows['Level'].values.tolist()[0], 'Source': 'Both', 'Value': matched_cell_val}
    #find matched col for both in nodes_df, then find rows with those IDS in links df and connect them
    update_idxs = links_df[links_df['Target'].isin(matched_cell_rows.index.values.tolist())].index.values.tolist()
    links_df.loc[update_idxs, 'Target'] = [nodes_df[nodes_df['Source']=='Both'].index.values] * len(update_idxs)
    if save:
        links_df.to_csv(savepath+'/sankey_links_df.csv')
    return links_df

def generate_sankey_diagram(seg_comp_df, barcode, save=True, savepath='/allen/programs/celltypes/workgroups/hct/emilyg/reseg_project/new_seg/'):
    """
    Creates sankey diagram comparing two segmentation results using segmentation comparison dataframe created by collect_mutual_match_and_doublets function
    INPUTS
        seg_comp_df: dataframe with unique index describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, mutual matches, and putative doublets. 
        barcode: unique identifier for section, used for creating the save name for the dataframe, string.
        save: whether to save the final results. Bool, default true
        savepath: path to which results should be saved, string   
    OUTPUTS
        fig: sankey diagram figure
    """
    nodes_df = create_node_df_sankey(seg_comp_df, barcode, save, savepath)
    links_df = create_link_df_sankey(nodes_df, barcode, save, savepath)
    # Sankey plot setup
    data_trace = dict(type='sankey',node = dict(line = dict(color = "black",width = 0),
                                                label =  nodes_df['Label'].dropna(axis=0, how='any'),
                                                color = nodes_df['Color']
                                               ),
                      link = dict(
                          source = links_df['Source'].dropna(axis=0, how='any'),
                          target = links_df['Target'].dropna(axis=0, how='any'),
                          value = links_df['Value'].dropna(axis=0, how='any'),
                          color = links_df['Link Color'].dropna(axis=0, how='any'),
                      )
                     )
    
    source_1 = seg_comp_df.source.unique().tolist()[0]
    source_2 = seg_comp_df.source.unique().tolist()[1]
    layout = dict(title = "Cell segmentation comparison human "+barcode+' '+source_1+' vs '+source_2, height = 772, font = dict(size = 10))
    fig = go.Figure(dict(data=[data_trace], layout=layout))
    fig.show()
    #save interactive file
    if save:
        fig.write_html(savepath+"/segmentation_comparison_"+barcode+"_"+source_1+'_'+source_2+".html")
    return fig