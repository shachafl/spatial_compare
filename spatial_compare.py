import anndata
import anndata as ad
import inspect
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def find_matched_groups(df0,df1,data_names=["data 0", "data 1"],
                        category="subclass",
                        n_top_groups=100, 
                        n_shared_groups=30,
                        min_n_cells=100,
                       exclude_group_string="NN",
                       plot_stuff=False, figsize=[10,10]):

    sorted_counts={}
    _other_columns = [[c for c in df0.columns if c!=category][0],[c for c in df1.columns if c!=category][0]]

    for ii,df in enumerate([df0,df1]):    
        sorted_counts[ii] = df.loc[:,[category,_other_columns[ii]]].groupby(category,observed=True).count().sort_values(_other_columns[ii], ascending=False)
        
        
    in_top_N_0 = sorted_counts[0].loc[[c for c in sorted_counts[0].index if exclude_group_string not in c],:].head(n_top_groups)
    in_top_N_1 = sorted_counts[1].loc[[c for c in sorted_counts[1].index if exclude_group_string not in c],:].head(n_top_groups)
    
    
    shared_top=list(in_top_N_0.loc[in_top_N_0.index.isin(in_top_N_1.index),:].index)
    
    if len(shared_top) > n_shared_groups:
        shared_top = shared_top[:n_shared_groups+1]
    
    prop_corr = np.corrcoef( in_top_N_0.loc[shared_top,:][_other_columns[0]],in_top_N_1.loc[shared_top,:][_other_columns[1]])[1][0]
    
    if plot_stuff:
            plt.figure(figsize=figsize)
            plt.loglog(in_top_N_0.loc[shared_top,:][_other_columns[0]],in_top_N_1.loc[shared_top,:][_other_columns[1]],'.')
            for g in shared_top:
                if len(g)>50:
                    continue
                plt.text(in_top_N_0.loc[g,:][_other_columns[0]],in_top_N_1.loc[g,:][_other_columns[1]], g)
            plt.text(np.mean(in_top_N_0.loc[shared_top,:][_other_columns[0]])*.5,
                     np.mean(in_top_N_1.loc[shared_top,:][_other_columns[1]])*1.5, "correlation = "+str(prop_corr)[:6])
            plt.xlabel(data_names[0])
            plt.ylabel(data_names[1])
            plt.title("cell type abundances")
            plt.plot([np.min(in_top_N_0.loc[shared_top,:][_other_columns[0]]),np.max(in_top_N_0.loc[shared_top,:][_other_columns[0]])],
                     [np.min(in_top_N_0.loc[shared_top,:][_other_columns[0]]),np.max(in_top_N_0.loc[shared_top,:][_other_columns[0]])],'--')
                     
            plt.axis('equal')
    return {"category_top_values":shared_top,
            "category":category,
            "proportion_correlation":prop_corr,
           "n0":in_top_N_0,
           "n1":in_top_N_1}
            

def compare_expression(ad_0,ad_1,data_names=["data 0", "data 1"],
                       category="cluster_name",
                       category_values=[], plot_stuff=False,
                      min_mean_expression=2,
                      min_genes_to_compare =5):
    # group cells
    if len(category_values)==0:
        raise ValueError("please supply a list of values for the category "+category)

    # identify the intersection of genes between the 2 datasets:
    shared_genes = list(set(ad_0.var.index)&set(ad_1.var.index))
    print(len(shared_genes))
    category_records=[]
    gene_ratio_dfs = {}
    for category_value in category_values:
        group_mask_0 = ad_0.obs[category]==category_value
        group_mask_1 = ad_1.obs[category]==category_value
        
    
        if np.sum(group_mask_0)==0 or np.sum(group_mask_1)==0:
            raise ValueError("at least 1 input has no cells in "+category+" == "+category_value)

        means_0 = np.array(np.mean(ad_0[group_mask_0, shared_genes].X,axis=0)).flatten()
        means_1 = np.array(np.mean(ad_1[group_mask_1, shared_genes].X,axis=0)).flatten()
        means_0_gt_min = np.nonzero(means_0>min_mean_expression)[0]
        means_1_gt_min = np.nonzero(means_1>min_mean_expression)[0]
        shared_above_mean = list(set(list(means_0_gt_min))& set(list(means_1_gt_min)))
        if len(shared_above_mean)<min_genes_to_compare:
            print(means_0)
            print(means_1)
            raise ValueError("less than "+str(min_genes_to_compare)+" shared genes above minimum mean = "+str(min_mean_expression))

        
        means_0 = means_0[shared_above_mean]
        means_1 = means_1[shared_above_mean]
        shared_genes = np.array(shared_genes)[shared_above_mean]
        p_coef = np.polynomial.Polynomial.fit(means_0,means_1,1).convert().coef
        category_records.append({category:category_value,
                                 "slope":p_coef[1],
                                 "mean_ratio": np.mean(means_1/means_0),
                                "correlation":np.corrcoef(means_0,means_1)[0][1],
                                "n_cells_0": np.sum(group_mask_0),
                                "n_cells_1":np.sum(group_mask_1),
                                "total_count_ratio": np.sum(ad_1[group_mask_1, shared_genes].X)/np.sum(ad_0[group_mask_0, shared_genes].X)
                                })
        
        gene_ratio_dfs[category_value] = pd.DataFrame(means_1/means_0,columns=["data1_data0_ratio"], index = shared_genes)
        
        if plot_stuff:
            plt.figure(figsize=[10,10])
            plt.title(category+": "+category_value+"\nmean counts per cell\ncorrelation: "+str(np.corrcoef(means_0,means_1)[0][1])[:4])
            low_expression = np.logical_and(means_0<1.0,means_1<1.0)
            plt.loglog(means_0[low_expression],means_1[low_expression],'.', color = [.5,0.5,0.5])
            plt.loglog(means_0[np.logical_not(low_expression)],means_1[np.logical_not(low_expression)],'.')

            plt.xlabel(data_names[0]+", N = "+str(np.sum(group_mask_0)))
            plt.ylabel(data_names[1]+", N = "+str(np.sum(group_mask_1)))

            
            for g in shared_genes:
                if (means_0[np.nonzero(np.array(shared_genes)==g)]==0 or 
                        means_1[np.nonzero(np.array(shared_genes)==g)]==0) or low_expression[np.array(shared_genes)==g]:
                    continue
                plt.text(means_0[np.nonzero(np.array(shared_genes)==g)]
                         ,means_1[np.nonzero(np.array(shared_genes)==g)], g, fontsize=10)
            plt.plot([np.min(means_0), np.max(means_0)], 
                     [np.min(means_0), np.max(means_0)],'--')
    gene_ratio_df = pd.concat(gene_ratio_dfs,axis=1)
    
    return {"category_results":pd.DataFrame.from_records(category_records),
            "gene_ratio_dataframe":gene_ratio_df}



def spatial_compare(ad_0:ad.AnnData, ad_1:ad.AnnData,
                    spatial_plot=False,obsm_key="spatial_cirro_grid", 
                    plot_legend=True,min_cells_to_plot = 10, **kwargs):

    fmr_kwargs = {key:kwargs[key] for key in inspect.signature(find_matched_groups).parameters.keys() if key in kwargs.keys()}
    ce_kwargs = {key:kwargs[key] for key in inspect.signature(compare_expression).parameters.keys()if key in kwargs.keys()}

    if spatial_plot:
        decimate_for_plot = 1
        plt.figure(figsize=[20,10])
        all_category_values = set(ad_0.obs[fmr_kwargs["category"]].unique()) | set(ad_1.obs[fmr_kwargs["category"]].unique()) 
        for c in all_category_values:
            plt.subplot(1,2,1)
            if np.sum(ad_0.obs[fmr_kwargs["category"]]==c)>min_cells_to_plot:
                label = c+": "+str(np.sum(ad_0.obs[fmr_kwargs["category"]]==c))
            else:
                label = None
            plt.plot(ad_0.obsm[obsm_key][ad_0.obs[fmr_kwargs["category"]]==c,0][::decimate_for_plot], 
                 ad_0.obsm[obsm_key][ad_0.obs[fmr_kwargs["category"]]==c,1][::decimate_for_plot],'.',
                     label = label, markersize=.5)
            plt.axis('equal')
            if plot_legend:
                plt.legend()
            plt.subplot(1,2,2)
            if np.sum(ad_1.obs[fmr_kwargs["category"]]==c)>min_cells_to_plot:
                label = c+": "+str(np.sum(ad_1.obs[fmr_kwargs["category"]]==c))
            else:
                label = None
            plt.plot(ad_1.obsm[obsm_key][ad_1.obs[fmr_kwargs["category"]]==c,0][::decimate_for_plot], 
                 ad_1.obsm[obsm_key][ad_1.obs[fmr_kwargs["category"]]==c,1][::decimate_for_plot],'.',
                     label = label, markersize=.5)
            plt.axis('equal')
            if plot_legend:
                plt.legend()
    match_results = find_matched_groups(ad_0.obs,ad_1.obs, 
                                        **fmr_kwargs)

    
    expression_results = compare_expression(ad_0,ad_1,
                       category_values=match_results["category_top_values"],
                       **ce_kwargs
                      )
    match_results["expression_results"] = expression_results
    return match_results



# function to run 3 rounds of Leiden clustering to an anndata object

def filter_and_cluster(input_ad, n_hvgs=2000, 
                       min_max_counts=3,min_cell_area=300,
                       min_transcript_counts=50,
                       n_pcs=[50,20,10],
                       n_iterations=[5,5,5]):
    
    low_detection_genes = input_ad.var.iloc[np.nonzero(np.max(input_ad.X, axis = 0)<=min_max_counts)].gene
    
    low_detection_genes
    
    # plt.figure(figsize=[20,20])
    # sns.pairplot(input_ad.obs.loc[np.logical_and(input_ad.obs.cell_area<min_cell_area,
    #                                                input_ad.obs.transcript_counts<min_transcript_counts),["transcript_counts", "cell_area","segmentation_method"]],
    #              hue="segmentation_method", markers='.')
    
    # throw out genes if no cells have more than 3 counts
    # and cells with less than 100 total counts (5k panel). 50 counts for 300 gene panel)
    # Normalizing to median total counts
    to_cluster = input_ad[input_ad.obs.transcript_counts>=min_transcript_counts,np.logical_not(input_ad.var.gene.isin(low_detection_genes))].copy()
    sc.pp.normalize_total(to_cluster)
    # Logarithmize the data
    sc.pp.log1p(to_cluster)
    
    sc.pp.highly_variable_genes(to_cluster, n_top_genes=n_hvgs)
    
    sc.tl.pca(to_cluster,n_comps=n_pcs[0])
    
    sc.pp.neighbors(to_cluster)
    
    sc.tl.umap(to_cluster)
    
    sc.pl.pca_variance_ratio(to_cluster, n_pcs[0], log=True)
    
    sc.pl.umap(
        to_cluster,
        color="subclass_name",
        # Setting a smaller point size to get prevent overlap
        size=2,
    )


    sc.tl.leiden(to_cluster, n_iterations=n_iterations[0])

    sc.pl.umap(to_cluster, color=["leiden","subclass_name"])
    
    # per cluster, repeat PCA and clustering...
    all_subs=[]
    for cl in to_cluster.obs.leiden.unique():
        print(cl)
        subcopy = to_cluster[to_cluster.obs.leiden==cl,:].copy()
        subcopy.obs.drop(columns=["leiden"], inplace=True)
        sc.tl.pca(subcopy,n_comps=n_pcs[1])    
        sc.pp.neighbors(subcopy)
        sc.tl.leiden(subcopy, n_iterations=n_iterations[1])
        sc.pl.umap(subcopy, color=["leiden","subclass_name", "supertype_name"])
        subcopy.obs["iterative_leiden"]=[cl+"_"+g for g in subcopy.obs.leiden]
        all_subs.append(subcopy)
    
    iterative_clusters_ad = ad.concat(all_subs)
    
    
    # one more round of clustering...
    # 
    all_subs=[]
    for cl in iterative_clusters_ad.obs.iterative_leiden.unique():
        subcopy = iterative_clusters_ad[iterative_clusters_ad.obs.iterative_leiden==cl,:].copy()
        subcopy.obs.drop(columns=["leiden"], inplace=True)
        sc.tl.pca(subcopy,n_comps=n_pcs[2])    
        sc.pp.neighbors(subcopy)
        sc.tl.leiden(subcopy, n_iterations=n_iterations[2])
        #sc.pl.umap(subcopy, color=["leiden","subclass_name", "supertype_name"])
        subcopy.obs["iterative_leiden_2"]=[cl+"_"+g for g in subcopy.obs.leiden]
        all_subs.append(subcopy)
    
    iterative_clusters_ad = ad.concat(all_subs)
    return iterative_clusters_ad






def generate_label_confusion(input_ad,
                             column_pairs=[["iterative_subclass","subclass_name"],
                                           ["iterative_leiden","supertype_name"],
                                           ["iterative_leiden","subclass_name"]]):
    
    
    #TODO make this function iterate over user-defined column pairs instead of repeating code
    #subclass confusion
    subclass_row_dfs = []
    for g in input_ad.obs[column_pairs[0][0]].unique():
        g_df = input_ad.obs.loc[input_ad.obs[column_pairs[0][0]]==g,column_pairs[0]].groupby(column_pairs[0][1],observed=False).count()
        g_df.rename(columns={column_pairs[0][0]:column_pairs[0][0]+"_"+g},inplace=True)
        subclass_row_dfs.append(g_df)
    tdf = pd.concat(subclass_row_dfs, axis=1)
    tdf_norm = tdf.copy().astype(float)
    for ii in tdf_norm.index:
        tdf_norm.loc[ii,:] = tdf_norm.loc[ii,:]/tdf_norm.loc[ii,:].sum()
    
    
    #supertype confusion
    supertype_row_dfs = []
    for g in input_ad.obs[column_pairs[1][0]].unique():
        g_df = input_ad.obs.loc[input_ad.obs[column_pairs[1][0]]==g,column_pairs[1]].groupby(column_pairs[1][1],observed=False).count()
        g_df.rename(columns={column_pairs[1][0]:column_pairs[1][0]+"_"+g},inplace=True)
        supertype_row_dfs.append(g_df)
    
    tdf2 = pd.concat(supertype_row_dfs, axis=1)
    tdf2_norm = tdf2.copy().astype(float)
    for ii in tdf2_norm.index:
        tdf2_norm.loc[ii,:] = tdf2_norm.loc[ii,:]/tdf2_norm.loc[ii,:].sum()
    
    
    
    #supertype subclass confusion
    supertype_subclass_row_dfs = []
    for g in input_ad.obs[column_pairs[2][0]].unique():
        g_df = input_ad.obs.loc[input_ad.obs[column_pairs[2][0]]==g,column_pairs[2]].groupby(column_pairs[2][1],observed=False).count()
        g_df.rename(columns={column_pairs[2][0]:column_pairs[2][0]+"_"+g},inplace=True)
        supertype_subclass_row_dfs.append(g_df)
    
    tdf3 = pd.concat(supertype_subclass_row_dfs, axis=1)
    tdf3_norm = tdf3.copy().astype(float)
    for ii in tdf3_norm.index:
        tdf3_norm.loc[ii,:] = tdf3_norm.loc[ii,:]/tdf3_norm.loc[ii,:].sum()

    return {"subclass_confusion":tdf_norm,
            "supertype_confusion":tdf2_norm,
            "subclass_iterative_supertype_confusion":tdf3_norm}




def get_column_ordering(df, ordered_rows):

    # function to get sorted columns matched to input ordered rows
    # runs through each entry in `ordered_rows` and adds first non-repeated column from the ranked list to that category
    # repeats until all the columns are done
    # ranking is based on the values in the array, in this case it is the fraction of cells that mapped to each input row
    flat_ordering = []
    empty_columns={r:[] for r in ordered_rows}
    while len(flat_ordering)<df.shape[1]:
        for r in ordered_rows:
            r_columns =[df.columns[rr] for rr in df.loc[r,:].argsort()[::-1] if  not (df.columns[rr] in flat_ordering)]
            if len(r_columns)>0 and len(flat_ordering)<df.shape[1]:
                empty_columns[r].append(r_columns[0])
                flat_ordering.append(r_columns[0])
                
    output=[]
    [output.extend(empty_columns[k]) for k in empty_columns]
    return output
