import anndata
import anndata as ad
import inspect
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

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