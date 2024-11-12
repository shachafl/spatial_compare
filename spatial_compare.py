import anndata as ad
import inspect
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from matplotlib.ticker import FormatStrFormatter
import warnings



DEFAULT_DATA_NAMES=["Data 0", "Data 1"]

           


class SpatialCompare():
    """
    A class for comparing spatial data between two AnnData objects.

    Parameters:
    ----------
    ad_0 : AnnData
        The first AnnData object to compare.
    ad_1 : AnnData
        The second AnnData object to compare.
    data_names : list, optional
        A list of names for the two datasets, by default "Data 0" and "Data 1"
    obsm_key : str, optional
        The key for the spatial data in the `obsm` attribute of the AnnData objects, by default "spatial_cirro_grid".
    category : str, optional
        The category to compare in the `obs` attribute of the AnnData objects, by default "subclass".

    Attributes:
    ----------
    ad_0 : AnnData
        The first AnnData object.
    ad_1 : AnnData
        The second AnnData object.
    obsm_key : str
        The key for the spatial data in the `obsm` attribute.
    category : str
        The category to compare in the `obs` attribute.
    shared_genes : list
        The list of genes that are shared between the two datasets.
    can_compare : bool
        A flag indicating if the datasets can be compared.
    data_names : list
        A list of names for the two datasets.
    spatial_compare_results : None or dict
        The results of the spatial comparison.

    Methods:
    -------
    set_category(category)
        Set the category to compare.
    spatial_plot(plot_legend=True, min_cells_to_plot=10, decimate_for_spatial_plot=1, figsize=[20,10], category_values=[])
        Plot the spatial data for the two datasets.
    de_novo_cluster(plot_stuff=False, correspondence_level="leiden_1")
        Perform de novo clustering on the two datasets.
    find_matched_groups(n_top_groups=100, n_shared_groups=30, min_n_cells=100, category_values=[], exclude_group_string="zzzzzzzzzzzzzzz", plot_stuff=False, figsize=[10,10])
        Find matched groups between the two datasets.
    compare_expression(category_values=[], plot_stuff=False, min_mean_expression=.2, min_genes_to_compare=5, min_cells=10)
        Compare gene expression between the two datasets.

    """

    def __init__(self, ad_0:ad.AnnData, ad_1:ad.AnnData, data_names=DEFAULT_DATA_NAMES, obsm_key="spatial_cirro_grid", category="subclass"):
        """
        Initialize the SpatialCompare object.

        Parameters:
        ----------
        ad_0 : AnnData
            The first AnnData object to compare.
        ad_1 : AnnData
            The second AnnData object to compare.
        data_names : list, optional
            A list of names for the two datasets, by default DEFAULT_DATA_NAMES.
        obsm_key : str, optional
            The key for the spatial data in the `obsm` attribute of the AnnData objects, by default "spatial_cirro_grid".
        category : str, optional
            The category to compare in the `obs` attribute of the AnnData objects, by default "subclass".
        """
        # implementation details...
class SpatialCompare():
    def __init__(self, ad_0:ad.AnnData, ad_1:ad.AnnData, 
                data_names=DEFAULT_DATA_NAMES,
                obsm_key="spatial_cirro_grid", 
                category="subclass",
                ):
        
        
        # validate self.ad_0 and self.ad_1:
        if not isinstance(ad_0, ad.AnnData) or not isinstance(ad_1, ad.AnnData):
            raise ValueError("ad_0 and self.ad_1 must be AnnData objects")
        
        self.ad_0 = ad_0
        self.ad_1 = ad_1
       
        # make sure both AnnData objects have raw counts in X:          


        #validate spatial key in obsm
        if obsm_key not in self.ad_0.obsm.keys() or obsm_key not in self.ad_1.obsm.keys():
            raise ValueError(" obsm key "+obsm_key+" is not in both input AnnData objects")    
        
        self.obsm_key = obsm_key

        #vaidate category in obs:

        if category not in self.ad_0.obs.columns or category not in self.ad_1.obs.columns:
            self.category = None
            self.can_compare = False
            warnings.warn("category "+category+" is not in both input AnnData objects")
        
        self.category = category 



        for adobj in [self.ad_0, self.ad_1]:
            if "gene" not in adobj.var.columns:
                adobj.var["gene"]=adobj.var.index


        # identify the intersection of genes between the 2 datasets:
        self.shared_genes = list(set(self.ad_0.var.index)&set(self.ad_1.var.index))
        print("input anndata objects have "+str(len(self.shared_genes))+ " shared genes")




        self.can_compare = True

        self.data_names = data_names

        self.spatial_compare_results = None


    def set_category(self, category):
        if category not in self.ad_0.obs.columns or category not in self.ad_1.obs.columns:
            self.category = None
            self.can_compare = False
            raise Warning("category "+category+" is not in both input AnnData objects") 
        
        self.category = category 
        self.can_compare = True

    def spatial_plot(self, plot_legend=True,min_cells_to_plot = 10, decimate_for_spatial_plot=1,figsize=[20,10], category_values=[]):

        plt.figure(figsize=figsize)
        all_category_values = set(self.ad_0.obs[self.category].unique()) | set(self.ad_1.obs[self.category].unique()) 
        if len(category_values)==0:
            category_values = all_category_values

        for c in category_values:
            plt.subplot(1,2,1)
            plt.title(self.data_names[0])
            if np.sum(self.ad_0.obs[self.category]==c)>min_cells_to_plot:
                label = c+": "+str(np.sum(self.ad_0.obs[self.category]==c))
            else:
                label = None

            plt.plot(self.ad_0.obsm[self.obsm_key][self.ad_0.obs[self.category]==c,0][::decimate_for_spatial_plot], 
                 self.ad_0.obsm[self.obsm_key][self.ad_0.obs[self.category]==c,1][::decimate_for_spatial_plot],'.',
                     label = label, markersize=.5)
            plt.axis('equal')
            if plot_legend:
                plt.legend(markerscale=5)
            plt.subplot(1,2,2)
            plt.title(self.data_names[1])
            if np.sum(self.ad_1.obs[self.category]==c)>min_cells_to_plot:
                label = c+": "+str(np.sum(self.ad_1.obs[self.category]==c))
            else:
                label = None
            plt.plot(self.ad_1.obsm[self.obsm_key][self.ad_1.obs[self.category]==c,0][::decimate_for_spatial_plot], 
                 self.ad_1.obsm[self.obsm_key][self.ad_1.obs[self.category]==c,1][::decimate_for_spatial_plot],'.',
                     label = label, markersize=.5)
            plt.axis('equal')
            if plot_legend:
                plt.legend(markerscale=5)

    def de_novo_cluster(self, plot_stuff=False, correspondence_level = "leiden_1"):
        self.ad_0 = filter_and_cluster_twice(self.ad_0, plot_stuff=plot_stuff)
        self.ad_1 = filter_and_cluster_twice(self.ad_1, plot_stuff=plot_stuff)
        find_best_match_groups(self.ad_0,self.ad_1,group_names=[correspondence_level,correspondence_level])
        self.can_compare = True
        self.category = "matched_leiden_clusters"
        return True



    def find_matched_groups(self,
                        n_top_groups=100, 
                        n_shared_groups=30,
                        min_n_cells=100,
                        category_values=[],
                        exclude_group_string="zzzzzzzzzzzzzzz",
                        plot_stuff=False, figsize=[10,10]):

        sorted_counts={}
        _other_columns = [[c for c in self.ad_0.obs.columns if c!=self.category][0],[c for c in self.ad_1.obs.columns if c!=self.category][0]]

        for ii,df in enumerate([self.ad_0.obs,self.ad_1.obs]):    
            sorted_counts[ii] = df.loc[:,[self.category,_other_columns[ii]]].groupby(self.category,observed=True).count().sort_values(_other_columns[ii], ascending=False)
            
        in_top_N_0 = sorted_counts[0].loc[[c for c in sorted_counts[0].index if exclude_group_string not in c],:].head(n_top_groups)
        in_top_N_1 = sorted_counts[1].loc[[c for c in sorted_counts[1].index if exclude_group_string not in c],:].head(n_top_groups)
        
        if len(category_values)>0:
            in_top_N_0 = in_top_N_0.loc[in_top_N_0.index.isin(category_values),:]
            in_top_N_1 = in_top_N_1.loc[in_top_N_1.index.isin(category_values),:]


        shared_top=list(in_top_N_0.loc[in_top_N_0.index.isin(in_top_N_1.index),:].index)
        
        if len(shared_top) > n_shared_groups:
            shared_top = shared_top[:n_shared_groups+1]
        
        prop_corr = np.corrcoef( in_top_N_0.loc[shared_top,:][_other_columns[0]],in_top_N_1.loc[shared_top,:][_other_columns[1]])[1][0]
        
        if plot_stuff:
                plt.figure(figsize=figsize)
                plt.loglog(in_top_N_0.loc[shared_top,:][_other_columns[0]],in_top_N_1.loc[shared_top,:][_other_columns[1]],'.')
                for g in shared_top:
                    # long names are annoying to read
                    if len(g)<50:
                        plt.text(in_top_N_0.loc[g,:][_other_columns[0]],in_top_N_1.loc[g,:][_other_columns[1]], g)
                plt.xlabel(self.data_names[0])
                plt.ylabel(self.data_names[1])
                plt.title("cell type abundances\n correlation = "+str(prop_corr)[:6])
                plt.plot([np.min(in_top_N_0.loc[shared_top,:][_other_columns[0]]),np.max(in_top_N_0.loc[shared_top,:][_other_columns[0]])],
                        [np.min(in_top_N_0.loc[shared_top,:][_other_columns[0]]),np.max(in_top_N_0.loc[shared_top,:][_other_columns[0]])],'--')
                        
                plt.axis('equal')
        return {"category_top_values":shared_top,
                "category":self.category,
                "proportion_correlation":prop_corr,
            "n0":in_top_N_0,
            "n1":in_top_N_1}
    
    def compare_expression(self, category_values=[], plot_stuff=False,
                      min_mean_expression=.2,
                      min_genes_to_compare =5, min_cells = 10):
        # group cells
        if len(category_values)==0:
            raise ValueError("please supply a list of values for the category "+self.category)

 
        category_records=[]
        gene_ratio_dfs = {}
        for category_value in category_values:
            group_mask_0 = self.ad_0.obs[self.category]==category_value
            group_mask_1 = self.ad_1.obs[self.category]==category_value
            
            
            if np.sum(group_mask_0)< min_cells or np.sum(group_mask_1)< min_cells:
                print("at least 1 input has less than "+str(min_cells)+" cells in "+self.category+" == "+category_value)
                continue

            means_0 = np.array(np.mean(self.ad_0[group_mask_0, self.ad_0.var.index.isin(self.shared_genes)].X,axis=0)).flatten()
            means_1 = np.array(np.mean(self.ad_1[group_mask_1, self.ad_1.var.index.isin(self.shared_genes)].X,axis=0)).flatten()
            means_0_gt_min = np.nonzero(means_0>min_mean_expression)[0]
            means_1_gt_min = np.nonzero(means_1>min_mean_expression)[0]
            above_means0 = self.ad_0.var[self.ad_0.var.index.isin(self.shared_genes)].iloc[means_0_gt_min]
            above_means1 = self.ad_1.var[self.ad_1.var.index.isin(self.shared_genes)].iloc[means_1_gt_min]
            shared_above_mean = [g for g in above_means1.index if g in above_means0.index]
            if len(shared_above_mean)<min_genes_to_compare:
                print( self.category+" "+category_value+" has less than "+str(min_genes_to_compare)+"\n shared genes above minimum mean = "+str(min_mean_expression))
                continue
            
            means_0 = np.array(np.mean(self.ad_0[group_mask_0, shared_above_mean].X,axis=0)).flatten()
            means_1 = np.array(np.mean(self.ad_1[group_mask_1, shared_above_mean].X,axis=0)).flatten()
            shared_genes = shared_above_mean
            p_coef = np.polynomial.Polynomial.fit(means_0,means_1,1).convert().coef
            category_records.append({self.category:category_value,
                                    "slope":p_coef[1],
                                    "mean_ratio": np.mean(means_1/means_0),
                                    "correlation":np.corrcoef(means_0,means_1)[0][1],
                                    "n_cells_0": np.sum(group_mask_0),
                                    "n_cells_1":np.sum(group_mask_1),
                                    "total_count_ratio": np.sum(self.ad_1[group_mask_1, shared_genes].X)/np.sum(self.ad_0[group_mask_0, shared_genes].X)
                                    })
            
            gene_ratio_dfs[category_value] = pd.DataFrame(means_1/means_0,columns=[self.data_names[1]+" / "+self.data_names[0]+" ratio"], index = shared_genes)
            
            if plot_stuff:
                plt.figure(figsize=[10,10])
                plt.title(self.category+": "+category_value+"\nmean counts per cell\ncorrelation: "+str(category_records[-1]["correlation"])[:4]+" mean ratio: "+str(category_records[-1]["mean_ratio"])[:4])
                low_expression = np.logical_and(means_0<1.0,means_1<1.0)
                plt.loglog(means_0[low_expression],means_1[low_expression],'.', color = [.5,0.5,0.5])
                plt.loglog(means_0[np.logical_not(low_expression)],means_1[np.logical_not(low_expression)],'.')

                plt.xlabel(self.data_names[0]+", N = "+str(np.sum(group_mask_0)))
                plt.ylabel(self.data_names[1]+", N = "+str(np.sum(group_mask_1)))

                
                for g in shared_genes:
                    if (means_0[np.nonzero(np.array(shared_genes)==g)]==0 or 
                            means_1[np.nonzero(np.array(shared_genes)==g)]==0) or low_expression[np.array(shared_genes)==g]:
                        continue
                    plt.text(means_0[np.nonzero(np.array(shared_genes)==g)]
                            ,means_1[np.nonzero(np.array(shared_genes)==g)], g, fontsize=10)
                plt.plot([np.min(means_0), np.max(means_0)], 
                        [np.min(means_0), np.max(means_0)],'--')
        print(gene_ratio_dfs.keys())
        if len(gene_ratio_dfs.keys())>0:
            
            gene_ratio_df = pd.concat(gene_ratio_dfs,axis=1)
        else:
            gene_ratio_df = None
        return {"data_names":self.data_names,
                "category_results":pd.DataFrame.from_records(category_records),
                "gene_ratio_dataframe":gene_ratio_df}

    def plot_detection_ratio(self,gene_ratio_dataframe, figsize=[15,15]):

        detection_ratio_plots(gene_ratio_dataframe,
                              data_names=self.data_names, figsize=figsize)

    def spatial_compare(self, **kwargs):

        if "category" in kwargs.keys():
            self.set_category(kwargs["category"])


        if not self.can_compare:
            raise ValueError("cannot compare AnnData objects without a category.\n run de novo clustering or supply a category with set_category()")
        

        fmr_kwargs = {key:kwargs[key] for key in inspect.signature(SpatialCompare.find_matched_groups).parameters.keys() if key in kwargs.keys()}
        ce_kwargs = {key:kwargs[key] for key in inspect.signature(SpatialCompare.compare_expression).parameters.keys()if key in kwargs.keys() and key != "category_values"}

        match_results = self.find_matched_groups(
                                        **fmr_kwargs)


        expression_results = self.compare_expression(
                       category_values=match_results["category_top_values"],
                       **ce_kwargs
                      )
        match_results["expression_results"] = expression_results
        return match_results

    def run_and_plot(self, **kwargs):
        if "category" in kwargs.keys():
            self.set_category(kwargs["category"])
    
        self.spatial_plot()
        self.spatial_compare_results = self.spatial_compare(plot_stuff=True, **kwargs)
        self.plot_detection_ratio(self.spatial_compare_results["expression_results"]["gene_ratio_dataframe"],figsize=[30,20])
        return True



def filter_and_cluster_twice(input_ad, n_hvgs=2000, 
                       min_max_counts=3,
                       min_transcript_counts=50,
                       n_pcs=[50,20],
                       n_iterations=[1,1], plot_stuff=True):
    """ 
    filter genes and cells, then run 2 rounds of Leiden clustering on an anndata object 
    resulting clusters are named "leiden_XX_YY" where "XX" and "YY" are the numbers of the first and second
    round clusters
    
    this function returns a modified copy of the anndata object with potentially fewer genes and cells!
    """
    
    # switch to array instead of sparse matrix:
    if isinstance(input_ad.X, sparse.csr.csr_matrix):
        input_ad.X = input_ad.X.toarray()

    print("converted to array  ")
    low_detection_genes = input_ad.var.iloc[np.nonzero(np.max(input_ad.X, axis = 0)<=min_max_counts)].gene
    
    # throw out genes if no cells have more than 3 counts
    to_cluster = input_ad[input_ad.obs.transcript_counts>=min_transcript_counts,[g for g in input_ad.var.gene if g not in low_detection_genes]]

    # this is to prevent multiple "leiden" columns from being added to the obs dataframe
    to_cluster.obs = to_cluster.obs.loc[:,[c for c in to_cluster.obs.columns if c not in ["leiden","leiden_0","leiden_1"]]]
    
    
    # Normalizing to median total counts
    print("normalizing total counts")
    sc.pp.normalize_total(to_cluster)
    # Logarithmize the data
    sc.pp.log1p(to_cluster)
    
    sc.pp.highly_variable_genes(to_cluster, n_top_genes=n_hvgs)
    
    sc.tl.pca(to_cluster,n_comps=n_pcs[0])
    
    sc.pp.neighbors(to_cluster)
    
    if plot_stuff:
        sc.tl.umap(to_cluster)
        sc.pl.pca_variance_ratio(to_cluster, n_pcs[0], log=True)
    

    sc.tl.leiden(to_cluster, n_iterations=n_iterations[0])
    to_cluster.obs.rename(columns={"leiden":"leiden_0"},inplace=True)
    
    if plot_stuff:
        sc.pl.umap(to_cluster, color=["leiden_0"],size=2)
    
    # per cluster, repeat PCA and clustering...
    # this duplicates the data! but it's necessary because the scanpy functions would create a copy of the subset anyway(?)
    print("2nd round of clustering")
    all_subs=[]
    for cl in to_cluster.obs.leiden_0.unique():
        subcopy = to_cluster[to_cluster.obs.leiden_0==cl,:].copy()
        print(subcopy.shape)
        sc.tl.pca(subcopy,n_comps=n_pcs[1])    
        sc.pp.neighbors(subcopy)
        sc.tl.leiden(subcopy, n_iterations=n_iterations[1])
        if plot_stuff:
            sc.pl.umap(subcopy, color=["leiden"])
        subcopy.obs["leiden_1"]=["leiden_"+str(cl).zfill(2)+"_"+str(g).zfill(2) for g in subcopy.obs.leiden]

        all_subs.append(subcopy)
    print("concatenating")
    iterative_clusters_ad = ad.concat(all_subs)

    return iterative_clusters_ad



def detection_ratio_plots(gene_ratio_df,data_names=DEFAULT_DATA_NAMES, figsize=[15,15]):
    
    sorted_genes=[str(s) for s in gene_ratio_df.mean(axis=1).sort_values().index.values]
    
    
    plt.figure(figsize=figsize)
    plt.subplot(3,1,1)
    p = sns.boxplot(gene_ratio_df.loc[sorted_genes,:].T, )
    p.set_yscale('log')
    p.set_xlabel("gene",fontsize=20)
    p.set_ylabel("detection ratio\n"+data_names[1]+" / "+data_names[0],fontsize=20)
    ax=plt.gca()
    ax.tick_params(axis='x',labelrotation=45, labelsize=10)
    ax.tick_params(axis='y',labelsize=20, which="major")
    ax.tick_params(axis='y',labelsize=10, which="minor")

    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    plt.plot(ax.get_xlim(), [1,1],'-')
    
    # sort columns for boxplot based on mean detection ratio:
    sorted_categories = gene_ratio_df.columns[np.argsort(gene_ratio_df.mean(axis=0))]

    plt.subplot(3,1,2)
    ax =sns.boxplot(gene_ratio_df.loc[:,sorted_categories])
    plt.ylabel("transcript detection ratio "+data_names[1]+" / "+data_names[0])
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
    plt.plot(ax.get_xlim(), [1,1])

    plt.subplot(3,1,3)
    h = plt.hist(np.log10(gene_ratio_df.values.ravel()), bins =100)
    median_ratio = np.nanmedian(gene_ratio_df.values.ravel())
    mean_ratio = np.nanmean(gene_ratio_df.values.ravel())
    plt.plot(np.log10(median_ratio)*np.ones([2]),[0,1.1*h[0].max()], label = "median ratio "+data_names[1]+" / "+data_names[0]+" : "+str(median_ratio )[:5])
    plt.plot(np.log10(mean_ratio)*np.ones([2]),[0,1.1*h[0].max()], label = "mean ratio "+data_names[1]+" / "+data_names[0]+" : "+str(mean_ratio )[:5])
    plt.legend()
    plt.title("ratio of cell segmented counts for shared genes\n "+data_names[1]+" / "+data_names[0], fontsize=20)
    plt.xlabel("log10(counts ratio)",fontsize=20)
    ax=plt.gca()
    ax.tick_params(axis='y',labelsize=12)
    plt.tight_layout()
    
def generate_label_confusion(input_ad,
                             column_pairs=[["iterative_subclass","subclass_name"],
                                           ["iterative_leiden","supertype_name"],
                                           ["iterative_leiden","subclass_name"]]):
    
    # identifies fraction correctly matching labels in pairs of columns supplied by user in `column_pairs`
    
    
    column_pair_results = []
    for column_pair in column_pairs:
        row_dfs = []
        for g in input_ad.obs[column_pairs[0][0]].unique():
            g_df = input_ad.obs.loc[input_ad.obs[column_pairs[0][0]]==g,column_pairs[0]].groupby(column_pairs[0][1],observed=False).count()
            g_df.rename(columns={column_pairs[0][0]:column_pairs[0][0]+"_"+g},inplace=True)
            row_dfs.append(g_df)
        tdf = pd.concat(row_dfs, axis=1)
        tdf_norm = tdf.copy().astype(float)
        for ii in tdf_norm.index:
            tdf_norm.loc[ii,:] = tdf_norm.loc[ii,:]/tdf_norm.loc[ii,:].sum()
        
        column_pair_results.append(tdf_norm)
        
    return {c_p[0]+"_"+c_p[1]:column_pair_results[ii] for ii, c_p in enumerate(column_pairs)}

def find_best_match_groups(ad0,ad1,group_names=["leiden_1","leiden_1"], in_place=True):

    g_o_m_0 = grouped_obs_mean(ad0, group_names[0])
    g_o_m_1 = grouped_obs_mean(ad1, group_names[1])
    
    # these can have different numbers of clusters and different sets of observed genes. 
    # subset to shared genes for this and deal with non-square correlation matrix
    
    shared_genes = list(g_o_m_0.loc[g_o_m_0.index.isin(g_o_m_1.index),:].index)
    correlation_array = np.corrcoef(pd.concat([g_o_m_0.loc[shared_genes,:], g_o_m_1.loc[shared_genes,:]], axis=1),rowvar=0)
    shape0 = g_o_m_0.loc[shared_genes,:].shape
    shape1 = g_o_m_1.loc[shared_genes,:].shape
    array_01= correlation_array[:shape0[1],shape0[1]:]
    array_10= array_01.T
    
    match_01={r:np.argmax(array_01[r,:]) for r in range(array_01.shape[0])}
    match_10={r:np.argmax(array_10[r,:]) for r in range(array_10.shape[0])}
    
    # find mutual matches:
    mutual_matches={}
    for k in match_01:
        if match_10[match_01[k]] == k:
            mutual_matches.update({k:match_01[k]})
    
    
    INPUT_CLUSTER_NAME = "leiden_0"


    if in_place:
        ad0.obs["matched_"+group_names[0]] = ""
        ad1.obs["matched_"+group_names[0]] = ""
    
        for mm in mutual_matches:
            
            match_mask0 = ad0.obs[INPUT_CLUSTER_NAME]==g_o_m_0.columns[mm]
            ad0.obs.loc[match_mask0,["matched_leiden_clusters"] ] = "matched_leiden_cluster_"+str(mm)
        
            match_mask1 = ad1.obs[INPUT_CLUSTER_NAME]==g_o_m_1.columns[mutual_matches[mm]]
            ad1.obs.loc[match_mask1,["matched_leiden_clusters"] ] = "matched_leiden_cluster_"+str(mm)
    else:
        ad0_out = ad0.copy()
        ad1_out = ad1.copy()
        
        ad0_out.obs["matched_"+group_names[0]] = ""
        ad1_out.obs["matched_"+group_names[0]] = ""
    
        for mm in mutual_matches:
            
            match_mask0 = ad0_out.obs[INPUT_CLUSTER_NAME]==g_o_m_0.columns[mm]
            ad0_out.obs.loc[match_mask0,["matched_"+group_names[0]] ] = "matched_leiden_cluster_"+str(mm)
        
            match_mask1 = ad1_out.obs[INPUT_CLUSTER_NAME]==g_o_m_1.columns[mutual_matches[mm]]
            ad1_out.obs.loc[match_mask1,["matched_"+group_names[0]] ] = "matched_leiden_cluster_"+str(mm)
        return (ad0,ad1)

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)).tolist()    
        
    return out


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




def spatial_detection_scores(reference: pd.DataFrame, query: pd.DataFrame,
                              plot_stuff=True,query_name: str="query data",comparison_column="transcript_counts",
                              category = "supercluster_name",
                              n_bins = 50, in_place=True, non_spatial = False):
    """
    Calculate and plot spatial detection scores for query data compared to reference data.
    
    Parameters:
        reference (pd.DataFrame): The reference data.
        query (pd.DataFrame): The query data.
        plot_stuff (bool, optional): Whether to plot the results. Defaults to True.
        query_name (str, optional): The name of the query data. Defaults to "query data".
        category (str, optional): The category column to compare. Defaults to "supercluster_name".
        n_bins (int, optional): The number of bins for spatial grouping. Defaults to 50.
        in_place (bool, optional): Whether to modify the query data in place. Defaults to True.
        non_spatial (bool, optional): Whether to compare to an ungrouped mean/std. Defaults to False.

    Returns:
        dict: A dictionary containing the bin image, extent, query data, and reference data (if in_place is False).
    """
    # code goes here

    if category not in reference.columns or category not in query.columns:
        raise ValueError("category "+category+" not in reference and query inputs")


    shared_category_values = list(set(reference[category].unique())& set(query[category].unique()))
    if (len(shared_category_values)<query[category].unique().shape[0] or 
        len(shared_category_values)<reference[category].unique().shape[0]):
        in_place=False
    
    
    if in_place:
        s2 = query.loc[query[category].isin(shared_category_values),:]
        s1 = reference.loc[reference[category].isin(shared_category_values),:]
    else:
        s2 = query.loc[query[category].isin(shared_category_values),:].copy()
        s1 = reference.loc[reference[category].isin(shared_category_values),:].copy()





    means=s1.groupby(category, observed=True)[comparison_column].mean()
    stds = s1.groupby(category, observed=True)[comparison_column].std()

    # if you want to compare to an ungrouped mean/std, try this:
    if non_spatial:
        means[:] = means.mean()
        stds[:] = stds.mean()
    

    s2["detection_relative_z_score"]=0.
    s2["detection_difference"]=0.
    s2["detection_ratio"]=0.

    for c,gb in s2.groupby(category, observed=True):
        if c not in shared_category_values:
            continue

        s2.loc[s2[category]==c,["detection_relative_z_score"]] = ( (s2.loc[s2[category]==c,[comparison_column]] - means[c])/stds[c]).values
        s2.loc[s2[category]==c,["detection_difference"]]=  (s2.loc[s2[category]==c,[comparison_column]] - means[c]).values
        s2.loc[s2[category]==c,["log_10_detection_ratio"]]=  np.log10((s2.loc[s2[category]==c,[comparison_column]] / means[c]).values)
    
    
    s2["xy_bucket"] = list(zip(pd.cut(s2.x_centroid, n_bins, labels=list(range(n_bins))),
                                   pd.cut(s2.y_centroid, n_bins, labels=list(range(n_bins)))))
    
    binx = s2.groupby("xy_bucket").x_centroid.mean()
    biny = s2.groupby("xy_bucket").y_centroid.mean()
    
    z_score = s2.groupby("xy_bucket").detection_relative_z_score.mean()
    difference = s2.groupby("xy_bucket").detection_difference.mean()
    log_ratio = s2.groupby("xy_bucket").log_10_detection_ratio.mean()
    n_cells = s2.groupby("xy_bucket").x_centroid.count()

    
    bin_image_z_score= np.zeros([n_bins,n_bins])
    bin_image_difference= np.zeros([n_bins,n_bins])
    bin_image_ratio= np.zeros([n_bins,n_bins])
    bin_image_counts = np.zeros([n_bins,n_bins])

    extent = [np.min(binx), np.max(binx),np.min(biny),np.max(biny)]
    for coord in binx.index:
        bin_image_z_score[coord[1],coord[0]] = z_score[coord]
        bin_image_difference[coord[1],coord[0]] = difference[coord]
        bin_image_ratio[coord[1],coord[0]] = log_ratio[coord]
        bin_image_counts[coord[1],coord[0]]= n_cells[coord]

    
    if plot_stuff:
        if non_spatial:
            title_string = "Non-spatial Detection Scores"
        else:
            title_string = "Spatial Detection Scores"
        min_maxes = {"detection z-score": [bin_image_z_score,[-1, 1]],
                     "total counts difference": [bin_image_difference,[-100, 100]],
                     "log10(detection ratio)": [bin_image_ratio,[-1, 1]]}
                     

        fig, axs = plt.subplots(1, 3,figsize=[15,5])
        if non_spatial:
            fig.suptitle(title_string+"\n"+query_name)
        else:
            fig.suptitle(title_string+"\n"+query_name+" grouped by "+category+" and spatially binned")
        for ii,plot_name in enumerate(min_maxes.keys()):
            ax = axs[ ii]
            pcm = ax.imshow( min_maxes[plot_name][0],extent = extent,
                                cmap='coolwarm_r', vmin=min_maxes[plot_name][1][0],vmax=min_maxes[plot_name][1][1])
            fig.colorbar(pcm, ax=ax,shrink=0.7)
            ax.set_title(query_name+"\n"+plot_name)

    if in_place:
        
        return dict(z_score_image=bin_image_z_score,
                difference_image=bin_image_difference,
                ratio_image=bin_image_ratio,
                extent = extent,
                count_image = bin_image_counts,
                query=True,reference=True)

    else:
        return dict(z_score_image=bin_image_z_score,
                difference_image=bin_image_difference,
                ratio_image=bin_image_ratio,
                extent = extent,
                count_image = bin_image_counts,
                query=s2,reference=s1)
                