import anndata as ad
import inspect
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.spatial import cKDTree
from pathlib import Path
import plotly.graph_objs as go
from matplotlib.ticker import FormatStrFormatter
import warnings


DEFAULT_DATA_NAMES = ["Data 0", "Data 1"]


class SpatialCompare:
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

    def __init__(
        self,
        ad_0: ad.AnnData,
        ad_1: ad.AnnData,
        data_names=DEFAULT_DATA_NAMES,
        obsm_key="spatial_cirro_grid",
        category="subclass",
    ):
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


class SpatialCompare:
    def __init__(
        self,
        ad_0: ad.AnnData,
        ad_1: ad.AnnData,
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

        # validate spatial key in obsm
        if (
            obsm_key not in self.ad_0.obsm.keys()
            or obsm_key not in self.ad_1.obsm.keys()
        ):
            raise ValueError(
                " obsm key " + obsm_key + " is not in both input AnnData objects"
            )

        self.obsm_key = obsm_key

        # vaidate category in obs:

        if (
            category not in self.ad_0.obs.columns
            or category not in self.ad_1.obs.columns
        ):
            self.category = None
            self.can_compare = False
            warnings.warn(
                "category " + category + " is not in both input AnnData objects"
            )

        self.category = category

        for adobj in [self.ad_0, self.ad_1]:
            if "gene" not in adobj.var.columns:
                adobj.var["gene"] = adobj.var.index

        # identify the intersection of genes between the 2 datasets:
        self.shared_genes = list(set(self.ad_0.var.index) & set(self.ad_1.var.index))
        print(
            "input anndata objects have "
            + str(len(self.shared_genes))
            + " shared genes"
        )

        self.can_compare = True

        self.data_names = data_names

        self.spatial_compare_results = None

    def set_category(self, category):
        if (
            category not in self.ad_0.obs.columns
            or category not in self.ad_1.obs.columns
        ):
            self.category = None
            self.can_compare = False
            raise Warning(
                "category " + category + " is not in both input AnnData objects"
            )

        self.category = category
        self.can_compare = True

    def spatial_plot(
        self,
        plot_legend=True,
        min_cells_to_plot=10,
        decimate_for_spatial_plot=1,
        figsize=[20, 10],
        category_values=[],
    ):

        plt.figure(figsize=figsize)
        all_category_values = set(self.ad_0.obs[self.category].unique()) | set(
            self.ad_1.obs[self.category].unique()
        )
        if len(category_values) == 0:
            category_values = all_category_values

        for c in category_values:
            plt.subplot(1, 2, 1)
            plt.title(self.data_names[0])
            if np.sum(self.ad_0.obs[self.category] == c) > min_cells_to_plot:
                label = c + ": " + str(np.sum(self.ad_0.obs[self.category] == c))
            else:
                label = None

            plt.plot(
                self.ad_0.obsm[self.obsm_key][self.ad_0.obs[self.category] == c, 0][
                    ::decimate_for_spatial_plot
                ],
                self.ad_0.obsm[self.obsm_key][self.ad_0.obs[self.category] == c, 1][
                    ::decimate_for_spatial_plot
                ],
                ".",
                label=label,
                markersize=0.5,
            )
            plt.axis("equal")
            if plot_legend:
                plt.legend(markerscale=5)
            plt.subplot(1, 2, 2)
            plt.title(self.data_names[1])
            if np.sum(self.ad_1.obs[self.category] == c) > min_cells_to_plot:
                label = c + ": " + str(np.sum(self.ad_1.obs[self.category] == c))
            else:
                label = None
            plt.plot(
                self.ad_1.obsm[self.obsm_key][self.ad_1.obs[self.category] == c, 0][
                    ::decimate_for_spatial_plot
                ],
                self.ad_1.obsm[self.obsm_key][self.ad_1.obs[self.category] == c, 1][
                    ::decimate_for_spatial_plot
                ],
                ".",
                label=label,
                markersize=0.5,
            )
            plt.axis("equal")
            if plot_legend:
                plt.legend(markerscale=5)

    def de_novo_cluster(self, plot_stuff=False, correspondence_level="leiden_1"):
        self.ad_0 = filter_and_cluster_twice(self.ad_0, plot_stuff=plot_stuff)
        self.ad_1 = filter_and_cluster_twice(self.ad_1, plot_stuff=plot_stuff)
        find_best_match_groups(
            self.ad_0,
            self.ad_1,
            group_names=[correspondence_level, correspondence_level],
        )
        self.can_compare = True
        self.category = "matched_leiden_clusters"
        return True

    def find_matched_groups(
        self,
        n_top_groups=100,
        n_shared_groups=30,
        min_n_cells=100,
        category_values=[],
        exclude_group_string="zzzzzzzzzzzzzzz",
        plot_stuff=False,
        figsize=[10, 10],
    ):

        sorted_counts = {}
        _other_columns = [
            [c for c in self.ad_0.obs.columns if c != self.category][0],
            [c for c in self.ad_1.obs.columns if c != self.category][0],
        ]

        for ii, df in enumerate([self.ad_0.obs, self.ad_1.obs]):
            sorted_counts[ii] = (
                df.loc[:, [self.category, _other_columns[ii]]]
                .groupby(self.category, observed=True)
                .count()
                .sort_values(_other_columns[ii], ascending=False)
            )

        in_top_N_0 = (
            sorted_counts[0]
            .loc[
                [c for c in sorted_counts[0].index if exclude_group_string not in c], :
            ]
            .head(n_top_groups)
        )
        in_top_N_1 = (
            sorted_counts[1]
            .loc[
                [c for c in sorted_counts[1].index if exclude_group_string not in c], :
            ]
            .head(n_top_groups)
        )

        if len(category_values) > 0:
            in_top_N_0 = in_top_N_0.loc[in_top_N_0.index.isin(category_values), :]
            in_top_N_1 = in_top_N_1.loc[in_top_N_1.index.isin(category_values), :]

        shared_top = list(
            in_top_N_0.loc[in_top_N_0.index.isin(in_top_N_1.index), :].index
        )

        if len(shared_top) > n_shared_groups:
            shared_top = shared_top[: n_shared_groups + 1]

        prop_corr = np.corrcoef(
            in_top_N_0.loc[shared_top, :][_other_columns[0]],
            in_top_N_1.loc[shared_top, :][_other_columns[1]],
        )[1][0]

        if plot_stuff:
            plt.figure(figsize=figsize)
            plt.loglog(
                in_top_N_0.loc[shared_top, :][_other_columns[0]],
                in_top_N_1.loc[shared_top, :][_other_columns[1]],
                ".",
            )
            for g in shared_top:
                # long names are annoying to read
                if len(g) < 50:
                    plt.text(
                        in_top_N_0.loc[g, :][_other_columns[0]],
                        in_top_N_1.loc[g, :][_other_columns[1]],
                        g,
                    )
            plt.xlabel(self.data_names[0])
            plt.ylabel(self.data_names[1])
            plt.title("cell type abundances\n correlation = " + str(prop_corr)[:6])
            plt.plot(
                [
                    np.min(in_top_N_0.loc[shared_top, :][_other_columns[0]]),
                    np.max(in_top_N_0.loc[shared_top, :][_other_columns[0]]),
                ],
                [
                    np.min(in_top_N_0.loc[shared_top, :][_other_columns[0]]),
                    np.max(in_top_N_0.loc[shared_top, :][_other_columns[0]]),
                ],
                "--",
            )

            plt.axis("equal")
        return {
            "category_top_values": shared_top,
            "category": self.category,
            "proportion_correlation": prop_corr,
            "n0": in_top_N_0,
            "n1": in_top_N_1,
        }

    def compare_expression(
        self,
        category_values=[],
        plot_stuff=False,
        min_mean_expression=0.2,
        min_genes_to_compare=5,
        min_cells=10,
    ):
        # group cells
        if len(category_values) == 0:
            raise ValueError(
                "please supply a list of values for the category " + self.category
            )

        category_records = []
        gene_ratio_dfs = {}
        for category_value in category_values:
            group_mask_0 = self.ad_0.obs[self.category] == category_value
            group_mask_1 = self.ad_1.obs[self.category] == category_value

            if np.sum(group_mask_0) < min_cells or np.sum(group_mask_1) < min_cells:
                print(
                    "at least 1 input has less than "
                    + str(min_cells)
                    + " cells in "
                    + self.category
                    + " == "
                    + category_value
                )
                continue

            means_0 = np.array(
                np.mean(
                    self.ad_0[
                        group_mask_0, self.ad_0.var.index.isin(self.shared_genes)
                    ].X,
                    axis=0,
                )
            ).flatten()
            means_1 = np.array(
                np.mean(
                    self.ad_1[
                        group_mask_1, self.ad_1.var.index.isin(self.shared_genes)
                    ].X,
                    axis=0,
                )
            ).flatten()
            means_0_gt_min = np.nonzero(means_0 > min_mean_expression)[0]
            means_1_gt_min = np.nonzero(means_1 > min_mean_expression)[0]
            above_means0 = self.ad_0.var[
                self.ad_0.var.index.isin(self.shared_genes)
            ].iloc[means_0_gt_min]
            above_means1 = self.ad_1.var[
                self.ad_1.var.index.isin(self.shared_genes)
            ].iloc[means_1_gt_min]
            shared_above_mean = [
                g for g in above_means1.index if g in above_means0.index
            ]
            if len(shared_above_mean) < min_genes_to_compare:
                print(
                    self.category
                    + " "
                    + category_value
                    + " has less than "
                    + str(min_genes_to_compare)
                    + "\n shared genes above minimum mean = "
                    + str(min_mean_expression)
                )
                continue

            means_0 = np.array(
                np.mean(self.ad_0[group_mask_0, shared_above_mean].X, axis=0)
            ).flatten()
            means_1 = np.array(
                np.mean(self.ad_1[group_mask_1, shared_above_mean].X, axis=0)
            ).flatten()
            shared_genes = shared_above_mean
            p_coef = np.polynomial.Polynomial.fit(means_0, means_1, 1).convert().coef
            category_records.append(
                {
                    self.category: category_value,
                    "slope": p_coef[1],
                    "mean_ratio": np.mean(means_1 / means_0),
                    "correlation": np.corrcoef(means_0, means_1)[0][1],
                    "n_cells_0": np.sum(group_mask_0),
                    "n_cells_1": np.sum(group_mask_1),
                    "total_count_ratio": np.sum(self.ad_1[group_mask_1, shared_genes].X)
                    / np.sum(self.ad_0[group_mask_0, shared_genes].X),
                }
            )

            gene_ratio_dfs[category_value] = pd.DataFrame(
                means_1 / means_0,
                columns=[self.data_names[1] + " / " + self.data_names[0] + " ratio"],
                index=shared_genes,
            )

            if plot_stuff:
                plt.figure(figsize=[10, 10])
                plt.title(
                    self.category
                    + ": "
                    + category_value
                    + "\nmean counts per cell\ncorrelation: "
                    + str(category_records[-1]["correlation"])[:4]
                    + " mean ratio: "
                    + str(category_records[-1]["mean_ratio"])[:4]
                )
                low_expression = np.logical_and(means_0 < 1.0, means_1 < 1.0)
                plt.loglog(
                    means_0[low_expression],
                    means_1[low_expression],
                    ".",
                    color=[0.5, 0.5, 0.5],
                )
                plt.loglog(
                    means_0[np.logical_not(low_expression)],
                    means_1[np.logical_not(low_expression)],
                    ".",
                )

                plt.xlabel(self.data_names[0] + ", N = " + str(np.sum(group_mask_0)))
                plt.ylabel(self.data_names[1] + ", N = " + str(np.sum(group_mask_1)))

                for g in shared_genes:
                    if (
                        means_0[np.nonzero(np.array(shared_genes) == g)] == 0
                        or means_1[np.nonzero(np.array(shared_genes) == g)] == 0
                    ) or low_expression[np.array(shared_genes) == g]:
                        continue
                    plt.text(
                        means_0[np.nonzero(np.array(shared_genes) == g)],
                        means_1[np.nonzero(np.array(shared_genes) == g)],
                        g,
                        fontsize=10,
                    )
                plt.plot(
                    [np.min(means_0), np.max(means_0)],
                    [np.min(means_0), np.max(means_0)],
                    "--",
                )
        print(gene_ratio_dfs.keys())
        if len(gene_ratio_dfs.keys()) > 0:

            gene_ratio_df = pd.concat(gene_ratio_dfs, axis=1)
        else:
            gene_ratio_df = None
        return {
            "data_names": self.data_names,
            "category_results": pd.DataFrame.from_records(category_records),
            "gene_ratio_dataframe": gene_ratio_df,
        }

    def plot_detection_ratio(self, gene_ratio_dataframe, figsize=[15, 15]):

        detection_ratio_plots(
            gene_ratio_dataframe, data_names=self.data_names, figsize=figsize
        )

    def spatial_compare(self, **kwargs):

        if "category" in kwargs.keys():
            self.set_category(kwargs["category"])

        if not self.can_compare:
            raise ValueError(
                "cannot compare AnnData objects without a category.\n run de novo clustering or supply a category with set_category()"
            )

        fmr_kwargs = {
            key: kwargs[key]
            for key in inspect.signature(
                SpatialCompare.find_matched_groups
            ).parameters.keys()
            if key in kwargs.keys()
        }
        ce_kwargs = {
            key: kwargs[key]
            for key in inspect.signature(
                SpatialCompare.compare_expression
            ).parameters.keys()
            if key in kwargs.keys() and key != "category_values"
        }

        match_results = self.find_matched_groups(**fmr_kwargs)

        expression_results = self.compare_expression(
            category_values=match_results["category_top_values"], **ce_kwargs
        )
        match_results["expression_results"] = expression_results
        return match_results

    def run_and_plot(self, **kwargs):
        if "category" in kwargs.keys():
            self.set_category(kwargs["category"])

        self.spatial_plot()
        self.spatial_compare_results = self.spatial_compare(plot_stuff=True, **kwargs)
        self.plot_detection_ratio(
            self.spatial_compare_results["expression_results"]["gene_ratio_dataframe"],
            figsize=[30, 20],
        )
        return True

    def collect_mutual_match_and_doublets(
        self,
        bc,
        save=True,
        nn_dist=2.5,
        reuse_saved=True,
        savepath=None,
        min_transcripts=40,
    ):
        """
        Runs all relevant functions in required order, generating dataframe summarizing comparisons between two segmentations for a single section.
        INPUTS
        bc: unique identifier for section segmentation comparison was run on. String.
        seg_name_a, seg_name_b: names of segmentations to be compared, string
        save: whether to save the intermediate and final results. Bool, default true
        nn_dist: distance cutoff for identifying likely same cell between segmentations. float, default 2.5 (um)
        min_transcripts: minimum number of transcripts needed to define a cell too low quality to be considered for mapping
        savepath: path to which you'd like to save your results. Required if using anndata objects and save = True, otherwise defaults to seg_b_path
        OUTPUTS:
        seg_comp_df: dataframe with unique index describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, mutual matches, and putative doublets.
        """

        if not savepath:
            savepath = seg_b_path
        # grab base comparison df
        seg_comp_df = get_segmentation_data(
            bc,
            self.ad_0,
            self.ad_1,
            self.data_names[0],
            self.data_names[1],
            save,
            savepath,
            reuse_saved,
            min_transcripts,
        )

        # mutual matches with both segmentations filtered
        dfa = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[0])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        dfb = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[1])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        col_name_both_filt = (
            "match_" + self.data_names[0] + "_filt_" + self.data_names[1] + "_filt"
        )
        seg_comp_df[col_name_both_filt] = seg_comp_df.index.map(
            get_mutual_matches(dfa, dfb, nn_dist)
        )  # filt both

        # mutual matches with seg a filtered only
        dfa = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[0])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        dfb = seg_comp_df[seg_comp_df["source"] == self.data_names[1]]
        col_name_afilt = (
            "match_" + self.data_names[0] + "_filt_" + self.data_names[1] + "_unfilt"
        )
        seg_comp_df[col_name_afilt] = seg_comp_df.index.map(
            get_mutual_matches(dfa, dfb, nn_dist)
        )  # a filt only

        # mutual matches with seg b filtered only
        dfa = seg_comp_df[seg_comp_df["source"] == self.data_names[0]]
        dfb = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[1])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        col_name_bfilt = (
            "match_" + self.data_names[0] + "_unfilt_" + self.data_names[1] + "_filt"
        )
        seg_comp_df[col_name_bfilt] = seg_comp_df.index.map(
            get_mutual_matches(dfa, dfb, nn_dist)
        )  # b filt only

        # save results
        if save:
            seg_comp_df.to_csv(
                savepath
                + bc
                + "_seg_comp_df_"
                + self.data_names[0]
                + "_and_"
                + self.data_names[1]
                + "_populated.csv",
                index=True,
            )
            print(
                "Saved to: "
                + savepath
                + bc
                + "_seg_comp_df_"
                + self.data_names[0]
                + "_and_"
                + self.data_names[1]
                + "_populated.csv"
            )
        return seg_comp_df

    def generate_sankey_diagram(
        self,
        seg_comp_df,
        bc,
        save=True,
        savepath="/allen/programs/celltypes/workgroups/hct/emilyg/",
    ):
        """
        Creates sankey diagram comparing two segmentation results using segmentation comparison dataframe created by collect_mutual_match_and_doublets function
        INPUTS
            seg_comp_df: dataframe with unique index describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, mutual matches, and putative doublets.
            bc: unique identifier for section, used for creating the save name for the dataframe, string.
            save: whether to save the final results. Bool, default true
            savepath: path to which results should be saved, string
        OUTPUTS
            fig: sankey diagram figure
        """
        nodes_df, unknown_unmatched_cells = create_node_df_sankey(
            seg_comp_df, bc, save, savepath
        )
        links_df = create_link_df_sankey(nodes_df, bc, save, savepath)
        # Sankey plot setup
        data_trace = dict(
            type="sankey",
            node=dict(
                line=dict(color="black", width=0),
                label=nodes_df["Label"].dropna(axis=0, how="any"),
                color=nodes_df["Color"],
            ),
            link=dict(
                source=links_df["Source"].dropna(axis=0, how="any"),
                target=links_df["Target"].dropna(axis=0, how="any"),
                value=links_df["Value"].dropna(axis=0, how="any"),
                color=links_df["Link Color"].dropna(axis=0, how="any"),
            ),
        )

        layout = dict(
            title="Cell segmentation comparison human "
            + bc
            + " "
            + self.data_names[0]
            + " vs "
            + self.data_names[1],
            height=772,
            font=dict(size=10),
        )
        fig = go.Figure(dict(data=[data_trace], layout=layout))
        # save interactive file
        if save:
            filepath = (
                savepath
                + "/segmentation_comparison_"
                + bc
                + "_"
                + self.data_names[0]
                + "_"
                + self.data_names[1]
                + ".html"
            )
            fig.write_html(filepath)
        fig.show()
        return fig, unknown_unmatched_cells, filepath

    def scaling_check(self, seg_comp_df):
        """
        View two sections plotted on same axes, to check for scaling issues
        INPUTS
            seg_comp_df:dataframe with unique index describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, mutual matches, and putative doublets.
        """
        # overview of both sections plotted on same axes
        fig, ax = plt.subplots()
        filtered_seg_a_df = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[0])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        filtered_seg_a_df.plot(
            "center_x",
            "center_y",
            kind="scatter",
            s=0.1,
            label="high quality cells " + self.data_names[0],
            ax=ax,
            color="blue",
            alpha=0.5,
        )
        filtered_seg_b_df = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[1])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        filtered_seg_b_df.plot(
            "center_x",
            "center_y",
            kind="scatter",
            s=0.3,
            label="high quality cells " + self.data_names[1],
            ax=ax,
            color="red",
            alpha=0.5,
        )
        plt.title("Both segmentations, overlaid")
        plt.axis("equal")

        # find center of one section to zoom in on
        arr = filtered_seg_a_df[["center_x", "center_y"]].values
        center = [np.mean(arr[:, 0]), np.mean(arr[:, 1])]
        x_range = [center[0] - center[0] * 0.05, center[0] + center[0] * 0.05]
        y_range = [center[1] - center[1] * 0.05, center[1] + center[1] * 0.05]

        # plot same area for both sections
        fig, ax = plt.subplots()
        filtered_seg_a_df = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[0])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        subs_a = filtered_seg_a_df[
            (filtered_seg_a_df["center_x"].between(x_range[0], x_range[1]))
            & (filtered_seg_a_df["center_y"].between(y_range[0], y_range[1]))
        ]
        subs_a.plot(
            "center_x",
            "center_y",
            kind="scatter",
            s=2,
            label=self.data_names[0],
            ax=ax,
            color="blue",
            alpha=0.5,
        )
        plt.title("Both segmentations, overlaid")
        filtered_seg_b_df = seg_comp_df[
            (seg_comp_df["source"] == self.data_names[1])
            & (seg_comp_df["low_quality_cells"] == False)
        ]
        subs_b = filtered_seg_b_df[
            (filtered_seg_b_df["center_x"].between(x_range[0], x_range[1]))
            & (filtered_seg_b_df["center_y"].between(y_range[0], y_range[1]))
        ]
        subs_b.plot(
            "center_x",
            "center_y",
            kind="scatter",
            s=10,
            label=self.data_names[1],
            ax=ax,
            facecolors="none",
            edgecolors="r",
            alpha=0.5,
        )
        plt.title("Both segmentations, overlaid, zoomed in to center")
        plt.axis("equal")
        return


def filter_and_cluster_twice(
    input_ad,
    n_hvgs=2000,
    min_max_counts=3,
    min_transcript_counts=50,
    n_pcs=[50, 20],
    n_iterations=[1, 1],
    plot_stuff=True,
):
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
    low_detection_genes = input_ad.var.iloc[
        np.nonzero(np.max(input_ad.X, axis=0) <= min_max_counts)
    ].gene

    # throw out genes if no cells have more than 3 counts
    to_cluster = input_ad[
        input_ad.obs.transcript_counts >= min_transcript_counts,
        [g for g in input_ad.var.gene if g not in low_detection_genes],
    ]

    # this is to prevent multiple "leiden" columns from being added to the obs dataframe
    to_cluster.obs = to_cluster.obs.loc[
        :,
        [
            c
            for c in to_cluster.obs.columns
            if c not in ["leiden", "leiden_0", "leiden_1"]
        ],
    ]

    # Normalizing to median total counts
    print("normalizing total counts")
    sc.pp.normalize_total(to_cluster)
    # Logarithmize the data
    sc.pp.log1p(to_cluster)

    sc.pp.highly_variable_genes(to_cluster, n_top_genes=n_hvgs)

    sc.tl.pca(to_cluster, n_comps=n_pcs[0])

    sc.pp.neighbors(to_cluster)

    if plot_stuff:
        sc.tl.umap(to_cluster)
        sc.pl.pca_variance_ratio(to_cluster, n_pcs[0], log=True)

    sc.tl.leiden(to_cluster, n_iterations=n_iterations[0])
    to_cluster.obs.rename(columns={"leiden": "leiden_0"}, inplace=True)

    if plot_stuff:
        sc.pl.umap(to_cluster, color=["leiden_0"], size=2)

    # per cluster, repeat PCA and clustering...
    # this duplicates the data! but it's necessary because the scanpy functions would create a copy of the subset anyway(?)
    print("2nd round of clustering")
    all_subs = []
    for cl in to_cluster.obs.leiden_0.unique():
        subcopy = to_cluster[to_cluster.obs.leiden_0 == cl, :].copy()
        print(subcopy.shape)
        sc.tl.pca(subcopy, n_comps=n_pcs[1])
        sc.pp.neighbors(subcopy)
        sc.tl.leiden(subcopy, n_iterations=n_iterations[1])
        if plot_stuff:
            sc.pl.umap(subcopy, color=["leiden"])
        subcopy.obs["leiden_1"] = [
            "leiden_" + str(cl).zfill(2) + "_" + str(g).zfill(2)
            for g in subcopy.obs.leiden
        ]

        all_subs.append(subcopy)
    print("concatenating")
    iterative_clusters_ad = ad.concat(all_subs)

    return iterative_clusters_ad


def detection_ratio_plots(
    gene_ratio_df, data_names=DEFAULT_DATA_NAMES, figsize=[15, 15]
):

    sorted_genes = [
        str(s) for s in gene_ratio_df.mean(axis=1).sort_values().index.values
    ]

    plt.figure(figsize=figsize)
    plt.subplot(3, 1, 1)
    p = sns.boxplot(
        gene_ratio_df.loc[sorted_genes, :].T,
    )
    p.set_yscale("log")
    p.set_xlabel("gene", fontsize=20)
    p.set_ylabel(
        "detection ratio\n" + data_names[1] + " / " + data_names[0], fontsize=20
    )
    ax = plt.gca()
    ax.tick_params(axis="x", labelrotation=45, labelsize=10)
    ax.tick_params(axis="y", labelsize=20, which="major")
    ax.tick_params(axis="y", labelsize=10, which="minor")

    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    plt.plot(ax.get_xlim(), [1, 1], "-")

    # sort columns for boxplot based on mean detection ratio:
    sorted_categories = gene_ratio_df.columns[np.argsort(gene_ratio_df.mean(axis=0))]

    plt.subplot(3, 1, 2)
    ax = sns.boxplot(gene_ratio_df.loc[:, sorted_categories])
    plt.ylabel("transcript detection ratio " + data_names[1] + " / " + data_names[0])
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha="right")
    plt.plot(ax.get_xlim(), [1, 1])

    plt.subplot(3, 1, 3)
    h = plt.hist(np.log10(gene_ratio_df.values.ravel()), bins=100)
    median_ratio = np.nanmedian(gene_ratio_df.values.ravel())
    mean_ratio = np.nanmean(gene_ratio_df.values.ravel())
    plt.plot(
        np.log10(median_ratio) * np.ones([2]),
        [0, 1.1 * h[0].max()],
        label="median ratio "
        + data_names[1]
        + " / "
        + data_names[0]
        + " : "
        + str(median_ratio)[:5],
    )
    plt.plot(
        np.log10(mean_ratio) * np.ones([2]),
        [0, 1.1 * h[0].max()],
        label="mean ratio "
        + data_names[1]
        + " / "
        + data_names[0]
        + " : "
        + str(mean_ratio)[:5],
    )
    plt.legend()
    plt.title(
        "ratio of cell segmented counts for shared genes\n "
        + data_names[1]
        + " / "
        + data_names[0],
        fontsize=20,
    )
    plt.xlabel("log10(counts ratio)", fontsize=20)
    ax = plt.gca()
    ax.tick_params(axis="y", labelsize=12)
    plt.tight_layout()


def generate_label_confusion(
    input_ad,
    column_pairs=[
        ["iterative_subclass", "subclass_name"],
        ["iterative_leiden", "supertype_name"],
        ["iterative_leiden", "subclass_name"],
    ],
):

    # identifies fraction correctly matching labels in pairs of columns supplied by user in `column_pairs`

    column_pair_results = []
    for column_pair in column_pairs:
        row_dfs = []
        for g in input_ad.obs[column_pairs[0][0]].unique():
            g_df = (
                input_ad.obs.loc[input_ad.obs[column_pairs[0][0]] == g, column_pairs[0]]
                .groupby(column_pairs[0][1], observed=False)
                .count()
            )
            g_df.rename(
                columns={column_pairs[0][0]: column_pairs[0][0] + "_" + g}, inplace=True
            )
            row_dfs.append(g_df)
        tdf = pd.concat(row_dfs, axis=1)
        tdf_norm = tdf.copy().astype(float)
        for ii in tdf_norm.index:
            tdf_norm.loc[ii, :] = tdf_norm.loc[ii, :] / tdf_norm.loc[ii, :].sum()

        column_pair_results.append(tdf_norm)

    return {
        c_p[0] + "_" + c_p[1]: column_pair_results[ii]
        for ii, c_p in enumerate(column_pairs)
    }


def find_best_match_groups(
    ad0, ad1, group_names=["leiden_1", "leiden_1"], in_place=True
):

    g_o_m_0 = grouped_obs_mean(ad0, group_names[0])
    g_o_m_1 = grouped_obs_mean(ad1, group_names[1])

    # these can have different numbers of clusters and different sets of observed genes.
    # subset to shared genes for this and deal with non-square correlation matrix

    shared_genes = list(g_o_m_0.loc[g_o_m_0.index.isin(g_o_m_1.index), :].index)
    correlation_array = np.corrcoef(
        pd.concat([g_o_m_0.loc[shared_genes, :], g_o_m_1.loc[shared_genes, :]], axis=1),
        rowvar=0,
    )
    shape0 = g_o_m_0.loc[shared_genes, :].shape
    shape1 = g_o_m_1.loc[shared_genes, :].shape
    array_01 = correlation_array[: shape0[1], shape0[1] :]
    array_10 = array_01.T

    match_01 = {r: np.argmax(array_01[r, :]) for r in range(array_01.shape[0])}
    match_10 = {r: np.argmax(array_10[r, :]) for r in range(array_10.shape[0])}

    # find mutual matches:
    mutual_matches = {}
    for k in match_01:
        if match_10[match_01[k]] == k:
            mutual_matches.update({k: match_01[k]})

    INPUT_CLUSTER_NAME = "leiden_0"

    if in_place:
        ad0.obs["matched_" + group_names[0]] = ""
        ad1.obs["matched_" + group_names[0]] = ""

        for mm in mutual_matches:

            match_mask0 = ad0.obs[INPUT_CLUSTER_NAME] == g_o_m_0.columns[mm]
            ad0.obs.loc[match_mask0, ["matched_leiden_clusters"]] = (
                "matched_leiden_cluster_" + str(mm)
            )

            match_mask1 = (
                ad1.obs[INPUT_CLUSTER_NAME] == g_o_m_1.columns[mutual_matches[mm]]
            )
            ad1.obs.loc[match_mask1, ["matched_leiden_clusters"]] = (
                "matched_leiden_cluster_" + str(mm)
            )
    else:
        ad0_out = ad0.copy()
        ad1_out = ad1.copy()

        ad0_out.obs["matched_" + group_names[0]] = ""
        ad1_out.obs["matched_" + group_names[0]] = ""

        for mm in mutual_matches:

            match_mask0 = ad0_out.obs[INPUT_CLUSTER_NAME] == g_o_m_0.columns[mm]
            ad0_out.obs.loc[match_mask0, ["matched_" + group_names[0]]] = (
                "matched_leiden_cluster_" + str(mm)
            )

            match_mask1 = (
                ad1_out.obs[INPUT_CLUSTER_NAME] == g_o_m_1.columns[mutual_matches[mm]]
            )
            ad1_out.obs.loc[match_mask1, ["matched_" + group_names[0]]] = (
                "matched_leiden_cluster_" + str(mm)
            )
        return (ad0, ad1)


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
        index=adata.var_names,
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
    empty_columns = {r: [] for r in ordered_rows}
    while len(flat_ordering) < df.shape[1]:
        for r in ordered_rows:
            r_columns = [
                df.columns[rr]
                for rr in df.loc[r, :].argsort()[::-1]
                if not (df.columns[rr] in flat_ordering)
            ]
            if len(r_columns) > 0 and len(flat_ordering) < df.shape[1]:
                empty_columns[r].append(r_columns[0])
                flat_ordering.append(r_columns[0])

    output = []
    [output.extend(empty_columns[k]) for k in empty_columns]
    return output


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
    cell_check = [x for x in seg_path.glob("*.csv") if "cellpose" in x.stem]
    if cell_check:
        cxg = pd.read_table(
            str(seg_path) + "/cellpose-cell-by-gene.csv", index_col=0, sep=","
        )
        metadata = pd.read_table(
            str(seg_path) + "/cellpose_metadata.csv", index_col=0, sep=","
        )
    else:
        cxg = pd.read_table(str(seg_path) + "/cell_by_gene.csv", index_col=0, sep=",")
        meta = ad.read_h5ad(str(seg_path) + "/metadata.h5ad")
        metadata = meta.obs

    # assemble seg df with filt cells col, xy cols, and seg name
    high_quality_cells = (
        cxg.index[np.where((transcripts_per_cell(cxg) >= min_transcripts))[0]]
        .astype(str)
        .values.tolist()
    )
    seg_df = metadata[["center_x", "center_y"]].copy()
    seg_df.index = seg_df.index.astype(str)
    seg_df.loc[:, "source"] = seg_name
    seg_df.loc[:, "low_quality_cells"] = np.where(
        seg_df.index.isin(high_quality_cells), False, True
    )
    seg_df.index = seg_name + "_" + seg_df.index
    return seg_df


def get_segmentation_data(
    bc,
    anndata_a,
    anndata_b,
    seg_name_a,
    seg_name_b,
    save,
    savepath,
    reuse_saved,
    min_transcripts,
):
    """
    Loads segmentation data from anndata object
    Inputs:
        bc: unique section identifier, string
        seg_name_a, seg_name_b: names of segmentations to be compared, string
        anndata_path_a, anndata_path_b: path to segmented anndata objects, string
        save: whether to save the intermediary dataframe (useful if running functions individually), bool
        savepath: path to which data should be saved. Defaults to seg_b_path if not specified.
        min_transcripts: minimum number of transcripts needed to define a cell too low quality to be considered for mapping
    RETURNS
        seg_comp_df: dataframe describing cell spatial locations (x and y), segmentation identifier, low quality cell identifier, with unique index
    """
    savepath = (
        savepath + bc + "_seg_comp_df_" + seg_name_a + "_and_" + seg_name_b + ".csv"
    )
    if Path(savepath).exists() and reuse_saved:
        seg_comp_df = pd.read_csv(savepath, index_col=0)
    else:
        seg_dfs = []
        for i in range(2):
            if i == 0:
                seg_h5ad = anndata_a
                name = seg_name_a
            else:
                seg_h5ad = anndata_b
                name = seg_name_b
            df_cols = ["center_x", "center_y", "source"]
            # sum across .X to get sum transcripts/cell > 40
            if "low_quality_cells" not in seg_h5ad.obs.columns.tolist():
                high_quality_cells = [
                    idx
                    for idx, x in enumerate(seg_h5ad.X.sum(axis=1))
                    if x >= min_transcripts
                ]
                seg_h5ad.obs["low_quality_cells"] = [True] * len(seg_h5ad.obs)
                seg_h5ad.obs.iloc[high_quality_cells, -1] = False
                df_cols.append("low_quality_cells")
            seg_h5ad.obs.loc[:, "source"] = name
            # save only necessary columns
            if "center_x" not in seg_h5ad.obs.columns.tolist():
                seg_h5ad.obs["center_x"] = seg_h5ad.obsm["spatial"][:, 0]
                seg_h5ad.obs["center_y"] = seg_h5ad.obsm["spatial"][:, 1]
            seg_df = seg_h5ad.obs[df_cols]
            seg_dfs.append(seg_df)
        seg_comp_df = pd.concat(seg_dfs, axis=0)
        seg_comp_df.index = seg_comp_df.index.astype(str)
    if save:
        seg_comp_df.to_csv(savepath)
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

    # get single nearest neighbor for each
    tree_a = cKDTree(dfa[["center_x", "center_y"]].copy())
    dists_a, inds_a = tree_a.query(
        dfb[["center_x", "center_y"]].copy(), 1
    )  # gives me index locs for dfa with neighbor dfb
    tree_b = cKDTree(dfb[["center_x", "center_y"]].copy())
    dists_b, inds_b = tree_b.query(
        dfa[["center_x", "center_y"]].copy(), 1
    )  # vice versa

    # get mutually matching pairs and save (index of seg_comp_df is str)
    match_to_dfb = pd.DataFrame(
        data=dfa.iloc[inds_a].index.values.tolist(),
        index=pd.Index(dfb.index, name="match_dfb_index"),
        columns=["match_dfa_index"],
    )
    match_to_dfb["same_cell"] = np.where(dists_a <= nn_dist, True, False)
    match_to_dfa = pd.DataFrame(
        data=dfa.index,
        index=pd.Index(dfb.iloc[inds_b].index.values.tolist(), name="match_dfb_index"),
        columns=["match_dfa_index"],
    )
    match_to_dfa["same_cell"] = np.where(dists_b <= nn_dist, True, False)
    mutual_matches = pd.merge(match_to_dfa, match_to_dfb, how="outer")
    mutual_matches = mutual_matches.set_index("match_dfa_index").join(
        match_to_dfa.reset_index().set_index("match_dfa_index"),
        how="left",
        rsuffix="_match",
    )
    mutual_match_dict = (
        mutual_matches[mutual_matches["same_cell"] == True]
        .drop(["same_cell", "same_cell_match"], axis=1)
        .to_dict()
    )
    inv_mutual_match_dict = {
        v: k for k, v in mutual_match_dict["match_dfb_index"].items()
    }
    # added union operwtor to dictionaries in 2020 for 3.9+ (pep 584), discussed https://stackoverflow.com/questions/38987/how-do-i-merge-two-dictionaries-in-a-single-expression-in-python
    mutual_matches_stacked = (
        mutual_match_dict["match_dfb_index"] | inv_mutual_match_dict
    )
    return mutual_matches_stacked


def create_node_df_sankey(
    seg_comp_df,
    barcode,
    save=True,
    savepath="/allen/programs/celltypes/workgroups/hct/emilyg/reseg_project/new_seg/",
):
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

    color_dict = {
        seg_comp_df.source.unique().tolist()[0]: "red",
        seg_comp_df.source.unique().tolist()[1]: "blue",
    }
    nodes_df = pd.DataFrame(columns=["Label", "Color", "Level", "Source", "Value"])
    unknown_unmatched_cells = {}

    for source, g in seg_comp_df.groupby("source"):
        new_rows = []
        low_q_and_match = g.groupby("low_quality_cells").agg(
            "count"
        )  # gives me low q t/f counts, and matched.
        total_cells = low_q_and_match.iloc[:, 0].sum()
        nodes_df.loc[len(nodes_df)] = {
            "Label": "total <br>" + str(total_cells),
            "Color": color_dict[source],
            "Level": 0,
            "Source": source,
            "Value": total_cells,
        }
        # low and normal quality cells
        nodes_df.loc[len(nodes_df)] = {
            "Label": "low quality cells <br>" + str(low_q_and_match.iloc[1, 1]),
            "Color": color_dict[source],
            "Level": 1,
            "Source": source,
            "Value": low_q_and_match.iloc[1, 1],
        }
        nodes_df.loc[len(nodes_df)] = {
            "Label": "normal quality cells <br>" + str(low_q_and_match.iloc[0, 1]),
            "Color": color_dict[source],
            "Level": 1,
            "Source": source,
            "Value": low_q_and_match.iloc[0, 1],
        }
        # matched and unmatched cells
        # because row 1 in agg is false only! so false x low (normal) and false x match (unmatch!)
        matched_cells = low_q_and_match.iloc[0, 1] - low_q_and_match.iloc[0, 3]
        nodes_df.loc[len(nodes_df)] = {
            "Label": "matched cells <br>" + str(matched_cells),
            "Color": color_dict[source],
            "Level": 2,
            "Source": source,
            "Value": matched_cells,
        }
        nodes_df.loc[len(nodes_df)] = {
            "Label": "unmatched cells <br>" + str(low_q_and_match.iloc[0, 3]),
            "Color": color_dict[source],
            "Level": 2,
            "Source": source,
            "Value": low_q_and_match.iloc[0, 3],
        }
        # raise flag if too few cells matched, may indicate scaling issue
        if matched_cells <= (0.1 * len(g)):
            raise ValueError(
                "The number of matched cells is less than 1% of total cells present. Please check your inputs for potential scaling issues."
            )

        rem_col_names = [
            x for x in g.columns[5:] if source + "_unfilt" not in x
        ]  # get remaining columns
        # pulling this to access unknown unmatched cells easily later
        high_q_df = g[g[g.columns.tolist()[3]] == False]
        unmatched_df = high_q_df[high_q_df[g.columns.tolist()[4]].isna() == False]
        rem_col_df = unmatched_df.loc[:, rem_col_names]
        rem_col_counts = rem_col_df.isna().value_counts().reset_index()
        name = " ".join(rem_col_names[0].split("_")) + "<br>"
        fcounts = rem_col_counts[rem_col_counts[rem_col_names[0]] == False][
            "count"
        ].values.tolist()[0]
        nodes_df.loc[len(nodes_df)] = {
            "Label": name + str(fcounts),
            "Color": color_dict[source],
            "Level": 3,
            "Source": source,
            "Value": fcounts,
        }
        # remaining unmatched cells (known or unknown)
        if len(rem_col_counts) > 1:
            rem_unmatched_cells = rem_col_counts[
                rem_col_counts[rem_col_names[0]] == True
            ]["count"].values.tolist()[0]
            unknown_unmatched_cells[source] = rem_col_df[
                rem_col_df[rem_col_names[0]].isna() == True
            ].index.values.tolist()
            name = " ".join(rem_col_names[-1].split("_"))
            nodes_df.loc[len(nodes_df)] = {
                "Label": name + str(rem_unmatched_cells),
                "Color": color_dict[source],
                "Level": 3,
                "Source": source,
                "Value": rem_unmatched_cells,
            }
    if save:
        nodes_df.to_csv(savepath + "/sankey_nodes_df" + barcode + ".csv")
    return nodes_df, unknown_unmatched_cells


def create_link_df_sankey(
    nodes_df,
    barcode,
    save=True,
    savepath="/allen/programs/celltypes/workgroups/hct/emilyg/reseg_project/new_seg/",
):
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

    link_color_dict = {"red": "rgb(205, 209, 228)", "blue": "rgb(205, 209, 228)"}
    links_df = pd.DataFrame(columns=["Source", "Target", "Value", "Link Color"])
    for i, row in nodes_df.iterrows():
        if "total" in row.Label:
            continue  # skip 0 since all comes from zero
        else:
            target = i
            value = row["Value"]
            source = nodes_df[
                (nodes_df["Level"] == row["Level"] - 1)
                & (nodes_df["Source"] == row["Source"])
            ].index.values.tolist()
            if len(source) > 1 and row["Level"] != 3:
                source = [
                    i
                    for i, x in nodes_df.iloc[source, :].iterrows()
                    if "un" not in x["Label"]
                    if "low" not in x["Label"]
                ][0]
            elif len(source) > 1 and row["Level"] == 3:
                source = [
                    i
                    for i, x in nodes_df.iloc[source, :].iterrows()
                    if "un" in x["Label"]
                ][0]
            else:
                source = source[0]
            prev_row = nodes_df.iloc[source, :]
            link_color = link_color_dict[prev_row["Color"]]
            links_df.loc[len(links_df)] = {
                "Source": source,
                "Target": target,
                "Value": value,
                "Link Color": link_color,
            }

    # special case: matched cells should have same name, new source, new color
    matched_cell_rows = nodes_df[
        (nodes_df["Label"].str.contains("matched cells"))
        & (~nodes_df["Label"].str.contains("un"))
    ]
    matched_cell_val = matched_cell_rows["Value"].values.tolist()[0]
    nodes_df.loc[len(nodes_df)] = {
        "Label": "matched cells <br> " + str(matched_cell_val),
        "Color": "violet",
        "Level": matched_cell_rows["Level"].values.tolist()[0],
        "Source": "Both",
        "Value": matched_cell_val,
    }
    # find matched col for both in nodes_df, then find rows with those IDS in links df and connect them
    update_idxs = links_df[
        links_df["Target"].isin(matched_cell_rows.index.values.tolist())
    ].index.values.tolist()
    links_df.loc[update_idxs, "Target"] = [
        nodes_df[nodes_df["Source"] == "Both"].index.values
    ] * len(update_idxs)
    if save:
        links_df.to_csv(savepath + "/sankey_links_df.csv")
    return links_df


def spatial_detection_scores(
    reference: pd.DataFrame,
    query: pd.DataFrame,
    plot_stuff=True,
    query_name: str = "query data",
    comparison_column="transcript_counts",
    category="supercluster_name",
    n_bins=50,
    in_place=True,
    non_spatial=False,
):
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
        raise ValueError("category " + category + " not in reference and query inputs")

    shared_category_values = list(
        set(reference[category].unique()) & set(query[category].unique())
    )
    if (
        len(shared_category_values) < query[category].unique().shape[0]
        or len(shared_category_values) < reference[category].unique().shape[0]
    ):
        in_place = False

    if in_place:
        s2 = query.loc[query[category].isin(shared_category_values), :]
        s1 = reference.loc[reference[category].isin(shared_category_values), :]
    else:
        s2 = query.loc[query[category].isin(shared_category_values), :].copy()
        s1 = reference.loc[reference[category].isin(shared_category_values), :].copy()

    means = s1.groupby(category, observed=True)[comparison_column].mean()
    stds = s1.groupby(category, observed=True)[comparison_column].std()

    # if you want to compare to an ungrouped mean/std, try this:
    if non_spatial:
        means[:] = means.mean()
        stds[:] = stds.mean()

    s2["detection_relative_z_score"] = 0.0
    s2["detection_difference"] = 0.0
    s2["detection_ratio"] = 0.0

    for c, gb in s2.groupby(category, observed=True):
        if c not in shared_category_values:
            continue

        s2.loc[s2[category] == c, ["detection_relative_z_score"]] = (
            (s2.loc[s2[category] == c, [comparison_column]] - means[c]) / stds[c]
        ).values
        s2.loc[s2[category] == c, ["detection_difference"]] = (
            s2.loc[s2[category] == c, [comparison_column]] - means[c]
        ).values
        s2.loc[s2[category] == c, ["log_10_detection_ratio"]] = np.log10(
            (s2.loc[s2[category] == c, [comparison_column]] / means[c]).values
        )

    s2["xy_bucket"] = list(
        zip(
            pd.cut(s2.x_centroid, n_bins, labels=list(range(n_bins))),
            pd.cut(s2.y_centroid, n_bins, labels=list(range(n_bins))),
        )
    )

    binx = s2.groupby("xy_bucket").x_centroid.mean()
    biny = s2.groupby("xy_bucket").y_centroid.mean()

    z_score = s2.groupby("xy_bucket").detection_relative_z_score.mean()
    difference = s2.groupby("xy_bucket").detection_difference.mean()
    log_ratio = s2.groupby("xy_bucket").log_10_detection_ratio.mean()
    n_cells = s2.groupby("xy_bucket").x_centroid.count()

    bin_image_z_score = np.zeros([n_bins, n_bins])
    bin_image_difference = np.zeros([n_bins, n_bins])
    bin_image_ratio = np.zeros([n_bins, n_bins])
    bin_image_counts = np.zeros([n_bins, n_bins])

    extent = [np.min(binx), np.max(binx), np.min(biny), np.max(biny)]
    for coord in binx.index:
        bin_image_z_score[coord[1], coord[0]] = z_score[coord]
        bin_image_difference[coord[1], coord[0]] = difference[coord]
        bin_image_ratio[coord[1], coord[0]] = log_ratio[coord]
        bin_image_counts[coord[1], coord[0]] = n_cells[coord]

    if plot_stuff:
        if non_spatial:
            title_string = "Non-spatial Detection Scores"
        else:
            title_string = "Spatial Detection Scores"
        min_maxes = {
            "detection z-score": [bin_image_z_score, [-1, 1]],
            "total counts difference": [bin_image_difference, [-100, 100]],
            "log10(detection ratio)": [bin_image_ratio, [-1, 1]],
        }

        fig, axs = plt.subplots(1, 3, figsize=[15, 5])
        if non_spatial:
            fig.suptitle(title_string + "\n" + query_name)
        else:
            fig.suptitle(
                title_string
                + "\n"
                + query_name
                + " grouped by "
                + category
                + " and spatially binned"
            )
        for ii, plot_name in enumerate(min_maxes.keys()):
            ax = axs[ii]
            pcm = ax.imshow(
                min_maxes[plot_name][0],
                extent=extent,
                cmap="coolwarm_r",
                vmin=min_maxes[plot_name][1][0],
                vmax=min_maxes[plot_name][1][1],
            )
            fig.colorbar(pcm, ax=ax, shrink=0.7)
            ax.set_title(query_name + "\n" + plot_name)

    if in_place:

        return dict(
            z_score_image=bin_image_z_score,
            difference_image=bin_image_difference,
            ratio_image=bin_image_ratio,
            extent=extent,
            count_image=bin_image_counts,
            query=True,
            reference=True,
        )

    else:
        return dict(
            z_score_image=bin_image_z_score,
            difference_image=bin_image_difference,
            ratio_image=bin_image_ratio,
            extent=extent,
            count_image=bin_image_counts,
            query=s2,
            reference=s1,
        )


def summarize_and_plot(
    spatial_density_results,
    min_cells_per_bin=50,
    z_score_limit=-0.5,
    area_frac_limit=0.1,
    figsize=[15, 15],
    plot_columns=9,
    plot_stuff=True,
    title_mapping={},
):
    """
    Summarizes and plots spatial density results.
    Parameters:
    - spatial_density_results (dict): A dictionary containing spatial density results.
    - min_cells_per_bin (int): The minimum number of cells per bin. Default is 50.
    - z_score_limit (float): The z-score limit. Default is -0.5.
    - area_frac_limit (float): The area fraction limit. Default is 0.1.
    - figsize (list): The figure size. Default is [15, 15].
    - plot_columns (int): The number of plot columns in figure subplot. Default is 9.
    - plot_stuff (bool): Whether to plot the results. Default is True.
    - title_mapping (dict): A dictionary mapping spatial_density_results keys to titles for the plots. Default is an empty dictionary.
    Returns:
    - results (list): A list of dictionaries containing the summarized results. Each dictionary contains the following keys:
        - key (float): same as original input key
        - bad_frac (float): The fraction of tissue bins with z-scores below the z-score limit.
        - area_frac_limit (float): The area fraction limit.
        - failed (bool): Whether the fraction of tissue bins with z-scores below the z-score limit is greater than or equal to the area fraction limit.
        - mean (float): The mean of the non-zero z-scores.
        - stdev (float): The standard deviation of the non-zero z-scores.
    """

    fails_vs_rnaseq = []
    results = []
    if plot_stuff:
        plt.figure(figsize=figsize)

    for ii, z in enumerate(sorted(list(spatial_density_results.keys()))):
        # estimate tissue bins:
        tissue_bins = spatial_density_results[z]["count_image"] > min_cells_per_bin
        bad_frac = np.sum(
            spatial_density_results[z]["z_score_image"][tissue_bins] <= z_score_limit
        ) / np.sum(tissue_bins)
        image_to_show = spatial_density_results[z]["z_score_image"].copy()
        image_to_show[np.logical_not(tissue_bins)] = 0

        if plot_stuff:
            if z in title_mapping:
                title_prefix = title_mapping[z]
            else:
                title_prefix = ""

            plt.subplot(
                1 + len(list(spatial_density_results.keys())) // plot_columns,
                plot_columns,
                ii + 1,
            )

            ax = plt.imshow(
                image_to_show,
                vmin=-1,
                vmax=1,
                extent=spatial_density_results[z]["extent"],
                cmap="coolwarm_r",
            )

            if bad_frac >= area_frac_limit:
                plt.title(title_prefix + " Fail", fontdict={"size": 12})
            else:
                plt.title(title_prefix, fontdict={"size": 12})
            plt.xticks([])
            plt.yticks([])

        if bad_frac >= area_frac_limit:
            fails_vs_rnaseq.append(z)

        results.append(
            dict(
                key=z,
                bad_frac=bad_frac,
                area_frac_limit=area_frac_limit,
                failed=bad_frac >= area_frac_limit,
                mean=np.mean(image_to_show[image_to_show != 0]),
                stdev=np.std(image_to_show[image_to_show != 0]),
            )
        )
    return results


def compare_reference_and_spatial(
    reference_anndata,
    spatial_anndata,
    category="MTG_subclass_name",
    layer_field=None,
    plot_stuff=True,
    target_obs_key="comparison_transcript_counts",
    ok_to_clobber=False,
):
    """
    Compare reference and spatial data based on a specified category.
    Parameters:
    - reference_anndata (AnnData): Annotated data matrix containing the reference data.
    - spatial_anndata (AnnData): Annotated data matrix containing the spatial data.
    - category (str): The category to compare the data on. Default is "MTG_subclass_name".
    - layer_field (str): The field to use for layer-based comparison. Default is None.
    - plot_stuff (bool): Whether to plot the comparison results. Default is True.
    - target_obs_key (str): The key to use for storing the comparison transcript counts. Default is "comparison_transcript_counts".
    - ok_to_clobber (bool): Whether it is okay to overwrite the target_obs_key in the reference_anndata. Default is False.
    Raises:
    - ValueError: If the target_obs_key already exists in the reference_anndata and ok_to_clobber is False.
    Returns:
    - None
    """
    if target_obs_key in reference_anndata.obs.columns:
        if not ok_to_clobber == True:
            raise ValueError(
                "obs key "
                + target_obs_key
                + " is already in the input reference .obs\n If desired, set ok_to_clobber to True"
            )
        else:
            print(
                "warning: modifying input reference anndata .obs field "
                + target_obs_key
            )

    if layer_field is not None:
        reference_anndata.obs[target_obs_key] = np.sum(
            reference_anndata.layers[layer_field], axis=1
        )
    else:
        reference_anndata.obs[target_obs_key] = np.sum(reference_anndata.X, axis=1)

    means = (
        reference_anndata.obs.loc[:, [category, target_obs_key]]
        .groupby(category, observed=True)
        .mean()
    )
    means["spatial_counts"] = 0.0

    for name in spatial_anndata.obs[category].unique():
        means.loc[name, ["spatial_counts"]] = (
            spatial_anndata.obs.loc[
                spatial_anndata.obs[category] == name, [target_obs_key]
            ]
            .mean()
            .values[0]
        )

    fit_values = np.polyfit(means.spatial_counts, means[target_obs_key], 1)

    scale_factor = 1 / fit_values[0]
    if plot_stuff:
        plt.figure()
        plt.plot(means.spatial_counts, means[target_obs_key], ".")
        plt.xlabel("spatial counts")
        plt.ylabel(target_obs_key + " in reference anndata")
        plt.figure()
        sns.regplot(x=means.spatial_counts, y=scale_factor * means[target_obs_key])
        plt.plot([0, 1000], [0, 1000], label="unity")
        # plt.plot(means.spatial_counts, means.spatial_counts*fit_values[0]+fit_values[1], label = "fit")
        plt.ylabel("scaled reference data")
        plt.title("\n scale = " + str(scale_factor)[:6])
        plt.legend()

    # scale transcript counts based on the linear fit.
    # could also do this per group.

    reference_anndata.obs[target_obs_key] = (
        scale_factor * reference_anndata.obs[target_obs_key]
    )
