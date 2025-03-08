import os
import sys
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co
import importlib

# import packages and print progress
def import_packages(packages):
    for i, package in enumerate(packages, start=1):
        globals()[package.split('.')[-1]] = importlib.import_module(package)

def main():
    # Print loading message
    print("Loading packages, this may take a few moments!")

    # List of packages to load
    packages = [
        'os', 'sys', 'matplotlib.colors', 'matplotlib.pyplot', 'numpy', 'pandas',
        'scanpy', 'seaborn', 'celloracle'
    ]
    
    # Import packages
    import_packages(packages)
    
    # Load Oracle and Links data produced earleir
    oracle = celloracle.load_hdf5("updated_Jens_subset_4000_with_psudotime.celloracle.oracle")
    links = celloracle.load_hdf5("updated_Jens_subset_4000_filtered.celloracle.links")
    links.filter_links()
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

    # Prompt user for "TF of interest"
    goi = input("Please enter the gene of interest (e.g., Sox9): ")

    try:
        # Neighbors and graph
        scanpy.pp.neighbors(oracle.adata)
        scanpy.tl.draw_graph(oracle.adata)
        scanpy.pl.draw_graph(oracle.adata, color=[goi, oracle.cluster_column_name], layer="imputed_count", use_raw=False, cmap="viridis")

        # Plot histogram
        scanpy.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").hist()
        plt.savefig(f"{goi}_expression_hist.pdf", format="pdf")

        # Get maximum expression level of the GOI
        max_expression = scanpy.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").max().values[0]
        print(f"The maximum expression level of {goi} is {max_expression}")

        # Prompt user for the type of simulation
        sim_type = input("Select simulation type (KO, overexpress, 2x max, 0.5x max): ").strip().lower()

        if sim_type == "ko":
            perturb_value = 0.0
        elif sim_type == "overexpress":
            perturb_value = max_expression
        elif sim_type == "2x max":
            perturb_value = 2 * max_expression
        elif sim_type == "0.5x max":
            perturb_value = 0.5 * max_expression
        else:
            print("Invalid selection. Exiting.")
            return

        oracle.simulate_shift(perturb_condition={goi: perturb_value}, n_propagation=3)

        oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
        oracle.calculate_embedding_shift(sigma_corr=0.05)
        
        n_grid = 40
        oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
        
        # Prompt user for min_mass
        min_mass = float(input("Please enter the minimum mass value (e.g., 2.6): "))
        oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

        fig, ax = plt.subplots(figsize=[6,6])
        scanpy.pl.embedding(adata=oracle.adata, basis=oracle.embedding_name, ax=ax, cmap="rainbow", color=["Pseudotime"])

        from celloracle.applications import Gradient_calculator
        gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
        
        gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
        gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
        gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly": 3}, plot=True)
        gradient.calculate_gradient()
        
        scale_dev = 40
        gradient.visualize_results(scale=scale_dev, s=5)
        
        fig, ax = plt.subplots(figsize=[6, 6])
        gradient.plot_dev_flow_on_grid(scale=scale_dev, ax=ax)
        
        from celloracle.applications import Oracle_development_module
        dev = Oracle_development_module()
        dev.load_differentiation_reference_data(gradient_object=gradient)
        dev.load_perturb_simulation_data(oracle_object=oracle)
        dev.calculate_inner_product()
        dev.calculate_digitized_ip(n_bins=10)
        dev.visualize_development_module_layout_0(s=5, scale_for_simulation=20, s_grid=30, scale_for_pseudotime=20, vm=1)

        # Save the final plot as PDF
        effect = sim_type.replace(" ", "_")
        final_plot_filename = f"{goi}_{effect}.pdf"
        plt.savefig(final_plot_filename, format="pdf")
        print(f"Plot saved as {final_plot_filename}")
    except Exception as e:
        print(f"Gene option invalid, please check you have spelled it right (case-sensitive). Other reasons could be:\n"
              f"i. Gene is not a transcription factor\n"
              f"ii. Gene is not in the top variable genes\n"
              f"iii. Transcription factor has a small number of links in the GRNs\nError: {e}")

    # Print final message with package versions
    print(f"These results were produced using CellOracle ({celloracle.__version__}), Scanpy ({scanpy.__version__}), "
          f"and Python ({sys.version.split()[0]}), Bayin Lab 2024")

if __name__ == "__main__":
    main()
