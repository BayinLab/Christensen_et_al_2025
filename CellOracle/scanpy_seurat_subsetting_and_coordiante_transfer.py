#import required libraries
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# Load Scanpy object
adata = sc.read_h5ad('all_cells_scanpy_top_4000_genes.h5ad')

# Load data extracted from Seurat object
seurat_data = pd.read_csv("cells_to_keep_and_annotations_UMAP.csv")

# Ensure barcodes are properly formatted
seurat_data['cell_barcodes'] = seurat_data['cell_barcodes'].str.replace('-\\d+$', '', regex=True)
adata.obs_names = adata.obs_names.str.replace('-\\d+$', '', regex=True)

# Subset Scanpy object to include only cells present in Seurat data
common_barcodes = adata.obs_names.intersection(seurat_data['cell_barcodes'])
print(f"Number of common barcodes: {len(common_barcodes)}")

adata_subset = adata[common_barcodes].copy()

# Merge the Seurat data with the Scanpy subset
merged_data = pd.merge(adata_subset.obs, seurat_data, left_index=True, right_on='cell_barcodes', how='left')
adata_subset.obs['annotations'] = merged_data['annotations'].values

# Inspect columns to identify UMAP coordinate names
print("Columns in merged_data:")
print(merged_data.columns)

# Update UMAP column names
umap_1 = 'umap_1'  
umap_2 = 'umap_2'  

# Check if UMAP columns exist in merged_data
if umap_1 in merged_data.columns and umap_2 in merged_data.columns:
    umap_coords = merged_data[[umap_1, umap_2]].values
    adata_subset.obsm['X_umap'] = umap_coords
else:
    print("UMAP columns not found in merged_data. Please update umap_1 and umap_2 with correct names.")

# Ensure annotations are categorical
adata_subset.obs['annotations'] = adata_subset.obs['annotations'].astype('category')

# Use the same colours that were previously used
colors = ["#e87d72", "#bd9c33", "#6db234", "#56bd97", "#51b3e6", "#a08bf8", "#e96bd2"]
adata_subset.uns['annotations_colors'] = colors

# Save the updated Scanpy object
adata_subset.write('Jens_subset_oracle_annotations_colours_umap_4000.h5ad')

# Plot to verify
sc.pl.umap(adata_subset, color="annotations")
plt.savefig("Jens_subset_oracle_annotations_colours_umap_4000.png")


# Optional, print UMAP coordinates for comparison
print("UMAP coordinates from Seurat data (first 5 cells):")
print(merged_data[[umap_1, umap_2]].head())

print("UMAP coordinates from Scanpy object (first 5 cells):")
print(adata_subset.obsm['X_umap'][:5])

