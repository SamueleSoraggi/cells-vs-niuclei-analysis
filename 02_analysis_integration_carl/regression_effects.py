import scanpy as sc
import numpy as np
import pandas as pd


adata = sc.read('./write/adata.dpt.lamian.h5ad')

from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import cross_val_score
from scipy.sparse import issparse

# Prepare features for the model
def prepare_features(adata):
    """Create feature matrix for batch effect modeling"""
    
    features = pd.DataFrame(index=adata.obs_names)
    
    # Batch indicator (this is what we want to remove)
    features['is_wholecell'] = (adata.obs['batch'] == 'Cells').astype(int)
    
    # Biological covariates (preserve these)
    le_celltype = LabelEncoder()
    features['celltype_encoded'] = le_celltype.fit_transform(adata.obs['clusters_single'])
    features['pseudotime'] = adata.obs['dpt_palantir_cellalign'].values
    
    # Interaction terms (batch might have different effect at different stages)
    features['batch_x_pseudotime'] = features['is_wholecell'] * features['pseudotime']
    features['batch_x_celltype'] = features['is_wholecell'] * features['celltype_encoded']
    
    return features, le_celltype

features_df, celltype_encoder = prepare_features(adata)
print(f"Feature matrix shape: {features_df.shape}")
features_df.head()

def correct_gene_expression(adata, gene_idx, features_df, model_type='rf'):
    """
    Build ML model to remove batch effect for a single gene
    
    Strategy:
    1. Train model to predict expression from all features (including batch)
    2. Make predictions twice: once with actual batch, once with batch=0 (nuclei)
    3. Corrected expression = original - (pred_actual - pred_batch0)
    
    This removes the batch-specific component while preserving biology
    """
    
    # Get expression for this gene
    if issparse(adata.X):
        y = adata.X[:, gene_idx].toarray().flatten()
    else:
        y = adata.X[:, gene_idx].flatten()
    
    X = features_df.values
    
    # Choose model
    if model_type == 'rf':
        model = RandomForestRegressor(
            n_estimators=100, 
            max_depth=8,
            min_samples_leaf=10,
            random_state=42,
            n_jobs=1
        )
    elif model_type == 'gbm':
        model = GradientBoostingRegressor(
            n_estimators=100,
            max_depth=5,
            learning_rate=0.1,
            random_state=42
        )
    
    # Train model
    model.fit(X, y)
    
    # Predict with actual features
    pred_actual = model.predict(X)
    
    # Predict with batch effect removed (set all cells to "nuclei" reference)
    X_corrected = X.copy()
    X_corrected[:, 0] = 0  # is_wholecell = 0
    X_corrected[:, 3] = 0  # batch_x_pseudotime = 0
    X_corrected[:, 4] = 0  # batch_x_celltype = 0
    pred_nobatch = model.predict(X_corrected)
    
    # Correct expression by removing batch-specific component
    batch_effect = pred_actual - pred_nobatch
    y_corrected = y - batch_effect
    
    return y_corrected, model, batch_effect

# Test on one gene
test_gene_idx = 10  # Change index to test different genes
test_gene_name = adata.var_names[test_gene_idx]
y_corrected, model, batch_effect = correct_gene_expression(adata, test_gene_idx, features_df)

print(f"Gene: {test_gene_name}")
print(f"Batch effect range: [{batch_effect.min():.3f}, {batch_effect.max():.3f}]")
print(f"Mean batch effect in Cells: {batch_effect[features_df['is_wholecell']==1].mean():.3f}")
print(f"Mean batch effect in Nuclei: {batch_effect[features_df['is_wholecell']==0].mean():.3f}")


from tqdm import tqdm
import multiprocessing as mp
from functools import partial

def correct_gene_wrapper(gene_idx, adata, features_df):
    """Wrapper for parallel processing"""
    try:
        y_corrected, _, _ = correct_gene_expression(adata, gene_idx, features_df, model_type='gbm')
        return gene_idx, y_corrected
    except Exception as e:
        print(f"Error processing gene {gene_idx}: {e}")
        return gene_idx, None

# Create corrected data matrix
n_genes = adata.n_vars
n_cells = adata.n_obs

print("Correcting batch effects for all genes using ML models...")
print(f"Processing {n_genes} genes across {n_cells} cells")
# Standard loop with progress bar (more compatible)
X_corrected = np.zeros((n_cells, n_genes))

for gene_idx in tqdm(range(n_genes)):
    try:
        y_corrected, _, _ = correct_gene_expression(adata, gene_idx, features_df, model_type='rf')
        X_corrected[:, gene_idx] = y_corrected
    except Exception as e:
        print(f"Error processing gene {gene_idx} ({adata.var_names[gene_idx]}): {e}")
        # Keep original expression if correction fails
        if issparse(adata.X):
            X_corrected[:, gene_idx] = adata.X[:, gene_idx].toarray().flatten()
        else:
            X_corrected[:, gene_idx] = adata.X[:, gene_idx].flatten()

print(f"Corrected expression matrix shape: {X_corrected.shape}")

# Store corrected data in new AnnData object
adata_corrected = adata.copy()
adata_corrected.layers['ml_corrected'] = X_corrected.copy()
adata_corrected.X = X_corrected.copy()

print("Batch-corrected data stored in adata_corrected.X and adata_corrected.layers['ml_corrected']")

adata_corrected.write('./write/adata.gbm_correction.h5ad')
