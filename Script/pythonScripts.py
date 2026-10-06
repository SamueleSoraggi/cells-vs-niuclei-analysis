### SDO outlyingness for scater PCA-based filtering (adjOutlyingness of the R robust package used in scater for scrna filtering)
###  (adjOutlyingness.R of the R `robust` package used in `scater` for scrna filtering in R)

import numpy as np
from robustbase import mad
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import scanpy as sc
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def directions(Xpca): #directions for points projection
    n, p = Xpca.shape
    P = Xpca[ np.random.choice(n, replace=False, size=p), :]
    G = np.zeros(p)
    E = np.ones(p)
    qrP, _ = np.linalg.qr(P)
    G2 = np.linalg.solve(qrP, E)
    return G2

def outlyingness(Xpca, N):
    Nx, p = Xpca.shape
    a = np.zeros([100,p])
    for i in range(N):
        a[i,:] = directions(Xpca)  
    all_values = np.zeros([Nx,N])
    for i in range(N):
        w = np.matmul(Xpca, a[i,:])
        all_values[:,i] = abs(w-np.median(w))/mad(w)
    return all_values

def PCAfiltering(adata, obs_subset=None, p=2, N=100, plot=True, random_seed=42):

    np.random.seed(42)
    print(f'--- PCA of dimension {p} on the following metadata:')
    if obs_subset==None:
        X = adata.obs
        print(list(adata.obs.columns))
    else:
        X = adata.obs[list(obs_subset)]
        print(obs_subset)
    print(f'on a sample with {X.shape[0]} datapoints')
    X2 = scale(X)
    Xpca = PCA(n_components=p).fit_transform(X2)
    
    print(f'--- Calculate outlyingness on {N} randomly sampled axae, (random seed {random_seed})')
    outlying_values = outlyingness(Xpca, N)
    sup_outlyingness = np.max(outlying_values, 1)
    mad_SDO = mad(sup_outlyingness)
    median_SDO = np.median(sup_outlyingness)
    outliers = sup_outlyingness > median_SDO + 5*mad_SDO
    adata.obs['SDO_outliers'] = pd.Categorical(outliers)
    adata.obs['SDO_outlyingness'] = sup_outlyingness
    
    print(f'--- Found {np.sum(outliers)} outliers')
    print('--- Returned annotated data object containing')
    print('--- * adata.obs["SDO_outliers"]: boolean variable identifying outliers')
    print('--- * adata.obs["SDO_outlyingness"]: outlyingness of each cell')
    
    if(plot):
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20,4), gridspec_kw={'wspace':0.9})
        x1 = sns.scatterplot(x=Xpca[:,0], y=Xpca[:,1], hue=sup_outlyingness, ax=ax1)
        x1.set_title('Outlyingness on\nmetadata PCA')
        x2 = sns.scatterplot(x=Xpca[:,0], y=Xpca[:,1], hue=outliers, ax=ax2)
        x2.set_title('Outliers on\nmetadata PCA')
        x3 = sns.scatterplot(x='total_counts', y='n_genes_by_counts', data=adata.obs, hue=outliers, ax=ax3)
        x=x3.set_title('Outliers on scatterplot\n#Transcripts vs #Genes')

    return adata

###PCA detection of highly dependent technical features

def dependentFeatures(adata, obs_subset=['total_counts'], n_pcs = 10):
    
    from sklearn.linear_model import LinearRegression as LR
    from sklearn.metrics import mean_squared_error, r2_score
    from matplotlib.lines import Line2D  
    
    L = len(obs_subset)
    X = pd.DataFrame(adata.obsm['X_pca'][:,:n_pcs], columns=[f'PC{i}' for i in range(n_pcs)])
    Y = adata.obs[obs_subset]
    scores = np.zeros([n_pcs, L])
    
    for n in range(n_pcs):
        for l in range(L):
            #reg = LR().fit(X[f'PC{n}'].reshape(-1,1), Y[obs_subset[l]].reshape(-1,1))
            reg = LR().fit(np.array(X[f'PC{n}']).reshape(-1,1), np.array(Y[obs_subset[l]]).reshape(-1,1))
            pred = reg.predict(np.array(X[f'PC{n}']).reshape(-1,1))          
            scores[n, l] = r2_score(Y[obs_subset[l]], pred)
            
    max_scores = np.max(scores, 0)
    argmax_scores = np.argmax(scores, 0)
 
    fig, ax = plt.subplots(L, 1, figsize=(6,6*L), gridspec_kw={'hspace':0.4})
    for l in range(L):
        x1 = sns.regplot(x=X[f'PC{argmax_scores[l]}'], y=Y[obs_subset[l]], ax=ax[l], scatter_kws={'alpha':0.3})
        x = x1.set_title(f'Best lin.regression:\n {obs_subset[l]} VS PC {argmax_scores[l]}\nR score {max_scores[l]}')
        x1.set(xlabel=f'{obs_subset[l]}', ylabel=f'PC {argmax_scores[l]}')
        
###calculate markers' scores from a dictionary
def marker_score(markers_dict, adata, N_samples=100, random_seed=42):
    np.random.seed(random_seed)
    markers_list = []
    N_genes = adata.shape[1]
    random_genes = np.unique( np.random.randint(low=0, high=N_genes, size=N_samples) )
    gene_names = adata.var_names[random_genes]
    for i in markers_dict:
        markers_list.append(f'{i}_score')
        adata.obs[f'{i}_score'] = np.array( np.mean(adata[:,markers_dict[i]].X,1) - np.mean(adata[:,gene_names].X,(0,1)) )
    return markers_list, adata

###function to rename clusters from a dictionary
def rename_clusters(names_dict, names_obs):
    clusters = pd.Categorical(names_obs)
    clusters=clusters.rename_categories(names_dict)
    cluster_array = np.array(clusters)
    split_array = [ i.split('.')[0] for i in cluster_array ]
    clusters = pd.Categorical(split_array)
    return clusters

###use scores instead of manual names (as above) to rename clusters
def clustersByScores(adata, markers_scores, leidenClusters):
    clusters = pd.Categorical(leidenClusters)
    scoresTable = adata.obs[markers_scores]
    clusterUnique = np.unique(leidenClusters)
    newNames = pd.Series(index=leidenClusters)
    for CLST in clusterUnique:
        meanScores = np.mean( scoresTable.loc[leidenClusters==CLST,:], 0)
        newId = meanScores.index[ np.argmax(meanScores) ].split('_')[0]
        newNames[CLST] = newId
    return(pd.Categorical(newNames))


def plot_genes_single(adata, 
                      N_grid=100, #bins 
                      GENE='MKI67', 
                      layer='MAGIC_counts', 
                      dataset_name='CELLS', 
                      plot_size = (10,4)):
    
    plt.rcParams['figure.figsize'] = plot_size
    
    from pygam import LinearGAM, s

    gam = LinearGAM()

    y = adata[:,GENE].layers[layer]
    
    y[y<0]=0
    X = np.empty( (len(y),1), dtype=object)
    X[:,0] = adata.obs['dpt_cellalign'] #pseudotimes

    y = y[adata.obs['dpt_cellalign']>=0]
    X = X[adata.obs['dpt_cellalign']>=0,]

    gam.gridsearch(X , y)

    XX = gam.generate_X_grid(term=0, n=N_grid)

    f=plt.figure();
    fig, axs = plt.subplots(1,1);
    f=axs.hlines(xmin=0, xmax=1, y=0, colors='Black', alpha=.25)
    f=sns.lineplot(x=XX[:, 0], y=gam.partial_dependence(term=0, X=XX), ax=axs, c='Red', linewidth=3 )
    f=sns.lineplot(x=XX[:, 0], y=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,0], ax=axs, c='Blue', alpha=.5)
    f=sns.lineplot(x=XX[:, 0], y=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,1], ax=axs, c='Blue', alpha=.5)

    f=axs.fill_between(x=XX[:, 0], y1=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,0], y2=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,1], alpha=.5)
    f=axs.set_title(f'Pseudotime curve for gene {GENE} - {dataset_name} data')

    


def plot_genes_sidebyside(adata_cells, 
                          adata_nuclei, 
                          N_grid=100,
                          GENE='MKI67', 
                          layer='MAGIC_counts', 
                          dataset_name='CELLS', 
                          plot_size=(12,4)):
    
    plt.rcParams['figure.figsize'] = plot_size
    
    from pygam import LinearGAM, s

    gam = LinearGAM()

    y = adata_cells[:,GENE].layers[layer]
    y[y<0]=0
    X = np.empty( (len(y),1), dtype=object)
    X[:,0] = adata_cells.obs['dpt_cellalign'] #pseudotimes
    y = y[adata_cells.obs['dpt_cellalign']>=0]
    X = X[adata_cells.obs['dpt_cellalign']>=0,]

    gam.gridsearch(X , y, )

    XX = gam.generate_X_grid(term=0, n=N_grid)

    f=plt.figure();
    fig, (ax1, ax2) = plt.subplots(1,2);
    
    f=ax1.hlines(xmin=0, xmax=1, y=0, colors='Black', alpha=.25)
    f=sns.lineplot(x=XX[:, 0], y=gam.partial_dependence(term=0, X=XX), ax=ax1, c='Red', linewidth=3 )
    f=sns.lineplot(x=XX[:, 0], y=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,0], ax=ax1, c='Blue', alpha=.5)
    f=sns.lineplot(x=XX[:, 0], y=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,1], ax=ax1, c='Blue', alpha=.5)

    f=ax1.fill_between(x=XX[:, 0], y1=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,0], y2=gam.partial_dependence(term=0, X=XX, width=.95)[1][:,1], alpha=.5)
    f=ax1.set_title(f'Pseudotime curve for gene {GENE} - CELLS data')


    gam2 = LinearGAM()

    y = adata_nuclei[:,GENE].layers[layer]
    y[y<0]=0
    X = np.empty( (len(y),1), dtype=object)
    X[:,0] = adata_nuclei.obs['dpt_cellalign']
    y = y[adata_nuclei.obs['dpt_cellalign']>=0]
    X = X[adata_nuclei.obs['dpt_cellalign']>=0,]

    gam2.gridsearch(X , y, )

    XX = gam2.generate_X_grid(term=0, n=100)

    f=ax2.hlines(xmin=0, xmax=1, y=0, colors='Black', alpha=.25)
    f=sns.lineplot(x=XX[:, 0], y=gam2.partial_dependence(term=0, X=XX), ax=ax2, c='Red', linewidth=3 )
    f=sns.lineplot(x=XX[:, 0], y=gam2.partial_dependence(term=0, X=XX, width=.95)[1][:,0], ax=ax2, c='Blue', alpha=.5)
    f=sns.lineplot(x=XX[:, 0], y=gam2.partial_dependence(term=0, X=XX, width=.95)[1][:,1], ax=ax2, c='Blue', alpha=.5)

    f=ax2.fill_between(x=XX[:, 0], y1=gam2.partial_dependence(term=0, X=XX, width=.95)[1][:,0], y2=gam2.partial_dependence(term=0, X=XX, width=.95)[1][:,1], alpha=.5)
    f=ax2.set_title(f'Pseudotime curve for gene {GENE} - NUCLEI data')   


def plot_genes_overlap(adata_cells, 
                       adata_nuclei, 
                       N_grid=100, 
                       GENE='MKI67', 
                       layer='MAGIC_counts', 
                       cluster_key='spermatogenesis_single', 
                       plot_size=(12,4)):
    
    plt.rcParams['figure.figsize'] = plot_size
    
    from pygam import LinearGAM, s
    from pygam import ExpectileGAM
    from scipy import stats

    gam1 = LinearGAM()
    #gam1 = ExpectileGAM(expectile=0.25)

    y1 = adata_cells[:,GENE].layers[layer]
    y1[y1<0]=0 #negative numbers from MAGIC to zero
    #y1 = stats.zscore(np.array(y1)) #scale to better compare the two datasets
    X1 = np.empty( (len(y1),1), dtype=object)
    X1[:,0] = adata_cells.obs['dpt_cellalign']
    y1 = y1[adata_cells.obs['dpt_cellalign']>=0]
    X1 = X1[adata_cells.obs['dpt_cellalign']>=0,]

    gam1.gridsearch(X1 , y1)

    XX1 = gam1.generate_X_grid(term=0, n=N_grid)
    
    NAMES1 = adata_cells.obs[cluster_key].cat.categories
    CC1 = list()
    CLST1 = np.empty(len(NAMES1))
    for i in range(len(NAMES1)):
        CLST1[i] = np.mean( adata_cells.obs['dpt_cellalign'][adata_cells.obs[cluster_key] == NAMES1[i]] )
    for i in range(len(XX1[:,0])):
        CC1.append( NAMES1[ np.argmin( np.abs( CLST1-XX1[i,0] ) ) ] )
    
    gam2 = LinearGAM()
    #gam2 = ExpectileGAM(expectile=0.25)
    
    y2 = adata_nuclei[:,GENE].layers[layer]
    y2[y2<0]=0
    #y2 = stats.zscore(np.array(y2))
    X2 = np.empty( (len(y2),1), dtype=object)
    X2[:,0] = adata_nuclei.obs['dpt_cellalign']
    y2 = y2[adata_nuclei.obs['dpt_cellalign']>=0]
    X2 = X2[adata_nuclei.obs['dpt_cellalign']>=0,]

    gam2.gridsearch(X2 , y2 )

    XX2 = gam2.generate_X_grid(term=0, n=N_grid)

    NAMES2 = adata_nuclei.obs[cluster_key].cat.categories
    CC2 = list()
    CLST2 = np.empty(len(NAMES2))
    for i in range(len(NAMES2)):
        CLST2[i] = np.mean( adata_nuclei.obs['dpt_cellalign'][adata_nuclei.obs[cluster_key] == NAMES2[i]] )
    for i in range(len(XX2[:,0])):
        CC2.append( NAMES2[ np.argmin( np.abs( CLST2-XX2[i,0] ) ) ] )
    
    hue1 = np.repeat('Cells', N_grid)
    hue2 = np.repeat('Nuclei', N_grid)
    
    yy1 = gam1.partial_dependence(term=0, X=XX1)
    yy2 = gam2.partial_dependence(term=0, X=XX2)

    ci1low= gam1.partial_dependence(term=0, X=XX1, width=.95)[1][:,0] 
    ci2low= gam2.partial_dependence(term=0, X=XX2, width=.95)[1][:,0] 
    cimin = [min([i,j]) for i,j in zip(ci1low,ci2low)]

    ci1high= gam1.partial_dependence(term=0, X=XX1, width=.95)[1][:,1]
    ci2high= gam2.partial_dependence(term=0, X=XX2, width=.95)[1][:,1]
    cimax = [max([i,j]) for i,j in zip(ci1high,ci2high)]
    
    plt.figure();
    fig, ax1 = plt.subplots(1,1);
    
    
    f=sns.lineplot(x=np.concatenate( (XX1[:, 0], XX2[:, 0]) ) , y=np.concatenate( (yy1, yy2) ) , 
                   hue = np.concatenate( (hue1, hue2) ) , ax=ax1, c='Red', linewidth=3, estimator=None)
    y_coords = [line.get_ydata() for line in f.lines]
    f=ax1.hlines(xmin=0, xmax=1, y=0, colors='Black', alpha=.25)
    #f=sns.lineplot(x=XX1[:, 0], y=cimin, ax=ax1, c='Blue', alpha=.5)
    #f=sns.lineplot(x=XX1[:, 0], y=cimax, ax=ax1, c='Blue', alpha=.5)
    f=sns.lineplot(x=XX1[:, 0], y=ci1low, ax=ax1, c='Blue', alpha=.5)
    f=sns.lineplot(x=XX1[:, 0], y=ci1high, ax=ax1, c='Blue', alpha=.5)
    f=sns.lineplot(x=XX2[:, 0], y=ci2low, ax=ax1, c='Red', alpha=.5)
    f=sns.lineplot(x=XX2[:, 0], y=ci2high, ax=ax1, c='Red', alpha=.5)
    f=ax1.fill_between(x=XX1[:, 0], y1=ci1low, y2=ci1high, alpha=.5)
    f=ax1.fill_between(x=XX2[:, 0], y1=ci2low, y2=ci2high, alpha=.5)


    
    #f=ax1.fill_between(x=XX1[:, 0], y1=cimin, y2=cimax, alpha=.5)
    f=ax1.vlines(ymin=min([min(y_coords[0]),min(y_coords[1])]),ymax=max( [max(y_coords[0]),max(y_coords[1])] ), x=0, colors='Black', alpha=.25, linestyles='dashed')
    for i in NAMES1:
        idx1=[j==i for j in CC1]
        idx2=[j==i for j in CC2]
        Xtime = max( [max(XX1[idx1,0]), max(XX2[idx2,0] )]  )
        f=ax1.vlines(ymin=min([min(y_coords[0]),min(y_coords[1])]),ymax=max([max(y_coords[0]),max(y_coords[1])]), x=Xtime, colors='Black', alpha=.25, linestyles='dashed')
        f=ax1.text(s=i, x=Xtime, y=min([min(y1),min(y2)]), rotation=90, )
    
    ax1.set_title(f'Pseudotime curves (normalized) of {GENE} in Nuclei and Cells')
    ax1.set_xlabel("Pseudotime")
    ax1.set_ylabel("Expression curve")
    plt.savefig(f'figures/genes_overlap_plot_{GENE}_{layer}.png', bbox_inches='tight', transparent=True )






def calc_genes_overlap(adata_cells, 
                       adata_nuclei,
                       N_grid=100, 
                       GENE='MKI67', 
                       layer='MAGIC_counts', 
                       cluster_key='spermatogenesis_single', 
                       plot_size=(12,4)):
    
    plt.rcParams['figure.figsize'] = plot_size

    N=N_grid
    from pygam import LinearGAM, s
    from scipy import stats

    gam1 = LinearGAM()

    y1 = adata_cells[:,GENE].layers[layer]
    y1[y1<0]=0
    #y1 = stats.zscore(np.array(y1))
    X1 = np.empty( (len(y1),1), dtype=object)
    X1[:,0] = adata_cells.obs['dpt_cellalign']

    gam1.gridsearch(X1 , y1)

    XX1 = gam1.generate_X_grid(term=0, n=N)
    
    NAMES1 = adata_cells.obs[cluster_key].cat.categories
    CC1 = list()
    CLST1 = np.empty(len(NAMES1))
    for i in range(len(NAMES1)):
        CLST1[i] = np.mean( adata_cells.obs['dpt_cellalign'][adata_cells.obs[cluster_key] == NAMES1[i]] )
    for i in range(len(XX1[:,0])):
        CC1.append( NAMES1[ np.argmin( np.abs( CLST1-XX1[i,0] ) ) ] )
    
    gam2 = LinearGAM()
    
    y2 = adata_nuclei[:,GENE].layers[layer]
    y2[y2<0]=0
    #y2 = stats.zscore(np.array(y2))
    X2 = np.empty( (len(y2),1), dtype=object)
    X2[:,0] = adata_nuclei.obs['dpt_cellalign']

    gam2.gridsearch(X2 , y2 )

    XX2 = gam2.generate_X_grid(term=0, n=N)

    NAMES2 = adata_nuclei.obs[cluster_key].cat.categories
    CC2 = list()
    CLST2 = np.empty(len(NAMES2))
    for i in range(len(NAMES2)):
        CLST2[i] = np.mean( adata_nuclei.obs['dpt_cellalign'][adata_nuclei.obs[cluster_key] == NAMES2[i]] )
    for i in range(len(XX2[:,0])):
        CC2.append( NAMES2[ np.argmin( np.abs( CLST2-XX2[i,0] ) ) ] )

    yy1 = gam1.partial_dependence(term=0, X=XX1)
    yy2 = gam2.partial_dependence(term=0, X=XX2)

    return yy1, yy2, XX1[:,0], XX2[:,0]
