try:
    import scanpy as sc
except ImportError:
    _has_scanpy = False
else:
    _has_scanpy = True


# https://www.nature.com/articles/nbt.3192
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTYMK", "RFC2", "RPA2", "NASP", "POLD3", "RPA3", "PRIM1", "UHRF1", "MSH2",
    "ATAD2", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD2", "MSH6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CENPF", "TACC3",
    "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF2C", "MAD2L1", "BUB1B",
    "NDC80", "SKA3", "CENPE", "TTK", "CDC20", "BIRC5", "KIF20B", "PSRC1", "ANLN", "KIF23",
    "AURKA", "PLK1", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
]

def cluster_cells(adata, targets="all", do_log=True, do_norm=True, resolution=1):
    """Cluster cells into subgroups via the Leiden algorithm

    Args:
      adata : Cytoprofiling anndata object

    Returns:
      adata with PCA, neighborhood graph, and leiden clusters
    """
    if not _has_scanpy:
        raise ImportError("scanpy is required for cluster_cells")
    
    if not targets == "all":
        adata = adata[:, targets]

    adata.layers["counts"] = adata.X.copy()

    if do_norm:
        sc.pp.normalize_total(adata)

    if do_log:
        sc.pp.log1p(adata)

    sc.tl.pca(adata)

    sc.pp.neighbors(adata)

    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2)
    return adata

def assign_cell_phase(adata, s_genes=S_GENES, g2m_genes=G2M_GENES, do_log=True, do_norm=True):
    """Assign cell phase based on known targets

    Args:
      adata : Cytoprofiling anndata object

    Returns:
      adata with cell phase scores and classifications
    """
    if not _has_scanpy:
        raise ImportError("scanpy is required for assign_cell_phase")

    s_genes = [x for x in s_genes if x in adata.var_names]
    g2m_genes = [x for x in g2m_genes if x in adata.var_names]

    adata.layers["counts"] = adata.X.copy()

    if do_norm:
        sc.pp.normalize_total(adata)

    if do_log:
        sc.pp.log1p(adata)

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    return adata