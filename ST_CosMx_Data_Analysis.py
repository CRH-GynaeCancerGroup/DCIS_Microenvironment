# Python Codes for Spatial Transcriptomics (CosMx) Data Analyses

# ==============================================================================
# STEP 1 — Quality Check of Supervised Cell Type Annotation
# ==============================================================================
# Before starting the analysis, prepare the following files:
# 1. metadata.csv (exported from AtoMx)
# 2. polygon.csv (exported from AtoMx)
# 3. exprMat (exported from AtoMx)
# 4. core-patient_info.csv (created per TMA by user)
# The "core-patient_info" file comprise FOV, Core, Patient ID, Clinical Group, Use_FOV (Y=1 or N=0)
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)

# STEP 1-A: Load Required Files
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from pathlib import Path
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from IPython.display import display, Markdown
import matplotlib.patches as mpatches

Tk().withdraw()

# 1. Load metadata.csv
meta_path = askopenfilename(title="Select METADATA CSV", filetypes=[("CSV files", "*.csv")])
cell_meta = pd.read_csv(meta_path)
display(Markdown("### ✅ Loaded metadata"))
display(cell_meta.head())
print("Meta shape:", cell_meta.shape)

# 2. Load polygon.csv
seg_path = askopenfilename(title="Select POLYGONS CSV", filetypes=[("CSV files", "*.csv")])
polygons = pd.read_csv(seg_path)
display(Markdown("### ✅ Loaded polygons"))
display(polygons.head())
print("Polygons shape:", polygons.shape)

# 3. Load exprMat (gene count matrix)
expr_path = askopenfilename(title="Select exprMat CSV", filetypes=[("CSV files", "*.csv")])
expr_raw = pd.read_csv(expr_path)
display(Markdown("### ✅ Loaded raw expression table"))
display(expr_raw.head())
print("Raw shape:", expr_raw.shape)

# STEP 1-B: Define Parameters
TMA_ID = "TMA1"
CELLTYPE_COL = "RNA_TMA1_20260409_TK_Cell.Typing.InSituType.1_1_clusters"
CELL_ID_COL = "cell"
FOV_COL = "fov"
X_COL = "CenterX_global_px"
Y_COL = "CenterY_global_px"
CORE_COL = "core"

BASE_MAIN_DIR = Path.home() / "Desktop" / f"{TMA_ID}_STEP1_results"
BASE_MAIN_DIR.mkdir(parents=True, exist_ok=True)
print(f"✅ Output directory created: {BASE_MAIN_DIR}")

# STEP 1-C: Build expr_mat
if not {"fov", "cell_ID"}.issubset(expr_raw.columns):
    raise ValueError("exprMat must contain fov and cell_ID")

expr_mat = expr_raw.copy()
expr_mat["fov"] = expr_mat["fov"].astype(str).str.strip()
expr_mat["cell_ID"] = expr_mat["cell_ID"].astype(str).str.strip()
expr_mat["cell_key"] = expr_mat["fov"] + "_" + expr_mat["cell_ID"]

gene_cols = [c for c in expr_mat.columns if c not in ["fov", "cell_ID", "cell_key"]]
expr_mat[gene_cols] = expr_mat[gene_cols].astype("int32")
expr_mat = expr_mat[["cell_key"] + gene_cols]

print("\n✅ Expression matrix prepared")
print("Shape:", expr_mat.shape)

# STEP 1-D: Fix metadata cell_key
cell_meta = cell_meta.copy()
cell_meta["cell_id"] = cell_meta["cell_id"].astype(str).str.strip()

def parse_cell_id(x):
    parts = str(x).split("_")
    return (parts[1], parts[2], parts[3]) if len(parts) == 4 else (None, None, None)

CELL_ID_SOURCE = "cell_id"
parsed = cell_meta[CELL_ID_SOURCE].astype(str).apply(parse_cell_id)

cell_meta["slide"] = parsed.str[0]
cell_meta["fov_parsed"] = parsed.str[1]
cell_meta["cell_ID_parsed"] = parsed.str[2]
cell_meta["cell_key"] = cell_meta["fov_parsed"].astype("string") + "_" + cell_meta["cell_ID_parsed"].astype("string")

cell_meta = cell_meta[cell_meta["fov_parsed"].notna() & cell_meta["cell_ID_parsed"].notna()]
cell_meta = cell_meta[(cell_meta["fov_parsed"] != "None") & (cell_meta["cell_ID_parsed"] != "None")]

print("\nSample cell_keys:", cell_meta["cell_key"].head())

if "cell_key" not in expr_mat.columns:
    if ("fov" in expr_mat.columns) and ("cell_ID" in expr_mat.columns):
        expr_mat["cell_key"] = expr_mat["fov"].astype(str) + "_" + expr_mat["cell_ID"].astype(str)
    else:
        expr_mat["cell_key"] = expr_mat.index.astype(str)

expr_mat["cell_key"] = expr_mat["cell_key"].astype(str).str.strip()
print("expr_mat cell_key sample:", list(expr_mat["cell_key"].head(5)))

# STEP 1-E: Data Alignment
expr_mat["cell_key"] = expr_mat["cell_key"].astype(str).str.strip()
cell_meta["cell_key"] = cell_meta["cell_key"].astype(str).str.strip()

print("\nUnique expr keys:", expr_mat["cell_key"].nunique())
print("Unique meta keys:", cell_meta["cell_key"].nunique())

common = sorted(set(expr_mat["cell_key"]) & set(cell_meta["cell_key"]))
print("Common keys:", len(common))

if len(common) == 0:
    raise ValueError("No overlap between exprMat and metadata — check ID format")

expr_mat = expr_mat.set_index("cell_key")
cell_meta = cell_meta.set_index("cell_key")
expr_mat = expr_mat.loc[common]
cell_meta = cell_meta.loc[common]

print("Aligned cells:", len(common))
print("Perfect alignment:", all(expr_mat.index == cell_meta.index))
print("Duplicate check - expr:", expr_mat.index.duplicated().sum(), "meta:", cell_meta.index.duplicated().sum())

# STEP 1-F: QC Features
cell_meta["nCount"] = expr_mat.sum(axis=1)
cell_meta["nFeature"] = (expr_mat > 0).sum(axis=1)
cell_meta["log_nCount"] = np.log1p(cell_meta["nCount"])
cell_meta["log_nFeature"] = np.log1p(cell_meta["nFeature"])
cell_meta["qc_score"] = cell_meta["log_nCount"] + cell_meta["log_nFeature"]

print(cell_meta[["nCount", "nFeature", "qc_score"]].head())

# STEP 1-G: Original Cell Type QC Sanity Check
# 1. Load core-patient_info
meta_path = askopenfilename(title="Select core_patient_info CSV", filetypes=[("CSV files", "*.csv")])
core_map = pd.read_csv(meta_path)
core_map["fov"] = core_map["fov"].astype(str).str.strip()
core_map["core"] = core_map["core"].astype(str).str.strip()
core_map = core_map[core_map["use_fov"] == 1].copy()

display(Markdown("### ✅ Loaded core_map"))
display(core_map.head())

# 2. Ensure consistent formatting
cell_meta["fov"] = cell_meta["fov"].astype(str).str.strip()
core_map["fov"] = pd.to_numeric(core_map["fov"], errors="coerce")
core_map["fov"] = core_map["fov"].round(0).astype("Int64").astype(str)

print("Overlap:", len(set(cell_meta["fov"]) & set(core_map["fov"])))

# 3. Attach core to metadata
cell_meta = cell_meta.merge(core_map[["fov", "core", "patient_id"]], on="fov", how="left", validate="m:1")

print("\nCore assignment in metadata:", cell_meta["core"].isna().sum(), "cells without core")
assert "core" in cell_meta.columns, "Core merge failed in metadata"
assert cell_meta["core"].notna().sum() > 0, "No valid core assignments found"

# 4. Define parameters & Colors
out_orig = BASE_MAIN_DIR / f"{TMA_ID}_CellTyping_Original"
out_orig.mkdir(exist_ok=True)
core_img_dir = out_orig / "core_images"
core_img_dir.mkdir(exist_ok=True)

COLOR_MAP = {
    "Normal.Epi": "#c5fc3d", "Tumour.Epi": "#058511", "FOLR2.Mac": "#ed0606", "IL4I1.Mac": "#ed0606",
    "SPP1.Mac": "#ed0606", "TREM2.Mac": "#ed0606", "C.Mo": "#ed0606", "cDC1": "#ff7f0e", "cDC2": "#ff7f0e",
    "CCR7.DC": "#ff7f0e", "Mast": "#f6e606", "CD4.T": "#52a1ce", "CD8.T": "#8bc6e7", "Treg": "#08517a",
    "Naive.B": "#c88df6", "Memory.B": "#c88df6", "Plasma": "#a65628", "NK": "#5b0568", "iCAF": "#09dceb",
    "myCAF": "#09dceb", "CD34.CAF": "#09dceb", "Blood.EC": "#fa91ca", "Lymph.EC": "#fa91ca",
    "Myoepithelial": "#8e17e9", "Pericyte": "#d2c98d", "Adipocyte": "#a7a7a7", "Unknown": "#FFFFFF"
}
BACKGROUND_COLOR = "#ffffff"

def map_colors(labels):
    return [COLOR_MAP.get(lbl, "#FFFFFF") for lbl in labels]

# 5. Compute proportions
prop = cell_meta[CELLTYPE_COL].value_counts(normalize=True)
prop.to_csv(out_orig / "celltype_proportions.csv")
print("✅ Exported cell type proportions")

# 6. Plot core images
if "core" not in cell_meta.columns:
    raise ValueError("core column missing in metadata after merge")

for core in sorted(cell_meta["core"].dropna().unique()):
    sub = cell_meta[cell_meta["core"] == core]
    if sub.empty: continue

    plt.figure(figsize=(6,6))
    ax = plt.gca()
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.set_aspect("equal")

    if X_COL not in sub.columns or Y_COL not in sub.columns:
        raise ValueError("Spatial coordinates missing")

    plt.scatter(sub[X_COL], sub[Y_COL], c=map_colors(sub[CELLTYPE_COL]), s=3, edgecolors="black", linewidths=0.15, alpha=0.9)

    unique_labels = sorted(sub[CELLTYPE_COL].dropna().unique())
    handles = [mpatches.Patch(color=COLOR_MAP.get(lbl, "#FFFFFF"), label=lbl) for lbl in unique_labels]
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=6, frameon=False)

    plt.title(f"Core {core} - Original", color="white")
    plt.axis("off")
    plt.savefig(core_img_dir / f"core_{core}_original.png", dpi=300, facecolor=BACKGROUND_COLOR, bbox_inches="tight")
    plt.close()

print("✅ Exported all core images")

# STEP 1-H: Check RNA Quality per Cell Type
# 1. Violin plots
def violin_plot(values, groups, title, outpath):
    values = np.asarray(values)
    groups = np.asarray(groups)
    mask = ~pd.isna(values) & ~pd.isna(groups)
    values, groups = values[mask], groups[mask]
    unique_groups = np.unique(groups)
    data_to_plot = [values[groups == g] for g in unique_groups]

    plt.figure(figsize=(14, 5))
    plt.violinplot(data_to_plot, showmeans=False, showmedians=True)
    plt.xticks(range(1, len(unique_groups) + 1), unique_groups, rotation=90)
    plt.title(title)
    plt.ylabel(title)
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

required_cols = ["nCount", "nFeature", CELLTYPE_COL]
missing = [c for c in required_cols if c not in cell_meta.columns]
if len(missing) > 0:
    raise ValueError(f"Missing required columns in cell_meta: {missing}")

violin_plot(cell_meta["nCount"].values, cell_meta[CELLTYPE_COL].values, "nCount vs Cell Type", out_orig / "nCount_violin.png")
violin_plot(cell_meta["nFeature"].values, cell_meta[CELLTYPE_COL].values, "nFeature vs Cell Type", out_orig / "nFeature_violin.png")
print("✅ STEP1-H: QC violin plots generated")

# 2. QC scatter
df = cell_meta[["nCount", "nFeature"]].copy().replace([np.inf, -np.inf], np.nan).dropna()

plt.figure(figsize=(6,5))
plt.scatter(df["nCount"], df["nFeature"], s=2, alpha=0.4)
plt.xlabel("nCount")
plt.ylabel("nFeature")
plt.title("QC structure")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.savefig(out_orig / "QC_scatter.png", dpi=300)
plt.close()
print("✅ STEP1-H: QC scatter plot generated")

# 3. Export QC summary table
qc_out = out_orig / "QC_summary"
qc_out.mkdir(exist_ok=True)

global_stats = {
    "nCount_p1": np.percentile(cell_meta["nCount"], 1), "nCount_p5": np.percentile(cell_meta["nCount"], 5),
    "nCount_p10": np.percentile(cell_meta["nCount"], 10), "nFeature_p1": np.percentile(cell_meta["nFeature"], 1),
    "nFeature_p5": np.percentile(cell_meta["nFeature"], 5), "nFeature_p10": np.percentile(cell_meta["nFeature"], 10),
}
global_df = pd.DataFrame([global_stats])
global_df.to_csv(qc_out / "global_thresholds.csv", index=False)
print("✅ Exported global QC thresholds")

summary_list = []
for ct in sorted(cell_meta[CELLTYPE_COL].dropna().unique()):
    sub = cell_meta[cell_meta[CELLTYPE_COL] == ct]
    row = {
        "cell_type": ct, "n_cells": len(sub),
        "median_nCount": sub["nCount"].median(), "median_nFeature": sub["nFeature"].median(),
        "p5_nCount": np.percentile(sub["nCount"], 5), "p5_nFeature": np.percentile(sub["nFeature"], 5),
        "p10_nCount": np.percentile(sub["nCount"], 10), "p10_nFeature": np.percentile(sub["nFeature"], 10),
    }
    summary_list.append(row)

summary_df = pd.DataFrame(summary_list)
nCount_cut = global_stats["nCount_p5"]
nFeature_cut = global_stats["nFeature_p5"]

summary_df["pct_below_nCount_cut"] = cell_meta.groupby(CELLTYPE_COL)["nCount"].apply(lambda x: (x < nCount_cut).mean()).values
summary_df["pct_below_nFeature_cut"] = cell_meta.groupby(CELLTYPE_COL)["nFeature"].apply(lambda x: (x < nFeature_cut).mean()).values

summary_df.to_csv(qc_out / "celltype_qc_summary.csv", index=False)
print("✅ Exported cell-type QC summary")

# ==============================================================================
# STEP 2 — Identify Stromal Cell Subsets in Spatially Distinct Niches
# ==============================================================================
# Define stromal cell localisation within the following spatially distinct niches:
# ID: stromal cells locate within ductal (luminal) cavities inside tumor nests
# IC: stromal cells localise INSIDE tumour nest and CONTACT with epithelial cells
# OC: stromal cells locate OUTSIDE tumor nests and CONTACT with epithelial cells
# ON: stromal cells locate OUTSIDE tumor nests within NEAR epithelium (4 cell-distance)
# OF: stromal cells locate OUTSIDE tumor nests FAR from epithelium

# Before starting the analysis, prepare the following files:
# 1. metadata.csv (exported from AtoMx)
# 2. polygon.csv (exported from AtoMx)
# 3. exprMat (exported from AtoMx)
# 4. core-patient_info.csv (created per TMA by user)
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)
 
# STEP 2-A: Load Required Files
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pickle
from pathlib import Path
from tkinter import Tk, filedialog
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from tqdm import tqdm

plt.ioff()
Tk().withdraw() 

# 1. Load Metadata
meta_path = filedialog.askopenfilename(title="Select METADATA CSV", filetypes=[("CSV files", "*.csv")])
if not meta_path: raise FileNotFoundError("Metadata file not selected.")
cell_meta = pd.read_csv(meta_path)

# 2. Load Polygons
seg_path = filedialog.askopenfilename(title="Select POLYGONS CSV", filetypes=[("CSV files", "*.csv")])
if not seg_path: raise FileNotFoundError("Polygons file not selected.")
seg = pd.read_csv(seg_path)

# 3. Load Core/Patient Info
core_path = filedialog.askopenfilename(title="Select CORE/PATIENT INFO CSV", filetypes=[("CSV files", "*.csv")])
if not core_path: raise FileNotFoundError("Core info file not selected.")
core_info = pd.read_csv(core_path)

print(f"✅ Loaded: Metadata ({cell_meta.shape[0]} cells), Core Info ({core_info.shape[0]} cores), Polygons ({seg.shape[0]} points)")

# STEP 2-B: Preprocessing & QC 
CELL_ID_COL = "cell"
FOV_COL = "fov"
CELLTYPE_COL = "RNA_TMA1_20260409_TK_Cell.Typing.InSituType.1_1_clusters"
X_COL = "CenterX_global_px"
Y_COL = "CenterY_global_px"
CORE_COL = "core"
PATIENT_COL = "patient_id"
TMA_ID = "TMA1"

# Define Epithelial vs Stromal
tumor_epi_types = ["Tumour.Epi"]
normal_epi_types = ["Normal.Epi"]
epi_types = tumor_epi_types + normal_epi_types

cell_meta["is_tumor_epi"] = cell_meta[CELLTYPE_COL].isin(tumor_epi_types)
cell_meta["is_normal_epi"] = cell_meta[CELLTYPE_COL].isin(normal_epi_types)
cell_meta["is_epithelial"] = cell_meta[CELLTYPE_COL].isin(epi_types)
cell_meta["is_stromal"] = (~cell_meta["is_epithelial"]).fillna(False)

# Assign Stromal Subclasses
cell_meta["stromal_subset"] = None
stromal_mask = cell_meta["is_stromal"]
cell_meta.loc[stromal_mask, "stromal_subset"] = cell_meta.loc[stromal_mask, CELLTYPE_COL]

# Merge Core Info & Filter FOVs
cell_meta = cell_meta.merge(core_info, on=FOV_COL, how="left")
if not cell_meta[CORE_COL].notna().all():
    raise ValueError("Some FOVs missing core annotation.")

cell_meta = cell_meta[cell_meta["use_fov"] == 1].copy()
print(f"   Filtered to analysable FOVs: {cell_meta.shape[0]} cells remaining.")

# Core Epithelial Content QC
core_epi_fraction = cell_meta.groupby(CORE_COL)["is_epithelial"].mean()

def core_sort_key(core):
    match = re.match(r"([A-Za-z]+)(\d+)", str(core))
    return (match.group(1), int(match.group(2))) if match else (str(core), 0)

core_epi_fraction = core_epi_fraction.sort_index(key=lambda idx: idx.map(core_sort_key))

bad_cores = core_epi_fraction[core_epi_fraction < 0.01].index.tolist()
if bad_cores:
    print(f"⚠️ Excluding {len(bad_cores)} cores with <1% epithelial content.")
    cell_meta_qc = cell_meta[~cell_meta[CORE_COL].isin(bad_cores)].copy()
else:
    cell_meta_qc = cell_meta.copy()

print(f"✅ QC Complete: {cell_meta_qc.shape[0]} cells in {cell_meta_qc[CORE_COL].nunique()} cores.")

# STEP 2-C: Build Cell Polygons
CELL_POLY_X = "x_global_px"
CELL_POLY_Y = "y_global_px"

cell_polygons = {}
for cid, df_poly in seg.groupby(CELL_ID_COL):
    coords = df_poly[[CELL_POLY_X, CELL_POLY_Y]].values
    if len(coords) >= 3:
        cell_polygons[cid] = Polygon(coords)

print(f"✅ Built {len(cell_polygons)} cell polygons.")

# Create per-core polygon dictionary
cell_polygons_by_core = {}
for core, df_core in cell_meta_qc.groupby(CORE_COL):
    valid_ids = [cid for cid in df_core[CELL_ID_COL] if cid in cell_polygons]
    cell_polygons_by_core[core] = {cid: cell_polygons[cid] for cid in valid_ids}

# STEP 2-D: Stromal Niche Definition (OC/ON/OF → IC/ID)
MEDIAN_EPI_AREA = 8192.5
MEDIAN_STROMAL_AREA = 6274.5
MEDIAN_GAP_PX = 185

CELL_DIST_PX = int(np.sqrt(MEDIAN_STROMAL_AREA / np.pi))
MAX_NEAR_DIST = 4 * CELL_DIST_PX
FRAGMENT_AREA_THRESH = 6.0 * MEDIAN_EPI_AREA
MIN_CELLS_PER_NEST = 6
OUTER_GAP_FILL_DIST = CELL_DIST_PX
SMALL_HOLE_FILL_DIST = MEDIAN_GAP_PX / 2
SMALL_HOLE_AREA_THRESH = 1.0 * MEDIAN_EPI_AREA

def detect_epi_clusters(epi_polys):
    if len(epi_polys) == 0: return []
    u = unary_union(epi_polys)
    return [u] if isinstance(u, Polygon) else list(u.geoms)

def classify_epi_clusters(clusters, epi_poly_map):
    tumour_clusters, fragments = [], []
    epi_centroids = [p.centroid for p in epi_poly_map.values()]
    for poly in clusters:
        n_cells = sum(poly.contains(c) for c in epi_centroids)
        if n_cells >= MIN_CELLS_PER_NEST and poly.area >= FRAGMENT_AREA_THRESH:
            tumour_clusters.append(poly)
        else:
            fragments.append(poly)
    return tumour_clusters, fragments

def assign_outer_niche(stromal_df, tumour_epi_polys, core_polys):
    outer = []
    for _, row in stromal_df.iterrows():
        cid, s_poly = row["cell"], core_polys.get(row["cell"])
        if s_poly is None or len(tumour_epi_polys) == 0:
            outer.append(np.nan); continue
        d = min(s_poly.distance(epi_poly) for epi_poly in tumour_epi_polys)
        if d <= CELL_DIST_PX: outer.append("OC")
        elif d <= MAX_NEAR_DIST: outer.append("ON")
        else: outer.append("OF")
    out = stromal_df.copy()
    out["outer_niche"] = outer
    return out

def unite_tumour_nest_epi(epi_polys):
    if len(epi_polys) == 0: return []
    u = unary_union(epi_polys).buffer(OUTER_GAP_FILL_DIST).buffer(-OUTER_GAP_FILL_DIST)
    return [u] if isinstance(u, Polygon) else list(u.geoms)

def fill_and_preserve_holes(nest_poly):
    filled = nest_poly.buffer(SMALL_HOLE_FILL_DIST).buffer(-SMALL_HOLE_FILL_DIST)
    polys = [filled] if isinstance(filled, Polygon) else list(filled.geoms)
    cleaned = []
    for poly in polys:
        keep_holes = [h for h in poly.interiors if Polygon(h).area >= SMALL_HOLE_AREA_THRESH]
        cleaned.append(Polygon(poly.exterior.coords, holes=[h.coords for h in keep_holes]))
    return unary_union(cleaned)

def build_final_nests(unioned_nests):
    final = []
    for nest in unioned_nests:
        filled = fill_and_preserve_holes(nest)
        polys = [filled] if isinstance(filled, Polygon) else list(filled.geoms)
        for p in polys:
            final.append({"poly": p, "shell": Polygon(p.exterior), "holes": [Polygon(h) for h in p.interiors]})
    return final

def assign_inner_niche(stromal_df, final_nests, core_polys):
    final = stromal_df.copy()
    final["final_niche"] = final["outer_niche"]
    for i, row in final.iterrows():
        cid, s_poly = row["cell"], core_polys.get(row["cell"])
        if s_poly is None: continue
        c = s_poly.centroid
        for nest in final_nests:
            for hole in nest["holes"]:
                if hole.covers(c):
                    final.at[i, "final_niche"] = "ID"; break
            else:
                if nest["shell"].covers(c):
                    final.at[i, "final_niche"] = "IC"
            if final.at[i, "final_niche"] in ["IC", "ID"]: break
    return final

# Main Loop
stromal_all, core_to_nests, core_to_fragments, core_to_polys = [], {}, {}, {}

for core, df_core in tqdm(cell_meta_qc.groupby(CORE_COL), desc="Processing Cores"):
    core_polys = cell_polygons_by_core.get(core, {})
    if len(core_polys) == 0: continue

    epi_ids = df_core.loc[df_core["is_epithelial"], "cell"]
    stromal_df = df_core[df_core["is_stromal"]].copy()
    epi_polys = {cid: core_polys[cid] for cid in epi_ids if cid in core_polys}
    if len(epi_polys) == 0: continue

    clusters = detect_epi_clusters(list(epi_polys.values()))
    tumour_clusters, fragments = classify_epi_clusters(clusters, epi_polys)
    
    tumour_epi_cells = [poly for cid, poly in epi_polys.items() if any(tc.covers(poly.centroid) for tc in tumour_clusters)]
    stromal_outer = assign_outer_niche(stromal_df, tumour_epi_cells, core_polys)
    
    tumour_union = unite_tumour_nest_epi(tumour_epi_cells)
    final_nests = build_final_nests(tumour_union)
    stromal_final = assign_inner_niche(stromal_outer, final_nests, core_polys)

    stromal_all.append(stromal_final)
    core_to_nests[core], core_to_fragments[core], core_to_polys[core] = final_nests, fragments, core_polys

stromal_niches_master = pd.concat(stromal_all, ignore_index=True)
print(f"✅ Niche Definition Complete: {stromal_niches_master.shape[0]} stromal cells assigned.")

# STEP 2-E: QC Plotting
BASE_MAIN_DIR = Path.home() / "Desktop" / f"{TMA_ID}_STEP2_results"
QC_DIR = BASE_MAIN_DIR / f"{TMA_ID}_niche_definition_QC"
QC_DIR.mkdir(parents=True, exist_ok=True)

QC_COLORS = {"background": "#1a1a1a", "cell_outline": "#000000", "tumor_nest": "#cbe833", "fragment": "#87CEEB", "hole": "#808080", "IC": "#a65628", "ID": "#984ea3", "OC": "#ed0606", "ON": "#ff7f0e", "OF": "#f6e606"}

def export_core_qc_image(core):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal', adjustable='box')
    ax.set_facecolor(QC_COLORS["background"])
    
    for poly in core_to_polys[core].values():
        if poly.is_empty: continue
        x, y = poly.exterior.xy
        ax.plot(x, y, color=QC_COLORS["cell_outline"], lw=0.25, alpha=0.5)
    
    for nest in core_to_nests[core]:
        x, y = nest["poly"].exterior.xy
        ax.fill(x, y, color=QC_COLORS["tumor_nest"], alpha=0.35)
        for h in nest["holes"]:
            xh, yh = h.exterior.xy
            ax.fill(xh, yh, color=QC_COLORS["hole"], alpha=1.0)
    
    for frag in core_to_fragments.get(core, []):
        if frag.is_empty: continue
        x, y = frag.exterior.xy
        ax.fill(x, y, color=QC_COLORS["fragment"], alpha=0.6)
    
    df = stromal_niches_master[stromal_niches_master["core"] == core]
    for niche in ["IC","ID","OC","ON","OF"]:
        sub = df[df["final_niche"] == niche]
        for cid in sub["cell"]:
            poly = core_to_polys[core].get(cid)
            if poly is None or poly.is_empty: continue
            x, y = poly.exterior.xy
            ax.fill(x, y, color=QC_COLORS[niche], alpha=0.85)
    
    ax.axis('off')
    
    all_x, all_y = [], []
    for poly in core_to_polys[core].values():
        if not poly.is_empty:
            all_x.extend(poly.exterior.xy[0])
            all_y.extend(poly.exterior.xy[1])
    
    if all_x and all_y:
        margin = 50
        ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
        ax.set_ylim(min(all_y) - margin, max(all_y) + margin)
    
    ax.set_title(f"Core {core}", color="white", pad=10)
    fig.savefig(QC_DIR / f"{core}_niche_QC.png", dpi=150, bbox_inches='tight')
    plt.close(fig)

cores_to_plot = sorted(core_to_nests.keys()) 
print(f"Generating QC plots for ALL {len(cores_to_plot)} cores...")
for core in tqdm(cores_to_plot, desc="Plotting All Cores"):
    export_core_qc_image(core)

print(f"✅ QC Plots saved to {QC_DIR}")

# STEP 2-F: Count Stromal Cell Subsets Per Niche
STROMAL_CELL_TYPES = ['FOLR2.Mac', 'IL4I1.Mac', 'SPP1.Mac', 'TREM2.Mac', 'C.Mo', 'cDC1', 'cDC2', 'CCR7.DC', 'Mast', 'CD4.T', 'CD8.T', 'Treg', 'Naive.B', 'Memory.B', 'Plasma', 'NK', 'iCAF', 'myCAF', 'CD34.CAF', 'Myoepithelial', 'Pericyte', 'Adipocyte', 'Blood.EC', 'Lymph.EC']

df = stromal_niches_master.copy()

if "patient_id" not in core_info.columns or "clinical_group" not in core_info.columns:
    raise KeyError(f"core_info missing required columns. Found: {list(core_info.columns)}")

core_info_unique = core_info[[CORE_COL, "patient_id", "clinical_group"]].drop_duplicates(subset=[CORE_COL])
print(f"✅ Core metadata prepared: {core_info_unique.shape[0]} cores")

cols_to_drop = [c for c in ["patient_id", "clinical_group"] if c in df.columns]
if cols_to_drop:
    print(f"⚠️ Dropping existing columns {cols_to_drop} from stromal data to avoid merge conflicts...")
    df = df.drop(columns=cols_to_drop)

df = df.merge(core_info_unique, on=CORE_COL, how='left')

missing_meta = df["patient_id"].isna().sum()
if missing_meta > 0:
    print(f"⚠️ Warning: {missing_meta} cells could not be matched to core metadata.")

df_stromal = df[df[CELLTYPE_COL].isin(STROMAL_CELL_TYPES)].copy()
print(f"✅ Filtered to {len(df_stromal)} stromal cells of interest")

core_niche_counts = df_stromal.groupby([CORE_COL, "final_niche", CELLTYPE_COL]).size().reset_index(name="count")

core_niche_counts_wide = core_niche_counts.pivot_table(
    index=[CORE_COL, "final_niche"], 
    columns=CELLTYPE_COL, 
    values="count", 
    fill_value=0
).reset_index()

for col in STROMAL_CELL_TYPES:
    if col not in core_niche_counts_wide.columns:
        core_niche_counts_wide[col] = 0

core_meta = df[[CORE_COL, "patient_id", "clinical_group"]].drop_duplicates()
final_cols = [CORE_COL, "final_niche"] + STROMAL_CELL_TYPES + ["patient_id", "clinical_group"]
core_niche_counts_wide = core_niche_counts_wide.merge(core_meta, on=CORE_COL, how='left')[final_cols]

output_counts = BASE_MAIN_DIR / f"{TMA_ID}_core_niche_cell_counts.csv"
core_niche_counts_wide.to_csv(output_counts, index=False)
print(f"✅ Counts saved: {output_counts}")

# STEP 2-G: Count Epithelial Subtypes in Nest & Fragment
EPITHELIAL_TYPES = ["Tumour.Epi", "Normal.Epi"]
results = []

for core, df_core in cell_meta_qc.groupby(CORE_COL):
    core_polys = cell_polygons_by_core.get(core, {})
    if len(core_polys) == 0: continue
    epi_ids = df_core.loc[df_core["is_epithelial"], "cell"]
    epi_polys = {cid: core_polys[cid] for cid in epi_ids if cid in core_polys}
    if len(epi_polys) == 0: continue

    clusters = detect_epi_clusters(list(epi_polys.values()))
    tumour_clusters, fragments = classify_epi_clusters(clusters, epi_polys)
    cluster_list = [("nest", i+1, p) for i, p in enumerate(tumour_clusters)] + [("fragment", i+1, p) for i, p in enumerate(fragments)]

    for ctype, cid, poly in cluster_list:
        assigned = [c for c, p in epi_polys.items() if poly.covers(p.centroid)]
        if not assigned: continue
        subtypes = df_core.set_index("cell").loc[assigned, CELLTYPE_COL]
        counts = {st: (subtypes == st).sum() for st in EPITHELIAL_TYPES}
        results.append({"core": core, "cluster_type": ctype, "cluster_id": cid, "n_cells_total": len(assigned), **counts})

epi_cluster_counts = pd.DataFrame(results)
epi_cluster_counts.to_csv(BASE_MAIN_DIR / f"{TMA_ID}_epithelial_counts_per_cluster.csv", index=False)
print(f"✅ Epithelial counts saved: {BASE_MAIN_DIR / f'{TMA_ID}_epithelial_counts_per_cluster.csv'}")

# ==============================================================================
# STEP 3 — Export Key Objects
# ==============================================================================
# Export data files required for subsequent analyses
# This includes: cell_meta_full.csv, stromal_niches_full.csv, core_to_nest.pkl, and cell_polygons_by_core.pkl

EXPORT_DIR = Path.home() / "Desktop" / f"{TMA_ID}_STEP3_results"
EXPORT_DIR.mkdir(exist_ok=True, parents=True)

full_export = cell_meta_qc.merge(stromal_niches_master[["cell", "outer_niche", "final_niche"]], on="cell", how="left")
full_export.to_csv(EXPORT_DIR / "cell_meta_full.csv", index=False)
stromal_niches_master.to_csv(EXPORT_DIR / "stromal_niches_full.csv", index=False)

with open(EXPORT_DIR / "core_to_nests.pkl", "wb") as f:
    pickle.dump(core_to_nests, f)
with open(EXPORT_DIR / "cell_polygons_by_core.pkl", "wb") as f:
    pickle.dump(cell_polygons_by_core, f)

print(f"🎉 Pipeline Complete! Results saved to {EXPORT_DIR}")

# ==============================================================================
# STEP 4 — Identify Stromal Cell Subsets within TREM2 Mac Milieu
# ==============================================================================
# Identify stromal cells localised within four-cell-diameters of TREM2+ Mac in OC/ON niches (milieu)
# To identify cells within other Mac milieu or other niches, change parameters REFERENCE_STROMAL_SUBSETS or TARGET_NICHES
# This analysis requires the following files:
# 1. stromal_niches_full.csv (from STEP 3)
# 2. cell_polygons_by_core.pkl (from STEP 3)
# 3. cell_meta_full.csv (from STEP 3, used in STEP4-G)
# 4. core-patient_info.csv (created by user)
# 5. patient_group.csv (created by user): This file comprises patient_id and clinical_group
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)

# STEP 4-A: Load Files
import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from tqdm import tqdm
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import re
from IPython.display import display
import matplotlib.pyplot as plt

Tk().withdraw()

stromal_path = askopenfilename(title="Select stromal_niches_full.csv")
stromal_df = pd.read_csv(stromal_path)
print(f"✅ Loaded stromal_niches_full: {stromal_df.shape[0]} rows")

cell_poly_path = askopenfilename(title="Select cell_polygons_by_core.pkl")
with open(cell_poly_path, "rb") as f:
    cell_polygons_by_core = pickle.load(f)
print("✅ Loaded cell_polygons_by_core.pkl")

core_info_path = askopenfilename(title="Select core_patient_info.csv")
core_info = pd.read_csv(core_info_path).drop_duplicates(subset=["fov"])

cols_to_drop = [c for c in stromal_df.columns if c in ["core", "patient_id", "clinical_group"] and c != "fov"]
if cols_to_drop:
    stromal_df = stromal_df.drop(columns=cols_to_drop)

stromal_df = stromal_df.merge(core_info[["fov", "core", "patient_id", "clinical_group"]], on="fov", how="left")
print(f"✅ Loaded core_patient_info: {core_info.shape}")

patient_group_path = askopenfilename(title="Select patient_group.csv")
temp_df = pd.read_csv(patient_group_path)
clean_cols = ~temp_df.columns.str.contains("^Unnamed")
patient_group = temp_df.loc[:, clean_cols]
print(f"✅ Loaded patient_group: {patient_group.shape}")

# STEP 4-B: Parameters
CORE_COL, CELL_COL = "core", "cell"
REFERENCE_STROMAL_SUBSETS = ["TREM2.Mac"]
TARGET_NICHES = ["OC", "ON"]
MEDIAN_STROMAL_AREA = 6274.5
CELL_DIST_PX = int(np.sqrt(MEDIAN_STROMAL_AREA / np.pi))
DIST_THRESHOLD = 4 * CELL_DIST_PX
TMA_ID = "TMA1"

OUTPUT_DIR = Path.home() / "Desktop" / f"{TMA_ID}_STEP4_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# STEP 4-C: Filter stromal cells in OC/ON & Define Milieu
stromal_df = stromal_df[stromal_df["final_niche"].isin(TARGET_NICHES)].copy()

required_cols = ["core", "patient_id", "clinical_group"]
missing = [c for c in required_cols if c not in stromal_df.columns]
if missing: raise ValueError(f"❌ Missing required columns: {missing}")

missing_cores = set(stromal_df["core"].unique()) - set(cell_polygons_by_core.keys())
print(f"⚠️ Missing cores in polygon dict: {len(missing_cores)}")

stromal_df["is_reference"] = stromal_df["stromal_subset"].isin(REFERENCE_STROMAL_SUBSETS)
stromal_df["inside_milieu"] = False

for core, df_core in tqdm(stromal_df.groupby(CORE_COL), desc="STEP10 distance"):
    poly_dict = cell_polygons_by_core.get(core, {})
    ref_ids = df_core[df_core["is_reference"]][CELL_COL].tolist()
    ref_polys = [poly_dict[cid] for cid in ref_ids if cid in poly_dict and not poly_dict[cid].is_empty]
    
    if len(ref_polys) == 0: continue

    for idx, row in df_core.iterrows():
        poly = poly_dict.get(row[CELL_COL])
        if poly is None or poly.is_empty: continue
        min_dist = min(poly.distance(rp) for rp in ref_polys)
        if min_dist <= DIST_THRESHOLD:
            stromal_df.at[idx, "inside_milieu"] = True

stromal_df["condition"] = np.where(stromal_df["inside_milieu"], "Inside_milieu", "Outside_milieu")

# STEP 4-D: Merge metadata
stromal_df = stromal_df.merge(patient_group, on="patient_id", how="left", suffixes=("", "_pg"))

if "clinical_group_pg" in stromal_df.columns:
    stromal_df["clinical_group"] = stromal_df["clinical_group"].fillna(stromal_df["clinical_group_pg"])
    stromal_df.drop(columns=["clinical_group_pg"], inplace=True)

if "patient_id_x" in stromal_df.columns: stromal_df["patient_id"] = stromal_df["patient_id_x"]
elif "patient_id_y" in stromal_df.columns: stromal_df["patient_id"] = stromal_df["patient_id_y"]

stromal_df = stromal_df.drop(columns=[c for c in stromal_df.columns if c.endswith("_x") or c.endswith("_y")], errors="ignore")

print(f"⚠️ Missing patient_id: {stromal_df['patient_id'].isna().sum()}")
print(f"⚠️ Missing clinical_group: {stromal_df['clinical_group'].isna().sum()}")
print(f"⚠️ Missing fov: {stromal_df['fov'].isna().sum()}")

# STEP 4-E: Summary table
counts_long = stromal_df.groupby(["core", "patient_id", "clinical_group", "stromal_subset", "condition"]).size().reset_index(name="count")
counts_wide = counts_long.pivot_table(index=["core", "patient_id", "clinical_group", "condition"], columns="stromal_subset", values="count", fill_value=0).reset_index()

inside_df = counts_wide[counts_wide["condition"] == "Inside_milieu"].drop(columns=["condition"]).copy()
outside_df = counts_wide[counts_wide["condition"] == "Outside_milieu"].drop(columns=["condition"]).copy()

valid_cores = inside_df["core"].unique()
inside_df = inside_df[inside_df["core"].isin(valid_cores)]
outside_df = outside_df[outside_df["core"].isin(valid_cores)]

all_cols = inside_df.columns.union(outside_df.columns)
inside_df = inside_df.reindex(columns=all_cols, fill_value=0)
outside_df = outside_df.reindex(columns=all_cols, fill_value=0)

meta_cols = ["core", "patient_id", "clinical_group"]
subset_order = ["TREM2.Mac", "FOLR2.Mac", "IL4I1.Mac", "SPP1.Mac", "C.Mo", "cDC1", "cDC2", "CCR7.DC", "Mast", "CD4.T", "CD8.T", "Treg", "Naive.B", "Memory.B", "Plasma", "NK", "Myoepithelial", "iCAF", "myCAF", "CD34.CAF", "Pericyte", "Blood.EC", "Lymph.EC", "Adipocyte"]

for col in subset_order:
    if col not in inside_df.columns: inside_df[col] = 0
    if col not in outside_df.columns: outside_df[col] = 0

inside_df = inside_df[meta_cols + subset_order]
outside_df = outside_df[meta_cols + subset_order]
subset_cols = subset_order.copy()

merge_keys = ["core", "patient_id", "clinical_group"]
total_df = pd.merge(inside_df, outside_df, on=merge_keys, how="inner", suffixes=("_in", "_out")).fillna(0)

for col in subset_cols:
    total_df[col] = total_df[f"{col}_in"] + total_df[f"{col}_out"]
total_df = total_df[merge_keys + subset_cols]

prop_df = pd.merge(inside_df, total_df, on=merge_keys, how="left", suffixes=("_in", "_total"))
for col in subset_cols:
    denom = prop_df[f"{col}_total"].replace(0, np.nan)
    prop_df[col] = prop_df[f"{col}_in"] / denom
prop_df = prop_df.fillna(0)[merge_keys + subset_cols]

inside_df["Milieu"], outside_df["Milieu"], total_df["Milieu"], prop_df["Milieu"] = "Inside", "Outside", "IN+OUT", "prop_inside"
final_df = pd.concat([inside_df, outside_df, total_df, prop_df], ignore_index=True)
final_df = final_df.rename(columns={"clinical_group": "Group", "patient_id": "Patient", "core": "Core"})

meta_order = ["Milieu", "Group", "Patient", "Core"]
other_cols = [c for c in final_df.columns if c not in meta_order]
final_df = final_df[meta_order + other_cols]

def core_sort_key(core):
    match = re.match(r"([A-Z]+)(\d+)", str(core))
    return (match.group(1), int(match.group(2))) if match else (core, 0)

milieu_order = ["Inside", "Outside", "IN+OUT", "prop_inside"]
final_df["Milieu"] = pd.Categorical(final_df["Milieu"], categories=milieu_order, ordered=True)
final_df = final_df.sort_values(by=["Milieu", "Core", "Group", "Patient"], key=lambda col: col.map(core_sort_key) if col.name == "Core" else col, kind="mergesort")

print("✅ STEP 4-E completed")
print("Shape:", final_df.shape, "Missing values:", final_df.isna().sum().sum())

# STEP 4-F: Save outputs
out_path = OUTPUT_DIR / f"{TMA_ID}_stromal_cell_TREM2_milieu.csv"
final_df.to_csv(out_path, index=False)

print("✅ STEP 4-F completed")
print("Saved file:", out_path)
display(final_df.head())

# STEP 4-G: QC Plot of ALL cells inside/Outside TREM2 Mac milieu
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from tqdm import tqdm
import numpy as np
import pandas as pd
from pathlib import Path

Tk().withdraw()

print("Please select cell_meta_full.csv (exported from STEP 3)")
cell_meta_path = askopenfilename(title="Select cell_meta_full CSV")
cell_meta_full = pd.read_csv(cell_meta_path)
print(f"✅ Loaded cell_meta_full: {cell_meta_full.shape[0]} cells")

if "core" not in cell_meta_full.columns or cell_meta_full["core"].isna().any():
    print("⚠️ Core info missing in cell_meta_full. Loading core_patient_info...")
    core_info_path = askopenfilename(title="Select core_patient_info CSV")
    core_info_qc = pd.read_csv(core_info_path).drop_duplicates(subset=["fov"])
    
    FOV_COL = "fov"
    cell_meta_full[FOV_COL] = pd.to_numeric(cell_meta_full[FOV_COL], errors="coerce").astype("Int64").astype(str)
    core_info_qc[FOV_COL] = pd.to_numeric(core_info_qc[FOV_COL], errors="coerce").astype("Int64").astype(str)
    
    cols_to_drop = [c for c in cell_meta_full.columns if c in core_info_qc.columns and c != FOV_COL]
    cell_meta_full = cell_meta_full.merge(core_info_qc.drop(columns=cols_to_drop, errors="ignore"), on=FOV_COL, how="left")
    print("✅ Metadata merged with core info")
else:
    print("✅ Core info already present in cell_meta_full")

CORE_COL = "core"
assert cell_meta_full[CORE_COL].isna().sum() == 0, "❌ Missing core assignments in cell_meta_full"

if "condition" not in stromal_df.columns:
    print("🔄 Re-calculating milieu condition (variable not found)...")
    TARGET_NICHES = ["OC", "ON"]
    REFERENCE_STROMAL_SUBSETS = ["TREM2.Mac"]
    MEDIAN_STROMAL_AREA = 6274.5
    CELL_DIST_PX = int(np.sqrt(MEDIAN_STROMAL_AREA / np.pi))
    DIST_THRESHOLD = 4 * CELL_DIST_PX
    
    stromal_df_qc = stromal_df[stromal_df["final_niche"].isin(TARGET_NICHES)].copy()
    stromal_df_qc["is_reference"] = stromal_df_qc["stromal_subset"].isin(REFERENCE_STROMAL_SUBSETS)
    stromal_df_qc["inside_milieu"] = False
    
    for core, df_core in tqdm(stromal_df_qc.groupby(CORE_COL), desc="Computing distances"):
        poly_dict = cell_polygons_by_core.get(core, {})
        ref_ids = df_core[df_core["is_reference"]]["cell"].tolist()
        ref_polys = [poly_dict[cid] for cid in ref_ids if cid in poly_dict and not poly_dict[cid].is_empty]
        if len(ref_polys) == 0: continue
        for idx, row in df_core.iterrows():
            poly = poly_dict.get(row["cell"])
            if poly is None or poly.is_empty: continue
            if min(poly.distance(rp) for rp in ref_polys) <= DIST_THRESHOLD:
                stromal_df_qc.at[idx, "inside_milieu"] = True
    stromal_df_qc["condition"] = np.where(stromal_df_qc["inside_milieu"], "Inside_milieu", "Outside_milieu")
else:
    print("✅ Using existing 'stromal_df' with pre-calculated conditions")
    stromal_df_qc = stromal_df.copy()

TMA_ID = "TMA1"
OUT_DIR = Path.home() / "Desktop" / f"{TMA_ID}_QC_ALL_STROMAL"
OUT_DIR.mkdir(exist_ok=True)

COLORS = {
    "outline": "#000000", 
    "tumor": "#90EE90", 
    "outside": "#87CEEB", 
    "inside": "#0000FF", 
    "reference": "#FF0000"
}

print(f"Generating QC plots for {stromal_df_qc[CORE_COL].nunique()} cores...")

for core in tqdm(stromal_df_qc[CORE_COL].unique(), desc="QC plotting"):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.axis("off")
    
    poly_dict = cell_polygons_by_core.get(core, {})
    
    for poly in poly_dict.values():
        if not poly.is_empty:
            x, y = poly.exterior.xy
            ax.plot(x, y, color=COLORS["outline"], lw=0.3, alpha=0.3, zorder=1)

    if core in core_to_nests:
        for nest in core_to_nests[core]:
            x, y = nest["poly"].exterior.xy
            ax.fill(x, y, color=COLORS["tumor"], alpha=0.2, zorder=2)
            for h in nest["holes"]:
                if not h.is_empty:
                    xh, yh = h.exterior.xy
                    ax.fill(xh, yh, color="#808080", alpha=1.0, zorder=4)

    df_core = stromal_df_qc[stromal_df_qc[CORE_COL] == core]
    for _, row in df_core.iterrows():
        if row["is_reference"]: continue
        poly = poly_dict.get(row["cell"])
        if poly is None or poly.is_empty: continue
        color = COLORS["inside"] if row["condition"] == "Inside_milieu" else COLORS["outside"]
        x, y = poly.exterior.xy
        ax.fill(x, y, color=color, alpha=0.7, zorder=3)

    ref_cells = df_core[df_core["is_reference"]]
    for cid in ref_cells["cell"]:
        poly = poly_dict.get(cid)
        if poly is not None and not poly.is_empty:
            x, y = poly.exterior.xy
            ax.fill(x, y, color=COLORS["reference"], alpha=1.0, zorder=5)

    plt.title(f"{core} — All Stromal Cells")
    plt.savefig(OUT_DIR / f"{core}.png", dpi=300, bbox_inches="tight")
    plt.close()

print(f"✅ QC plots saved to: {OUT_DIR}")

# ==============================================================================
# STEP 5 — Determine Distances of Rec+ Epi from Closest Lig+ TREM2 Mac
# ==============================================================================
# This analysis requires the following files:
# 1. cell_meta_full.csv (from STEP3)
# 2. stromal_niche_full.csv (from STEP3)
# 3. core_to_nests.pkl (from STEP3)
# 4. cell_polygons_by_core.pkl (from STEP3) 
# 5. exprMat.csv (from AtoMx)

# Candidate ligand-receptor pairs were predefined by analysing snRNA-seq data (R code)
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)
# The output data are used to plot distribution densities via another script (R code) 

# STEP 5-A: Load required data files
import pandas as pd
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from tqdm import tqdm
from pathlib import Path
import pickle
from shapely.ops import nearest_points
import os

Tk().withdraw()

cell_meta_path = askopenfilename(title="Select cell_meta_full CSV file")
cell_meta_full = pd.read_csv(cell_meta_path)
print(f"✅ Loaded cell_meta_full: {cell_meta_full.shape[0]} cells")

stromal_niches_path = askopenfilename(title="Select stromal_niches_full CSV file")
stromal_niches_full = pd.read_csv(stromal_niches_path)
print(f"✅ Loaded stromal_niches_full: {stromal_niches_full.shape[0]} rows")

cell_polygons_path = askopenfilename(title="Select cell_polygons_by_core.pkl file")
with open(cell_polygons_path, "rb") as f:
    cell_polygons_by_core = pickle.load(f)
print("✅ Loaded cell_polygons_by_core.pkl")

core_to_nests_path = askopenfilename(title="Select core_to_nests.pkl file")
with open(core_to_nests_path, "rb") as f:
    core_to_nests = pickle.load(f)
print("✅ Loaded core_to_nests.pkl")

expr_path = askopenfilename(title="Select COUNT MATRIX (exprMat) CSV file", filetypes=[("CSV files", "*.csv")])
expr_mat = pd.read_csv(expr_path)

if "cell_ID" not in expr_mat.columns or "fov" not in expr_mat.columns:
    raise ValueError("exprMat must contain 'cell_ID' and 'fov'")

expr_mat["cell_ID"] = expr_mat["cell_ID"].astype(str)
expr_mat["fov"] = expr_mat["fov"].astype(str)
expr_mat["cell_key"] = expr_mat["fov"] + "_" + expr_mat["cell_ID"]
expr_mat = expr_mat.set_index("cell_key").drop(columns=["cell_ID", "fov"], errors="ignore")
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce").fillna(0).astype("int32")

print(f"✅ Expression matrix: {expr_mat.shape}")
print("✅ STEP5-A completed")

# STEP 5-B: Rebuild metadata & Alignment
CORE_COL = "core_x" if "core_x" in cell_meta_full.columns else "core"
cell_meta_full["core"] = cell_meta_full[CORE_COL]

def parse_cell_id(x):
    parts = str(x).split("_")
    return (parts[1], parts[2], parts[3]) if len(parts) == 4 else (None, None, None)

parsed = cell_meta_full["cell_id"].astype(str).apply(parse_cell_id)
cell_meta_full["slide"] = parsed.str[0]
cell_meta_full["fov_parsed"] = parsed.str[1]
cell_meta_full["cell_ID_parsed"] = parsed.str[2]
cell_meta_full["cell_key"] = cell_meta_full["fov_parsed"].astype(str) + "_" + cell_meta_full["cell_ID_parsed"].astype(str)
cell_meta_full["cell"] = cell_meta_full["cell"].astype(str)

print("Using core column:", CORE_COL)
print("Unique cores:", cell_meta_full["core"].nunique())

matched = cell_meta_full["cell_key"].isin(expr_mat.index).sum()
print(f"Match rate: {matched/len(cell_meta_full):.2%}")
print("✅ STEP5-B completed")

# STEP 5-C: Define parameters
TMA_ID = "TMA1"
CELL_DIST_PX = 44
DIST_THRESHOLD = 4 * CELL_DIST_PX
PX_TO_UM = 0.12  
LR_PAIRS = [("SPP1", "CD44"), ("LGALS9", "CD44"), ("GRN", "SORT1"), ("SEMA4C", "PLXNB2"), ("F11R", "F11R"), ("CD99", "CD99")]

print("✅ STEP5-C completed")

# STEP 5-D: Define ligand positive TREM2 Mac
trem2_mac = cell_meta_full.loc[cell_meta_full["stromal_subset"].eq("TREM2.Mac")].copy()
trem2_mac["cell_key"] = trem2_mac["cell_key"].astype(str)
trem2_mac = trem2_mac[trem2_mac["cell_key"].isin(expr_mat.index)].copy().set_index("cell_key")
print(f"Total TREM2.Mac (aligned): {trem2_mac.shape[0]}")

expr_mat.columns = expr_mat.columns.str.strip()
expr_mat.index = expr_mat.index.astype(str)
expr_bin = expr_mat.astype(bool)

lig_cols = []
for lig, rec in LR_PAIRS:
    col_name = f"lig_pos_{lig}"
    lig_cols.append(col_name)
    if lig not in expr_bin.columns:
        trem2_mac[col_name] = False
    else:
        trem2_mac[col_name] = expr_bin.loc[trem2_mac.index, lig].to_numpy()

print("\n===== STEP5-D SANITY CHECK =====")
for lig, rec in LR_PAIRS:
    col = f"lig_pos_{lig}"
    if col in trem2_mac.columns:
        pos = trem2_mac[col].sum()
        pct = pos / len(trem2_mac) if len(trem2_mac) > 0 else 0
        print(f"{col}: {pos} ({pct:.2%})")
print("NaN check:", trem2_mac[lig_cols].isna().sum().sum())
print("===== END =====\n")
print("✅ STEP5-D completed")

# STEP 5-E: Define receptor positive epithelial cells in tumour nests
epi_cells = cell_meta_full.loc[cell_meta_full["is_epithelial"].eq(True)].copy()
epi_cells["cell_key"] = epi_cells["cell_key"].astype(str)
epi_cells = epi_cells[epi_cells["cell_key"].isin(expr_mat.index)].copy().set_index("cell_key")
print(f"Epithelial cells (aligned): {len(epi_cells)}")

epi_cells["epi_in_tumour_nest"] = False
for core, df_core in epi_cells.groupby(CORE_COL):
    if core not in core_to_nests: continue
    nests = core_to_nests[core]
    poly_dict = cell_polygons_by_core.get(core, {})
    for cell_key, row in df_core.iterrows():
        poly = poly_dict.get(row["cell"])
        if poly is None or poly.is_empty: continue
        centroid = poly.centroid
        for nest in nests:
            if nest["shell"].covers(centroid):
                epi_cells.at[cell_key, "epi_in_tumour_nest"] = True
                break

epi_nest = epi_cells[epi_cells["epi_in_tumour_nest"]].copy()
print("Epithelial in tumor nests:", len(epi_nest))

epi_nest["epi_state"] = np.select(
    [epi_nest["is_normal_epi"].eq(True), epi_nest["is_tumor_epi"].eq(True)],
    ["Normal.Epi", "Tumour.Epi"],
    default="Other"
)
epi_nest = epi_nest[epi_nest["epi_state"] != "Other"]

expr_bin = expr_mat.astype(bool)
rec_cols = []

for lig, rec in LR_PAIRS:
    col_name = f"rec_pos_{lig}"
    rec_cols.append(col_name)
    genes = rec.split("_") if isinstance(rec, str) else [rec]
    genes = [g for g in genes if g in expr_bin.columns]
    
    if len(genes) == 0:
        epi_nest[col_name] = False
        continue
    
    valid_cells = epi_nest.index.intersection(expr_bin.index)
    expr_sub = expr_bin.loc[valid_cells, genes]
    result = expr_sub.all(axis=1) if len(genes) > 1 else expr_sub.iloc[:, 0]
    
    epi_nest[col_name] = False
    epi_nest.loc[valid_cells, col_name] = result.to_numpy()

print("\n===== STEP5-E SANITY CHECK =====")
normal_nest = epi_nest.loc[epi_nest["epi_state"].eq("Normal.Epi")]
tumor_nest  = epi_nest.loc[epi_nest["epi_state"].eq("Tumour.Epi")]
print("Normal.Epi in nest:", len(normal_nest))
print("Tumour.Epi in nest:", len(tumor_nest))

for lig, rec in LR_PAIRS:
    col = f"rec_pos_{lig}"
    rec_label = rec.replace("_", "–") if isinstance(rec, str) else str(rec)
    if col in epi_nest.columns:
        print(f"{rec_label} ({lig}) - Normal: {normal_nest[col].sum()}, Tumor: {tumor_nest[col].sum()}")
print("===== END =====\n")
print("✅ STEP5-E completed")

# STEP 5-F: Compute distance of Rec+ Epi from closest Lig+ TREM2.Mac
results = []

for core, epi_df in epi_nest.groupby(CORE_COL):
    if core not in cell_polygons_by_core: continue
    poly_dict = cell_polygons_by_core[core]
    
    trem2_core = trem2_mac[trem2_mac[CORE_COL] == core].copy()
    trem2_core["poly"] = trem2_core["cell"].map(poly_dict)
    trem2_core = trem2_core.dropna(subset=["poly"])

    for lig, rec in LR_PAIRS:
        lig_cells = trem2_core[trem2_core[f"lig_pos_{lig}"]]
        if lig_cells.empty: continue
        
        lig_polys = [p for p in lig_cells["poly"] if p is not None and not p.is_empty]
        if len(lig_polys) == 0: continue

        epi_subset = epi_df[epi_df[f"rec_pos_{lig}"]]
        for _, row in epi_subset.iterrows():
            epi_poly = poly_dict.get(row["cell"])
            if epi_poly is None or epi_poly.is_empty: continue
            
            min_dist = min(epi_poly.distance(lp) for lp in lig_polys) * PX_TO_UM
            results.append({
                "TMA_id": TMA_ID, "core": core, "cell_id": row["cell"],
                "epi_state": row["epi_state"], "LR_pair": f"{lig}-{rec}",
                "min_distance_um": min_dist
            })

distance_df = pd.DataFrame(results)
print("Distance table:", distance_df.shape)
print(distance_df["min_distance_um"].describe())
print("✅ STEP5-F completed")

# STEP 5-G: Export data file
outdir = Path.home() / "Desktop" / f"{TMA_ID}_STEP5_results"
os.makedirs(outdir, exist_ok=True)

outpath = os.path.join(outdir, f"{TMA_ID}_LR_distance.csv")
distance_df.to_csv(outpath, index=False)

print("✅ Saved:", outpath)
print("✅ STEP5 completed")

# ==============================================================================
# STEP 6 — Compute Percent of Rec+ Epi within 40 µm of Lig+ TREM2 Mac
# ==============================================================================
# This analysis requires the following files:
# 1. LR_distance.csv (from STEP5 for all TMAs)
# 2. core_patient_group.csv (created by user)
# The core_patient_group includes TMA ID, Core ID, Patient ID, Clinical Group across all TMAs

# STEP 6-A: Load required data files
import pandas as pd
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from pathlib import Path
import os

Tk().withdraw()

NUM_FILES = 6
THRESHOLD_UM = 40

df_list = []
for i in range(NUM_FILES):
    file_path = askopenfilename(title=f"Select LR_distance CSV {i+1}/{NUM_FILES}", filetypes=[("CSV files", "*.csv")])
    if not file_path: raise ValueError("Missing input file")
    df = pd.read_csv(file_path)
    df["source_file"] = Path(file_path).name
    df_list.append(df)
    print(f"Loaded: {file_path}")

all_df = pd.concat(df_list, ignore_index=True)
print("\nCombined shape:", all_df.shape)

meta_path = askopenfilename(title="Select core_patient_group metadata CSV", filetypes=[("CSV files", "*.csv")])
meta_df = pd.read_csv(meta_path)

required = ["TMA_id", "core", "patient_id", "clinical_group"]
missing = [c for c in required if c not in meta_df.columns]
if missing: raise ValueError(f"Missing columns: {missing}")

all_df = all_df.merge(meta_df, on=["TMA_id", "core"], how="left")
print("After merge:", all_df.shape)

# STEP 6-B: Define Epi type and metric
all_df["epi_state"] = "Total.Epi"
all_df["within_40um"] = all_df["min_distance_um"] <= THRESHOLD_UM

# STEP 6-C: Per-core analysis
core_stats = all_df.groupby(["TMA_id", "core", "LR_pair"]).apply(lambda x: pd.Series({"N_cells": len(x), "Pct_within_40um": x["within_40um"].mean() * 100, "Mean_distance_um": x["min_distance_um"].mean(), "Median_distance_um": x["min_distance_um"].median()})).reset_index()
print("\nSTEP6-C (per core):")
print(core_stats.head())

# STEP 6-D: Per-patient analysis
patient_stats = all_df.groupby(["patient_id", "clinical_group", "LR_pair"]).apply(lambda x: pd.Series({"N_cells": len(x), "Pct_within_40um": x["within_40um"].mean() * 100, "Mean_distance_um": x["min_distance_um"].mean(), "Median_distance_um": x["min_distance_um"].median()})).reset_index()
print("\nSTEP6-D (per patient):")
print(patient_stats.head())

# STEP 6-E: Export data
outdir = Path.home() / "Desktop" / "STEP6_Pct40_results"
os.makedirs(outdir, exist_ok=True)

core_stats.to_csv(outdir / "STEP6_per_core_Pct40.csv", index=False)
patient_stats.to_csv(outdir / "STEP6_per_patient_Pct40.csv", index=False)

print("\nSaved to:")
print(outdir)
print("DONE.")

# ==============================================================================
# STEP 7 — Visualise Ligand/Receptor Gene Expression per Core
# ==============================================================================
# This analysis requires the following files:
# 1. cell_meta_full.csv (from STEP3)
# 2. stromal_niche_full.csv (from STEP3)
# 3. cell_polygons_by_core (from STEP3)
# 4. core_to_nests.pkl (from STEP3)
# 5. exprMat.csv (from AtoMx)
# Candidate ligand-receptor pairs were predefined by analysing snRNA-seq data (R code)
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)

# STEP 7-A: Load required data files
import pandas as pd
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from tqdm import tqdm
from pathlib import Path
import pickle
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

Tk().withdraw()

cell_meta_path = askopenfilename(title="Select cell_meta_full CSV file")
cell_meta_full = pd.read_csv(cell_meta_path)
print(f"✅ Loaded cell_meta_full: {cell_meta_full.shape[0]} cells")

stromal_niches_path = askopenfilename(title="Select stromal_niches_full CSV file")
stromal_niches_full = pd.read_csv(stromal_niches_path)
print(f"✅ Loaded stromal_niches_full: {stromal_niches_full.shape[0]} rows")

poly_path = askopenfilename(title="Select cell_polygons_by_core.pkl file")
with open(poly_path, "rb") as f: cell_polygons_by_core = pickle.load(f)
print("✅ Loaded cell_polygons_by_core.pkl")

nests_path = askopenfilename(title="Select core_to_nests.pkl file")
with open(nests_path, "rb") as f: core_to_nests = pickle.load(f)
print("✅ Loaded core_to_nests.pkl")

expr_path = askopenfilename(title="Select COUNT MATRIX (exprMat) CSV file", filetypes=[("CSV files", "*.csv")])
expr_mat = pd.read_csv(expr_path)
if "cell_ID" not in expr_mat.columns or "fov" not in expr_mat.columns: raise ValueError("exprMat must contain 'cell_ID' and 'fov'")

expr_mat["cell_ID"] = expr_mat["cell_ID"].astype(str)
expr_mat["fov"] = expr_mat["fov"].astype(str)
expr_mat["cell_key"] = expr_mat["fov"] + "_" + expr_mat["cell_ID"]
expr_mat = expr_mat.set_index("cell_key").drop(columns=["cell_ID", "fov"], errors="ignore")
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce").fillna(0).astype("int32")

print(f"✅ Expression matrix: {expr_mat.shape}")
print("✅ STEP7-A completed")

# STEP 7-B: Develop common key system
CORE_COL = "core_x" if "core_x" in cell_meta_full.columns else "core"
cell_meta_full["core"] = cell_meta_full[CORE_COL]

def parse_cell_id(x):
    parts = str(x).split("_")
    return (parts[1], parts[2], parts[3]) if len(parts) == 4 else (None, None, None)

parsed = cell_meta_full["cell_id"].astype(str).apply(parse_cell_id)
cell_meta_full["slide"] = parsed.str[0]
cell_meta_full["fov_parsed"] = parsed.str[1]
cell_meta_full["cell_ID_parsed"] = parsed.str[2]
cell_meta_full["cell_key"] = cell_meta_full["fov_parsed"].astype(str) + "_" + cell_meta_full["cell_ID_parsed"].astype(str)
cell_meta_full["cell"] = cell_meta_full["cell"].astype(str)

print("Matched metadata cells:", cell_meta_full["cell_key"].isin(expr_mat.index).sum())

meta_keys = set(cell_meta_full["cell_key"].dropna().astype(str))
expr_keys = set(expr_mat.index.astype(str))
meta_expr_overlap = len(meta_keys & expr_keys)

print(f"[Sanity] Overlap rate: {meta_expr_overlap / max(len(meta_keys), 1)}")
if meta_expr_overlap == 0: raise ValueError("❌ Metadata does NOT match expr_mat → key parsing incorrect")

print("✅ STEP7-B completed")

# STEP 7-C: Define parameters
TMA_ID = "TMA1"
genes_to_plot = ["SPP1", "CD44", "GRN", "SORT1", "LGALS9", "F11R", "CD99"]
CELL_DIST_PX = 44
DIST_THRESHOLD = 4 * CELL_DIST_PX

OUT_DIR = Path.home() / "Desktop" / f"{TMA_ID}_STEP7_results"
(OUT_DIR / "Panel1").mkdir(parents=True, exist_ok=True)
(OUT_DIR / "Panel2").mkdir(parents=True, exist_ok=True)

global_min = 0
global_max = np.percentile(expr_mat.values.flatten(), 99)
norm_expr = Normalize(vmin=global_min, vmax=global_max)
cmap_expr = plt.cm.Blues

print("✅ STEP7-C completed")

# STEP 7-D: Define epithelial cell states
CLUSTER_COL = "RNA_TMA1_20260409_TK_Cell.Typing.InSituType.1_1_clusters"
meta_lookup = cell_meta_full.set_index("cell")
stromal_lookup = stromal_niches_full.set_index("cell")

epi_cells = cell_meta_full[cell_meta_full[CLUSTER_COL].isin(["Tumour.Epi", "Normal.Epi"])].copy()
epi_lookup = {}

print("Reconstructing epithelial spatial states...")
for core in tqdm(cell_meta_full["core"].unique()):
    poly_dict = cell_polygons_by_core.get(core, {})
    stromal_core = stromal_niches_full[stromal_niches_full["core"] == core]
    trem2_oc_on = stromal_core[(stromal_core["stromal_subset"] == "TREM2.Mac") & (stromal_core["final_niche"].isin(["OC", "ON"]))]["cell"].tolist()
    trem2_polys = [poly_dict[c] for c in trem2_oc_on if c in poly_dict]
    epi_core = epi_cells[epi_cells["core"] == core]
    
    for _, row in epi_core.iterrows():
        cid, poly = row["cell"], poly_dict.get(row["cell"])
        if poly is None or poly.is_empty: continue
        
        in_nest = any(nest["shell"].covers(poly.centroid) for nest in core_to_nests.get(core, []))
        inside_milieu = False
        if len(trem2_polys) > 0 and in_nest:
            min_dist = min(poly.distance(tp) for tp in trem2_polys)
            inside_milieu = min_dist <= DIST_THRESHOLD
        
        final_label = "Inside_milieu" if (in_nest and inside_milieu) else ("Outside_milieu" if in_nest else "Other_epithelial")
        epi_lookup[cid] = {"in_nest": in_nest, "inside_milieu": inside_milieu, "final_label": final_label}

print("✅ STEP7-D completed")

# STEP 7-E: Plot helper functions
def panel1_fill(cell_type, epi_label, niche):
    if cell_type in ["Tumour.Epi", "Normal.Epi"]: return "#00A651" if epi_label == "Inside_milieu" else ("#7CFC00" if epi_label == "Outside_milieu" else "none")
    if cell_type == "TREM2.Mac": return "#FD5800" if niche in ["OC", "ON"] else "#FF69B4"
    return "none"

def panel2_outline(cell_type, epi_label, niche):
    if cell_type in ["Tumour.Epi", "Normal.Epi"]: return "#00A651" if epi_label == "Inside_milieu" else "#90EE90"
    if cell_type == "TREM2.Mac": return "#FD5800" if niche in ["OC", "ON"] else "#FF69B4"
    return "#B0B0B0"

print("✅ STEP7-E completed")

# STEP 7-F: Plot per core
for core in tqdm(cell_meta_full["core"].unique(), desc="Plotting cores"):
    poly_dict = cell_polygons_by_core.get(core, {})
    if len(poly_dict) == 0: continue

    core_polys = {cid: poly.exterior.xy for cid, poly in poly_dict.items() if poly is not None and not poly.is_empty}
    
    fig1, ax1 = plt.subplots(figsize=(8, 8))
    for cid, (x, y) in core_polys.items():
        if cid not in meta_lookup.index: continue
        cell_type = meta_lookup.at[cid, CLUSTER_COL]
        epi_label = epi_lookup.get(cid, {}).get("final_label", None)
        niche = stromal_lookup.at[cid, "final_niche"] if cid in stromal_lookup.index else None
        ax1.fill(x, y, facecolor=panel1_fill(cell_type, epi_label, niche), edgecolor="black", linewidth=0.2)
    
    ax1.set_title(f"{core} — Cell type / niche map")
    ax1.set_aspect("equal")
    ax1.axis("off")
    fig1.savefig(OUT_DIR / "Panel1" / f"{core}_panel1.png", dpi=300, bbox_inches="tight")
    plt.close(fig1)

    valid_cells = [cid for cid in core_polys.keys() if cid in meta_lookup.index]
    if len(valid_cells) == 0: continue
    
    for gene in genes_to_plot:
        if gene not in expr_mat.columns: continue
        
        fig2, ax2 = plt.subplots(figsize=(8, 8))
        draw_list = []
        
        for cid in valid_cells:
            cell_key = meta_lookup.at[cid, "cell_key"]
            if cell_key not in expr_mat.index: continue
            
            x, y = core_polys[cid]
            cell_type = meta_lookup.at[cid, CLUSTER_COL]
            epi_label = epi_lookup.get(cid, {}).get("final_label", None)
            niche = stromal_lookup.at[cid, "final_niche"] if cid in stromal_lookup.index else None
            
            expr_val = expr_mat.at[cell_key, gene]
            fill_color = cmap_expr(norm_expr(expr_val))
            edge = panel2_outline(cell_type, epi_label, niche)
            
            if cell_type == "TREM2.Mac": z = 4 if niche in ["OC", "ON"] else 3
            elif cell_type in ["Tumour.Epi", "Normal.Epi"]: z = 2 if epi_label == "Inside_milieu" else 1
            else: z = 0
            
            draw_list.append((z, x, y, fill_color, edge))
        
        draw_list.sort(key=lambda t: t[0])
        
        for z, x, y, fill_color, edge in draw_list:
            ax2.fill(x, y, facecolor=fill_color, edgecolor="none", antialiased=False, zorder=z)
            ax2.fill(x, y, facecolor="none", edgecolor=edge, linewidth=0.7, antialiased=False, joinstyle="miter", zorder=z + 0.1)
        
        ax2.set_aspect("equal", adjustable="box")
        ax2.autoscale()
        ax2.set_title(f"{core} — {gene}")
        ax2.axis("off")
        
        sm = ScalarMappable(norm=norm_expr, cmap=cmap_expr)
        sm.set_array([])
        fig2.colorbar(sm, ax=ax2, fraction=0.046, pad=0.04)
        
        if len(draw_list) > 0: fig2.savefig(OUT_DIR / "Panel2" / f"{core}_{gene}.png", dpi=300, bbox_inches="tight")
        else: print(f"[Panel2] WARNING: no cells plotted for {core} | {gene}")
        
        plt.close(fig2)

print("✅ STEP7 completed successfully")

# ==============================================================================
# STEP 8 — Count Number of Lig+ TREM2 Mac in OC/ON paired with Rec+ Epi
# ==============================================================================
# This analysis requires the following files:
# 1. cell_meta_full.csv (from STEP3)
# 2. stromal_niche_full.csv (from STEP3)
# 3. cell_polygons_by_core (from STEP3)
# 4. core_to_nests.pkl (from STEP3)
# 5. exprMat.csv (from AtoMx)
# Candidate ligand-receptor pairs were predefined by analysing snRNA-seq data (R code)
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)

# STEP 8-A: Load required data files
import pandas as pd
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from tqdm import tqdm
from pathlib import Path
import pickle
from shapely.ops import nearest_points
import os
from IPython.display import display

Tk().withdraw()

cell_meta_path = askopenfilename(title="Select cell_meta_full CSV file")
cell_meta_full = pd.read_csv(cell_meta_path)
print(f"✅ Loaded cell_meta_full: {cell_meta_full.shape[0]} cells")

stromal_niches_path = askopenfilename(title="Select stromal_niches_full CSV file")
stromal_niches_full = pd.read_csv(stromal_niches_path)
print(f"✅ Loaded stromal_niches_full: {stromal_niches_full.shape[0]} rows")

cell_polygons_path = askopenfilename(title="Select cell_polygons_by_core.pkl file")
with open(cell_polygons_path, "rb") as f: cell_polygons_by_core = pickle.load(f)
print("✅ Loaded cell_polygons_by_core.pkl")

core_to_nests_path = askopenfilename(title="Select core_to_nests.pkl file")
with open(core_to_nests_path, "rb") as f: core_to_nests = pickle.load(f)
print("✅ Loaded core_to_nests.pkl")

expr_path = askopenfilename(title="Select COUNT MATRIX (exprMat) CSV file", filetypes=[("CSV files", "*.csv")])
expr_mat = pd.read_csv(expr_path)

if "cell_ID" not in expr_mat.columns or "fov" not in expr_mat.columns:
    raise ValueError("exprMat must contain 'cell_ID' and 'fov'")

expr_mat["cell_ID"] = expr_mat["cell_ID"].astype(str)
expr_mat["fov"] = expr_mat["fov"].astype(str)
expr_mat["cell_key"] = expr_mat["fov"] + "_" + expr_mat["cell_ID"]
expr_mat = expr_mat.set_index("cell_key").drop(columns=["cell_ID", "fov"], errors="ignore")
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce").fillna(0).astype("int32")

print(f"✅ Expression matrix: {expr_mat.shape}")
print("✅ STEP8-A completed")

# STEP 8-B: Rebuild metadata
CORE_COL = "core_x" if "core_x" in cell_meta_full.columns else "core"
cell_meta_full["core"] = cell_meta_full[CORE_COL]

def parse_cell_id(x):
    parts = str(x).split("_")
    return (parts[1], parts[2], parts[3]) if len(parts) == 4 else (None, None, None)

parsed = cell_meta_full["cell_id"].astype(str).apply(parse_cell_id)
cell_meta_full["slide"] = parsed.str[0]
cell_meta_full["fov_parsed"] = parsed.str[1]
cell_meta_full["cell_ID_parsed"] = parsed.str[2]
cell_meta_full["cell_key"] = cell_meta_full["fov_parsed"].astype(str) + "_" + cell_meta_full["cell_ID_parsed"].astype(str)
cell_meta_full["cell"] = cell_meta_full["cell"].astype(str)

cell_meta_full = cell_meta_full.merge(stromal_niches_full[["cell", "final_niche"]], on="cell", how="left")

matched = cell_meta_full["cell_key"].isin(expr_mat.index).sum()
print(f"Match rate: {matched/len(cell_meta_full):.2%}")
print("✅ STEP8-B completed")

# STEP 8-C: Define parameters
TMA_ID = "TMA1"
CELL_COL = "cell"
REFERENCE_STROMAL_SUBSETS = ["TREM2.Mac"]
TARGET_NICHES = ["OC", "ON"]
CELL_DIST_PX = 44
DIST_THRESHOLD = 4 * CELL_DIST_PX
PX_TO_UM = 0.12  

LR_PAIRS = [("SPP1", "CD44"), ("LGALS9", "CD44"), ("GRN", "SORT1"), ("F11R", "F11R"), ("CD99", "CD99")]
print("✅ STEP8-C completed")

# STEP 8-D: Filter TREM2 Mac in OC/ON
stromal_df = stromal_niches_full[stromal_niches_full["final_niche"].isin(TARGET_NICHES)].copy()
required_cols = ["core", "patient_id", "clinical_group"]
missing = [c for c in required_cols if c not in stromal_df.columns]
if missing: raise ValueError(f"❌ Missing required columns: {missing}")

missing_cores = set(stromal_df["core"].unique()) - set(cell_polygons_by_core.keys())
print(f"⚠️ Missing cores in polygon dict: {len(missing_cores)}")

stromal_df["is_reference"] = stromal_df["stromal_subset"].isin(REFERENCE_STROMAL_SUBSETS)
stromal_df["inside_milieu"] = False

# STEP 8-E: Define ligand positive TREM2.Mac
trem2_mac = cell_meta_full.loc[cell_meta_full["stromal_subset"].eq("TREM2.Mac")].copy()
trem2_mac["cell_key"] = trem2_mac["cell_key"].astype(str)
trem2_mac = trem2_mac[trem2_mac["cell_key"].isin(expr_mat.index)].copy().set_index("cell_key")
print(f"Total TREM2.Mac (aligned): {trem2_mac.shape[0]}")

expr_mat.columns = expr_mat.columns.str.strip()
expr_mat.index = expr_mat.index.astype(str)
expr_bin = expr_mat.astype(bool)

lig_cols = []
for lig, rec in LR_PAIRS:
    col_name = f"lig_pos_{lig}"
    lig_cols.append(col_name)
    if lig not in expr_bin.columns:
        trem2_mac[col_name] = False
    else:
        trem2_mac[col_name] = expr_bin.loc[trem2_mac.index, lig].to_numpy()

print("\n===== STEP8-E SANITY CHECK =====")
for lig, rec in LR_PAIRS:
    col = f"lig_pos_{lig}"
    if col in trem2_mac.columns:
        pos = trem2_mac[col].sum()
        pct = pos / len(trem2_mac) if len(trem2_mac) > 0 else 0
        print(f"{col}: {pos} ({pct:.2%})")
print("NaN check:", trem2_mac[lig_cols].isna().sum().sum())
print("===== END =====\n")
print("✅ STEP8-E completed")

# STEP 8-F: Define receptor positive epithelial cells in tumour nests
epi_cells = cell_meta_full.loc[cell_meta_full["is_epithelial"].eq(True)].copy()
epi_cells["cell_key"] = epi_cells["cell_key"].astype(str)
epi_cells = epi_cells[epi_cells["cell_key"].isin(expr_mat.index)].copy().set_index("cell_key")
print(f"Epithelial cells (aligned): {len(epi_cells)}")

epi_cells["epi_in_tumour_nest"] = False
for core, df_core in epi_cells.groupby(CORE_COL):
    if core not in core_to_nests: continue
    nests = core_to_nests[core]
    poly_dict = cell_polygons_by_core.get(core, {})
    for cell_key, row in df_core.iterrows():
        poly = poly_dict.get(row["cell"])
        if poly is None or poly.is_empty: continue
        centroid = poly.centroid
        for nest in nests:
            if nest["shell"].covers(centroid):
                epi_cells.at[cell_key, "epi_in_tumour_nest"] = True
                break

epi_nest = epi_cells[epi_cells["epi_in_tumour_nest"]].copy()
print("Epithelial in tumor nests:", len(epi_nest))

epi_nest["epi_state"] = np.select(
    [epi_nest["is_normal_epi"].eq(True), epi_nest["is_tumor_epi"].eq(True)],
    ["Normal.Epi", "Tumour.Epi"],
    default="Other"
)
epi_nest = epi_nest[epi_nest["epi_state"] != "Other"]

expr_bin = expr_mat.astype(bool)
rec_cols = []

for lig, rec in LR_PAIRS:
    col_name = f"rec_pos_{lig}"
    rec_cols.append(col_name)
    genes = rec.split("_") if isinstance(rec, str) else [rec]
    genes = [g for g in genes if g in expr_bin.columns]
    
    if len(genes) == 0:
        epi_nest[col_name] = False
        continue
    
    valid_cells = epi_nest.index.intersection(expr_bin.index)
    expr_sub = expr_bin.loc[valid_cells, genes]
    result = expr_sub.all(axis=1) if len(genes) > 1 else expr_sub.iloc[:, 0]
    
    epi_nest[col_name] = False
    epi_nest.loc[valid_cells, col_name] = result.to_numpy()

print("\n===== STEP8-F SANITY CHECK =====")
normal_nest = epi_nest.loc[epi_nest["epi_state"].eq("Normal.Epi")]
tumor_nest  = epi_nest.loc[epi_nest["epi_state"].eq("Tumour.Epi")]
print("Normal.Epi in nest:", len(normal_nest))
print("Tumour.Epi in nest:", len(tumor_nest))

for lig, rec in LR_PAIRS:
    col = f"rec_pos_{lig}"
    rec_label = rec.replace("_", "–") if isinstance(rec, str) else str(rec)
    if col in epi_nest.columns:
        print(f"{rec_label} ({lig}) - Normal: {normal_nest[col].sum()}, Tumor: {tumor_nest[col].sum()}")
print("===== END =====\n")
print("✅ STEP8-F completed")

# STEP 8-G: Define Lig+ TREM2.Mac associated with Rec+ Epi within 40 µm
PAIR_THRESHOLD_UM = 40
pair_results = []

if "final_niche" not in epi_nest.columns:
    raise ValueError("final_niche missing in epi_nest — merge step not applied")

epi_ocon = epi_nest.copy()

for core in cell_polygons_by_core.keys():
    poly_dict = cell_polygons_by_core[core]
    
    trem2_core = trem2_mac[(trem2_mac[CORE_COL] == core) & (trem2_mac["final_niche"].isin(TARGET_NICHES))].copy()
    if trem2_core.empty: continue
    
    trem2_core["poly"] = trem2_core["cell"].map(poly_dict)
    trem2_core = trem2_core.dropna(subset=["poly"])
    
    epi_core = epi_ocon[epi_ocon[CORE_COL] == core].copy()
    if epi_core.empty: continue
    
    epi_core["poly"] = epi_core["cell"].map(poly_dict)
    epi_core = epi_core.dropna(subset=["poly"])
    
    epi_polys = {row["cell"]: row["poly"] for _, row in epi_core.iterrows()}
    
    for lig, rec in LR_PAIRS:
        lig_subset = trem2_core[trem2_core[f"lig_pos_{lig}"]]
        rec_subset = epi_core[epi_core[f"rec_pos_{lig}"]]
        if lig_subset.empty or rec_subset.empty: continue
        
        rec_polys = {row["cell"]: epi_polys.get(row["cell"]) for _, row in rec_subset.iterrows()}
        
        for _, mac in lig_subset.iterrows():
            mac_poly = mac["poly"]
            if mac_poly is None: continue
            
            min_dist = np.inf
            for epi_cell, epi_poly in rec_polys.items():
                if epi_poly is None: continue
                d = mac_poly.distance(epi_poly) * PX_TO_UM
                if d < min_dist: min_dist = d
            
            pair_results.append({
                "patient_id": mac["patient_id"], "core": core, "cell": mac["cell"],
                "ligand": lig, "receptor": rec, "paired": min_dist <= PAIR_THRESHOLD_UM,
                "min_distance_um": min_dist
            })

pair_df = pd.DataFrame(pair_results)
print("✅ STEP8-G completed")
print(pair_df.head())

# STEP 8-H: Final counting table
if "pair_df" not in globals(): raise ValueError("Run STEP 8-G first")

summary_rows = []
group_cols = ["patient_id", CORE_COL]
trem2_ocon = trem2_mac[trem2_mac["final_niche"].isin(TARGET_NICHES)].copy()
trem2_all = trem2_mac.copy()

for (patient_id, core), df in trem2_ocon.groupby(group_cols):
    row = {
        "patient_id": patient_id, "TMA_id": TMA_ID, "core": core,
        "Total TREM2 Mac in OC+ON niche": len(df),
        "Total TREM2 Mac across all niches": len(trem2_all[(trem2_all["patient_id"] == patient_id) & (trem2_all[CORE_COL] == core)])
    }
    
    for lig, rec in LR_PAIRS:
        lig_col = f"lig_pos_{lig}"
        row[f"{lig}+ TREM2 Mac in OC+ON niche"] = df[lig_col].sum()
        
        paired_df = pair_df[
            (pair_df["patient_id"] == patient_id) & (pair_df["core"] == core) &
            (pair_df["ligand"] == lig) & (pair_df["receptor"] == rec) & (pair_df["paired"] == True)
        ]
        row[f"{lig}+ TREM2 Mac in OC+ON niche associate with {rec}+ Epi"] = paired_df["cell"].nunique()
    
    summary_rows.append(row)

final_df = pd.DataFrame(summary_rows)

ordered_cols = ["patient_id", "TMA_id", "core", "Total TREM2 Mac in OC+ON niche", "Total TREM2 Mac across all niches"]
for lig, rec in LR_PAIRS:
    ordered_cols.append(f"{lig}+ TREM2 Mac in OC+ON niche")
    ordered_cols.append(f"{lig}+ TREM2 Mac in OC+ON niche associate with {rec}+ Epi")

final_df = final_df[ordered_cols]

print("\n===== FINAL QC =====")
print(final_df.head())
print("====================\n")

# STEP 8-I: Save outputs
OUTPUT_DIR = Path.home() / "Desktop" / f"{TMA_ID}_STEP8_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

out_path = OUTPUT_DIR / f"{TMA_ID}_Ligand+_TREM2_count.csv"
final_df.to_csv(out_path, index=False)

print("✅ STEP8 completed")
print("Saved file:", out_path)
display(final_df.head())

# ==============================================================================
# STEP 9 — Analyse Correlations between Macrophage Markers
# ==============================================================================
# Investigate correlations between % of Lig+ TREM2 Mac and % TREM2 Mac
# Both proportions were calculated among all macrophages across niches using data from STEP 3 and STEP 8
# Analysis was conducted at patient level per candidate ligand (e.g., GRN, SPP1, LGALS9)
# Before analysis, prepare TREM2 biomarker data file per ligand that include: 
# Sample_ID, Group, TREM2_Mac_percent, and Lpos_TREM2Mac_RposEpi_percent

# STEP 9-A: Check/install packages
%matplotlib inline
import os, sys, warnings, importlib.util, subprocess
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

packages = {"pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib", "openpyxl": "openpyxl", "scipy": "scipy"}
for module, package in packages.items():
    if importlib.util.find_spec(module) is None:
        print(f"Installing {package}..."); subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# STEP 9-B: Select input file
root = tk.Tk(); root.withdraw()
file_path = filedialog.askopenfilename(title="Select TREM2 biomarker data file", filetypes=[("Excel files", "*.xlsx"), ("CSV files", "*.csv"), ("All files", "*.*")])
root.destroy()
if not file_path: messagebox.showinfo("Cancelled", "No input file was selected."); sys.exit()
print(f"\nInput file: {file_path}")

# STEP 9-C: Create output folder
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
output_folder = os.path.join(desktop, "STEP9_Correlation_Plots_AllGroups_GRN")
os.makedirs(output_folder, exist_ok=True)
print(f"Output folder: {output_folder}")

# STEP 9-D: Load data
if file_path.lower().endswith(".xlsx"): df = pd.read_excel(file_path)
elif file_path.lower().endswith(".csv"): df = pd.read_csv(file_path)
else: raise ValueError("Please select an Excel (.xlsx) or CSV (.csv) file.")

# STEP 9-E: Required columns
required_columns = ["Sample_ID", "Group", "TREM2_Mac_percent", "Lpos_TREM2Mac_RposEpi_percent"]
missing = [col for col in required_columns if col not in df.columns]
if missing: messagebox.showerror("Missing columns", "These required columns are missing:\n\n" + "\n".join(missing)); raise ValueError(f"Missing columns: {missing}")

# STEP 9-F: Clean data
df["Group"] = df["Group"].astype(str).str.strip().str.upper()
df["TREM2_Mac_percent"] = pd.to_numeric(df["TREM2_Mac_percent"], errors="coerce")
df["Lpos_TREM2Mac_RposEpi_percent"] = pd.to_numeric(df["Lpos_TREM2Mac_RposEpi_percent"], errors="coerce")
df = df.dropna(subset=["Group", "TREM2_Mac_percent", "Lpos_TREM2Mac_RposEpi_percent"])

valid_groups = ["NR", "RD", "RI"]
df = df[df["Group"].isin(valid_groups)]
if df.empty: messagebox.showerror("Error", "No data found for groups NR, RD, or RI after cleaning."); sys.exit()

color_map = {"NR": "#5B8CFF", "RD": "#FFB300", "RI": "#E91E63"}
x_col, y_col = "TREM2_Mac_percent", "Lpos_TREM2Mac_RposEpi_percent"

rho_global, p_global = spearmanr(df[x_col], df[y_col])

group_stats = {}
for group in valid_groups:
    subset = df[df["Group"] == group]
    if len(subset) > 2:
        rho, p = spearmanr(subset[x_col], subset[y_col])
        group_stats[group] = {"rho": rho, "p": p, "n": len(subset)}
    else:
        group_stats[group] = {"rho": np.nan, "p": np.nan, "n": len(subset)}

# STEP 9-G: Generate dot plots
plt.figure(figsize=(10, 7))
plt.scatter(df[x_col], df[y_col], color="purple", s=60, edgecolors='black', linewidth=0.5, alpha=0.7, label="All Groups")
plt.xlabel("% TREM2+ Mac", fontsize=12)
plt.ylabel("% L+ TREM2+ Mac associated with R+ Epi", fontsize=12)
plt.title(f"Biomarker Correlation (All Groups)\nSpearman rho = {rho_global:.2f}, p = {p_global:.2g}", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=True, fancybox=True, shadow=True)
plt.subplots_adjust(right=0.75)
plot1_path = os.path.join(output_folder, "Correlation_all_no_distinction.png")
plt.savefig(plot1_path, dpi=300, bbox_inches='tight'); plt.close()
print(f"Saved: {plot1_path}")

plt.figure(figsize=(12, 7))
for group in valid_groups:
    subset = df[df["Group"] == group]
    stats = group_stats[group]
    label = f"{group} (n={stats['n']}, ρ={stats['rho']:.2f}, p={stats['p']:.2g})" if not np.isnan(stats['rho']) else f"{group} (n={stats['n']})"
    plt.scatter(subset[x_col], subset[y_col], color=color_map[group], s=60, edgecolors='black', linewidth=0.5, alpha=0.7, label=label)

plt.xlabel("% TREM2+ Mac", fontsize=12)
plt.ylabel("% L+ TREM2+ Mac associated with R+ Epi", fontsize=12)
plt.title(f"Biomarker Correlation by Group", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=True, fancybox=True, shadow=True)
plt.subplots_adjust(right=0.70)
plot2_path = os.path.join(output_folder, "Correlation_all_colored_by_group.png")
plt.savefig(plot2_path, dpi=300, bbox_inches='tight'); plt.close()
print(f"Saved: {plot2_path}")

print("\n==========================================")
print("PLOT GENERATION COMPLETE")
print("==========================================")
print(f"Results saved to: {output_folder}")

# ==============================================================================
# STEP 10 — Determine Distances of R+ myCAF from Closest L+ TREM2 Mac
# ==============================================================================
# This analysis requires the following files:
# 1. cell_meta_full.csv (from STEP3)
# 2. stromal_niche_full.csv (from STEP3)
# 3. core_to_nests.pkl (from STEP3)
# 4. cell_polygons_by_core.pkl (from STEP3) 
# 5. exprMat.csv (from AtoMx)

# Candidate ligand-receptor pairs were predefined by analysing snRNA-seq data (R code)
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)
# The output data are used to plot distribution densities via another script (R code) 

# STEP 10-A: Load required files 
import pandas as pd
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from tqdm import tqdm
from pathlib import Path
import pickle
from shapely.ops import nearest_points
import os

Tk().withdraw()

cell_meta_path = askopenfilename(title="Select cell_meta_full CSV file")
cell_meta_full = pd.read_csv(cell_meta_path)
print(f"✅ Loaded cell_meta_full: {cell_meta_full.shape[0]} cells")

stromal_niches_path = askopenfilename(title="Select stromal_niches_full CSV file")
stromal_niches_full = pd.read_csv(stromal_niches_path)
print(f"✅ Loaded stromal_niches_full: {stromal_niches_full.shape[0]} rows")

cell_polygons_path = askopenfilename(title="Select cell_polygons_by_core.pkl file")
with open(cell_polygons_path, "rb") as f: cell_polygons_by_core = pickle.load(f)
print("✅ Loaded cell_polygons_by_core.pkl")

core_to_nests_path = askopenfilename(title="Select core_to_nests.pkl file")
with open(core_to_nests_path, "rb") as f: core_to_nests = pickle.load(f)
print("✅ Loaded core_to_nests.pkl")

expr_path = askopenfilename(title="Select COUNT MATRIX (exprMat) CSV file", filetypes=[("CSV files", "*.csv")])
expr_mat = pd.read_csv(expr_path)

if "cell_ID" not in expr_mat.columns or "fov" not in expr_mat.columns:
    raise ValueError("exprMat must contain 'cell_ID' and 'fov'")

expr_mat["cell_ID"] = expr_mat["cell_ID"].astype(str)
expr_mat["fov"] = expr_mat["fov"].astype(str)
expr_mat["cell_key"] = expr_mat["fov"] + "_" + expr_mat["cell_ID"]
expr_mat = expr_mat.set_index("cell_key").drop(columns=["cell_ID", "fov"], errors="ignore")
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce").fillna(0).astype("int32")

print(f"✅ Expression matrix: {expr_mat.shape}")
print("✅ STEP10-A completed")

# STEP 10-B: Rebuild metadata
CORE_COL = "core_x" if "core_x" in cell_meta_full.columns else "core"
cell_meta_full["core"] = cell_meta_full[CORE_COL]

def parse_cell_id(x):
    parts = str(x).split("_")
    return (parts[1], parts[2], parts[3]) if len(parts) == 4 else (None, None, None)

parsed = cell_meta_full["cell_id"].astype(str).apply(parse_cell_id)
cell_meta_full["slide"] = parsed.str[0]
cell_meta_full["fov_parsed"] = parsed.str[1]
cell_meta_full["cell_ID_parsed"] = parsed.str[2]
cell_meta_full["cell_key"] = cell_meta_full["fov_parsed"].astype(str) + "_" + cell_meta_full["cell_ID_parsed"].astype(str)
cell_meta_full["cell"] = cell_meta_full["cell"].astype(str)

print("Using core column:", CORE_COL)
matched = cell_meta_full["cell_key"].isin(expr_mat.index).sum()
print(f"Match rate: {matched/len(cell_meta_full):.2%}")
print("✅ STEP10-B completed")

# STEP 10-C: Define parameters
TMA_ID = "TMA1"
CELL_DIST_PX = 44
DIST_THRESHOLD = 4 * CELL_DIST_PX
PX_TO_UM = 0.12  

LR_PAIRS = [("SPP1", "CD44"), ("SPP1", "ITGAV_ITGB1"), ("SPP1", "ITGAV_ITGB5"), ("SPP1", "ITGA5_ITGB1"), ("LGALS9", "CD44"), ("GRN", "SORT1"), ("CD99", "CD99")]
print("✅ STEP10-C completed")

# STEP 10-D: Define ligand positive TREM2 Mac
trem2_mac = cell_meta_full.loc[cell_meta_full["stromal_subset"].eq("TREM2.Mac")].copy()
trem2_mac["cell_key"] = trem2_mac["cell_key"].astype(str)
trem2_mac = trem2_mac[trem2_mac["cell_key"].isin(expr_mat.index)].copy().set_index("cell_key")
print(f"Total TREM2.Mac (aligned): {trem2_mac.shape[0]}")

expr_mat.columns = expr_mat.columns.str.strip()
expr_mat.index = expr_mat.index.astype(str)
expr_bin = expr_mat.astype(bool)

lig_cols = []
for lig, rec in LR_PAIRS:
    col_name = f"lig_pos_{lig}"
    lig_cols.append(col_name)
    if lig not in expr_bin.columns:
        trem2_mac[col_name] = False
    else:
        trem2_mac[col_name] = expr_bin.loc[trem2_mac.index, lig].to_numpy()

print("\n===== STEP10-D SANITY CHECK =====")
for lig, rec in LR_PAIRS:
    col = f"lig_pos_{lig}"
    if col in trem2_mac.columns:
        pos = trem2_mac[col].sum()
        pct = pos / len(trem2_mac) if len(trem2_mac) > 0 else 0
        print(f"{col}: {pos} ({pct:.2%})")
print("NaN check:", trem2_mac[lig_cols].isna().sum().sum())
print("===== END =====\n")
print("✅ STEP10-D completed")

# STEP 10-E: Define receptor positive myCAFs
mycaf_cells = cell_meta_full.loc[cell_meta_full["stromal_subset"].eq("myCAF")].copy()
mycaf_cells["cell_key"] = mycaf_cells["cell_key"].astype(str)
mycaf_cells = mycaf_cells[mycaf_cells["cell_key"].isin(expr_mat.index)].copy().set_index("cell_key")
print(f"myCAF cells (aligned): {len(mycaf_cells)}")

mycaf_cells["mycaf_state"] = "myCAF"
expr_bin = expr_mat.astype(bool)
rec_cols = []

for lig, rec in LR_PAIRS:
    col_name = f"rec_pos_{lig}"
    rec_cols.append(col_name)
    genes = rec.split("_") if isinstance(rec, str) else [rec]
    genes = [g for g in genes if g in expr_bin.columns]
    
    if len(genes) == 0:
        mycaf_cells[col_name] = False
        continue
    
    valid_cells = mycaf_cells.index.intersection(expr_bin.index)
    expr_sub = expr_bin.loc[valid_cells, genes]
    result = expr_sub.all(axis=1) if len(genes) > 1 else expr_sub.iloc[:, 0]
    
    mycaf_cells[col_name] = False
    mycaf_cells.loc[valid_cells, col_name] = result.to_numpy()

print("\n===== STEP10-E SANITY CHECK =====")
print("Total myCAF cells:", len(mycaf_cells))
for lig, rec in LR_PAIRS:
    col = f"rec_pos_{lig}"
    rec_label = rec.replace("_", "–") if isinstance(rec, str) else str(rec)
    if col in mycaf_cells.columns:
        print(f"myCAF {rec_label}+: {mycaf_cells[col].sum()}")
print("===== END =====\n")
print("✅ STEP10-E completed")

# STEP 10-F: Compute distance of Rec+ myCAF from closest Lig+ TREM2.Mac
results = []

for core, mycaf_df in mycaf_cells.groupby(CORE_COL):
    if core not in cell_polygons_by_core: continue
    poly_dict = cell_polygons_by_core[core]
    
    trem2_core = trem2_mac[trem2_mac[CORE_COL] == core].copy()
    trem2_core["poly"] = trem2_core["cell"].map(poly_dict)
    trem2_core = trem2_core.dropna(subset=["poly"])
    
    for lig, rec in LR_PAIRS:
        lig_cells = trem2_core[trem2_core[f"lig_pos_{lig}"]]
        if lig_cells.empty: continue
        
        lig_polys = [p for p in lig_cells["poly"] if p is not None and not p.is_empty]
        if len(lig_polys) == 0: continue
        
        mycaf_subset = mycaf_df[mycaf_df[f"rec_pos_{lig}"]]
        
        for _, row in mycaf_subset.iterrows():
            mycaf_poly = poly_dict.get(row["cell"])
            if mycaf_poly is None or mycaf_poly.is_empty: continue
            
            min_dist = min(mycaf_poly.distance(lp) for lp in lig_polys) * PX_TO_UM
            results.append({
                "TMA_id": TMA_ID, "core": core, "cell_id": row["cell"],
                "cell_state": row["mycaf_state"], "LR_pair": f"{lig}-{rec}",
                "min_distance_um": min_dist
            })

distance_df = pd.DataFrame(results)
print("Distance table:", distance_df.shape)
if not distance_df.empty: print(distance_df["min_distance_um"].describe())
else: print("WARNING: No distance results were generated.")
print("✅ STEP10-F completed")

# STEP 10-G: Export data file
outdir = Path.home() / "Desktop" / f"{TMA_ID}_STEP10_myCAF_results"
os.makedirs(outdir, exist_ok=True)

outpath = os.path.join(outdir, f"{TMA_ID}_LR_distance_myCAF.csv")
distance_df.to_csv(outpath, index=False)

print("✅ Saved:", outpath)
print("✅ STEP10 completed")

# ==============================================================================
# STEP 11 — Compute Percent of R+ myCAF within 40 µm of L+ TREM2 Mac
# ==============================================================================
# This analysis requires the following files:
# 1. LR_distance.csv (from STEP10 for all TMAs)
# 2. core_patient_group.csv (created by user)
# The core_patient_group includes TMA ID, Core ID, Patient ID, Clinical Group across all TMAs

# STEP 11-A: Load required data files
import pandas as pd
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from pathlib import Path
import os

Tk().withdraw()

NUM_FILES = 6
THRESHOLD_UM = 40

df_list = []
for i in range(NUM_FILES):
    file_path = askopenfilename(title=f"Select LR_distance CSV {i+1}/{NUM_FILES}", filetypes=[("CSV files", "*.csv")])
    if not file_path: raise ValueError("Missing input file")
    df = pd.read_csv(file_path)
    df["source_file"] = Path(file_path).name
    df_list.append(df)
    print(f"Loaded: {file_path}")

all_df = pd.concat(df_list, ignore_index=True)
print("\nCombined shape:", all_df.shape)

meta_path = askopenfilename(title="Select core_patient_group metadata CSV", filetypes=[("CSV files", "*.csv")])
meta_df = pd.read_csv(meta_path)

required = ["TMA_id", "core", "patient_id", "clinical_group"]
missing = [c for c in required if c not in meta_df.columns]
if missing: raise ValueError(f"Missing columns: {missing}")

all_df = all_df.merge(meta_df, on=["TMA_id", "core"], how="left")
print("After merge:", all_df.shape)

# STEP 11-B: Define myCAF and metric
all_df["cell_state"] = "myCAF"
all_df["within_40um"] = all_df["min_distance_um"] <= THRESHOLD_UM

# STEP 11-C: Per-core analysis
core_stats = all_df.groupby(["TMA_id", "core", "LR_pair"]).apply(lambda x: pd.Series({"N_cells": len(x), "Pct_within_40um": x["within_40um"].mean() * 100, "Mean_distance_um": x["min_distance_um"].mean(), "Median_distance_um": x["min_distance_um"].median()})).reset_index()
print("\nSTEP11 (per core):")
print(core_stats.head())

# STEP 11-D: Per-patient analysis
patient_stats = all_df.groupby(["patient_id", "clinical_group", "LR_pair"]).apply(lambda x: pd.Series({"N_cells": len(x), "Pct_within_40um": x["within_40um"].mean() * 100, "Mean_distance_um": x["min_distance_um"].mean(), "Median_distance_um": x["min_distance_um"].median()})).reset_index()
print("\nSTEP11 (per patient):")
print(patient_stats.head())

# STEP 11-E: Export data
outdir = Path.home() / "Desktop" / "STEP11_Pct40_results"
os.makedirs(outdir, exist_ok=True)

core_stats.to_csv(outdir / "STEP11_per_core_Pct40.csv", index=False)
patient_stats.to_csv(outdir / "STEP11_per_patient_Pct40.csv", index=False)

print("\nSaved to:")
print(outdir)
print("DONE.")

# ==============================================================================
# STEP 12 — Count Number of L+ TREM2 Mac in OC/ON paired with R+ myCAF
# ==============================================================================
# This analysis requires the following files:
# 1. cell_meta_full.csv (from STEP3)
# 2. stromal_niche_full.csv (from STEP3)
# 3. cell_polygons_by_core (from STEP3)
# 4. core_to_nests.pkl (from STEP3)
# 5. exprMat.csv (from AtoMx)
# Candidate ligand-receptor pairs were predefined by analysing snRNA-seq data (R code)
# Conduct the analysis per TMA (i.e., repeat same analysis for 6 TMAs)

# STEP 12-A: Load required data files
import pandas as pd
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from tqdm import tqdm
from pathlib import Path
import pickle
from shapely.ops import nearest_points
import os
from IPython.display import display

Tk().withdraw()

cell_meta_path = askopenfilename(title="Select cell_meta_full CSV file")
cell_meta_full = pd.read_csv(cell_meta_path)
print(f"✅ Loaded cell_meta_full: {cell_meta_full.shape[0]} cells")

stromal_niches_path = askopenfilename(title="Select stromal_niches_full CSV file")
stromal_niches_full = pd.read_csv(stromal_niches_path)
print(f"✅ Loaded stromal_niches_full: {stromal_niches_full.shape[0]} rows")

cell_polygons_path = askopenfilename(title="Select cell_polygons_by_core.pkl file")
with open(cell_polygons_path, "rb") as f: cell_polygons_by_core = pickle.load(f)
print("✅ Loaded cell_polygons_by_core.pkl")

core_to_nests_path = askopenfilename(title="Select core_to_nests.pkl file")
with open(core_to_nests_path, "rb") as f: core_to_nests = pickle.load(f)
print("✅ Loaded core_to_nests.pkl")

expr_path = askopenfilename(title="Select COUNT MATRIX (exprMat) CSV file", filetypes=[("CSV files", "*.csv")])
expr_mat = pd.read_csv(expr_path)

if "cell_ID" not in expr_mat.columns or "fov" not in expr_mat.columns:
    raise ValueError("exprMat must contain 'cell_ID' and 'fov'")

expr_mat["cell_ID"] = expr_mat["cell_ID"].astype(str)
expr_mat["fov"] = expr_mat["fov"].astype(str)
expr_mat["cell_key"] = expr_mat["fov"] + "_" + expr_mat["cell_ID"]
expr_mat = expr_mat.set_index("cell_key").drop(columns=["cell_ID", "fov"], errors="ignore")
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce").fillna(0).astype("int32")

print(f"✅ Expression matrix: {expr_mat.shape}")
print("✅ STEP12-A completed")

# STEP 12-B: Rebuild metadata
CORE_COL = "core_x" if "core_x" in cell_meta_full.columns else "core"
cell_meta_full["core"] = cell_meta_full[CORE_COL]

def parse_cell_id(x):
    parts = str(x).split("_")
    return (parts[1], parts[2], parts[3]) if len(parts) == 4 else (None, None, None)

parsed = cell_meta_full["cell_id"].astype(str).apply(parse_cell_id)
cell_meta_full["slide"] = parsed.str[0]
cell_meta_full["fov_parsed"] = parsed.str[1]
cell_meta_full["cell_ID_parsed"] = parsed.str[2]
cell_meta_full["cell_key"] = cell_meta_full["fov_parsed"].astype(str) + "_" + cell_meta_full["cell_ID_parsed"].astype(str)
cell_meta_full["cell"] = cell_meta_full["cell"].astype(str)

cell_meta_full = cell_meta_full.merge(stromal_niches_full[["cell", "final_niche"]], on="cell", how="left")

matched = cell_meta_full["cell_key"].isin(expr_mat.index).sum()
print(f"Match rate: {matched/len(cell_meta_full):.2%}")
print("✅ STEP12-B completed")

# STEP 12-C: Define parameters
TMA_ID = "TMA1"
CELL_COL = "cell"
REFERENCE_STROMAL_SUBSETS = ["TREM2.Mac"]
TARGET_NICHES = ["OC", "ON"]
CELL_DIST_PX = 44
DIST_THRESHOLD = 4 * CELL_DIST_PX
PX_TO_UM = 0.12  

LR_PAIRS = [("SPP1", "CD44"), ("LGALS9", "CD44"), ("GRN", "SORT1"), ("CD99", "CD99")]
print("✅ STEP12-C completed")

# STEP 12-D: Filter TREM2 Mac in OC/ON
stromal_df = stromal_niches_full[stromal_niches_full["final_niche"].isin(TARGET_NICHES)].copy()
required_cols = ["core", "patient_id", "clinical_group"]
missing = [c for c in required_cols if c not in stromal_df.columns]
if missing: raise ValueError(f"❌ Missing required columns: {missing}")

missing_cores = set(stromal_df["core"].unique()) - set(cell_polygons_by_core.keys())
print(f"⚠️ Missing cores in polygon dict: {len(missing_cores)}")

stromal_df["is_reference"] = stromal_df["stromal_subset"].isin(REFERENCE_STROMAL_SUBSETS)
stromal_df["inside_milieu"] = False

# STEP 12-E: Define ligand positive TREM2.Mac
trem2_mac = cell_meta_full.loc[cell_meta_full["stromal_subset"].eq("TREM2.Mac")].copy()
trem2_mac["cell_key"] = trem2_mac["cell_key"].astype(str)
trem2_mac = trem2_mac[trem2_mac["cell_key"].isin(expr_mat.index)].copy().set_index("cell_key")
print(f"Total TREM2.Mac (aligned): {trem2_mac.shape[0]}")

expr_mat.columns = expr_mat.columns.str.strip()
expr_mat.index = expr_mat.index.astype(str)
expr_bin = expr_mat.astype(bool)

lig_cols = []
for lig, rec in LR_PAIRS:
    col_name = f"lig_pos_{lig}"
    lig_cols.append(col_name)
    if lig not in expr_bin.columns:
        trem2_mac[col_name] = False
    else:
        trem2_mac[col_name] = expr_bin.loc[trem2_mac.index, lig].to_numpy()

print("\n===== STEP12-E SANITY CHECK =====")
for lig, rec in LR_PAIRS:
    col = f"lig_pos_{lig}"
    if col in trem2_mac.columns:
        pos = trem2_mac[col].sum()
        pct = pos / len(trem2_mac) if len(trem2_mac) > 0 else 0
        print(f"{col}: {pos} ({pct:.2%})")
print("NaN check:", trem2_mac[lig_cols].isna().sum().sum())
print("===== END =====\n")
print("✅ STEP12-E completed")

# STEP 12-F: Define receptor positive myCAFs
if "stromal_subset" not in cell_meta_full.columns: raise ValueError("Column 'stromal_subset' not found")

mycaf_cells = cell_meta_full.loc[cell_meta_full["stromal_subset"].eq("myCAF")].copy()
mycaf_cells["cell_key"] = mycaf_cells["cell_key"].astype(str)
mycaf_cells = mycaf_cells[mycaf_cells["cell_key"].isin(expr_mat.index)].copy()
if mycaf_cells.empty: raise ValueError("No myCAF cells found aligned with expression matrix.")

mycaf_cells = mycaf_cells.set_index("cell_key")
print(f"✅ myCAF cells (aligned): {len(mycaf_cells)}")

if "final_niche" not in mycaf_cells.columns:
    print("⚠️ Merging 'final_niche' into mycaf_cells...")
    mycaf_reset = mycaf_cells.reset_index().merge(stromal_niches_full[["cell", "final_niche"]].drop_duplicates(), on="cell", how="left")
    mycaf_cells = mycaf_reset.set_index("cell_key")
    if "final_niche" not in mycaf_cells.columns: raise ValueError("Failed to merge 'final_niche'.")
    print("✅ 'final_niche' successfully added.")

mycaf_cells["mycaf_state"] = "myCAF"
expr_bin = expr_mat.astype(bool)
rec_cols = []

for lig, rec in LR_PAIRS:
    col_name = f"rec_pos_{lig}"
    rec_cols.append(col_name)
    genes = rec.split("_") if isinstance(rec, str) else [rec]
    genes = [g for g in genes if g in expr_bin.columns]
    if len(genes) == 0: mycaf_cells[col_name] = False; continue
    valid_cells = mycaf_cells.index.intersection(expr_bin.index)
    if len(valid_cells) == 0: mycaf_cells[col_name] = False; continue
    expr_sub = expr_bin.loc[valid_cells, genes]
    result = expr_sub.all(axis=1) if len(genes) > 1 else expr_sub.iloc[:, 0]
    mycaf_cells[col_name] = False
    mycaf_cells.loc[valid_cells, col_name] = result.to_numpy()

print("\n===== STEP12-F SANITY CHECK =====")
print(f"Total myCAF cells: {len(mycaf_cells)}, 'final_niche' present: {'final_niche' in mycaf_cells.columns}")
for lig, rec in LR_PAIRS:
    col = f"rec_pos_{lig}"
    if col in mycaf_cells.columns:
        rec_label = rec.replace("_", "–") if isinstance(rec, str) else str(rec)
        print(f"myCAF {rec_label}+ ({lig}): {mycaf_cells[col].sum()}")
print("===== END =====\n")
print("✅ STEP12-F completed")

# STEP 12-G: Spatial Pairing (Lig+ TREM2.Mac <-> Rec+ myCAF)
PAIR_THRESHOLD_UM = 40
pair_results = []

if "final_niche" not in mycaf_cells.columns: raise ValueError("final_niche missing in mycaf_cells")

mycaf_ocon = mycaf_cells.reset_index()[lambda x: x["final_niche"].isin(TARGET_NICHES)].copy()
print(f"myCAFs in OC/ON niches: {len(mycaf_ocon)}")

if "final_niche" not in trem2_mac.columns:
    print("⚠️ Merging 'final_niche' into trem2_mac...")
    trem2_mac = trem2_mac.reset_index().merge(stromal_niches_full[["cell", "final_niche"]].drop_duplicates(), on="cell", how="left").set_index("cell_key")

print(f"🔄Processing {len(cell_polygons_by_core.keys())} cores...")

for core in tqdm(cell_polygons_by_core.keys(), desc="Calculating Distances"):
    poly_dict = cell_polygons_by_core[core]
    
    t2_temp = trem2_mac.reset_index() if trem2_mac.index.name == "cell_key" else trem2_mac.copy()
    trem2_core = t2_temp[(t2_temp[CORE_COL] == core) & (t2_temp["final_niche"].isin(TARGET_NICHES))].copy()
    if trem2_core.empty: continue
    if "cell" not in trem2_core.columns: trem2_core["cell"] = trem2_core.index
    trem2_core["poly"] = trem2_core["cell"].map(poly_dict)
    trem2_core = trem2_core.dropna(subset=["poly"])
    
    mycaf_core = mycaf_ocon[mycaf_ocon[CORE_COL] == core].copy()
    if mycaf_core.empty: continue
    if "cell" not in mycaf_core.columns: mycaf_core["cell"] = mycaf_core.index
    mycaf_core["poly"] = mycaf_core["cell"].map(poly_dict)
    mycaf_core = mycaf_core.dropna(subset=["poly"])
    
    mycaf_polys = {row["cell"]: row["poly"] for _, row in mycaf_core.iterrows()}
    
    for lig, rec in LR_PAIRS:
        lig_col, rec_col = f"lig_pos_{lig}", f"rec_pos_{lig}"
        if lig_col not in trem2_core.columns or rec_col not in mycaf_core.columns: continue
        
        lig_subset = trem2_core[trem2_core[lig_col]]
        rec_subset = mycaf_core[mycaf_core[rec_col]]
        if lig_subset.empty or rec_subset.empty: continue
        
        rec_polys = {row["cell"]: mycaf_polys.get(row["cell"]) for _, row in rec_subset.iterrows()}
        
        for _, mac in lig_subset.iterrows():
            mac_poly = mac["poly"]
            if mac_poly is None: continue
            min_dist = min((mac_poly.distance(m_poly) * PX_TO_UM for m_poly in rec_polys.values() if m_poly is not None), default=np.inf)
            
            p_id = mac.get("patient_id", np.nan)
            if pd.isna(p_id) and "patient_id" in mac.index: p_id = mac["patient_id"]
            
            pair_results.append({
                "patient_id": p_id, "core": core, "cell": mac["cell"],
                "ligand": lig, "receptor": rec, "paired": min_dist <= PAIR_THRESHOLD_UM,
                "min_distance_um": min_dist
            })

pair_df = pd.DataFrame(pair_results)
print("✅ STEP12-G completed")
print(pair_df.head() if not pair_df.empty else "⚠️ No pairs found within threshold.")

# STEP 12-H: Final counting table
if "pair_df" not in globals(): raise ValueError("Run STEP 12-G first")

summary_rows = []
group_cols = ["patient_id", CORE_COL]
trem2_ocon = trem2_mac[trem2_mac["final_niche"].isin(TARGET_NICHES)].copy()
trem2_all = trem2_mac.copy()

for (patient_id, core), df in trem2_ocon.groupby(group_cols):
    row = {
        "patient_id": patient_id, "TMA_id": TMA_ID, "core": core,
        "Total TREM2 Mac in OC+ON niche": len(df),
        "Total TREM2 Mac across all niches": len(trem2_all[(trem2_all["patient_id"] == patient_id) & (trem2_all[CORE_COL] == core)])
    }
    
    for lig, rec in LR_PAIRS:
        lig_col = f"lig_pos_{lig}"
        row[f"{lig}+ TREM2 Mac in OC+ON niche"] = df[lig_col].sum()
        
        paired_df = pair_df[
            (pair_df["patient_id"] == patient_id) & (pair_df["core"] == core) &
            (pair_df["ligand"] == lig) & (pair_df["receptor"] == rec) & (pair_df["paired"] == True)
        ]
        row[f"{lig}+ TREM2 Mac in OC+ON niche associate with {rec}+ myCAF"] = paired_df["cell"].nunique()
    
    summary_rows.append(row)

final_df = pd.DataFrame(summary_rows)

ordered_cols = ["patient_id", "TMA_id", "core", "Total TREM2 Mac in OC+ON niche", "Total TREM2 Mac across all niches"]
for lig, rec in LR_PAIRS:
    ordered_cols.append(f"{lig}+ TREM2 Mac in OC+ON niche")
    ordered_cols.append(f"{lig}+ TREM2 Mac in OC+ON niche associate with {rec}+ myCAF")

final_df = final_df[ordered_cols]

print("\n===== FINAL QC =====")
print(final_df.head())
print("====================\n")

# STEP 12-I: Save outputs
OUTPUT_DIR = Path.home() / "Desktop" / f"{TMA_ID}_STEP12_myCAF_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

out_path = OUTPUT_DIR / f"{TMA_ID}_Ligand+_TREM2_count_myCAF.csv"
final_df.to_csv(out_path, index=False)

print("✅ STEP12 completed")
print("Saved file:", out_path)
display(final_df.head())
