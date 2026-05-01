import os
import pandas as pd
from cobra.io import read_sbml_model

# -----------------------------
# SETTINGS
# -----------------------------
models_dir = "models"

GROWTH_THRESHOLD = 0.001
SECRETION_UB = 1000.0
FE_BOUND = -0.0066   # iron

# -----------------------------
# FIND ALL XML FILES
# -----------------------------
xml_files = []
for root, dirs, files in os.walk("models"):
    for f in files:
        if f.endswith(".xml"):
            xml_files.append(os.path.join(root, f))

print("Total strain models found:", len(xml_files))

# -----------------------------
# SCFM2 CONSTRAINTS
# -----------------------------
SCFM2_BOUNDS = {
    # environment
    "EX_C00001[e]": -1000.0,   # H2O
    "EX_C00011[e]": -1000.0,   # CO2
    "EX_C00007[e]": -1.0,      # O2

    # carbon
    "EX_C00031[e]": -3.0,      # glucose
    "EX_C00186[e]": -9.3,      # lactate

    # amino acids
    "EX_C00133[e]": -1.78,     # alanine
    "EX_C00062[e]": -0.31,     # arginine
    "EX_C00049[e]": -0.83,     # aspartate
    "EX_C00097[e]": -0.16,     # cysteine
    "EX_C00025[e]": -1.55,     # glutamate
    "EX_C00037[e]": -1.20,     # glycine
    "EX_C00135[e]": -0.52,     # histidine
    "EX_C00407[e]": -1.12,     # isoleucine
    "EX_C00123[e]": -1.61,     # leucine
    "EX_C00047[e]": -2.13,     # lysine
    "EX_C00073[e]": -0.63,     # methionine
    "EX_C00079[e]": -0.53,     # phenylalanine
    "EX_C00148[e]": -1.66,     # proline
    "EX_C00065[e]": -1.45,     # serine
    "EX_C00188[e]": -1.07,     # threonine
    "EX_C00078[e]": -0.01,     # tryptophan
    "EX_C00082[e]": -0.80,     # tyrosine
    "EX_C00183[e]": -1.12,     # valine
    "EX_C00077[e]": -0.68,     # ornithine

    # others
    "EX_C00009[e]": -2.55,         # phosphate
    "EX_C00244[e]": -0.35,         # nitrate
    "EX_C00059[e]": -0.27,         # sulfate
    "EX_C01342[e]": -2.28,         # NH3
    "EX_CHEBI29033[e]": FE_BOUND,  # Fe2+
    "EX_C00140[e]": -0.30          # N-acetylglucosamine
}


def apply_scfm2_constraints(model, bounds, secretion_ub=1000.0):
    for rxn in model.exchanges:
        rxn.lower_bound = 0.0
        rxn.upper_bound = secretion_ub

    applied = {}
    for rxn_id, lb in bounds.items():
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = float(lb)
            rxn.upper_bound = float(secretion_ub)
            applied[rxn_id] = float(lb)

    return model, applied

# -----------------------------
# APPLY CONSTRAINTS + CALCULATE GROWTH
# -----------------------------
wt_results = []

for i, model_file in enumerate(xml_files, start=1):
    try:
        model = read_sbml_model(model_file)
        model, applied = apply_scfm2_constraints(model, SCFM2_BOUNDS, SECRETION_UB)

        sol = model.optimize()

        growth = float(sol.objective_value) if (
            sol.status == "optimal" and sol.objective_value is not None
        ) else 0.0

        wt_results.append({
            "strain": model_file,
            "status": sol.status,
            "wt_growth": growth,
            "growing": growth > GROWTH_THRESHOLD,
            "n_constraints_applied": len(applied)
        })

        if i % 10 == 0:
            print(f"{i}/{len(xml_files)} strains processed")

    except Exception as e:
        wt_results.append({
            "strain": model_file,
            "status": "error",
            "wt_growth": None,
            "growing": False,
            "n_constraints_applied": None,
            "error": str(e)
        })

# -----------------------------
# SAVE RESULTS
# -----------------------------
wt_df = pd.DataFrame(wt_results)
wt_df.to_csv("all_strains_wt_growth_scfm2.csv", index=False)

print("Saved: all_strains_wt_growth_scfm2.csv")
print("Growing strains:", wt_df["growing"].sum())
print("Non-growing strains:", (~wt_df["growing"]).sum())

import pandas as pd
import numpy as np
from cobra.io import read_sbml_model

GROWTH_THRESHOLD = 0.001

# -----------------------------
# Target single genes
# -----------------------------
single_targets = {
    "AFUA_1G13740": "Shikimate pathway",
    "AFUA_6G04820": "Anthranilate / tryptophan biosynthesis",
    "AFUA_1G06940": "Chorismate synthase",
    "AFUA_1G09050": "Phospholipid metabolism",
    "AFUA_7G01220": "Ergosterol biosynthesis",
    "AFUA_3G04220": "Fatty acid synthesis",
    "AFUA_5G07750": "Heme biosynthesis",
    "AFUA_6G08440": "Heme biosynthesis",
    "AFUA_4G08340": "Sterol/isoprenoid/respiratory metabolism",
    "AFUA_6G03750": "Purine biosynthesis"
}

# -----------------------------
# Double + triple combinations
# -----------------------------
combined_targets = {
    "AFUA_1G13740 + AFUA_7G01220": {
        "genes": ["AFUA_1G13740", "AFUA_7G01220"],
        "pathway": "Shikimate pathway + ergosterol biosynthesis"
    },
    "AFUA_5G07750 + AFUA_3G04220": {
        "genes": ["AFUA_5G07750", "AFUA_3G04220"],
        "pathway": "Heme biosynthesis + fatty acid synthesis"
    },
    "AFUA_6G03750 + AFUA_1G09050": {
        "genes": ["AFUA_6G03750", "AFUA_1G09050"],
        "pathway": "Purine biosynthesis + phospholipid metabolism"
    },
    "AFUA_1G13740 + AFUA_6G04820 + AFUA_1G06940": {
        "genes": ["AFUA_1G13740", "AFUA_6G04820", "AFUA_1G06940"],
        "pathway": "Combined shikimate/aromatic amino acid biosynthesis"
    }
}

# -----------------------------
# Helper function
# -----------------------------
def classify_ko(model, genes, wt_growth, threshold=GROWTH_THRESHOLD):
    gene_ids = {g.id for g in model.genes}

    missing = [g for g in genes if g not in gene_ids]
    if missing:
        return "gene_missing"

    with model:
        for gene_id in genes:
            model.genes.get_by_id(gene_id).knock_out()

        sol = model.optimize()

        if sol.status != "optimal" or sol.objective_value is None:
            return "lethal"

        ko_growth = float(sol.objective_value)

        if ko_growth < threshold:
            return "lethal"
        else:
            return "grew"

# -----------------------------
# Run across all strain models
# -----------------------------
summary_results = []

for model_file in xml_files:
    try:
        model = read_sbml_model(model_file)
        model, applied = apply_scfm2_constraints(model, SCFM2_BOUNDS, SECRETION_UB)

        wt_sol = model.optimize()

        if wt_sol.status != "optimal" or wt_sol.objective_value is None:
            continue

        wt_growth = float(wt_sol.objective_value)

        if wt_growth < GROWTH_THRESHOLD:
            continue

        # single knockouts
        for gene_id, pathway in single_targets.items():
            result = classify_ko(model, [gene_id], wt_growth)

            summary_results.append({
                "test": gene_id,
                "type": "single",
                "pathway": pathway,
                "result": result
            })

        # double/triple knockouts
        for combo_name, info in combined_targets.items():
            ko_type = "triple" if len(info["genes"]) == 3 else "double"
            result = classify_ko(model, info["genes"], wt_growth)

            summary_results.append({
                "test": combo_name,
                "type": ko_type,
                "pathway": info["pathway"],
                "result": result
            })

    except Exception as e:
        print(f"Error with {model_file}: {e}")

# -----------------------------
# Summarise counts only
# -----------------------------
results_df = pd.DataFrame(summary_results)

summary_table = (
    results_df
    .groupby(["test", "type", "pathway", "result"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

for col in ["lethal", "grew", "gene_missing"]:
    if col not in summary_table.columns:
        summary_table[col] = 0

summary_table = summary_table[
    ["test", "type", "pathway", "lethal", "grew", "gene_missing"]
]

summary_table.to_csv("strain_specific_ko_summary_counts.csv", index=False)

summary_table

control_genes = [
    "AFUA_8G05440",
    "AFUA_2G01950"
]

for gene_id in control_genes:
    print(f"\nTesting control {gene_id}...")

    if gene_id not in [g.id for g in model.genes]:
        print("Gene not present")
        continue

    with model:
        model.genes.get_by_id(gene_id).knock_out()
        ko_sol = model.optimize()
        ko_growth = float(ko_sol.objective_value) if (
            ko_sol.status == "optimal" and ko_sol.objective_value is not None
        ) else 0.0

        print("KO status:", ko_sol.status)
        print("KO growth:", ko_growth)
        print("Relative growth:", ko_growth / wt_growth if wt_growth > 0 else None)

def classify(reaction):
    r = reaction.lower()

    if any(k in r for k in ["imp", "amp", "gmp", "ump", "phosphoribosyl"]):
        return "Purine / pyrimidine metabolism"

    if any(k in r for k in ["acyl", "lipid", "phosphatidyl", "fatty"]):
        return "Fatty acid / lipid metabolism"

    if any(k in r for k in ["thiamine", "nad", "biotin", "pantothenate"]):
        return "Cofactor / vitamin metabolism"

    if any(k in r for k in ["sterol", "mevalonate", "squalene"]):
        return "Sterol / isoprenoid metabolism"

    if any(k in r for k in ["heme", "porphyrin"]):
        return "Heme metabolism"

    if any(k in r for k in ["folate", "tetrahydro"]):
        return "Folate metabolism"

    return "Central carbon / energy metabolism"

import pandas as pd
import matplotlib.pyplot as plt

# load adjusted results
df = pd.read_csv("all_strains_wt_growth_scfm2.csv")

# replace these with your original counts if known
orig_growing = 20
orig_nongrowing = 232

adj_growing = int(df["growing"].sum())
adj_nongrowing = len(df) - adj_growing

fig, axes = plt.subplots(1,2, figsize=(12,5))

# Panel A
labels = ["Original\nSCFM2", "Adjusted\nFe²⁺"]
grow = [orig_growing, adj_growing]
nogrow = [orig_nongrowing, adj_nongrowing]

axes[0].bar(labels, grow, label="Growing")
axes[0].bar(labels, nogrow, bottom=grow, label="Non-growing")
axes[0].set_ylabel("Number of strains")
axes[0].set_title("Growth feasibility across strain models")
axes[0].legend()

# Panel B
axes[1].hist(df["wt_growth"].dropna(), bins=20)
axes[1].axvline(0.001, linestyle="--")
axes[1].set_xlabel("Growth rate (mmol/gDW/h)")
axes[1].set_ylabel("Number of strains")
axes[1].set_title("Growth rate distribution (adjusted Fe²⁺)")

plt.tight_layout()
plt.savefig("figure4_strain_growth.png", dpi=300)
plt.show()

import os
print(os.listdir())

import cobra
print(cobra.__version__)

import pandas as pd
import numpy as np

# Baseline / wild-type growth under current SCFM2 constraints
baseline_solution = model.optimize()

if baseline_solution.status != "optimal":
    raise ValueError(f"Baseline model is not optimal. Status: {baseline_solution.status}")

baseline_growth = baseline_solution.objective_value
print("Baseline growth:", baseline_growth)

target_genes = {
    "AFUA_1G13740": "Shikimate pathway",
    "AFUA_6G04820": "Anthranilate / tryptophan biosynthesis",
    "AFUA_1G06940": "Chorismate synthase",
    "AFUA_1G09050": "Phospholipid metabolism",
    "AFUA_7G01220": "Ergosterol biosynthesis",
    "AFUA_3G04220": "Fatty acid synthesis",
    "AFUA_5G07750": "Heme biosynthesis",
    "AFUA_6G08440": "Heme biosynthesis",
    "AFUA_4G08340": "Sterol/isoprenoid/respiratory metabolism",
    "AFUA_6G03750": "Purine biosynthesis"
}

single_ko_results = []

gene_ids_in_model = [g.id for g in model.genes]

for gene_id, pathway in target_genes.items():
    with model:
        if gene_id not in gene_ids_in_model:
            single_ko_results.append({
                "Gene": gene_id,
                "Pathway": pathway,
                "KO growth": np.nan,
                "Relative growth": np.nan,
                "Growth reduction (%)": np.nan,
                "Phenotype": "Gene not found",
                "Solver status": "Gene not found"
            })
            continue

        model.genes.get_by_id(gene_id).knock_out()
        solution = model.optimize()

        if solution.status == "optimal":
            ko_growth = solution.objective_value
            relative_growth = ko_growth / baseline_growth if baseline_growth > 0 else np.nan
            growth_reduction = (1 - relative_growth) * 100

            if ko_growth < 0.001:
                phenotype = "Lethal"
            elif relative_growth < 0.5:
                phenotype = "Strong growth defect"
            elif relative_growth < 0.9:
                phenotype = "Moderate growth defect"
            else:
                phenotype = "No major growth defect"

        else:
            ko_growth = 0
            relative_growth = 0
            growth_reduction = 100
            phenotype = "Lethal / infeasible"

        single_ko_results.append({
            "Gene": gene_id,
            "Pathway": pathway,
            "Growth reduction (%)": growth_reduction,
            "Phenotype": phenotype,
        })

single_ko_df = pd.DataFrame(single_ko_results)

# Round for results table
single_ko_df["Growth reduction (%)"] = single_ko_df["Growth reduction (%)"].round(1)

single_ko_df

double_combinations = [
    {
        "Combination": "AFUA_1G13740 + AFUA_7G01220",
        "Genes": ["AFUA_1G13740", "AFUA_7G01220"],
        "Pathways": "Shikimate pathway + ergosterol biosynthesis"
    },
    {
        "Combination": "AFUA_5G07750 + AFUA_3G04220",
        "Genes": ["AFUA_5G07750", "AFUA_3G04220"],
        "Pathways": "Heme biosynthesis + fatty acid synthesis"
    },
    {
        "Combination": "AFUA_6G03750 + AFUA_1G09050",
        "Genes": ["AFUA_6G03750", "AFUA_1G09050"],
        "Pathways": "Purine biosynthesis + phospholipid metabolism"
    }
]

double_ko_results = []

for combo in double_combinations:
    with model:
        missing_genes = []

        for gene_id in combo["Genes"]:
            if gene_id in gene_ids_in_model:
                model.genes.get_by_id(gene_id).knock_out()
            else:
                missing_genes.append(gene_id)

        if missing_genes:
            double_ko_results.append({
                "Combination": combo["Combination"],
                "Pathways": combo["Pathways"],
                "KO growth": np.nan,
                "Relative growth": np.nan,
                "Growth reduction (%)": np.nan,
                "Phenotype": "Gene not found",
                "Solver status": "Missing: " + ", ".join(missing_genes)
            })
            continue

        solution = model.optimize()
        growth, relative_growth, reduction, phenotype = classify_growth(solution, baseline_growth)

        double_ko_results.append({
            "Combination": combo["Combination"],
            "Pathways": combo["Pathways"],
            "Growth reduction (%)": reduction,
            "Phenotype": phenotype,
        })

double_ko_df = pd.DataFrame(double_ko_results)

double_ko_df["Growth reduction (%)"] = double_ko_df["Growth reduction (%)"].round(1)

double_ko_df

triple_shikimate_genes = ["AFUA_1G13740", "AFUA_6G04820", "AFUA_1G06940"]

with model:
    missing_genes = []

    for gene_id in triple_shikimate_genes:
        if gene_id in gene_ids_in_model:
            model.genes.get_by_id(gene_id).knock_out()
        else:
            missing_genes.append(gene_id)

    if missing_genes:
        triple_ko_df = pd.DataFrame([{
            "Combination": " + ".join(triple_shikimate_genes),
            "Pathways": "Combined shikimate / aromatic amino acid biosynthesis",
            "KO growth": np.nan,
            "Relative growth": np.nan,
            "Growth reduction (%)": np.nan,
            "Phenotype": "Gene not found",
            "Solver status": "Missing: " + ", ".join(missing_genes)
        }])
    else:
        solution = model.optimize()
        growth, relative_growth, reduction, phenotype = classify_growth(solution, baseline_growth)

        triple_ko_df = pd.DataFrame([{
            "Combination": " + ".join(triple_shikimate_genes),
            "Pathways": "Combined shikimate/aromatic amino acid biosyn.",
            "Growth reduction (%)": reduction,
            "Phenotype": phenotype,
        }])

triple_ko_df["Growth reduction (%)"] = triple_ko_df["Growth reduction (%)"].round(1)

triple_ko_df

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------
def add_box(ax, x, y, text, width=2.7, height=0.7):
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.08",
        linewidth=1.2,
        facecolor="white",
        edgecolor="black"
    )
    ax.add_patch(box)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=8
    )

def add_arrow(ax, x1, y1, x2, y2):
    arrow = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="->",
        mutation_scale=12,
        linewidth=1.1,
        color="black"
    )
    ax.add_patch(arrow)

# -----------------------------
# FIGURE SETUP
# -----------------------------
fig, ax = plt.subplots(figsize=(13, 8))
ax.set_xlim(0, 13)
ax.set_ylim(0, 10)
ax.axis("off")

# -----------------------------
# SHIKIMATE PATHWAY
# -----------------------------
ax.text(0.5, 9.2, "Aromatic amino acid biosynthesis", fontsize=10, weight="bold")

add_box(ax, 0.5, 8.2, "Shikimate pathway\nAFUA_1G13740")
add_arrow(ax, 3.2, 8.55, 4.0, 8.55)

add_box(ax, 4.0, 8.2, "Chorismate synthase\nAFUA_1G06940")
add_arrow(ax, 6.7, 8.55, 7.5, 8.55)

add_box(ax, 7.5, 8.2, "Tryptophan branch\nAFUA_6G04820")

# -----------------------------
# STEROL / ISOPRENOID
# -----------------------------
ax.text(0.5, 7.2, "Sterol / isoprenoid metabolism", fontsize=10, weight="bold")

add_box(ax, 0.5, 6.2, "Mevalonate / isoprenoid\nprecursors")
add_arrow(ax, 3.2, 6.55, 4.0, 6.55)

add_box(ax, 4.0, 6.2, "Squalene synthase\nAFUA_7G01220")
add_arrow(ax, 6.7, 6.55, 7.5, 6.55)

add_box(ax, 7.5, 6.2, "Ergosterol\nbiosynthesis")
add_arrow(ax, 10.2, 6.55, 10.6, 6.55)

add_box(
    ax,
    10.6,
    6.2,
    "Respiratory / sterol-linked\nmetabolism\nAFUA_4G08340",
    width=2.2
)

# -----------------------------
# FATTY ACID / PHOSPHOLIPID
# -----------------------------
ax.text(0.5, 5.2, "Fatty acid and phospholipid metabolism", fontsize=10, weight="bold")

add_box(ax, 0.5, 4.2, "Fatty acid synthesis\nAFUA_3G04220")
add_arrow(ax, 3.2, 4.55, 4.0, 4.55)

add_box(ax, 4.0, 4.2, "Membrane lipid\nprecursors")
add_arrow(ax, 6.7, 4.55, 7.5, 4.55)

add_box(ax, 7.5, 4.2, "PE → PC methylation\nAFUA_1G09050")

# -----------------------------
# HEME
# -----------------------------
ax.text(0.5, 3.2, "Heme biosynthesis", fontsize=10, weight="bold")

add_box(ax, 0.5, 2.2, "Porphyrin pathway\nAFUA_6G08440")
add_arrow(ax, 3.2, 2.55, 4.0, 2.55)

add_box(ax, 4.0, 2.2, "Ferrochelatase\nAFUA_5G07750")
add_arrow(ax, 6.7, 2.55, 7.5, 2.55)

add_box(ax, 7.5, 2.2, "Heme-dependent\nmetabolism")

# -----------------------------
# PURINE
# -----------------------------
ax.text(0.5, 1.5, "Purine biosynthesis", fontsize=10, weight="bold")

add_box(ax, 0.5, 0.5, "PRPP / purine\nbiosynthesis\nAFUA_6G03750")
add_arrow(ax, 3.2, 0.85, 4.0, 0.85)

add_box(ax, 4.0, 0.5, "Nucleotide\nproduction")

# -----------------------------
# TITLE + SAVE
# -----------------------------
plt.title(
    "Metabolic Distribution of Prioritised Essential Targets",
    fontsize=13,
    weight="bold"
)

plt.tight_layout()
plt.savefig("prioritised_gene_pathway_map.png", dpi=300, bbox_inches="tight")
plt.show()
