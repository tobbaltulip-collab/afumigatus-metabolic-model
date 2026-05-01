#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import glob
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import single_gene_deletion
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import cobra
print('numpy version:',np.__version__)
print('pandas version:',pd.__version__)
print('matplotlib version:',matplotlib.__version__)
print('cobra version:',cobra.__version__)


# In[4]:


from cobra.io import read_sbml_model

model = read_sbml_model("Pan_Aspergillus_fumigatus.xml")


# In[5]:


print('current model objective function:','\n',model.objective)


# In[6]:


model.reactions.get_by_id('R3625')


# In[7]:


solution = model.optimize()
print('growth of the cell when optimised for cell growth:',solution.objective_value,
      'mmol/gDW/hour, i.e.',1/(solution.objective_value),'cell divisions per hour')


# In[8]:


solution_df = solution.to_frame()
reactions = []
for n in range(len(solution_df)): 
    reactions.append(model.reactions.get_by_id((solution_df.index.to_list())[n]).reaction)
solution_df['reactions'] = reactions
original_solution_df = solution_df 
original_solution_df


# In[9]:


model.summary()


# In[10]:


summary = model.summary()
original_uptakes = summary.uptake_flux 
original_secretions = summary.secretion_flux 
print(original_secretions)


# In[11]:


from cobra.io import read_sbml_model

model = read_sbml_model("Pan_Aspergillus_fumigatus.xml")

print("Reactions:", len(model.reactions))
print("Metabolites:", len(model.metabolites))
print("Genes:", len(model.genes))
print("Objective:", model.objective)


# In[12]:


model.medium


# In[13]:


from cobra.medium import minimal_medium

baseline_growth = model.slim_optimize()

mm = minimal_medium(model, baseline_growth)

print(mm)


# In[14]:


print(model.objective.expression)


# In[15]:


for rxn in model.exchanges:
    print(rxn.id, "|", rxn.name)


# In[16]:


print(model.objective.direction)


# In[17]:


ex = [(r.id, r.name) for r in model.exchanges]
ex[:50]


# In[18]:


baseline_medium = model.medium.copy()
baseline_growth = model.optimize().objective_value
print("Baseline growth:", baseline_growth)


# In[19]:


candidates = [r for r in model.reactions if "maint" in (r.name or "").lower()]
if candidates:
    print (candidates[0])


# In[20]:


def find_exchange_ids(keywords):
    keywords = [k.lower() for k in keywords]
    hits = []
    for rxn in model.exchanges:
        hay = f"{rxn.id} {rxn.name}".lower()
        if any(k in hay for k in keywords):
            mets = list(rxn.metabolites.keys())
            met = mets[0] if mets else None
            hits.append((rxn.id, rxn.name, met.id if met else None, met.name if met else None))
    return hits

for label, keys in {
    "TRYPTOPHAN": ["tryptophan", "trp"],
    "PHENYLALANINE": ["phenylalanine", "phe"],
    "TYROSINE": ["tyrosine", "tyr"],
    "IRON": ["iron", "fe2", "fe3", "ferric", "ferrous"],
    "ZINC": ["zinc", "zn2"],
}.items():
    hits = find_exchange_ids(keys)
    print("\n", label)
    if not hits:
        print("  No matches. (Try KEGG IDs instead.)")
    else:
        for h in hits[:20]:
            print(" ", h)


# In[21]:


def apply_medium(model, updates: dict, base=None):
    m = (base if base is not None else model.medium).copy()
    m.update(updates)
    model.medium = m
    return model.medium


# In[22]:


baseline_growth = model.slim_optimize()
target_growth = 0.5 * baseline_growth
print("Baseline:", baseline_growth, "Target:", target_growth)


# In[23]:


baseline_growth = model.slim_optimize()
mm_dict = {k: float(v) for k, v in mm.items()}

with model:
    model.medium = mm_dict
    g = model.slim_optimize()
    print("Baseline:", baseline_growth)
    print("Growth on mm:", g)


# In[24]:


for ex_id in mm.index:
    r = model.reactions.get_by_id(ex_id)
    met = list(r.metabolites.keys())[0]
    print(ex_id, "|", r.name, "|", met.id, "|", met.name)


# In[25]:


import pandas as pd

def medium_table(model):
    rows = []
    for ex_id, ub in model.medium.items():
        r = model.reactions.get_by_id(ex_id)
        met = list(r.metabolites.keys())[0] if r.metabolites else None
        rows.append({
            "ex_id": ex_id,
            "ub": float(ub),
            "rxn_name": r.name,
            "met_id": getattr(met, "id", None),
            "met_name": getattr(met, "name", None),
            "rxn_bounds": r.bounds
        })
    return pd.DataFrame(rows).sort_values("ub", ascending=False)

medium_table(model)


# In[26]:


closed = []
for r in model.exchanges:
    lb, ub = r.bounds
    if lb == 0 and ub != 0:
        met = list(r.metabolites.keys())[0] if r.metabolites else None
        closed.append((r.id, r.name, getattr(met, "id", None), getattr(met, "name", None), r.bounds))

closed_df = pd.DataFrame(closed, columns=["ex_id","rxn_name","met_id","met_name","bounds"])
closed_df.head(50)


# In[27]:


model.optimize()


# In[29]:


[z.id for z in model.exchanges if "zn" in z.name.lower() or "zinc" in z.name.lower()]


# In[30]:


import cobra
from cobra.io import read_sbml_model

def build_closed_medium(model, default_ub=1000.0):
    for rxn in model.exchanges:
        rxn.lower_bound = 0.0
        rxn.upper_bound = default_ub


def apply_allowlist_bounds(model, allow_lb_by_rxn, default_ub=1000.0):
    for rxn_id, lb in allow_lb_by_rxn.items():
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = float(lb)
            rxn.upper_bound = float(default_ub)


def set_lung_scfm2_mmol(
    model_file,
    biomass_id,
    o2_lb=-1.0,         
    co2_lb=-1000.0,
    water_lb=-1000.0,
    default_ub=1000.0,
):
    model = read_sbml_model(model_file)

    build_closed_medium(model, default_ub=default_ub)

    allow = {}

    if "EX_C00001[e]" in model.reactions:   # H2O
        allow["EX_C00001[e]"] = water_lb

    if "EX_C00011[e]" in model.reactions:   # CO2
        allow["EX_C00011[e]"] = co2_lb

    if "EX_C00007[e]" in model.reactions:   # O2
        allow["EX_C00007[e]"] = float(o2_lb)

    scfm2_mmol = {
        "EX_C00031[e]": -3.00,      # D-glucose 
        "EX_C00186[e]": -9.30,      # (S)-lactate
        "EX_C00133[e]": -1.78,      # L-alanine
        "EX_C00062[e]": -0.31,      # L-arginine
        "EX_C00049[e]": -0.83,      # L-aspartate
        "EX_C00097[e]": -0.16,      # L-cysteine
        "EX_C00025[e]": -1.55,      # L-glutamate
        "EX_C00037[e]": -1.20,      # glycine
        "EX_C00135[e]": -0.52,      # L-histidine
        "EX_C00407[e]": -1.12,      # L-isoleucine
        "EX_C00123[e]": -1.61,      # L-leucine
        "EX_C00047[e]": -2.13,      # L-lysine
        "EX_C00073[e]": -0.63,      # L-methionine
        "EX_C00079[e]": -0.53,      # L-phenylalanine
        "EX_C00148[e]": -1.66,      # L-proline
        "EX_C00065[e]": -1.45,      # L-serine
        "EX_C00188[e]": -1.07,      # L-threonine
        "EX_C00078[e]": -0.01,      # L-tryptophan
        "EX_C00082[e]": -0.80,      # L-tyrosine
        "EX_C00183[e]": -1.12,      # L-valine
        "EX_C00077[e]": -0.68,      # L-ornithine
        "EX_C00009[e]": -2.55,      # phosphate total = 1.30 + 1.25 mM
        "EX_C00244[e]": -0.35,      # nitrate (KNO3)
        "EX_C00059[e]": -0.27,      # sulfate (K2SO4)
        "EX_C01342[e]": -2.28,      # NH3 
        "EX_CHEBI29033[e]": -0.0036,# Fe2+ = 3.60 uM = 0.0036 mM
        "EX_C00140[e]": -0.30,      # N-acetylglucosamine = 0.30 mM
    }

    for rxn_id, lb in scfm2_mmol.items():
        if rxn_id in model.reactions:
            allow[rxn_id] = lb

    apply_allowlist_bounds(model, allow, default_ub=default_ub)

    sol = model.optimize()
    growth = float(sol.fluxes.get(biomass_id, sol.objective_value))

    print("SCFM2-like growth:", growth)
    print("\nOpened reactions:")
    for rxn_id in sorted(allow.keys()):
        print(rxn_id, allow[rxn_id])

    print("\nActive uptake fluxes:")
    for rxn_id in sorted(allow.keys()):
        if rxn_id in sol.fluxes and sol.fluxes[rxn_id] < -1e-9:
            print(rxn_id, float(sol.fluxes[rxn_id]))

    return model, allow, growth, sol


# In[31]:


model, allow, growth, sol = set_lung_scfm2_mmol(
    model_file="Pan_Aspergillus_fumigatus.xml",
    biomass_id="R3625",
    o2_lb=-1.0
)


# In[32]:


model.summary()


# In[33]:


model.reactions.get_by_id("R3625").bounds


# In[34]:


print("=== ACTIVE UPTAKE FLUXES ===")
for rxn_id, flux in sol.fluxes.items():
    if flux < -1e-6 and rxn_id.startswith("EX_"):
        rxn = model.reactions.get_by_id(rxn_id)
        print(f"  {rxn_id:<25} {rxn.name:<30} {flux:.6f}")

print("\n=== UNEXPECTED UPTAKE FLUXES (not in your allow list) ===")
for rxn_id, flux in sol.fluxes.items():
    if flux < -1e-6 and rxn_id.startswith("EX_") and rxn_id not in allow:
        rxn = model.reactions.get_by_id(rxn_id)
        print(f"  !! {rxn_id:<25} {rxn.name:<30} {flux:.6f}")

print("\n=== OBJECTIVE FUNCTION ===")
print(model.objective)

critical = {
    "EX_CHEBI29033[e]": "Fe2+ (must be very low -0.0036)",
    "EX_C00007[e]":     "O2  (must be -2.0 for hypoxic)",
    "EX_C00078[e]":     "Tryptophan (must be very low -0.01)",
    "EX_C00079[e]":     "Phenylalanine",
    "EX_C00031[e]":     "Glucose",
}
print("\n=== CRITICAL CONSTRAINT CHECK ===")
for rxn_id, label in critical.items():
    if rxn_id in model.reactions:
        lb = model.reactions.get_by_id(rxn_id).lower_bound
        print(f"  {label:<45} lb = {lb}")
    else:
        print(f"  {rxn_id} NOT FOUND")


# In[35]:


print("Lactate flux:", sol.fluxes.get("EX_C00186[e]", "not in solution"))
print("Glutamine flux:", sol.fluxes.get("EX_C00064[e]", "not in solution"))


# In[36]:


import pandas as pd
from cobra.flux_analysis import single_gene_deletion

# ── Baseline growth ──────────────────────────────────────────────────────────
baseline_growth = sol.objective_value
threshold = 0.001  # from Mirhakkak et al. 2023
print(f"Baseline growth (SCFM2 lung): {baseline_growth:.6f}")

# ── Run single gene deletions ────────────────────────────────────────────────
print(f"\nRunning knockouts for {len(model.genes)} genes...")
print("This may take a few minutes...\n")

deletion_results = single_gene_deletion(model)

# ── Process results ──────────────────────────────────────────────────────────
deletion_results["growth"]              = deletion_results["growth"].fillna(0)
deletion_results["gene"]                = deletion_results["ids"].apply(lambda x: list(x)[0])
deletion_results["essential"]           = deletion_results["growth"] < threshold
deletion_results["percent_of_baseline"] = (deletion_results["growth"] / baseline_growth * 100).round(2)

# ── Split results ────────────────────────────────────────────────────────────
essential    = deletion_results[deletion_results["essential"]].sort_values("growth")
nonessential = deletion_results[~deletion_results["essential"]].sort_values("growth", ascending=False)

print(f"Total genes tested:   {len(deletion_results)}")
print(f"Essential genes:      {len(essential)}")
print(f"Non-essential genes:  {len(nonessential)}")

# ── Print top essential genes ─────────────────────────────────────────────────
print("\n=== Essential genes (antifungal target candidates) ===")
print(essential[["gene", "growth", "percent_of_baseline"]].to_string())

# ── Save results ──────────────────────────────────────────────────────────────
deletion_results.to_csv("knockouts_SCFM2.csv", index=False)
essential.to_csv("essential_genes_SCFM2.csv",  index=False)
print("\nSaved: knockouts_SCFM2.csv")
print("Saved: essential_genes_SCFM2.csv")


# In[37]:


print("=== Essential genes with annotations ===")
for _, row in essential.iterrows():
    gene_id = row["gene"]
    if gene_id in [g.id for g in model.genes]:
        gene = model.genes.get_by_id(gene_id)
        rxn_names = [r.name for r in gene.reactions][:3] 
        print(f"  {gene_id}  →  {rxn_names}")


# In[38]:


# ── Single gene knockout of AFUA_1G09050 ────────────────────────────────────

target_gene = "AFUA_1G09050"

print(f"Baseline growth (SCFM2 lung): {baseline_growth:.6f} mmol/gDW/h")

with model:
    gene = model.genes.get_by_id(target_gene)
    
    print(f"\nGene {target_gene} is involved in {len(gene.reactions)} reaction(s):")
    for rxn in gene.reactions:
        print(f"  {rxn.id:<15} {rxn.name}")
        print(f"    Current bounds: lb={rxn.lower_bound}, ub={rxn.upper_bound}")
    
    gene.knock_out()
    
    ko_solution = model.optimize()
    ko_growth   = ko_solution.objective_value if ko_solution.objective_value else 0.0
    
    growth_change = ko_growth - baseline_growth
    percent_change = (growth_change / baseline_growth) * 100

    print(f"\n=== Knockout Result: {target_gene} ===")
    print(f"  Gene name/function : PE methyltransferase (phospholipid biosynthesis)")
    print(f"  Baseline growth    : {baseline_growth:.6f} mmol/gDW/h")
    print(f"  Knockout growth    : {ko_growth:.6f} mmol/gDW/h")
    print(f"  Growth change      : {growth_change:.6f} mmol/gDW/h")
    print(f"  Percent change     : {percent_change:.2f}%")
    print(f"  Status             : {'ESSENTIAL - growth abolished' if ko_growth < 0.001 else 'Non-essential - growth maintained'}")
    print(f"  Solution status    : {ko_solution.status}")

restored_growth = model.optimize().objective_value
print(f"\nModel restored — baseline growth confirmed: {restored_growth:.6f}")


# In[39]:


import pandas as pd

# ── Define the three shikimate / AAA biosynthesis targets ───────────────────
aaa_targets = {
    "AFUA_1G13740": "Shikimate pathway (shikimate dehydrogenase + EPSP synthase)",
    "AFUA_6G04820": "Chorismate:L-glutamine aminotransferase (anthranilate synthase)",
    "AFUA_1G06940": "Chorismate synthase",
}

print(f"Baseline growth (SCFM2 lung): {baseline_growth:.6f} mmol/gDW/h")
print(f"Essentiality threshold:        0.001000 mmol/gDW/h")
print("=" * 65)

# ── KNOCKOUT 1: Each gene individually ──────────────────────────────────────
print("\n── Individual Knockouts ──────────────────────────────────────")

individual_results = {}

for gene_id, function in aaa_targets.items():
    with model:
        gene = model.genes.get_by_id(gene_id)
        gene.knock_out()
        sol = model.optimize()
        ko_growth = max(sol.objective_value, 0.0) if sol.objective_value else 0.0

        essential = ko_growth < 0.001
        pct_change = ((ko_growth - baseline_growth) / baseline_growth * 100)

        individual_results[gene_id] = {
            "growth":     ko_growth,
            "essential":  essential,
            "pct_change": pct_change
        }

        print(f"\n  Gene     : {gene_id}")
        print(f"  Function : {function}")
        print(f"  Growth   : {ko_growth:.6f} mmol/gDW/h  ({pct_change:.1f}%)")
        print(f"  Status   : {'ESSENTIAL ✓' if essential else 'NON-ESSENTIAL — alternative route exists!'}")

        # If non-essential, check what the model is using instead
        if not essential:
            print(f"  → Model maintained growth without this gene")
            print(f"  → Checking active uptake fluxes for alternative routes:")
            for rxn_id, flux in sol.fluxes.items():
                if flux < -1e-6 and rxn_id.startswith("EX_"):
                    rxn = model.reactions.get_by_id(rxn_id)
                    # Highlight aromatic amino acid uptake specifically
                    if any(aa in rxn.name.lower() for aa in 
                           ["phenyl", "tryptophan", "tyrosine", "chorismate",
                            "anthranilate", "shikimate", "aromatic"]):
                        print(f"    !! AAA UPTAKE DETECTED: {rxn_id} {rxn.name} flux={flux:.6f}")
                    else:
                        print(f"       {rxn_id:<25} {rxn.name:<35} flux={flux:.6f}")

# ── KNOCKOUT 2: All three genes simultaneously ───────────────────────────────
print("\n── Combined Knockout (all 3 shikimate genes simultaneously) ──")

with model:
    # Knock out all three at once
    for gene_id in aaa_targets.keys():
        model.genes.get_by_id(gene_id).knock_out()

    sol_combined = model.optimize()
    combined_growth = max(sol_combined.objective_value, 0.0) if sol_combined.objective_value else 0.0
    combined_essential = combined_growth < 0.001
    combined_pct = ((combined_growth - baseline_growth) / baseline_growth * 100)

    print(f"\n  Genes knocked out : {list(aaa_targets.keys())}")
    print(f"  Combined growth   : {combined_growth:.6f} mmol/gDW/h  ({combined_pct:.1f}%)")
    print(f"  Status            : {'ESSENTIAL ✓ — no alternative route' if combined_essential else 'NON-ESSENTIAL — alternative route found!'}")

    if not combined_essential:
        print(f"\n  !! IMPORTANT: Model found alternative route to obtain AAA")
        print(f"  Active uptake fluxes after combined knockout:")
        for rxn_id, flux in sol_combined.fluxes.items():
            if flux < -1e-6 and rxn_id.startswith("EX_"):
                rxn = model.reactions.get_by_id(rxn_id)
                print(f"    {rxn_id:<25} {rxn.name:<40} flux={flux:.6f}")
                if any(aa in rxn.name.lower() for aa in
                       ["phenyl", "tryptophan", "tyrosine", "chorismate",
                        "anthranilate", "shikimate"]):
                    print(f"    ^^ AROMATIC AMINO ACID UPTAKE — this is the alternative route!")
    else:
        print(f"\n  Combined knockout is lethal — no alternative route available")
        print(f"  These three genes together are essential for AAA supply")
        print(f"  This validates them as antifungal targets per Mirhakkak et al. 2023")

# ── KNOCKOUT 3: Block biosynthesis AND close AAA uptake ──────────────────────
print("\n── Combined Knockout + AAA uptake blocked ─────────────────────")
print("   (simulating drug that blocks both synthesis AND transport)")

with model:
    # Knock out all three biosynthesis genes
    for gene_id in aaa_targets.keys():
        model.genes.get_by_id(gene_id).knock_out()

    # Also block external aromatic amino acid uptake
    aaa_exchanges = {
        "EX_C00079[e]": "L-phenylalanine",
        "EX_C00078[e]": "L-tryptophan",
        "EX_C00082[e]": "L-tyrosine",
        "EX_C00493[e]": "shikimate",
        "EX_C00108[e]": "anthranilate",
    }

    for rxn_id, name in aaa_exchanges.items():
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = 0  # block uptake
            print(f"  Blocked uptake: {rxn_id} ({name})")

    sol_full = model.optimize()
    full_growth = max(sol_full.objective_value, 0.0) if sol_full.objective_value else 0.0
    full_essential = full_growth < 0.001

    print(f"\n  Growth after biosynthesis KO + uptake blocked: {full_growth:.6f} mmol/gDW/h")
    print(f"  Status: {'LETHAL ✓ — confirms AAA are essential and no other route exists' if full_essential else 'Still growing — another alternative route exists!'}")

# ── SUMMARY TABLE ────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)

summary_data = []
for gene_id, res in individual_results.items():
    summary_data.append({
        "condition":  gene_id,
        "function":   aaa_targets[gene_id],
        "growth":     round(res["growth"], 6),
        "essential":  res["essential"],
        "pct_change": round(res["pct_change"], 2)
    })

summary_data.append({
    "condition":  "All 3 combined KO",
    "function":   "Full shikimate pathway block",
    "growth":     round(combined_growth, 6),
    "essential":  combined_essential,
    "pct_change": round(combined_pct, 2)
})

summary_data.append({
    "condition":  "Combined KO + uptake blocked",
    "function":   "Biosynthesis + transport block",
    "growth":     round(full_growth, 6),
    "essential":  full_essential,
    "pct_change": round(((full_growth - baseline_growth) / baseline_growth * 100), 2)
})

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

# Save
summary_df.to_csv("AAA_biosynthesis_knockouts.csv", index=False)
print("\nSaved: AAA_biosynthesis_knockouts.csv")


# In[40]:


import pandas as pd
from collections import Counter

def categorize_reduction(reduction_pct):
    if reduction_pct <= 0:
        return "no reduction"
    elif reduction_pct < 10:
        return "0-10%"
    elif reduction_pct < 20:
        return "10-20%"
    elif reduction_pct < 50:
        return "20-50%"
    elif reduction_pct < 80:
        return "50-80%"
    elif reduction_pct < 99:
        return "80-99%"
    else:
        return "lethal"


def single_gene_knockout_screen(model, biomass_id="R3625", lethal_threshold=1e-6):
    """
    Run single-gene knockouts on the CURRENT model state.
    Assumes SCFM2 constraints are already applied to model.
    """

    # Wild-type growth
    wt_sol = model.optimize()
    wt_growth = float(wt_sol.fluxes.get(biomass_id, wt_sol.objective_value))

    print("=" * 50)
    print("WILD-TYPE GROWTH UNDER CURRENT MEDIUM")
    print("=" * 50)
    print("WT growth:", wt_growth)

    if wt_growth <= lethal_threshold:
        raise ValueError("Wild-type growth is zero or too close to zero. Check medium constraints first.")

    results = []

    for gene in model.genes:
        with model:
            gene.knock_out()
            sol = model.optimize()

            if sol.status != "optimal" or sol.objective_value is None:
                ko_growth = 0.0
            else:
                ko_growth = float(sol.fluxes.get(biomass_id, sol.objective_value))

            reduction_pct = max(0.0, (wt_growth - ko_growth) / wt_growth * 100)

            if ko_growth <= lethal_threshold:
                reduction_pct = 100.0

            results.append({
                "gene_id": gene.id,
                "gene_name": getattr(gene, "name", ""),
                "wt_growth": wt_growth,
                "ko_growth": ko_growth,
                "reduction_pct": reduction_pct,
                "category": categorize_reduction(reduction_pct)
            })

    df = pd.DataFrame(results)
    df = df.sort_values(by=["reduction_pct", "ko_growth"], ascending=[False, True]).reset_index(drop=True)

    summary = df["category"].value_counts().reindex(
        ["no reduction", "0-10%", "10-20%", "20-50%", "50-80%", "80-99%", "lethal"],
        fill_value=0
    )

    print("\n" + "=" * 50)
    print("KNOCKOUT SUMMARY")
    print("=" * 50)
    for cat, count in summary.items():
        print(f"{cat}: {count}")

    print("\nTop 20 most severe knockouts:")
    print(df[["gene_id", "gene_name", "ko_growth", "reduction_pct", "category"]].head(20).to_string(index=False))

    return df, summary


# In[41]:


model, allow, growth, sol = set_lung_scfm2_mmol(
    model_file="Pan_Aspergillus_fumigatus.xml",
    biomass_id="R3625",
    o2_lb=-1.0
)

df_ko, summary_ko = single_gene_knockout_screen(model, biomass_id="R3625")


# In[42]:


df_ko.to_csv("single_gene_knockout_results_scfm2.csv", index=False)
print("Saved knockout table to single_gene_knockout_results_scfm2.csv")


# In[43]:


lethal_genes = df_ko[df_ko["category"] == "lethal"]
print(lethal_genes[["gene_id", "gene_name", "ko_growth"]].to_string(index=False))


# In[44]:


important_genes = df_ko[df_ko["reduction_pct"] >= 50]
print(important_genes[["gene_id", "gene_name", "ko_growth", "reduction_pct"]].to_string(index=False))


# In[45]:


def knockout_gene_set(model, gene_ids, biomass_id="R3625"):
    with model:
        for gid in gene_ids:
            if gid in [g.id for g in model.genes]:
                model.genes.get_by_id(gid).knock_out()
            else:
                print(f"Gene not found: {gid}")

        sol = model.optimize()
        growth = float(sol.fluxes.get(biomass_id, sol.objective_value))

        print("\nAAA biosynthesis knockout result")
        print("Knocked out genes:", gene_ids)
        print("Growth:", growth)

        print("\nActive uptake fluxes after knockout:")
        for rxn in model.exchanges:
            flux = sol.fluxes[rxn.id]
            if flux < -1e-6:
                print(rxn.id, float(flux))

        return sol, growth


# In[46]:


def knockout_aaa_pathway(model, gene_list, biomass_id="R3625"):
    with model:
        for gene_id in gene_list:
            if gene_id in [g.id for g in model.genes]:
                model.genes.get_by_id(gene_id).knock_out()
            else:
                print(f"Missing gene: {gene_id}")

        sol = model.optimize()

        if sol.status != "optimal" or sol.objective_value is None:
            growth = 0.0
        else:
            growth = float(sol.fluxes.get(biomass_id, sol.objective_value))

        if abs(growth) < 1e-6:
            growth = 0.0

        print("\nAAA BIOSYNTHESIS KNOCKOUT")
        print("Growth:", growth)

        print("\nAAA uptake after knockout:")
        for rxn in ["EX_C00078[e]", "EX_C00079[e]", "EX_C00082[e]"]:
            if rxn in sol.fluxes.index:
                print(rxn, float(sol.fluxes[rxn]))
            else:
                print(rxn, "not found")

        return sol, growth


aaa_genes = [
    "AFUA_1G13740",
    "AFUA_3G14800",
    "AFUA_1G11590",
    "AFUA_3G14850",
    "AFUA_1G11610",
    "AFUA_4G06460",
    "AFUA_1G06940",
    "AFUA_2G09110",
    "AFUA_6G12110",
    "AFUA_5G05690",
    "AFUA_2G10450",
    "AFUA_2G13630",
    "AFUA_6G12580",
    "AFUA_1G13090",
    "AFUA_4G11980",
    "AFUA_6G13520",
    "AFUA_2G13250"
]

sol_aaa, growth_aaa = knockout_aaa_pathway(model, aaa_genes, biomass_id="R3625")


# In[47]:


def assign_process_label(reaction_names):
    text = reaction_names.lower()

    if any(x in text for x in ["shikimate", "chorismate", "tryptophan", "phenylalanine", "tyrosine", "arogenate", "prephenate", "anthranilate"]):
        return "aromatic amino acid metabolism"
    elif any(x in text for x in ["fatty acid", "acyl", "stearoyl", "phosphatidyl", "lipid"]):
        return "lipid metabolism"
    elif any(x in text for x in ["atpase", "nad", "nadp", "oxidoreductase", "dehydrogenase", "respiratory", "oxygen"]):
        return "energy/redox metabolism"
    elif any(x in text for x in ["ergosterol", "sterol", "lanosterol"]):
        return "sterol metabolism"
    elif any(x in text for x in ["ribose", "purine", "pyrimidine", "adenine", "guanine", "xanthine"]):
        return "nucleotide metabolism"
    elif any(x in text for x in ["phosphate", "sulfate", "ammonia", "glutamine", "glutamate", "nitrogen"]):
        return "nitrogen / inorganic metabolism"
    elif any(x in text for x in ["transport", "exchange"]):
        return "transport"
    else:
        return "other / unclear"


# In[48]:


df_lethal_process["process"] = df_lethal_process["reaction_names"].apply(assign_process_label)
print(df_lethal_process[["gene_id", "process", "reaction_names"]].to_string(index=False))


# In[ ]:


# Get genes with small (0–10%) reduction
mild_genes = df_ko[
    (df_ko["reduction_pct"] > 0) &
    (df_ko["reduction_pct"] <= 10)
]

# Print them clearly
print("\nGENES WITH 0–10% GROWTH REDUCTION:\n")
print(mild_genes[["gene_id", "gene_name", "ko_growth", "reduction_pct"]].to_string(index=False))

# Optional: check count (should be 18)
print("\nCount:", len(mild_genes))


# In[ ]:


print(df_ko[["gene_id", "ko_growth", "reduction_pct", "category"]]
      .sort_values("reduction_pct", ascending=False)
      .head(40)
      .to_string(index=False))


# In[ ]:


print(df_ko["category"].value_counts(dropna=False))


# In[ ]:


fresh_summary = df_ko["category"].value_counts().reindex(
    ["no reduction", "0-10%", "10-20%", "20-50%", "50-80%", "80-99%", "lethal"],
    fill_value=0
)

print(fresh_summary)


# In[ ]:


target_label = [x for x in df_ko["category"].unique() if "10%" in str(x)][0]

print("Target label is:", repr(target_label))

print(
    df_ko[df_ko["category"] == target_label]
    [["gene_id", "gene_name", "ko_growth", "reduction_pct", "category"]]
    .to_string(index=False)
)

print("\nCount:", len(df_ko[df_ko["category"] == target_label]))


# In[ ]:


def classify_reduction(x):
    if x >= 99.999:
        return "lethal"
    elif x == 0:
        return "no reduction"
    elif 0 < x <= 10:
        return "0-10%"
    elif 10 < x <= 20:
        return "10-20%"
    elif 20 < x <= 50:
        return "20-50%"
    elif 50 < x <= 80:
        return "50-80%"
    elif 80 < x < 99.999:
        return "80-99%"
    else:
        return "check"

df_ko["category"] = df_ko["reduction_pct"].apply(classify_reduction)

print(df_ko["category"].value_counts())


# In[ ]:


import pandas as pd
import re

# ---------- 1) get lethal genes ----------
lethal_genes = df_ko[df_ko["category"] == "lethal"]["gene_id"].tolist()

# ---------- 2) keyword-based broad pathway classifier ----------
PATHWAY_RULES = {
    "Shikimate / aromatic amino acid metabolism": [
        "shikimate", "chorismate", "anthranilate", "arogenate",
        "prephenate", "phenylalanine", "tyrosine", "tryptophan",
        "kynurenine", "indole"
    ],
    "Folate / PABA / one-carbon metabolism": [
        "aminodeoxychorismate", "4-amino-4-deoxychorismate",
        "p-aminobenzo", "folate", "tetrahydrofolate", "thymidylate"
    ],
    "Purine / pyrimidine metabolism": [
        "purine", "pyrimidine", "orotate", "omp", "saicar",
        "adenylosuccinate", "aicar", "guanylate", "ctp",
        "nucleoside", "nucleotide", "phosphoribosyl"
    ],
    "Fatty acid / lipid metabolism": [
        "fatty acid", "acyl", "acyl-carrier", "malonyl", "stearoyl",
        "beta-ketoacyl", "lipase", "phosphatid", "cardiolipin",
        "inositol", "triacylglycerol"
    ],
    "Sterol / ergosterol / isoprenoid metabolism": [
        "sterol", "ergosterol", "squalene", "mevalonate",
        "farnesyl", "isopentenyl", "diphosphate"
    ],
    "Heme / porphyrin metabolism": [
        "heme", "porphyrin", "uroporphyrin", "coproporphyrin",
        "protoporphyrin", "ferrochelatase", "aminolevulinate",
        "porphobilinogen"
    ],
    "Cofactor / vitamin metabolism": [
        "riboflavin", "lumazine", "pantothenate", "coa",
        "thiamine", "nicotinate", "nad", "nadp", "vitamin"
    ],
    "Cell wall / chitin metabolism": [
        "chitin", "glcnac", "udp-glcnac", "cell wall"
    ],
    "Central carbon / energy metabolism": [
        "glycolysis", "glyceraldehyde", "pyruvate", "tca",
        "citrate", "succinate", "malate", "atp", "oxidoreductase"
    ]
}

def classify_reaction_name(name: str):
    n = (name or "").lower()
    hits = []
    for pathway, words in PATHWAY_RULES.items():
        if any(w in n for w in words):
            hits.append(pathway)
    if not hits:
        return "Unclassified / needs manual check"
    return "; ".join(sorted(set(hits)))

# ---------- 3) build gene -> reactions table ----------
rows = []

for gene_id in lethal_genes:
    if gene_id not in [g.id for g in model.genes]:
        rows.append({
            "gene_id": gene_id,
            "reaction_id": None,
            "reaction_name": None,
            "reaction_equation": None,
            "broad_pathway": "Missing from model / check ID"
        })
        continue

    gene = model.genes.get_by_id(gene_id)
    rxns = list(gene.reactions)

    if not rxns:
        rows.append({
            "gene_id": gene_id,
            "reaction_id": None,
            "reaction_name": None,
            "reaction_equation": None,
            "broad_pathway": "No associated reactions"
        })
        continue

    for rxn in rxns:
        rows.append({
            "gene_id": gene_id,
            "reaction_id": rxn.id,
            "reaction_name": rxn.name,
            "reaction_equation": rxn.reaction,
            "broad_pathway": classify_reaction_name(rxn.name)
        })

gene_rxn_table = pd.DataFrame(rows)

# ---------- 4) make a cleaner one-row-per-gene summary ----------
summary_rows = []

for gene_id, sub in gene_rxn_table.groupby("gene_id"):
    pathways = sorted(set(sub["broad_pathway"].dropna()))
    rxn_names = sorted(set([x for x in sub["reaction_name"].dropna() if x]))
    summary_rows.append({
        "gene_id": gene_id,
        "n_reactions": len(sub),
        "broad_pathways": " | ".join(pathways),
        "reaction_names": " || ".join(rxn_names[:5])  # first 5 names to keep it readable
    })

gene_summary = pd.DataFrame(summary_rows).sort_values("gene_id")

# ---------- 5) pathway counts ----------
pathway_count_rows = []
for pathway, sub in gene_rxn_table.groupby("broad_pathway"):
    pathway_count_rows.append({
        "broad_pathway": pathway,
        "n_unique_genes": sub["gene_id"].nunique()
    })

pathway_counts = pd.DataFrame(pathway_count_rows).sort_values(
    "n_unique_genes", ascending=False
)

print("\n=== ONE ROW PER LETHAL GENE ===")
print(gene_summary.to_string(index=False))

print("\n=== PATHWAY COUNTS ===")
print(pathway_counts.to_string(index=False))

# optional save
gene_rxn_table.to_csv("lethal_gene_reaction_table.csv", index=False)
gene_summary.to_csv("lethal_gene_summary.csv", index=False)
pathway_counts.to_csv("lethal_gene_pathway_counts.csv", index=False)


# In[49]:


import matplotlib.pyplot as plt

pathways = [
    "Nucleotide\nbiosynthesis",
    "Cofactor\nbiosynthesis",
    "Ergosterol\nbiosynthesis",
    "Fatty acid\nbiosynthesis",
    "Heme/porphyrin\nbiosynthesis",
    "Phospholipid\nbiosynthesis",
    "Shikimate/AAA\nbiosynthesis",
    "Cell wall/chitin",
    "Other"
]
counts  = [21, 18, 9, 7, 7, 7, 3, 3, 5]
colours = ["#1565c0","#6a1b9a","#2e7d32","#e65100",
           "#b71c1c","#00838f","#f9a825","#37474f","#78909c"]

plt.figure(figsize=(9, 7))
wedges, texts, autotexts = plt.pie(
    counts,
    labels=pathways,
    colors=colours,
    autopct="%1.1f%%",
    startangle=140,
    pctdistance=0.8
)
plt.title("Distribution of 80 Essential Genes by Metabolic Pathway\nunder SCFM2 Lung Constraints",
          fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig("figure2_pathway_distribution.png", dpi=300)
plt.show()

