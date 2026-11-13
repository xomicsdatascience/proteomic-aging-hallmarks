import pandas as pd

nes = pd.read_csv("results/nes_matrix_all_drugs.csv", index_col=0)
fdr = pd.read_csv("results/fdr_q_matrix_all_drugs.csv", index_col=0)
meta = pd.read_csv("results/de_summary_with_metadata.csv")

known_aging_kws = [
    "lonafarnib", "nordihydroguaiaretic", "nicotinamide", "lipoic", "methylene blue",
    "baicalein", "dimethyl fumarate", "monomethyl fumarate", "simvastatin", "atorvastatin",
    "eicosapentaenoic", "dfo", "mercaptopurine", "cabazitaxel", "mifepristone",
    "tucidinostat", "valemetostat", "bufexamac", "amodiaquine", "amiodarone",
    "binimetinib", "erlotinib", "lapatinib", "abemaciclib", "vemurafenib",
    "ferulic", "oryzanol", "danthron", "diacerein", "rilpivirine", "efavirenz"
]

print("=== Known / Literature Aging-Related Drugs in Dataset ===")
matched = [d for d in nes.index if any(k in d.lower() for k in known_aging_kws)]
for d in matched:
    cat = meta.loc[meta["drug"] == d, "Category"].values[0]
    sig_cols = fdr.columns[fdr.loc[d] < 0.05].tolist()
    sig_str = ", ".join([f"{c} ({nes.loc[d, c]:+.2f}, q={fdr.loc[d, c]:.3f})" for c in sig_cols])
    print(f"- {d} [{cat}]: {sig_str}")

print("\n=== DEMAGALHAES_AGING_UP / DN strongest hits (q < 0.05) ===")
for col in ["DEMAGALHAES_AGING_UP", "DEMAGALHAES_AGING_DN"]:
    sub = pd.DataFrame({"NES": nes[col], "q": fdr[col]}).query("q < 0.05").sort_values("NES")
    print(f"\n{col} (n_sig={len(sub)}):")
    print("  Most negative:", sub.head(5).to_dict("index"))
    print("  Most positive:", sub.tail(5).to_dict("index"))
