import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# -----------------------------
# 1. CONFIGURATION
# -----------------------------
CONFIG = dict(
    reconstructed_tsv="output_reconstruction/ancestor_gene_sets_by_node.tsv",
    true_root_karyotype="output_simulator/karyotypes_true_root.txt",
    output_png="output_reconstruction/DotPlot_ReconstructedRoot_vs_TrueRoot.png",
    root_label="Root",
)


# -----------------------------
# 2. DATA LOADING
# -----------------------------
def _load_true_root_karyotype(path):
    if not path or not os.path.exists(path):
        return None

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        engine="python",
        comment="-",
        dtype=str,
        names=["raw"],
    ).fillna("")

    karyo = {}
    for raw in df["raw"].tolist():
        if not raw:
            continue
        if raw.startswith("Species"):
            continue
        parts = str(raw).split()
        if len(parts) < 3:
            continue
        sp_name = parts[0]
        if sp_name != "TrueRoot":
            continue
        cid = parts[1]
        start_idx = 3 if (len(parts) >= 4 and parts[2].isdigit()) else 2
        genes = [g.replace("-", "") for g in parts[start_idx:]]
        karyo[cid] = genes

    return karyo


# -----------------------------
# 3. DOT PLOT GENERATION
# -----------------------------
def generate_root_vs_trueroot_dotplot(reconstructed_tsv_path, true_root_karyo_path, out_path, root_label="Root"):
    if not os.path.exists(reconstructed_tsv_path):
        print(f"Reconstructed TSV not found: {reconstructed_tsv_path}")
        return False
    true_root = _load_true_root_karyotype(true_root_karyo_path)
    if not true_root:
        print(f"True root karyotype not found or empty: {true_root_karyo_path}")
        return False

    df = pd.read_csv(reconstructed_tsv_path, sep="\t", dtype=str).fillna("")
    if "Node" not in df.columns or "ChrID" not in df.columns or "Genes" not in df.columns:
        print("Reconstructed TSV missing required columns (Node, ChrID, Genes)")
        return False
    df = df[df["Node"] == root_label]

    rec_karyo = {}
    rec_order = []
    for row in df.itertuples(index=False):
        cid = getattr(row, "ChrID", "")
        genes_str = getattr(row, "Genes", "")
        if not cid:
            continue
        genes = [g for g in str(genes_str).split(" ") if g]
        rec_karyo[cid] = genes
        rec_order.append(cid)

    if not rec_karyo:
        print("No reconstructed karyotype data found for root.")
        return False

    def _flatten(karyo, order):
        genes = []
        for cid in order:
            genes.extend(karyo.get(cid, []))
        return genes

    true_order = sorted(true_root.keys(), key=str)
    true_flat = _flatten(true_root, true_order)
    true_pos = {g: i for i, g in enumerate(true_flat)}

    gene_to_true_chr = {}
    true_chr_boundaries = []
    current_x = 0
    cmap_true = plt.get_cmap("tab20")
    cmap_rec = plt.get_cmap("tab20b")
    true_colors = {cid: cmap_true(i % 20) for i, cid in enumerate(true_order)}
    for cid in true_order:
        start_x = current_x
        for g in true_root.get(cid, []):
            gene_to_true_chr[g] = cid
            current_x += 1
        end_x = current_x
        true_chr_boundaries.append((cid, start_x, end_x, true_colors[cid]))

    rec_genes_by_cid = rec_karyo
    rec_colors = {cid: cmap_rec(i % 20) for i, cid in enumerate(rec_order)}

    rec_flat = []
    rec_chr_boundaries = []
    current_y = 0
    for cid in rec_order:
        start_y = current_y
        for g in rec_genes_by_cid.get(cid, []):
            rec_flat.append(g)
            current_y += 1
        end_y = current_y
        rec_chr_boundaries.append((cid, start_y, end_y, rec_colors[cid]))
    rec_pos = {g: i for i, g in enumerate(rec_flat)}

    xs = []
    ys = []
    cs = []
    for g, x in true_pos.items():
        y = rec_pos.get(g)
        if y is None:
            continue
        xs.append(x)
        ys.append(y)
        cs.append(true_colors.get(gene_to_true_chr.get(g), "gray"))

    total_genes_x = len(true_flat)
    total_genes_y = len(rec_flat)

    if not xs:
        print("No matching genes between true root and reconstructed root.")
        return False

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    font_family = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.scatter(xs, ys, c=cs, s=8, marker=".", edgecolors="none", alpha=0.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_title("TrueRoot", fontsize=16, pad=10, y=1.04, fontweight="bold", fontfamily=font_family)
    ax.set_ylabel("ReconstructedRoot", fontsize=16, labelpad=80, fontweight="bold", fontfamily=font_family)

    ax.set_xlim(0, total_genes_x)
    ax.set_ylim(0, total_genes_y)
    ax.invert_yaxis()
    for (cid, sx, ex, col) in true_chr_boundaries:
        ax.axvline(x=sx, color="black", linestyle="-", linewidth=0.5)
        ax.axvline(x=ex, color="black", linestyle="-", linewidth=0.5)
        bar_height = max(1, int(total_genes_y * 0.02))
        rect = mpatches.Rectangle((sx, -bar_height), ex - sx, bar_height, color=col, clip_on=False)
        ax.add_patch(rect)
        ax.text((sx + ex) / 2, -bar_height * 1.5, cid, ha="center", va="bottom", fontsize=10, fontweight="bold", fontfamily=font_family)

    for (cid, sy, ey, col) in rec_chr_boundaries:
        ax.axhline(y=sy, color="black", linestyle="-", linewidth=0.5)
        ax.axhline(y=ey, color="black", linestyle="-", linewidth=0.5)
        ax.text(-0.02, (sy + ey) / 2, cid, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=10, fontweight="bold", fontfamily=font_family, clip_on=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Dot plot saved to {out_path}")
    return True


# -----------------------------
# 4. MAIN
# -----------------------------
def main():
    print("Generating root comparison dot plot...")
    success = generate_root_vs_trueroot_dotplot(
        reconstructed_tsv_path=CONFIG["reconstructed_tsv"],
        true_root_karyo_path=CONFIG["true_root_karyotype"],
        out_path=CONFIG["output_png"],
        root_label=CONFIG["root_label"],
    )
    if success:
        print("Done.")
    else:
        print("Failed to generate dot plot. Ensure reconstruction has been run first.")


if __name__ == "__main__":
    main()
