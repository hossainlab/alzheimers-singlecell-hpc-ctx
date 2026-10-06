import csv
import re

with open('tables/consensus_target_prioritization_full.csv', 'r') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

print("=== TOP 10 OVERALL CANDIDATES ===")
print(f"{'Gene':<15} | {'Score':<8} | {'Cell Breadth':<12} | {'Log2FC':<8} | {'-Log10 Padj':<12} | {'PPI Deg':<8}")
print("-" * 75)
for r in rows[:10]:
    print(f"{r['gene']:<15} | {float(r['Consensus_Score']):<8.3f} | {r['cell_type_breadth']:<12} | {float(r['mean_abs_log2FC']):<8.3f} | {float(r['mean_neg_log10_padj']):<12.1f} | {r['ppi_degree']:<8}")

canonical = [r for r in rows if not re.match(r'^(AC\d|AL\d|AP\d|LINC\d|MIR|FP\d|BX\d|CR\d|CTD\d)', r['gene'])]

print("\n=== TOP 10 CANONICAL / PROTEIN-CODING TARGETS ===")
print(f"{'Gene':<15} | {'Score':<8} | {'Cell Breadth':<12} | {'Log2FC':<8} | {'-Log10 Padj':<12} | {'PPI Deg':<8}")
print("-" * 75)
for r in canonical[:10]:
    print(f"{r['gene']:<15} | {float(r['Consensus_Score']):<8.3f} | {r['cell_type_breadth']:<12} | {float(r['mean_abs_log2FC']):<8.3f} | {float(r['mean_neg_log10_padj']):<12.1f} | {r['ppi_degree']:<8}")

with open('tables/top5_canonical_protein_targets.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=reader.fieldnames)
    writer.writeheader()
    writer.writerows(canonical[:5])
