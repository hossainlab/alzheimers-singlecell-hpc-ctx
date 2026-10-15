import os
import glob
import json

with open('tables/geo_metadata_summary.json', 'r') as f:
    geo_data = json.load(f)

gsm_map = {}
for uid in geo_data['result']['uids']:
    item = geo_data['result'][uid]
    acc = item.get('accession', '')
    if acc.startswith('GSM'):
        gsm_map[acc] = item.get('title', '')

raw_pattern = 'data/GSE163577_RAW/*.tar.gz' if os.path.exists('data/GSE163577_RAW') else 'GSE163577_RAW/*.tar.gz'
files = sorted(glob.glob(raw_pattern))
parsed_samples = []

for fpath in files:
    fname = os.path.basename(fpath)
    gsm = fname.split('_')[0]
    title = gsm_map.get(gsm, 'Unknown')
    
    # Check region: Ctx vs Hpc
    if 'Ctx' in title or 'Ctx' in fname:
        region = 'Cortex (Ctx)'
    elif 'Hpc' in title or 'Hpc' in fname:
        region = 'Hippocampus (Hpc)'
    else:
        region = 'Unknown'
        
    # Check condition: AD vs Control
    # Look at title first
    if title.startswith('AD') or ' AD ' in title:
        cond = 'AD'
    elif title.startswith('Control') or ' Control ' in title:
        cond = 'Control'
    else:
        if '_AD_' in fname or 'AD_Ctx' in fname:
            cond = 'AD'
        elif '_C_' in fname or '_C_Ctx' in fname or '1OC_Ctx' in fname or 'O3C_Ctx' in fname or 'O9C_Ctx' in fname:
            cond = 'Control'
        else:
            cond = 'Unknown'

    parsed_samples.append({
        'gsm': gsm,
        'title': title,
        'filename': fname,
        'region': region,
        'condition': cond,
        'size_mb': round(os.path.getsize(fpath) / (1024 * 1024), 2)
    })

print(f"Total samples: {len(parsed_samples)}")
print("-" * 80)
print(f"{'GSM':<12} | {'Title':<16} | {'Region':<18} | {'Condition':<10} | {'Size (MB)':<10}")
print("-" * 80)
for s in parsed_samples:
    print(f"{s['gsm']:<12} | {s['title']:<16} | {s['region']:<18} | {s['condition']:<10} | {s['size_mb']:<10}")

with open('tables/sample_metadata.json', 'w') as f:
    json.dump(parsed_samples, f, indent=2)

# Summaries
import collections
region_counts = collections.Counter(s['region'] for s in parsed_samples)
cond_counts = collections.Counter((s['region'], s['condition']) for s in parsed_samples)

print("\n=== Sample Counts by Region and Condition ===")
for (reg, c), count in sorted(cond_counts.items()):
    print(f"{reg} - {c}: {count} samples")
