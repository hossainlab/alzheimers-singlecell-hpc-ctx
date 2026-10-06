import urllib.request
import json
import os

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=200163577,304982107,304982106,304982105,304982104,304982103,304982102,304982101,304982100,304982099,304982098,304982097,304982096,304982095,304982094,304982093,304982092,304982091,304982090,304982089,304982088,304982087,304982086,304982085,304982084,304982083&retmode=json'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as resp:
    data = json.loads(resp.read().decode('utf-8'))

samples = []
for uid in data['result']['uids']:
    item = data['result'][uid]
    acc = item.get('accession', '')
    title = item.get('title', '')
    summary = item.get('summary', '')
    sample_info = {
        'uid': uid,
        'accession': acc,
        'title': title,
        'summary': summary
    }
    samples.append(sample_info)
    print(f"{acc} | {title} | {summary[:60]}")

with open('tables/geo_metadata_summary.json', 'w') as f:
    json.dump(data, f, indent=2)
