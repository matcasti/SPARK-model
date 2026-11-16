# ==============================================================================
# tools/prep/blockade_maneuvers.py -- hyperventilation and breath-hold times of the saline period of each session, from the free-text experiment logs.
#   python3 -I tools/prep/blockade_maneuvers.py <BlockerData.zip or logs folder> <output csv>
# For each subject and day: begin of hyperventilation (hv0), begin of breath-hold (bh0) and its end (bh1) in the FIRST cycle after the saline baselines, plus the time of the first
# drug administration (drug0), which closes the saline period. Reviewed by hand against the logs (the printed table).
# ==============================================================================
import re, sys, csv, os, zipfile
folder, out = sys.argv[1], sys.argv[2]
def read_log(src, s, day):
    """Raw bytes of the free-text log of subject s, day `day`: from the zip (member 'Experiment Logs/Subject<s>_day<d>.doc') or from a folder."""
    name = f'Subject{s}_day{day}.doc'
    if src.lower().endswith('.zip'):
        with zipfile.ZipFile(src) as z: return z.read('Experiment Logs/' + name)
    return open(os.path.join(src, name), 'rb').read()
EV = re.compile(r'^\s*(\d+)\*?\s+(-?\d*\.\d+)\s+.*?(\d\d-\d\d-\d\d)\s+(\d\d:\d\d:\d\d)\s*(.*)$')
rows = []
for day in (1, 2):
    for s in range(1, 8):
        raw = read_log(folder, s, day); ev = {}
        for m in re.findall(rb'[\x20-\x7e]{3,}', raw):
            mm = EV.match(m.decode('latin-1'))
            if mm: ev[(int(mm.group(1)), float(mm.group(2)))] = mm.group(5).strip().lower()
        events = sorted((t, txt) for (i, t), txt in ev.items())
        drug0 = next((t for t, x in events if re.search(r'atrop|propranolol|prop\b|prorpanolol|first dose|1st dos|dose now', x) and t > 600), None)
        hv0 = next((t for t, x in events if 'hyperv' in x and re.search(r'begin|start|hyperv', x) and 'end' not in x and t > 400), None)
        bh0 = next((t for t, x in events if re.search(r'breath ?hold|breathold|breahthold', x) and re.search(r'begin|start', x) and hv0 is not None and t > hv0), None)
        if bh0 is None and hv0 is not None: bh0 = next((t for t, x in events if re.search(r'breath ?hold|breathold|breahthold', x) and 'end' not in x and t > hv0), None)
        bh1 = next((t for t, x in events if (re.search(r'end|ned', x) and re.search(r'hold|breathold', x) or 'breathing again' in x) and bh0 is not None and t > bh0 and (drug0 is None or t < drug0)), None)   # 'ned' = typo of 'end' in subject 6
        rows.append(dict(subject=s, day=day, hv0=hv0, bh0=bh0, bh1=bh1, drug0=drug0))
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['subject', 'day', 'hv0', 'bh0', 'bh1', 'drug0']); w.writeheader(); [w.writerow(r) for r in rows]
for r in rows: print(r, '| breath-hold s:', None if r['bh0'] is None or r['bh1'] is None else round(r['bh1'] - r['bh0']))
