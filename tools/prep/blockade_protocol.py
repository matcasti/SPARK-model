# ==============================================================================
# tools/prep/blockade_protocol.py -- segment table (subject, day, drug condition, breathing, start, end) from the experiment logs of the Lancaster
# "Blocker data" (VCU, 1998). The logs are free-text Word 6 files: only their printable text is read (nothing is executed or parsed as a document).
#   python3 -I tools/prep/blockade_protocol.py <BlockerData.zip or folder with the 14 logs> <output csv>
# Rules: a segment starts at an event "begin/start ... baseline | uncontrolled | 0.25 Hz | 0.1 Hz" and lasts 300 s (the logged ends are 300-305 s later).
# Condition: written in the event text (saline, low dose atropine, high dose atropine / atropine, propranolol) or, if absent, the last drug state:
#   day 1 saline -> low-dose atropine -> high-dose atropine -> propranolol on top (= double blockade); day 2 saline -> propranolol alone.
# ==============================================================================
import re, sys, csv, os, zipfile, glob
folder, out = sys.argv[1], sys.argv[2]
def read_log(src, s, day):
    """Raw bytes of the free-text log of subject s, day `day`: from the zip (member 'Experiment Logs/Subject<s>_day<d>.doc') or from a folder."""
    name = f'Subject{s}_day{day}.doc'
    if src.lower().endswith('.zip'):
        with zipfile.ZipFile(src) as z: return z.read('Experiment Logs/' + name)
    return open(os.path.join(src, name), 'rb').read()
EV = re.compile(r'^\s*(\d+)\*?\s+(-?\d*\.\d+)\s+.*?(\d\d-\d\d-\d\d)\s+(\d\d:\d\d:\d\d)\s*(.*)$')
rows = []; problems = []
for day in (1, 2):
    for s in range(1, 8):
        raw = read_log(folder, s, day)
        lines = [m.decode('latin-1') for m in re.findall(rb'[\x20-\x7e]{3,}', raw)]
        ev = {}
        for ln in lines:
            m = EV.match(ln)
            if m: ev[(int(m.group(1)), float(m.group(2)))] = m.group(5).strip()
        events = sorted(((t, i, txt) for (i, t), txt in ev.items()))
        state = 'saline'
        for t, i, txt in events:
            low = re.sub(r'bas(e)?line|basl?ine|basline', 'baseline', txt.lower()).replace('salline', 'saline').replace('prorpanolol', 'propranolol').replace('propranalol', 'propranolol')
            if day == 1:
                if 'low dose' in low and ('in now' in low or low.strip().startswith('low dose') or 'atropine' in low) and 'begin' not in low and 'baseline' not in low and '0.25' not in low and '.25' not in low:
                    state = 'atropine_low'
                elif ('first dose atropine' in low or 'first atropine' in low or 'atropine infusion' in low or low.startswith('1 ml atropine') or low.startswith('atropine 1ml') or '1ml atropine' in low
                      or 'ml atropine' in low or 'dose of atropine' in low or 'dose atropine' in low) and 'low dose' not in low and state in ('atropine_low', 'saline'):
                    state = 'atropine_high'
                if re.search(r'propranolol|prorpanolol|propranalol|2 cc prop|first propr|prop\b', low) and state in ('atropine_high', 'atropine_low') and not re.search(r'begin|start|baseline|0\.25|\.25|hyperv|breath', low):
                    state = 'double'
            else:
                if re.search(r'propranolol|prop\b|dose now|first dose|1st dos|rest of prop|all in now|2\.5 ml|1st dose|2nd dose|3rd dose|4th ddose', low) and not re.search(r'begin|start|baseline|\.25|hyperv|breath', low) and t > 1000:
                    state = 'propranolol'
            is_start = re.search(r'\b(begin|start)\b', low) is not None or low.startswith('baseline-s') or (low.startswith('uncontrolled') and 'end' not in low) or re.match(r'^0\.1 ?hz', low) is not None or re.match(r'^(\.25|0\.25) ?hz', low) is not None or low.startswith('baseline')
            is_end = re.search(r'\b(end|stop)\b', low) is not None
            if not is_start or is_end and not re.search(r'\b(begin|start)\b', low): continue
            if re.search(r'hyperv|breath ?hold|breathold|breahthold|valsalva|cal\b|doppler|admin', low): continue
            if re.search(r'0\.25|\.25 ?hz', low): breathing = 'paced_0.25'
            elif re.search(r'0\.1 ?hz', low): breathing = 'paced_0.1'
            elif 'baseline' in low or 'uncontrolled' in low or low.startswith('baseline-s'):
                breathing = 'spontaneous'
            else: continue
            cond = None
            if 'saline' in low or 'salline' in low or '-s ' in low + ' ' or low.endswith('-s begin'): cond = 'saline'
            elif 'low dose' in low: cond = 'atropine_low'
            elif 'high dose' in low and 'atrop' in low: cond = 'atropine_high'
            elif 'prop' in low: cond = 'double' if day == 1 else 'propranolol'
            elif 'atrop' in low: cond = 'atropine_high'
            else: cond = state if state != 'saline' or day == 2 else state
            if cond == 'saline' and state not in ('saline',) and not ('saline' in low or 'salline' in low): cond = state
            rows.append(dict(subject=s, day=day, cond=cond, breathing=breathing, t0=t, t1=t + 300.0, text=txt[:60]))
        # no check here: the table is printed and reviewed
# drops: calibration or aborted starts that the log itself supersedes (reviewed by hand against the logs)
DROP = {(3, 2, 673.8), (5, 1, 1241.0), (7, 1, 158.2)}
rows = [r for r in rows if r['t0'] >= 90 and (r['subject'], r['day'], round(r['t0'], 1)) not in DROP]
rows.sort(key=lambda r: (r['subject'], r['day'], r['t0']))
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['subject', 'day', 'cond', 'breathing', 't0', 't1', 'text']); w.writeheader(); [w.writerow(r) for r in rows]
print(len(rows), 'segments written to', out)
