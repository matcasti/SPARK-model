# ==============================================================================
# tools/prep/blockade_rpeaks.py -- R-peak times from the ECG of the Lancaster "Blocker data", read straight from the zip (the other signals are never extracted).
#   python3 -I tools/prep/blockade_rpeaks.py <BlockerData.zip> <output folder>
# ECG.dat is two tab-separated ASCII columns (time s, value) at 500 Hz. Detector: 5-20 Hz band-pass, squared and integrated (120 ms), adaptive threshold,
# refined to the largest deflection within +-60 ms (works for either polarity). Writes rpeaks_Subject<s>_day<d>.csv and qc.csv (beats, median RR, RR range).
# ==============================================================================
import sys, io, os, zipfile, csv
import numpy as np, pandas as pd
from scipy import signal
zpath, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
z = zipfile.ZipFile(zpath); qc = []
for info in sorted(z.infolist(), key=lambda i: i.filename):
    if not info.filename.lower().endswith('/ecg.dat'): continue
    name = info.filename.split('/')[0]
    with z.open(info) as f: raw = f.read()                                   # CRC is verified by zipfile at the end of the read
    raw = raw.replace(b'\x00', b'')                                         # some files end in NUL padding
    df = pd.read_csv(io.BytesIO(raw), sep='\t', header=None, names=['t', 'v']).apply(pd.to_numeric, errors='coerce').dropna(); del raw
    t = df.t.values; x = df.v.values; fs = 1 / np.median(np.diff(t))
    b, a = signal.butter(2, [5 / (fs / 2), 20 / (fs / 2)], btype='band'); y = signal.filtfilt(b, a, x)
    e = np.convolve(y ** 2, np.ones(int(0.12 * fs)) / int(0.12 * fs), mode='same')
    thr = np.percentile(e, 99) * 0.15
    pk, _ = signal.find_peaks(e, height=thr, distance=int(0.30 * fs)); w = int(0.06 * fs)
    r = np.array(sorted({max(p - w, 0) + int(np.argmax(np.abs(y[max(p - w, 0):min(p + w, len(y) - 1)]))) for p in pk}))
    tr = t[r]; rr = np.diff(tr)
    np.savetxt(os.path.join(out, f'rpeaks_{name}.csv'), tr, fmt='%.4f')
    qc.append(dict(recording=name, seconds=round(t[-1]), fs=round(fs, 1), beats=len(tr), rr_median=round(float(np.median(rr)), 3), rr_min=round(float(rr.min()), 3), rr_max=round(float(rr.max()), 3),
                   rr_lt_0p3=int((rr < 0.3).sum()), rr_gt_2p5=int((rr > 2.5).sum()), ecg_median=round(float(np.median(x)), 1), ecg_p1=round(float(np.percentile(x, 1)), 1), ecg_p99=round(float(np.percentile(x, 99)), 1)))
    print(qc[-1], flush=True)
with open(os.path.join(out, 'qc.csv'), 'w', newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=list(qc[0].keys())); w_.writeheader(); [w_.writerow(q) for q in qc]
