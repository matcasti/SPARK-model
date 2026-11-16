# ==============================================================================
# tools/prep/blockade_msna.py -- the integrated nerve signal (muscle sympathetic nerve activity, MSNA) of the Lancaster data, read straight from the zip and
# reduced to 1 s means (the 500 Hz signal is never written to disk).   python3 -I tools/prep/blockade_msna.py <BlockerData.zip> <output folder>
# ==============================================================================
import sys, io, os, zipfile
import numpy as np, pandas as pd
zpath, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True); z = zipfile.ZipFile(zpath)
for info in sorted(z.infolist(), key=lambda i: i.filename):
    if not info.filename.endswith('/Integrated Nerve Signal.dat'): continue
    name = info.filename.split('/')[0]
    with z.open(info) as f: raw = f.read()
    raw = raw.replace(b'\x00', b'')                                         # some files end in NUL padding
    df = pd.read_csv(io.BytesIO(raw), sep='\t', header=None, names=['t', 'v']).apply(pd.to_numeric, errors='coerce').dropna(); del raw
    sec = np.floor(df.t.values).astype(int); g = df.groupby(sec).v.agg(['mean', 'std']).reset_index(); g.columns = ['t', 'mean', 'sd']
    g.to_csv(os.path.join(out, f'msna_{name}.csv'), index=False, float_format='%.5g'); print(name, len(g), 'seconds; median level', round(float(g['mean'].median()), 2), flush=True)
