#!/usr/bin/env python3
"""Analyze how simulation cost depends on number of cells."""

import sys

# Parse timing file
tlines = open('monitor/timing').readlines()
theader = tlines[0].split()
tdata = {}
for line in tlines[1:]:
    vals = line.split()
    if len(vals) == len(theader):
        row = {h: float(v) for h, v in zip(theader, vals)}
        step = int(row['Timestep'])
        if step >= 11:  # skip init+startup
            tdata[step] = row

# Parse grid file
glines = open('monitor/grid').readlines()
# Header spans 2 lines for grid file
gheader_line1 = glines[0].split()
# Columns: Timestep, Time, Nlvl, Nbox, Ncell, Compression, Maximum RSS, Minimum RSS, Average RSS
gheader = ['Timestep', 'Time', 'Nlvl', 'Nbox', 'Ncell', 'Compression', 'RSS_max', 'RSS_min', 'RSS_avg']
gdata = {}
for line in glines[1:]:
    vals = line.split()
    if len(vals) >= 6:
        step = int(float(vals[0]))
        if step >= 11:
            gdata[step] = {gheader[i]: float(vals[i]) for i in range(min(len(vals), len(gheader)))}

# Match: for each grid step, find the closest timing step
# Grid is written every N steps; use the grid's Ncell for surrounding timing steps
grid_steps = sorted(gdata.keys())

# Interpolate Ncell for each timing step
def get_ncell(tstep):
    """Get Ncell for a timing step by interpolating from grid data."""
    if tstep in gdata:
        return gdata[tstep]['Ncell']
    # Find bracketing grid steps
    lo = max((g for g in grid_steps if g <= tstep), default=None)
    hi = min((g for g in grid_steps if g >= tstep), default=None)
    if lo is None and hi is None:
        return None
    if lo is None:
        return gdata[hi]['Ncell']
    if hi is None:
        return gdata[lo]['Ncell']
    if lo == hi:
        return gdata[lo]['Ncell']
    # Linear interpolation
    frac = (tstep - lo) / (hi - lo)
    return gdata[lo]['Ncell'] * (1-frac) + gdata[hi]['Ncell'] * frac

# Build matched dataset
steps = sorted(tdata.keys())
matched = []
for s in steps:
    nc = get_ncell(s)
    if nc is not None:
        matched.append((nc, tdata[s]))

print(f'Matched {len(matched)} timesteps with cell count data')
print(f'Cell count range: {min(m[0] for m in matched):.2e} to {max(m[0] for m in matched):.2e}')
print()

# Bin by cell count
ncells = [m[0] for m in matched]
nc_min, nc_max = min(ncells), max(ncells)
nbins = 8
bin_edges = [nc_min + i*(nc_max-nc_min)/nbins for i in range(nbins+1)]

def mean(v): return sum(v)/len(v) if v else 0

bins = [[] for _ in range(nbins)]
for nc, row in matched:
    for i in range(nbins):
        if bin_edges[i] <= nc < bin_edges[i+1] or (i == nbins-1 and nc == bin_edges[i+1]):
            bins[i].append((nc, row))
            break

print('='*110)
print(f'{"Ncell_avg":>12s} {"count":>6s} | {"dQdt":>8s} {"prim":>8s} {"SL":>8s} {"FV":>8s} {"div":>8s} {"plic":>8s} {"visc":>8s} | {"total":>8s}')
print('='*110)
for i, b in enumerate(bins):
    if not b:
        continue
    nc_avg = mean([x[0] for x in b])
    dqdt = mean([x[1]['dQdt_max'] for x in b])
    prim = mean([x[1]['prim_max'] for x in b])
    sl   = mean([x[1]['sl_max'] for x in b])
    fv   = mean([x[1]['fv_max'] for x in b])
    div  = mean([x[1]['div_max'] for x in b])
    plic = mean([x[1]['plic_max'] for x in b])
    visc = mean([x[1]['visc_max'] for x in b])
    total = dqdt + plic + visc
    print(f'{nc_avg:12.0f} {len(b):6d} | {dqdt:8.3f} {prim:8.3f} {sl:8.3f} {fv:8.3f} {div:8.3f} {plic:8.3f} {visc:8.3f} | {total:8.3f}')

# Correlation analysis
print()
print('='*110)
print('CORRELATION: Pearson r of each phase with Ncell')
print('='*110)

def pearson(x, y):
    n = len(x)
    if n < 2: return 0
    mx, my = mean(x), mean(y)
    sx = (sum((xi-mx)**2 for xi in x)/(n-1))**0.5
    sy = (sum((yi-my)**2 for yi in y)/(n-1))**0.5
    if sx == 0 or sy == 0: return 0
    cov = sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))/(n-1)
    return cov/(sx*sy)

nc_list = [m[0] for m in matched]
for label, key in [('dQdt','dQdt_max'), ('prim','prim_max'), ('SL','sl_max'),
                    ('FV','fv_max'), ('div','div_max'), ('plic','plic_max'), ('visc','visc_max')]:
    vals = [m[1][key] for m in matched]
    r = pearson(nc_list, vals)
    print(f'  {label:>8s} vs Ncell: r = {r:.3f}')

# Per-cell cost
print()
print('='*110)
print('PER-CELL COST (microseconds per cell per timestep)')
print('='*110)
print(f'{"Ncell_avg":>12s} {"count":>6s} | {"dQdt/cell":>10s} {"prim/cell":>10s} {"SL/cell":>10s} {"FV/cell":>10s} {"div/cell":>10s}')
print('-'*110)
for i, b in enumerate(bins):
    if not b:
        continue
    nc_avg = mean([x[0] for x in b])
    dqdt = mean([x[1]['dQdt_max']/x[0]*1e6 for x in b])
    prim = mean([x[1]['prim_max']/x[0]*1e6 for x in b])
    sl   = mean([x[1]['sl_max']/x[0]*1e6 for x in b])
    fv   = mean([x[1]['fv_max']/x[0]*1e6 for x in b])
    div  = mean([x[1]['div_max']/x[0]*1e6 for x in b])
    print(f'{nc_avg:12.0f} {len(b):6d} | {dqdt:10.3f} {prim:10.3f} {sl:10.3f} {fv:10.3f} {div:10.3f}')

# Mixed cell dependence
if 'mixed_max' in theader:
    print()
    print('='*110)
    print('MIXED CELLS vs SL COST')
    print('='*110)
    mixed_list = [m[1]['mixed_max'] for m in matched]
    sl_list = [m[1]['sl_max'] for m in matched]
    r = pearson(mixed_list, sl_list)
    print(f'  Pearson r (mixed_max vs sl_max) = {r:.3f}')
    # Bin by mixed cells
    mx_min, mx_max = min(mixed_list), max(mixed_list)
    mbins = 6
    mbin_edges = [mx_min + i*(mx_max-mx_min)/mbins for i in range(mbins+1)]
    mbinned = [[] for _ in range(mbins)]
    for m in matched:
        mx = m[1]['mixed_max']
        for i in range(mbins):
            if mbin_edges[i] <= mx < mbin_edges[i+1] or (i == mbins-1 and mx == mbin_edges[i+1]):
                mbinned[i].append(m)
                break
    print(f'{"mixed_avg":>12s} {"count":>6s} | {"SL_max":>8s} {"FV_max":>8s} {"dQdt":>8s}')
    for i, b in enumerate(mbinned):
        if not b:
            continue
        mx_avg = mean([x[1]['mixed_max'] for x in b])
        sl = mean([x[1]['sl_max'] for x in b])
        fv = mean([x[1]['fv_max'] for x in b])
        dqdt = mean([x[1]['dQdt_max'] for x in b])
        print(f'{mx_avg:12.0f} {len(b):6d} | {sl:8.3f} {fv:8.3f} {dqdt:8.3f}')
