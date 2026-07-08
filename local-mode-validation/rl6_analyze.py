#!/usr/bin/env python3
"""Brief 26_0610-068: compare --ckpt vs stock Stockholm alignments.
Degapped-residue identity (accuracy-neutral test) + avg-PP delta per sequence."""
import sys, re

def parse_sto(path):
    """Return {name: (aligned_seq, pp_string)} concatenating multi-block."""
    seqs, pps = {}, {}
    order = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line == '//' or line.startswith('# STOCKHOLM'):
                continue
            if line.startswith('#=GR '):
                # #=GR <seqname> PP <ppstring>
                parts = line.split(None, 3)
                if len(parts) == 4 and parts[2] == 'PP':
                    pps[parts[1]] = pps.get(parts[1], '') + parts[3]
                continue
            if line.startswith('#'):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            name, aseq = parts
            if name not in seqs:
                order.append(name)
            seqs[name] = seqs.get(name, '') + aseq
    return order, seqs, pps

def degap(s):
    return re.sub(r'[-.]', '', s).upper()

# PP char -> approximate numeric (infernal scale: 0-9 -> 0.0-0.9 buckets, * = 0.95+)
def pp_avg(pp, aseq):
    vals = []
    for c, a in zip(pp, aseq):
        if a in '-.':
            continue
        if c == '*':
            vals.append(0.975)
        elif c.isdigit():
            vals.append(int(c)/10.0 + 0.05)
        # '.' in PP for gap positions skipped
    return sum(vals)/len(vals) if vals else float('nan')

def main(ck, st, label):
    ok, so, sp = parse_sto(ck)  # order, seqs, pps
    _,  to, tp = parse_sto(st)
    all_names = list(dict.fromkeys(list(so.keys())+list(to.keys())))
    n_deg_mismatch = 0
    maxppd = 0.0
    print(f"=== {label} ===")
    print(f"{'seq':<26} {'degap':<10} {'avgPP_ck':>9} {'avgPP_st':>9} {'dPP':>7}")
    for nm in all_names:
        a_ck = so.get(nm, ''); a_st = to.get(nm, '')
        d_ck = degap(a_ck); d_st = degap(a_st)
        deg = 'IDENT' if d_ck == d_st else 'DIFFER!!'
        if d_ck != d_st: n_deg_mismatch += 1
        p_ck = pp_avg(sp.get(nm,''), a_ck)
        p_st = pp_avg(tp.get(nm,''), a_st)
        dpp = p_ck - p_st
        if abs(dpp) == abs(dpp) and abs(dpp) > maxppd: maxppd = abs(dpp)
        print(f"{nm:<26} {deg:<10} {p_ck:>9.4f} {p_st:>9.4f} {dpp:>+7.4f}")
    print(f"  --> degapped mismatches: {n_deg_mismatch}/{len(all_names)}   max|avgPP delta|: {maxppd:.4f}")
    print()
    return n_deg_mismatch

if __name__ == '__main__':
    tot = 0
    for label, ck, st in [
        ('LOCAL non-trunc (rung-3)',  'rl6-data/cmp_local_notrunc_ck.sto',  'rl6-data/cmp_local_notrunc_st.sto'),
        ('LOCAL trunc (rung-4)',      'rl6-data/cmp_local_trunc_ck.sto',    'rl6-data/cmp_local_trunc_st.sto'),
        ('GLOBAL non-trunc (rung-3)', 'rl6-data/cmp_global_notrunc_ck.sto', 'rl6-data/cmp_global_notrunc_st.sto'),
        ('GLOBAL trunc (rung-4)',     'rl6-data/cmp_global_trunc_ck.sto',   'rl6-data/cmp_global_trunc_st.sto'),
    ]:
        tot += main(ck, st, label)
    print(f"TOTAL degapped-residue mismatches across all configs: {tot} (0 = accuracy-neutral, residues preserved)")
