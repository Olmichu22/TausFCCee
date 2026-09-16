#!/usr/bin/env python3
"""Compara bin a bin todos los histogramas de dos ficheros ROOT (recursivo)."""
import sys
import ROOT

ROOT.gROOT.SetBatch(True)


def walk(d, path=""):
    for k in d.GetListOfKeys():
        obj = k.ReadObj()
        name = f"{path}/{k.GetName()}"
        if obj.InheritsFrom("TDirectory"):
            yield from walk(obj, name)
        elif obj.InheritsFrom("TH1"):
            obj.SetDirectory(0)
            yield name, obj


def main(a, b):
    fa, fb = ROOT.TFile.Open(a), ROOT.TFile.Open(b)
    ha, hb = dict(walk(fa)), dict(walk(fb))
    print(f"A: {len(ha)} hists, B: {len(hb)} hists")
    only = set(ha) ^ set(hb)
    if only:
        print("Solo en uno:", sorted(only)[:20])
    ndiff = 0
    maxrel = 0.0
    for name in sorted(set(ha) & set(hb)):
        x, y = ha[name], hb[name]
        n = x.GetNcells()
        if n != y.GetNcells():
            print("NCELLS difiere:", name); ndiff += 1; continue
        for i in range(n):
            va, vb = x.GetBinContent(i), y.GetBinContent(i)
            ea, eb = x.GetBinError(i), y.GetBinError(i)
            if va != vb or ea != eb:
                rel = abs(va - vb) / max(abs(va), abs(vb), 1e-300)
                maxrel = max(maxrel, rel)
                if rel > 1e-9 or abs(ea - eb) > 1e-9 * max(ea, eb, 1e-300):
                    ndiff += 1
                    if ndiff <= 10:
                        print(f"DIFF {name} bin {i}: {va} vs {vb} (err {ea} vs {eb})")
                    break
    print(f"Histogramas con diferencias (>1e-9 rel): {ndiff}; max rel diff: {maxrel:.3g}")
    return 1 if (ndiff or only) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
