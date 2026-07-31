#!/usr/bin/env python3
"""obj_merge.py OUT.obj IN1.obj IN2.obj ...

Concatenate several world-coordinate OBJs into one (offsetting face indices), for
visualising a cube + its neighbours together. Verts may carry trailing colours;
they're passed through (mesh_render ignores them). Viz helper only.
"""
import sys

def main():
    if len(sys.argv) < 3:
        print("usage: obj_merge.py OUT.obj IN1.obj [IN2.obj ...]"); return 1
    out, ins = sys.argv[1], sys.argv[2:]
    vlines, flines, voff = [], [], 0
    for fn in ins:
        n = 0
        try:
            fh = open(fn)
        except OSError:
            sys.stderr.write("  skip (missing): %s\n" % fn); continue
        for line in fh:
            if line.startswith('v '):
                vlines.append(line); n += 1
            elif line.startswith('f '):
                idx = [int(p.split('/')[0]) for p in line.split()[1:]]
                flines.append('f ' + ' '.join(str(i + voff) for i in idx) + '\n')
        fh.close()
        sys.stderr.write("  +%-55s %d verts\n" % (fn.split('/')[-1], n))
        voff += n
    with open(out, 'w') as o:
        o.writelines(vlines); o.writelines(flines)
    sys.stderr.write("merged -> %s  (%d verts, %d faces)\n" % (out, voff, len(flines)))
    return 0

if __name__ == "__main__":
    sys.exit(main())
