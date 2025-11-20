
import sys

# set of SNPs in file named f to a dictionary d[(chr,tr,pos)] = base
def to_dict(f) :

    d = {}
    for line in open(f,'r') :
        c, p, b = line.strip().split()

        c, t = c.split('.')
        c = int(c.lstrip('g'))
        t = int(t.lstrip('t'))
        p = int(p)
        assert b in ['A', 'C', 'G', 'T']

        d[(c,t,p)] = b

    return d

# Main
#----------------------------------------------------------------------

# get universal set ss of SRRs and corresponding set of dictionaries
ss = []
ds = {}
for f in sys.argv[1:] :
    s, _ = f.split('.',1)
    _, s = s.rsplit('/',1)
    ss.append(s)

    ds[s] = to_dict(f)

# build universal set ctps of (chr,tr,pos) tuples
ctps = set()
for s in ss :
    for ctp in ds[s] :
        ctps.add(ctp)

# dump CSV with ctps as rows and ss as columns
print('chr.tr.pos\s', *ss, sep = ',')
for ctp in sorted(ctps) :

    # a row of the CSV
    xs = []
    for s in ss :
        x = '-'

        if ctp in ds[s] :
            x = ds[s][ctp]

        xs.append(x)

    if len(set(xs)) > 1 : # a row with some differentiation
        c,t,p = ctp

        print('{}.{}.{}'.format(c,t,p), *xs, sep = ',')
