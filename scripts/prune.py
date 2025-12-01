
import sys
from csv import DictReader

# Main
#----------------------------------------------------------------------

rows = DictReader(sys.stdin)
cols = rows.fieldnames

print(*cols, sep = ',')
for row in rows :
    s = set(row[x].strip() for x in cols[1:])

    if len(s) > 1 : # there is some variation
        print(*(row[x].strip() for x in cols), sep = ',')
