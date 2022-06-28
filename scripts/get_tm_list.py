import re
import requests

lines  = requests.get("https://zhanggroup.org/TM-align/benchmark/").text.split("\n")

pdblist = []
for l in lines:
    m = re.search(".pdb", l)
    if m:
        start, end = m.span()
        pdbid = l[start-5:end]
        pdblist.append(f"https://zhanggroup.org/TM-align/benchmark/{pdbid}")

with open("tm_pdblist.txt", 'w') as tm:
    for p in pdblist:
        tm.write(p + "\n")
