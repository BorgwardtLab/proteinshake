import os
import re
import requests

PATH = os.path.realpath(os.path.dirname(__file__))

lines  = requests.get("https://zhanggroup.org/TM-align/benchmark/").text.split("\n")

pdblist = []
for l in lines:
    m = re.search(".pdb", l)
    if m:
        start, end = m.span()
        pdbid = l[start-5:end]
        pdblist.append(f"https://zhanggroup.org/TM-align/benchmark/{pdbid}")

with open(os.path.join(PATH, "..", "torch_pdb", "pkg_data", "tm_pdblist.txt"), 'w') as tm:
    for p in pdblist:
        tm.write(p + "\n")
