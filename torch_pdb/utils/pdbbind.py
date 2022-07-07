import re

import torch

def affinity_parse(s):
    """ Parse the affinity string. e.g. `Kd=30uM`.
    Parameters
    ----------
    s: str
        Affinity measurement string to parse.
    Returns
    -------
    dict
        Dictionary containing parsed affinity information. `value` key stores
        the float value of the measurement. `operator` is the logical operator
        (e.g. `=`, `>`) applied to the value, `unit` is `uM, nM, pM` and
        `measure` is the type experimental measurement (e.g. `Kd, Ki, IC50`)
    """
    operator = "".join(re.findall(r"[=|<|>|~]", s))
    measures = ['Kd', 'Ki', 'IC50']
    for m in measures:
        if s.startswith(m):
            measure = m
            break
    value = float(re.search(r"\d+[.,]?\d*", s).group())
    unit = re.search(r"[m|u|n|f|p]M", s).group()

    return {'operator': operator,
            'measure': measure,
            'value': value,
            'unit': unit
            }

def parse_pdbbind_PL_index(index_path):
    """
    # ==============================================================================
    # List of protein-ligand complexes with known binding data in PDBbind v.2020
    # 19443 protein-ligand complexes in total, sorted by their release year.
    # Latest update: July 2021
    # PDB code, resolution, release year, binding data, reference, ligand name
    # ==============================================================================
    2tpi  2.10  1982  Kd=49uM       // 2tpi.pdf (2-mer)
    5tln  2.30  1982  Ki=0.43uM     // 5tln.pdf (BAN) incomplete ligand structure
    4tln  2.30  1982  Ki=190uM      // 4tln.pdf (LNO)
    4cts  2.90  1984  Kd<10uM       // 4cts.pdf (OAA)
    6rsa   NMR  1986  Ki=40uM       // 6rsa.pdf (UVC)
    1rnt  1.90  1987  Kd=6.5uM      // 1rnt.pdf (2GP)
    """
    data = {}
    with open(index_path, 'r') as ind_file:
        for line in ind_file:
            if line.startswith("#"):
                continue
            pre, post = line.split("//")
            pdbid, res, date, kd = pre.split()
            kd = affinity_parse(kd)

            lig_id = post.split("(")[1].rstrip(")")
            data[pdbid] = {'resolution': float(res),
                           'date': int(date),
                           'kd': kd,
                           'ligand_id': lig_id
                           }

    return data
