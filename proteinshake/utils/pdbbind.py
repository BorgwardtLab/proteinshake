import re

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
    > INDEX_refined_data.2020
    # ==============================================================================
    # List of the protein-ligand complexes in the PDBbind refined set v.2020
    # 5316 protein-ligand complexes in total, which are ranked by binding data
    # Latest update: July 2021
    # PDB code, resolution, release year, -logKd/Ki, Kd/Ki, reference, ligand name
    # ==============================================================================
    2r58  2.00  2007   2.00  Kd=10mM       // 2r58.pdf (MLY)
    3c2f  2.35  2008   2.00  Kd=10.1mM     // 3c2f.pdf (PRP)
    3g2y  1.31  2009   2.00  Ki=10mM       // 3g2y.pdf (GF4)
    3pce  2.06  1998   2.00  Ki=10mM       // 3pce.pdf (3HP)
    4qsu  1.90  2014   2.00  Kd=10mM       // 4qsu.pdf (TDR)
    4qsv  1.90  2014   2.00  Kd=10mM       // 4qsv.pdf (THM)
    """
    data = {}
    with open(index_path, 'r') as ind_file:
        for line in ind_file:
            if line.startswith("#"):
                continue
            pre, post = line.split("//")
            pdbid, res, date, neglog, kd = pre.split()
            kd = affinity_parse(kd)

            lig_id = post.split("(")[1].rstrip(")")

            # remove peptide ligands
            if lig_id.endswith('-mer'):
                continue
            data[pdbid] = {'resolution': float(res),
                           'date': int(date),
                           'kd': kd,
                           'neglog_aff': float(neglog),
                           'ligand_id': lig_id
                           }

    return data
