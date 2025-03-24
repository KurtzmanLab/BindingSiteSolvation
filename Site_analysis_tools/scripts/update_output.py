import pandas as pd

pd.set_option('future.no_silent_downcasting', True)

def update_totals(results):

    df_water = pd.DataFrame(columns=["Solute to Water", "Water to Solute", "Water Neighbors"])
    df_conf = pd.DataFrame(columns=["Rigid", "Flexible"])

    total_accept = {}
    total_donate = {}
    total_neighbors = {}
    rigid = {}
    flexible = {}

    for result in results:
        don = result['don']
        acc = result['acc']
        neighbors = result['neighbors']
        conformation = result['conf']

        for atom, ndon in don.items():
            df_water.loc[atom, "Solute to Water"] = total_donate.get(atom, 0) + ndon
            total_donate[atom] = total_donate.get(atom, 0) + ndon
        for atom, nacc in acc.items():
            df_water.loc[atom, "Water to Solute"] = total_accept.get(atom, 0) + nacc
            total_accept[atom] = total_accept.get(atom, 0) + nacc
        for atom, neigh in neighbors.items():
            df_water.loc[atom, "Water Neighbors"] = total_neighbors.get(atom, 0) + neigh
            total_neighbors[atom] = total_neighbors.get(atom, 0) + neigh

        for res, conf in conformation.items():
            res = str(res)
            ahr, bbr = 0, 0
            if conf == "ahr":
                ahr += 1
            elif conf == "bbr":
                bbr += 1

            if res not in rigid.keys():
                rigid[res] = ahr
            elif res in rigid.keys():
                rigid[res] += ahr

            if res not in flexible.keys():
                flexible[res] = bbr
            elif res in flexible.keys():
                flexible[res] += bbr


    for res in rigid.keys():
        df_conf.loc[res, "Rigid"] = rigid[res]
        df_conf.loc[res, "Flexible"] = flexible[res]


    df_water = df_water.fillna(0)

    return df_water, df_conf