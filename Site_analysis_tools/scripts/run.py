import time
import multiprocessing
from scripts.update_output import update_totals
from src.hbond_analysis import HbondAnalysis
from src.side_chain_analysis import SideChainAnalysis

target = ""				# target name
resi_list = [242, 245]			# residue indices
num_frames = range(80000, 100000)	# range(start_frame, stop_frame)


ahr_avg_file = target + "_ahr_avg.pdb"	# Rigid time-average file (.pdb)	
bbr_avg_file = target + "_bbr_avg.pdb"	# Flexible time-average file (.pdb)

path = ""
top_file = path + target + ".prmtop"	# target topology file
traj_file = path + target + ".nc"	# target trajectory file

Conf_Analysis = SideChainAnalysis(top_file, traj_file, ahr_avg_file, bbr_avg_file, resi_list)
Hbond = HbondAnalysis(top_file, traj_file, resi_list)


def process_frame(i):
    don = Hbond.num_of_donate(i)
    acc = Hbond.num_of_accept(i)
    neighbors = Hbond.water_neighbors(i)
    conf = Conf_Analysis.conformation(i)

    result = {'don': don, 'acc': acc, 'neighbors': neighbors, 'conf': conf}

    return result


def main():
    start_time = time.time()
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        output = pool.map(process_frame, num_frames)

        df, df_conf = update_totals(output)

        df.to_csv("results/"+target+"_water.csv")
        df_conf.to_csv("results/"+target+"_confs.csv")
        # print(df)
        # print(f"Rigid: {ahr}, Flexible: {bbr}")
        # print(df_conf)
    end_time = time.time()
    print(end_time - start_time)



if __name__ == "__main__":
    main()







