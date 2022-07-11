# platform: win
# env: base
# error: command 'babel' not found in windows cli
import multiprocessing
from config import *


def batchproc(f):
    def mol2tosdf(path, file):
        try:
            sdf_file = os.path.join(path, file.split('\\')[-1].replace('.mol2', '.sdf'))
            cmdline = "babel -imol2 {} -osdf {} -h".format(file, sdf_file)
            os.system(cmdline)
        except:
            print('converting error for %s' % file)
    pass
    p = os.path.join(ign_example, "sdf_files")
    return mol2tosdf(p, f)


if __name__ == '__main__':
    mol2path = os.path.join(ign_example, "mol2_files")
    sdfpath = os.path.join(ign_example, "sdf_files")
    num_process = 12
    if not os.path.exists(sdfpath):
        os.makedirs(sdfpath)

    entries = os.listdir(mol2path)
    files = [os.path.join(mol2path, entry) for entry in entries]
    pool = multiprocessing.Pool(24)
    pool.starmap(batchproc, zip(files))
    pool.close()
    pool.join()
