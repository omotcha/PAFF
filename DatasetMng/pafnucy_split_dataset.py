from sklearn.utils import shuffle
from config import *
import h5py

import os

import argparse

dataset_dir = base_dir

parser = argparse.ArgumentParser(description='Split dataset into training,'
                                 ' validation and test sets.')
parser.add_argument('--input_path', '-i', default=dataset_dir,
                    help='directory with pdbbind dataset')
parser.add_argument('--output_path', '-o', default=os.path.join(dataset_dir, 'hdfs'),
                    help='directory to store output files')
parser.add_argument('--size_val', '-s', type=int, default=1000,
                    help='number of samples in the validation set')
args = parser.parse_args()

# create files with the training and validation sets
with h5py.File('%s\\training_set.hdf' % args.output_path, 'w') as g, \
     h5py.File('%s\\validation_set.hdf' % args.output_path, 'w') as h:
    with h5py.File('%s\\refined.hdf' % args.input_path, 'r') as f:
        refined_shuffled = shuffle(list(f.keys()), random_state=123)
        for pdb_id in refined_shuffled[:args.size_val]:
            ds = h.create_dataset(pdb_id, data=f[pdb_id])
            ds.attrs['affinity'] = f[pdb_id].attrs['affinity']
        for pdb_id in refined_shuffled[args.size_val:]:
            ds = g.create_dataset(pdb_id, data=f[pdb_id])
            ds.attrs['affinity'] = f[pdb_id].attrs['affinity']
    with h5py.File('%s\\general.hdf' % args.input_path, 'r') as f:
        for pdb_id in f:
            ds = g.create_dataset(pdb_id, data=f[pdb_id])
            ds.attrs['affinity'] = f[pdb_id].attrs['affinity']

# # create a symlink for the test set
# os.symlink(os.path.abspath('%s\\core.hdf' % args.input_path),
#            os.path.abspath('%s\\test_set.hdf' % args.output_path))
# just manually copy /core.hdf to /hdfs/test_set.hdf
