import os 
import numpy as np 

default_email = '40268323@ads.qub.ac.uk'
default_partition = 'k2-math-physics'
default_time_limit = '1:00:00'
default_mem_per_cpu = 6

default_symmetry_path = ''
default_directory_of_data = '..'

default_n0 = 0.21
default_n = 10.2
default_dn = 0.2

default_rewrite_override = True

class Input:
    def __init__(self,
                 override=default_rewrite_override,
                 symmetry_path=default_symmetry_path,
                 data_directory=default_directory_of_data,
                 email=default_email,
                 partition=default_partition,
                 time_limit = default_time_limit,
                 mem_per_cpu = default_mem_per_cpu,
                 n0=default_n0,
                 nf=default_n,
                 dn=default_dn):
        self.override = override
        self.symmetry_path = symmetry_path 
        self.data_directory = data_directory
        self.email = email 
        self.partition = partition
        self.time_limit = time_limit
        self.mem_per_cpu = mem_per_cpu
        self.n0 = default_n0
        self.nf = nf 
        self.dn = dn


def read_in_symmetries(path_to_symmetries):
    #reads in the symmetries file
    symmetries = np.loadtxt(path_to_symmetries,dtype=str)
    return symmetries

def decode_symmetries(symmetries):
    #converts the symmetries file into the desired directories.
    desired_symmetries = []
    
    for row in symmetries:
        string = 'symdir_2j'+row[0] + '_p' + row[1]
        desired_symmetries.append(string)

    return desired_symmetries

def make_symmetry_directories(desired_directories,symmetries):

    #makes the directories - ignores ones that already exists.
    #ouputs the directories that didnt already exist - these are the ones that we need to consider.
    #directories = desired_directories.copy()
    directories_kept = []
    symmetries_kept = []
    for (index,directory) in enumerate(desired_directories):
        if os.path.isdir(directory):
            print(directory,' already exists - skipping')
        else:
            os.mkdir(directory)
            symmetries_kept.append(symmetries[index])
            directories_kept.append(desired_directories[index])

            #desired_directories.remove(directory)

    return directories_kept,symmetries_kept

def write_dstgb(symmetry,n0,nf,delta_n):
    #writes dstgb.


    symmetry_string = ' {} {} {}\n'.format(0,symmetry[0],symmetry[1])
    n_string = ' {} {} {}\n'.format(n0,nf,delta_n)

    dstgb = open('dstgb','w')
    dstgb.write(' &STGB IOPT2=1 IRAD=1 IPERT=1 IPRINT=0 &END\n')
    dstgb.write(symmetry_string)
    dstgb.write(n_string)
    dstgb.write(' -1 -1 -1')
    dstgb.close()

    return 0

def write_stgb_job_script(partition,time_limit,mem_gb,email,title):

    stgb_job = open('stgb.job','w')

    stgb_job.write('#!/bin/bash\n')

    stgb_job.write('#SBATCH --job-name={}\n'.format(title))

    stgb_job.write('#SBATCH --nodes=1\n')
    stgb_job.write('#SBATCH --ntasks-per-node=1\n')
    stgb_job.write('#SBATCH --mem={}G\n'.format(mem_gb) )
    stgb_job.write('#SBATCH --time={}\n'.format(time_limit))
    stgb_job.write('#SBATCH --partition={} \n'.format(partition))
    stgb_job.write('#SBATCH --error="error_stgb.out"\n')
    stgb_job.write('#SBATCH --mail-type=BEGIN,FAIL,END\n')
    stgb_job.write('#SBATCH --mail-user={}\n'.format(email))
    stgb_job.write('SCRATCH_DIRECTORY=/mnt/scratch2/users/${USER}/${SLURM_JOBID}.pex  ## scratch2 has 2PB size; no quota^M\n')

    stgb_job.write('module purge\n')
    stgb_job.write('module load mpi/openmpi/4.0.0/gcc-5.1.0\n')
    stgb_job.write('module load compilers/gcc/9.3.0\n')

    stgb_job.write('./stgb_blas.proto.x < dstgb')

    stgb_job.close()


def setup_many_dstgb_directories(directories,symmetries,n0,nf,delta_n,partition,time_limit,mem_gb,email,title,data_directory):

    assert(len(directories) == len(symmetries))

    print('linking data from directory:',data_directory)

    for ii in range(0,len(directories)):
        os.chdir(directories[ii])
        write_dstgb(symmetries[ii],n0,nf,delta_n)
        #print(symmetries[ii])
        title = ''.join(symmetries[ii])+'stgb'
        print(title)
        write_stgb_job_script(partition,time_limit,mem_gb,email,title)
        os.system('cp ../stgb_blas.proto.x .')
        os.system('sbatch stgb.job')
        os.system('ln -s ../{}/D* .'.format(data_directory))
        os.system('ln -s ../{}/H.DAT H.DAT'.format(data_directory))

        os.chdir('..')


def main(input:Input):

    n_0 = input.n0
    n_f = input.nf
    delta_n = input.dn

    email = input.email
    partition = input.partition
    time_limit = input.time_limit
    title = 'stgb'
    mem_gb = input.mem_per_cpu

    directory_of_data_files = input.data_directory

    symmetries = read_in_symmetries(input.symmetry_path)
    desired_directories = decode_symmetries(symmetries)

    print('the desired directories are:',desired_directories)

    over_ride = input.override
    directories_to_be_worked_with,symmetries_kept = make_symmetry_directories(desired_directories,symmetries)

    #if the user chooses to override, the directories will be entered anyway and the various stgb's written.
    if over_ride:
        directories_to_be_worked_with = desired_directories
        symmetries_kept = symmetries

    print('working in directories: ',directories_to_be_worked_with)
    print('working with symmetries: ',symmetries_kept)

    setup_many_dstgb_directories(directories_to_be_worked_with,symmetries_kept,n_0,n_f,delta_n,partition,time_limit,mem_gb,email,title,directory_of_data_files)

import json 
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-j', '--json',  help='path of json')
args = parser.parse_args()

if not args.json:
    input_default = Input()
    default = json.dumps(input_default.__dict__,indent=1)
    print(default) 

else: 
    with open(args.json, 'r') as j:
        contents = json.loads(j.read())

    input = Input(**contents)
    main(input)

