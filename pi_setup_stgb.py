import os 
import json 
import argparse

#sets up and (optionally) runs your bound state calculations for photoionisation.

#todo:
#refactor the functions so they actually use the object orientation.
#rn the OOP is only for easy transmission of the input. most everything else is still entirely functionally 
#programmed which seems wasteful when there is OOP in here anyway.  

#refactor on 15.04.24
#stgb only needs the Hamiltonian. For this reason I'll move the sym linking of the D files to the
#pstgb set up code, and there change it so the user gives the path to the B files which are then 
#sym linked.

#additionally - change a separate symmetry deck to a list of symmetries in the input deck.
#makes the input deck entirely self contained.

#also - you need to send the code to Connor after you do this 
#- you from the past


#this is hard coded right now. i might change that.
stgb_exec_path = '../stgb_blas.proto.x'

#what does the input deck do, Leo? 

#overwrite_existing_arrays - boolean - whether or not this code overwrites directories it finds with the same name as it calcualtes.
#symmetries - list of strings in format "2j parity".
#hamiltonian_path -str- - real (or relative) path of the Hamiltonian.dat
#email - str -email for slurm - leaving it blank should in theory stop it. (i havent tested that)
#partition - str -  slurm partition
#time_limit - str -slurm job time limit
#mem_per_cpu_GB - int -memory per cpu in GB
#run_jobs - boolean - whether or not the code sbatch's the jobs made. 
#for testing really, afterwards one might type sbatch symdir*/stgb.job  or the specific ones they want.

#n0,nf,dn refer to the quantum defect paramters for stgb. see the code.


#default stuff, so that the user easily knows how to change the input. in theory.

#the user might consider changing the defaults to suit their needs.
default_email = '40268323@ads.qub.ac.uk'
default_partition = 'k2-math-physics-debug'
default_time_limit = '1:00:00'
default_mem_per_cpu = 6

default_symmetries = ['1 0','1 1']
default_hamiltonian_path = '../H.DAT'

default_n0 = 0.21
default_n = 4.2
default_dn = 0.2

default_rewrite_override = True

default_run_jobs = True

class Input:
    #this object oriented business allows the input deck to be read in in the json format
    #which makes my life a bit easier.
    def __init__(self,
                 overwrite_existing_arrays=default_rewrite_override,
                 symmetries=default_symmetries,
                 hamiltonian_path=default_hamiltonian_path,
                 email=default_email,
                 partition=default_partition,
                 time_limit = default_time_limit,
                 mem_per_cpu_GB = default_mem_per_cpu,
                 run_jobs = default_run_jobs,
                 n0=default_n0,
                 nf=default_n,
                 dn=default_dn):
        self.overwrite_existing_arrays = overwrite_existing_arrays
        self.symmetries = symmetries 
        self.hamiltonian_path = hamiltonian_path
        self.email = email 
        self.partition = partition
        self.time_limit = time_limit
        self.mem_per_cpu_GB = mem_per_cpu_GB
        self.run_jobs = run_jobs
        self.n0 = n0
        self.nf = nf 
        self.dn = dn


def read_in_symmetries(input:Input):
    #reads in the symmetries in in the input deck, returns aa list[list[str]]
    symmetries_strings_list = input.symmetries
    symmetries = []
    for sym in symmetries_strings_list:

        split = sym.split()
        twoj = (split[0])
        pari = (split[1])

        symmetries.append([twoj,pari])
    return symmetries

def decode_symmetries(symmetries):
    #converts the symmetries file into the desired directories.
    desired_symmetries = []
    
    for row in symmetries:
        string = 'symdir_2j'+str(row[0]) + '_p' + str(row[1])
        desired_symmetries.append(string)

    return desired_symmetries

def make_symmetry_directories(desired_directories,symmetries):

    #makes the directories - ignores ones that already exists - UNLESS override is on.
    #ouputs the directories that didnt already exist - these are the ones that we need to consider.
    #directories = desired_directories.copy()
    directories_kept = []
    symmetries_kept = []
    for (index,directory) in enumerate(desired_directories):
        if os.path.isdir(directory):
            print(directory,' already exists WARNING')
        else:
            os.mkdir(directory)
            symmetries_kept.append(symmetries[index])
            directories_kept.append(desired_directories[index])

            #desired_directories.remove(directory)

    return directories_kept,symmetries_kept

def write_dstgb(symmetry,n0,nf,delta_n):
    
    #writes dstgb, based on the parameters you give it. 
    #you, the user, might consider changing some of it.


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

    current_directory = os.getcwd()

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
    stgb_job.write('SCRATCH_DIRECTORY=/mnt/scratch2/users/${USER}/${SLURM_JOBID}.pex  \n')
    stgb_job.write('export TMPDIR={}\n'.format(current_directory))
    stgb_job.write('cd {}\n'.format(current_directory))

    stgb_job.write('module purge\n')
    stgb_job.write('module load mpi/openmpi/4.0.0/gcc-5.1.0\n')
    stgb_job.write('module load compilers/gcc/9.3.0\n')

    stgb_job.write('./stgb_blas.proto.x < dstgb')

    stgb_job.close()


def setup_many_dstgb_directories(directories,symmetries,n0,nf,delta_n,partition,time_limit,mem_gb,email,title_end,hamiltonian_path,run_jobs):

    assert(len(directories) == len(symmetries))


    for ii in range(0,len(directories)):
        os.chdir(directories[ii])
        write_dstgb(symmetries[ii],n0,nf,delta_n)
        title = ''.join(symmetries[ii])+title_end

        print('writing job script with title ',title)
        #write out a stgb job script
        write_stgb_job_script(partition,time_limit,mem_gb,email,title)

        os.system('cp {}  .'.format(stgb_exec_path))
        os.system('ln -s {} H.DAT'.format(hamiltonian_path))
        print('linking hamiltonian from directory:',hamiltonian_path)

        if run_jobs:
            os.system('sbatch stgb.job')

        os.chdir('..')


def main(input:Input):

    n_0 = input.n0
    n_f = input.nf
    delta_n = input.dn

    email = input.email
    partition = input.partition
    time_limit = input.time_limit
    title = 'stgb'
    mem_gb = input.mem_per_cpu_GB

    directory_of_data_files = input.hamiltonian_path

    symmetries = read_in_symmetries(input)
    desired_directories = decode_symmetries(symmetries)

    print('the desired directories are:',desired_directories)

    over_ride = input.overwrite_existing_arrays
    directories_to_be_worked_with,symmetries_kept = make_symmetry_directories(desired_directories,symmetries)

    #if the user chooses to override, the directories will be entered anyway and the various stgb's written.
    if over_ride:
        print('override is on - so i will overwrite the directories that already exist')
        directories_to_be_worked_with = desired_directories
        symmetries_kept = symmetries

    print('working in directories: ',directories_to_be_worked_with)
    print('working with symmetries: ',symmetries_kept)
    run_jobs = input.run_jobs
    setup_many_dstgb_directories(directories_to_be_worked_with,symmetries_kept,n_0,n_f,delta_n,partition,time_limit,mem_gb,email,title,directory_of_data_files,run_jobs)



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

