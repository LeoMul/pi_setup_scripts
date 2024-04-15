import sys
import os 
import json 
import argparse

#sets up and (optionally) runs your pstgbf0 jobs.

#it 

#i am assuming that only one symmetry per elev file for now.

#todo - make better use of the oop. edge cases 

#example - the code right now produces all the possible allowed final states from the initial bound state.
#I should ask Connor if this is sound - or if the possible allowed final states should be limited to those in 
#the DSTG2.INP or the sizeBP files.



default_pstgbf0_exec_path = 'pstgbf0.x'


#default stuff so the user knows how to easily change the input.
default_email = '40268323@ads.qub.ac.uk'
default_partition = 'k2-math-physics'
default_time_limit = '1:00:00'

default_run_jobs = True

default_hamiltonian_path = 'H.DAT'
default_D_files_path = '../D*'
default_bound_state_directory_abs_path = ''

default_directory_list = ['symdir_2j1_p0','symdir_2j0_p1']

default_mem_per_cpu = 6
default_nodes = 1
default_num_cpu_per_node = 128

#these are in Rydbergs. NOT scaled Rydbergs. THe code does that for you.
default_E0 = 0.01
default_NE = 128
default_DE = 0.1

#this assumes you correctly set the guage in stgd. 
#I'm looking at you, past Leo, you rapscallion.
default_guage = 'len'

class Input:
    def __init__(self,
                 email=default_email,
                 partition=default_partition,
                 time_limit = default_time_limit,
                 mem_per_cpu_GB = default_mem_per_cpu,
                 nodes=default_nodes,
                 cpu_per_node = default_num_cpu_per_node,
                 run_jobs = default_run_jobs,
                 e0_ryd=default_E0,
                 num_e_points=default_NE,
                 delta_e_ryd = default_DE,
                 gauge = default_guage,
                 directory_list = default_directory_list,
                 hamiltonian_absolute_path = default_hamiltonian_path,
                 D_files_absolute_path = default_D_files_path,
                 bound_state_directory_abs_path = default_bound_state_directory_abs_path,
                 pstgbfx_absolute_path = default_pstgbf0_exec_path):
 
        self.email = email 
        self.partition = partition
        self.time_limit = time_limit
        self.cpu_per_node = cpu_per_node
        self.mem_per_cpu_GB = mem_per_cpu_GB
        self.nodes= nodes
        self.run_jobs = run_jobs
        self.e0_ryd = e0_ryd
        self.num_e_points = num_e_points
        self.delta_e_ryd = delta_e_ryd
        self.gauge = gauge
        self.pstgbfx_absolute_path = pstgbfx_absolute_path
        self.hamiltonian_absolute_path = hamiltonian_absolute_path
        self.D_files_absolute_path = D_files_absolute_path
        self.bound_state_directory_abs_path = bound_state_directory_abs_path
        self.directory_list = directory_list



def read_elev(elev_path):
    #reads Elev and gets the number of bound states.

    elev = open(elev_path,'r')
    elev_read = elev.readlines() 
    
    header = elev_read[0].replace('\n','')
    N = int(header.split()[-1])
    Z = int(header.split()[-2])

    initial_symmetry = (elev_read[4]).split()[0:3]    
    num_bound_states = int(elev_read[5].split()[0])
    elev.close()
    return initial_symmetry,num_bound_states,Z,N

def determine_final_bound_states(initial_bound_state_string_split):
    #initial_bound_state = initial_bound_state_string.split()

    assert(len(initial_bound_state_string_split)) == 3 ,'bound state invalid'

    final_bound_states = []
    new_2j = []
    old_parity = int(initial_bound_state_string_split[-1])
    old_2j = int(initial_bound_state_string_split[1])
    if old_parity == 1: 
        new_parity = 0
    elif old_parity == 0:
        new_parity = 1
    else:
        print('invalid parity in one of the initial boundstates:',initial_bound_state_string_split)
        sys.exit()
    
    if old_2j == 0: 
        new_2j.append(2)
    elif old_2j == 1:
        new_2j.append(1)
        new_2j.append(3)
    elif old_2j > 1:
        new_2j.append(old_2j - 2)
        new_2j.append(old_2j )
        new_2j.append(old_2j + 2)

    for twoj in new_2j:
        state = '0 {} {}'.format(twoj,new_parity)
        final_bound_states.append(state)



    return final_bound_states


def write_dstgbf0damp(final_bound_states,num_bound_states,num_points,e0,de,guage):

    input_deck = open('dstgbf0damp','w')
    input_deck.write('&STGF\n')
    input_deck.write('       NODAMP=1 IPERT=1 IPRINT=-2 IQDT=0 IMESH=1 IOPT1=2\n') 
    input_deck.write('       IAUGER=0 NTYP1=0 NTYP2I=0 NMIN=4 IGAUGE={} \n'.format(guage))
    input_deck.write('       IPHOTO=400 NPIEB={} /\n'.format(num_bound_states))
    input_deck.write('&MESH1 MXE={} E0={} EINCR={} &END\n'.format(int(num_points),e0,de))

    for bound_state in final_bound_states:
        input_deck.write('{}\n'.format(bound_state))

    input_deck.write('-1 -1 -1')

    input_deck.close()
    return 0

def write_pstgbf0_job_script(partition,nodes,cpu_per_node,time_limit,mem_gb,email,title):

    current_directory = os.getcwd()

    total_cpu = nodes * cpu_per_node

    pstgbf0_job = open('pstgbf.job','w')
    pstgbf0_job.write('#!/bin/bash\n')
    pstgbf0_job.write('#SBATCH --job-name={}\n'.format(title))
    pstgbf0_job.write('#SBATCH --nodes={}\n'.format(nodes))
    pstgbf0_job.write('#SBATCH --ntasks-per-node={}\n'.format(cpu_per_node))
    pstgbf0_job.write('#SBATCH --mem-per-cpu={}G\n'.format(mem_gb) )
    pstgbf0_job.write('#SBATCH --time={}\n'.format(time_limit))
    pstgbf0_job.write('#SBATCH --partition={} \n'.format(partition))
    pstgbf0_job.write('#SBATCH --error="error_stgb.out"\n')
    pstgbf0_job.write('#SBATCH --mail-type=BEGIN,FAIL,END\n')
    pstgbf0_job.write('#SBATCH --mail-user={}\n'.format(email))
    pstgbf0_job.write('SCRATCH_DIRECTORY=/mnt/scratch2/users/${USER}/${SLURM_JOBID}.pex  ## scratch2 has 2PB size; no quota^M\n')
    pstgbf0_job.write('export TMPDIR={}\n'.format(current_directory))
    pstgbf0_job.write('cd {}\n'.format(current_directory))
    pstgbf0_job.write('module purge\n')
    pstgbf0_job.write('module load mpi/openmpi/4.0.0/gcc-5.1.0\n')
    pstgbf0_job.write('module load compilers/gcc/9.3.0\n')
    pstgbf0_job.write('mpirun -np {} ./pstgbf0.x'.format(total_cpu))
    pstgbf0_job.close()

def run_many_pstgf(input:Input):
    directories = input.directory_list
    num_points = input.num_e_points
    e0=input.e0_ryd
    de=input.delta_e_ryd
    partition=input.partition
    nodes=input.nodes
    cpu_per_node=input.cpu_per_node
    time_limit=input.time_limit
    mem_gb=input.mem_per_cpu_GB
    email = input.email
    guage = input.gauge

    bound_path = input.bound_state_directory_abs_path

    d_path = input.D_files_absolute_path
    ham_path = input.hamiltonian_absolute_path
    run_jobs = input.run_jobs
    e0_new = e0 
    de_new = de

    if guage =='len':
        guage =0
    elif guage == 'vel':
        guage =1

    print(directories)

    for directory in directories:

        if not os.path.exists(directory):
            os.mkdir(directory)


        print('changing directory to,',directory)
        os.chdir(directory)
        os.system('cp {} pstgbf0.x'.format(input.pstgbfx_absolute_path))
        #os.system('cd {}'.format(directory))
        bound_symmetry_path = bound_path+'/'+directory

        elev_path = bound_symmetry_path+'/ELEV'

        if os.path.exists(elev_path):
            

            initial_symmetry,num_bound_states,Z,N = read_elev(elev_path)

            if num_bound_states != 0:

                eff = Z - N 
                print('found nuclear charge {} and num electrons {}'.format(Z,N))
                if eff != 0:
                    e0_new = e0 / eff**2 
                    de_new = de /eff ** 2 


                new_symmetries = determine_final_bound_states(initial_symmetry)
                write_dstgbf0damp(new_symmetries,num_bound_states,num_points,e0_new,de_new,guage)
                title = str(Z)+str(N)+str(''.join(initial_symmetry))+'bf'
                print('making job script with title',title)
                write_pstgbf0_job_script(partition,nodes,cpu_per_node,time_limit,mem_gb,email,title)
                os.system('ln -s {} H.DAT'.format(ham_path))
                os.system('ln -s {}/B* .'.format(bound_symmetry_path))

                os.system('ln -s {} .'.format(d_path))
                if run_jobs:
                    os.system('sbatch pstgbf.job')
            else:
                print(' i found no bound states, skipping sym',directory)
        else:
            print('i did not find an ELEV in',directory,' so i am skipping it.')
        print('*************')
        os.chdir('..')

def main(input:Input):

    run_many_pstgf(input)

    
    return 0 


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
