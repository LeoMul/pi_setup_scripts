import numpy as np 
import sys
import os 
#notes
#i am assuming that only one symmetry per elev file for now.

default_email = '40268323@ads.qub.ac.uk'
default_partition = 'k2-math-physics'
default_time_limit = '1:00:00'
default_directory_list_location = ''
default_mem_per_cpu = 6
default_nodes = 1
default_num_cpu_per_node = 128

default_E0 = 0.01
default_NE = 128
default_DE = 0.1
default_guage = 'len'

class Input:
    def __init__(self,
                 email=default_email,
                 partition=default_partition,
                 time_limit = default_time_limit,
                 mem_per_cpu_GB = default_mem_per_cpu,
                 nodes=default_nodes,
                 cpu_per_node = default_num_cpu_per_node,
                 e0_ryd=default_E0,
                 num_e_points=default_NE,
                 delta_e_ryd = default_DE,
                 gauge = default_guage,
                 location_of_directory_list = default_directory_list_location):
 
        self.email = email 
        self.partition = partition
        self.time_limit = time_limit
        self.cpu_per_node = cpu_per_node
        self.mem_per_cpu_GB = mem_per_cpu_GB
        self.nodes= nodes
        self.e0_ryd = e0_ryd
        self.num_e_points = num_e_points
        self.delta_e_ryd = delta_e_ryd
        self.gauge = gauge
        self.location_of_directory_list = location_of_directory_list



def read_elev(elev_path):
    #reads Elev and gets the number of bound states.

    elev = open(elev_path,'r')
    elev_read = elev.readlines() 
    
    header = elev_read[0].replace('\n','')
    N = int(header[-1])
    Z = int(header[-2])

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

    total_cpu = nodes * cpu_per_node

    stgb_job = open('pstgbf.job','w')

    stgb_job.write('#!/bin/bash\n')

    stgb_job.write('#SBATCH --job-name={}\n'.format(title))

    stgb_job.write('#SBATCH --nodes={}\n'.format(nodes))
    stgb_job.write('#SBATCH --ntasks-per-node={}\n'.format(cpu_per_node))
    stgb_job.write('#SBATCH --mem-per-cpu={}G\n'.format(mem_gb) )
    stgb_job.write('#SBATCH --time={}\n'.format(time_limit))
    stgb_job.write('#SBATCH --partition={} \n'.format(partition))
    stgb_job.write('#SBATCH --error="error_stgb.out"\n')
    stgb_job.write('#SBATCH --mail-type=BEGIN,FAIL,END\n')
    stgb_job.write('#SBATCH --mail-user={}\n'.format(email))
    stgb_job.write('SCRATCH_DIRECTORY=/mnt/scratch2/users/${USER}/${SLURM_JOBID}.pex  ## scratch2 has 2PB size; no quota^M\n')

    stgb_job.write('module purge\n')
    stgb_job.write('module load mpi/openmpi/4.0.0/gcc-5.1.0\n')
    stgb_job.write('module load compilers/gcc/9.3.0\n')

    stgb_job.write('mpirun -np {} ./pstgbf0.x'.format(total_cpu))

    stgb_job.close()

def run_many_pstgf(directories,eff_charge,num_points,e0,de,partition,nodes,cpu_per_node,time_limit,mem_gb,email,guage):

    e0_new = e0 
    de_new = de

    print(directories)
    for directory in directories:
        print('changing directory to,',directory)
        os.chdir(directory)
        os.system('cp ../pstgbf0.x .')
        #os.system('cd {}'.format(directory))
        
        initial_symmetry,num_bound_states,Z,N = read_elev('ELEV')
        eff = Z - N 

        if eff != 0:
            e0_new = e0 / eff**2 
            de_new = de /eff ** 2 
        

        new_symmetries = determine_final_bound_states(initial_symmetry)
        write_dstgbf0damp(new_symmetries,num_bound_states,num_points,e0_new,de_new,guage)
        title = str(''.join(initial_symmetry))+'bf'
        print(title)
        write_pstgbf0_job_script(partition,nodes,cpu_per_node,time_limit,mem_gb,email,title)
        os.system('sbatch pstgbf.job')
        os.chdir('..')

def main(input:Input):

    num_points = input.num_e_points
    e0 = input.e0_ryd
    d0 = input.delta_e_ryd

    guage = input.gauge #options: 'len','vel'

    partition = input.partition
    time_limit = input.time_limit
    email = input.email
    nodes = input.nodes
    num_cpu_per_node = input.cpu_per_node
    mem_per_cpu = input.mem_per_cpu_GB

    if guage =='len':
        guage_key =0
    elif guage == 'vel':
        guage_key =1

    eff = 1

    f = open(input.location_of_directory_list,'r')
    f_read = f.readlines()
    f.close()
    directories_raw = []
    for line in f_read:
        linelist = list(line)
        if not('#' in linelist):
            directories_raw.append(line.replace('\n',''))

    #directories_raw = [(np.loadtxt('directories_for_pstgbf',dtype=str))]
    
    print('found directories: ',directories_raw)
    run_many_pstgf(directories_raw,eff,num_points,e0,d0,partition,nodes,num_cpu_per_node,time_limit,mem_per_cpu,email,guage_key)

    
    return 0 


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
