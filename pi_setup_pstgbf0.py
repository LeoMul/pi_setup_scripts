import numpy as np 
import sys
import os 
#notes
#i am assuming that only one symmetry per elev file for now.


def read_elev(elev_path):
    #reads Elev and gets the number of bound states.

    elev = open(elev_path,'r')
    elev_read = elev.readlines() 
    initial_symmetry = (elev_read[4]).split()[0:3]    
    num_bound_states = int(elev_read[5].split()[0])
    elev.close()
    return initial_symmetry,num_bound_states

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

    e0 = e0 / eff_charge**2 
    de = de /eff_charge ** 2
    print(directories)
    for directory in directories:
        print('changing directory to,',directory)
        os.chdir(directory)
        os.system('cp ../pstgbf0.x .')
        #os.system('cd {}'.format(directory))
        initial_symmetry,num_bound_states = read_elev('ELEV')
        new_symmetries = determine_final_bound_states(initial_symmetry)
        write_dstgbf0damp(new_symmetries,num_bound_states,num_points,e0,de,guage)
        title = str(''.join(initial_symmetry))+'bf'
        print(title)
        write_pstgbf0_job_script(partition,nodes,cpu_per_node,time_limit,mem_gb,email,title)
        os.system('sbatch pstgbf.job')
        os.chdir('..')

def main():

    Z = 52 
    N = 50 #in target...
    num_points = 1
    e0 = 1.0 
    d0 = 1.0

    guage = 'len' #options: 'len','vel'

    partition = 'k2-math-physics-debug'
    time_limit = '0:10:00'
    email = '40268323@ads.qub.ac.uk'
    nodes = 1
    num_cpu_per_node = 128
    mem_per_cpu = 6
    eff = Z - N 

    if guage =='len':
        guage_key =0
    elif guage == 'vel':
        guage_key =1
    

    f = open('directories_for_pstgbf','r')
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

main()