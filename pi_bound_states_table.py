import os 

from glob import glob

#scans directories all child directories in the current directory and makes a table of ionisation potentials
#from each found level of the N+1 ion.

directories = glob(str(os.curdir)+'/*/')
#print(directories)

def main(directories):

    initial_symmeties = []
    ionisations = []
    defects = []
    for directory in directories:
        elev_path = directory+'/ELEV' 
        if os.path.exists(elev_path):
            initial_symmetry,ionisation_energies,n_princ = read_elev_for_table(elev_path)
            num = len(ionisation_energies)

            for ii in range(0,num):
                initial_symmeties.append(initial_symmetry)
            ionisations.extend(ionisation_energies)
            defects.extend(n_princ)
        else:
            print('no elev found in ',directory)

    ground = min(ionisations)
    ion_rel_to_ground = []
    for ion in ionisations:
        ion_rel_to_ground.append(ion-ground)


    sorted_defects = [x for _, x in sorted(zip( ion_rel_to_ground,defects))]
    initial_symmeties = [x for _, x in sorted(zip( ion_rel_to_ground,initial_symmeties))]
    #print(initial_symmeties)
    ionisations.sort()
    ion_rel_to_ground.sort()
    print('--------------------------------------------------')
    header = 'Index   IP (Ry)  Level (Ry)      QDQN   2J  parity'
    
    print(header)
    print('--------------------------------------------------')

    output = '{:5}  {:8.5f}    {:8.5f}  {:8.5f}  {:3}  {:4}'
    for jj in range(0,len(ionisations)):
        twoj = int(initial_symmeties[jj][1])
        p = int(initial_symmeties[jj][2])
        if p == 0:
            p = 'even'
        elif p==1:
            p = 'odd'
        else:
            p = 'error'
        print(output.format(jj+1,ionisations[jj],ion_rel_to_ground[jj],sorted_defects[jj],twoj,p))

    print('--------------------------------------------------')



def read_elev_for_table(elev_path):
    #reads Elev and gets the number of bound states.

    elev = open(elev_path,'r')
    elev_read = elev.readlines() 
    
    header = elev_read[0].replace('\n','')
    N = int(header.split()[-1])
    Z = int(header.split()[-2])

    scale = (Z-N)**2

    initial_symmetry = (elev_read[4]).split()[0:3]    
    num_bound_states = int(elev_read[5].split()[0])

    n_princ = []
    ionisation_energies = []

    position = 6
    for jj in range(0,num_bound_states):
        current_line_split = elev_read[jj + position].split()

        n_princ_quantum = current_line_split[-1]
        ionisation = current_line_split[-2]
        n_princ.append(float(n_princ_quantum))
        ionisation_energies.append(float(ionisation)*scale)

    elev.close()
    return initial_symmetry,ionisation_energies,n_princ

main(directories)