#!/share/apps/python-2.7.2/bin/python
#
# Either include pyrosetta in your path or run something like 
# source /home/ateufel/Pyrosetta/PyRosetta.monolith.ubuntu.release-80/SetPyRosettaEnvironment.sh
# so that pyrosetta is in the path
# run this program like: python forwardsim.py tiny_gene
# where tiny_gene is the name of the pdb file you want to use without the ".pdb" part
# note that your pdb files must start a residue 1, this may mean you have to renumber your pbd file 

# Imports
import sys
import argparse
import random, math, os 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, change_cys_state, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover

from toolbox import mutate_residue
from toolbox import cleanATOM
from time import time


def main():
    #takes name of pdb file without the extention
    args =  sys.argv	
    pdb_file = args[1]
    out_file = args[2]
    #set up timer to figure out how long the code took to run
    t0=time()

    # Initialize Rosetta.
    init(extra_options='-mute basic -mute core')

    # Constants
    PACK_RADIUS = 0.0
    #Amino acids, notice there is no C
    AAs = ("A","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    #Number of mutations to accept
    max_accept_mut = 1500
    #Population size
    N = 100
    #Beta (temp term)
    beta = .6

    #Prepare data headers
    data = ['Variant,Rosetta Score,"delta-delta-G",Probability,Generation\n']

    #Load and clean up pdb file  
    name=pdb_file+".pdb"
    cleanATOM(name)
    clean_name=pdb_file+".clean.pdb"
    initial_pose = pose_from_pdb(clean_name)

    #Set up ScoreFunction
    sf = get_fa_scorefxn()
       
    #Set up MoveMap This is where you turn the bb and side chain flexibility on and off
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)

    #Get the init score of the struct to calc the threshold
    pre_pre_packing_score = sf(initial_pose)

    min_mover = MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(sf)
    min_mover.min_type('dfpmin_armijo_nonmonotone')

    #Set threshold for selection 
    threshold = pre_pre_packing_score+10

    data.append('WT,' + str(pre_pre_packing_score) + ',0.0 ,0.0,0\n')

    #number of residues to select from
    n_res = initial_pose.total_residue()

    #start sim
    i=0
    gen=0
    while i < max_accept_mut:
            #update the number of generations that have pased
            gen+=1

	    print 'accepts:', i 

	    #pick a place to mutate
	    mut_location = random.randint(1, n_res)

	    #get the amino acid at that position
	    res = initial_pose.residue(mut_location)

	    #don't mess with C, just choose again
	    while(res.name1() == 'C'):
		mut_location = random.randint(1, n_res)
	    	#get the amino acid at that position
	    	res = initial_pose.residue(mut_location)


	    #choose the amino acid to mutate to
            toname = res.name1()
	    new_mut_key = random.randint(0,len(AAs)-1)
	    proposed_res = AAs[new_mut_key]
	  
	    #don't bother mutating to the same amino acid it just takes more time
	    while(proposed_res == res.name1()):
		new_mut_key = random.randint(0,len(AAs)-1)
	        proposed_res = AAs[new_mut_key]

	    #make the mutation
	    mutant_pose = mutate_residue(initial_pose, mut_location, proposed_res, PACK_RADIUS, sf)

	    #repack the mutant
	    task1 = standard_packer_task(mutant_pose)
	    task1.restrict_to_repacking()
            task1.or_include_current(True)
            packer_rotamers_mover1 = RotamerTrialsMover(sf,task1)
	    packer_rotamers_mover1.apply(mutant_pose)

	    #repack the init
            task2 = standard_packer_task(initial_pose)
	    task2.restrict_to_repacking()
	    task2.or_include_current(True)
	    pack_rotamers_mover2 = RotamerTrialsMover(sf, task2)
	    pack_rotamers_mover2.apply(initial_pose)

	    #apply minimization to both
	    min_mover.apply(mutant_pose)
	    min_mover.apply(initial_pose)
 
	    #score both structures
	    variant_score = sf(mutant_pose)
            post_pre_packing_score = sf(initial_pose)

	    #get the probability that the mutation will be accepted
	    probability = calc_prob_mh(variant_score, post_pre_packing_score, N, beta, threshold)
		
	    print "prob", probability
	    rand = random.random()
	    #print 'rand', rand

	    #test to see if mutation is accepted
	    if float(rand) < float(probability):
		print "accepted" 	
		
		#make a name for the new mutant
		variant_name = str(toname) + str(initial_pose.pdb_info().number(mut_location)) + str(proposed_res)

		# Assuming some burn in phase, make this zero if you want to store everything
		if i>=0:
		   #save name and energy change
	    	   data.append(variant_name + "," + str(variant_score) + "," + str(variant_score - post_pre_packing_score) + "," + str(probability) + "," + str(gen) + "\n")

		   pdb_name=str(i)+".pdb"	
		   mutant_pose.dump_pdb(pdb_name)

		#update the wildtype 
		initial_pose = mutant_pose

		#update number of accepts
	    	i+=1

    #end of sim
    print '\nMutations and scoring complete.'
    t1 = time()

    # Output results
    data_filename = pdb_file[:-5] + args[2]  + '.csv'
    with open(data_filename, "w") as f:
        f.writelines(data)

    print 'Data written to:', data_filename
    print 'program takes %f' %(t1-t0)


###assorted functions that have to do with scoring and prob of acceptance ####

#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):
 
  xi = calc_x(stab_org, beta, thresholds)
  xj = calc_x(stab_mut, beta, thresholds)

  if xj >= xi:
    return(float(1.0))
  else:
    exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))


def calc_x(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total += -math.log(safe_calc(exponent) + 1)
  return(total)


def safe_calc(exponent):
  if exponent > 700:
    print("system maxed")
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))


#Run main program
if __name__ == '__main__':
   main()


