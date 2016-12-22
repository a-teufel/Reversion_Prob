#!/share/apps/python-2.7.2/bin/python
# Calculates the probability of accepting a reversion to the ancesteral state over the previous 15 and next 15 accepted mutations
# Either include pyrosetta in your path or run something like 
# source /home/ateufel/Pyrosetta/PyRosetta.monolith.ubuntu.release-80/SetPyRosettaEnvironment.sh
# so that pyrosetta is in the path
# note that your pdb files must start a residue 1, this may mean you have to renumber your pbd file 
# VERY IMPORTANT, THE NUMBER OF PDB FILES IN THE DIRECTORY HAS TO BE RIGHT! It can only contain PDB files with numerical names made by the forward simulation! If you have other pdb files in the dir with numbers in there names it will throw the whole thing off. 

# Imports
import sys
import argparse
import random, math, os, glob 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, change_cys_state, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover

from toolbox import mutate_residue
from toolbox import cleanATOM
from time import time
import re

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


def main():
    #read in the file made by the forward sim
    args = sys.argv
    inputfile = args[1]
    data = open(inputfile)
    first_line = data.readlines()[1]
    var_line=first_line.split(',')
    start_stab=var_line[1]

    #the first entry in the file is the wild type structure, calc the threshold using this
    threshold=float(start_stab)+10
    print(threshold)
    
    # Initialize Rosetta.
    init(extra_options='-mute basic -mute core')

    # Constants
    PACK_RADIUS = 0
    #Population size
    N = 100
    #Beta (temp term)
    beta = .6
  
    #Set up ScoreFunction
    sf = get_fa_scorefxn()

    #Set up MoveMap.
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)

    min_mover = MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(sf)
    min_mover.min_type('dfpmin_armijo_nonmonotone')

    #Prepare data headers
    data = ['pdbfile_target,pdbfile_used,step,RevertTo,Change,Pos,From,OrgScore,RevScore,Change,Prob\n']

    # Get the reversions file, the output file the score_mutant_pdb has made
    variant_scores=open(inputfile)

    #get just the mutation we want to revert to
    lines= variant_scores.readlines()
    var_line=lines[500] #gets the Nth line how ever long you want the burn to be
    print "staring here", var_line  
    var_line=var_line.split(',')[0]
  
    var_loc=int(filter(str.isdigit, var_line))
    var_rev=var_line[:1]

    gen=1
    #get all the pdb files
    sort_list=sorted(glob.glob('*[0-9].pdb'), key=numericalSort)
    sort_list=sort_list[-1016:] #include the last 1000 and some pdbs, the 16 is because we want the ones that happened before the 500th mutation too. 

  
    for i in range(1,len(sort_list)-30):
      step=-15
      #calc reversion for next 15 moves
      for infile in sort_list[i:i+31]:

	#for each mutation	
        var_line=lines[gen+500] #gets the Nth line how ever long you want the burn to be
        var_line=var_line.split(',')[0]
	print(var_line)
        var_loc=int(filter(str.isdigit, var_line))
	var_rev=""
	old=""
	if(step<0):
        	var_rev=var_line[len(var_line)-1:len(var_line)]
		old=var_line[:1]
	
	else:
		var_rev=var_line[:1]
		old=var_line[len(var_line)-1:len(var_line)]

      	print "Current File Being Processed is: " + infile
        print "revering to:", var_rev
        print "at:", var_loc

	#get the pdb you want to revert and make the reversion
        initial_pose = pose_from_pdb(infile)
        mutant_pose = mutate_residue(initial_pose, var_loc , var_rev, PACK_RADIUS, sf)

	#repack mut
        task1 = standard_packer_task(mutant_pose)
	task1.restrict_to_repacking()
        task1.or_include_current(True)
        packer_rotamers_mover1 = RotamerTrialsMover(sf,task1)
	packer_rotamers_mover1.apply(mutant_pose)

	#repack init
        task2 = standard_packer_task(initial_pose)
	task2.restrict_to_repacking()
	task2.or_include_current(True)
	pack_rotamers_mover2 = RotamerTrialsMover(sf, task2)
	pack_rotamers_mover2.apply(initial_pose)

	#apply min mover
	min_mover.apply(mutant_pose)
	min_mover.apply(initial_pose)
	
	#get scores    
	variant_score = sf(mutant_pose)
        initial_score = sf(initial_pose)

	#get prob
        probability = calc_prob_mh(variant_score, initial_score, N, beta, threshold)

	print(str(gen+499)+".pdb"+","+str(infile)+","+str(step)+","+ str(var_line) + ","+str(var_rev)+","+str(var_loc)+","+str(old)+"," +str(initial_score) + "," + str(variant
_score) + "," + str(variant_score - initial_score)+ ","+ str(probability)+ "\n")
      	data.append(str(gen+499)+".pdb"+","+str(infile)+","+str(step)+","+ str(var_line) + ","+str(var_rev)+","+str(var_loc)+","+str(old)+"," +str(initial_score) + "," + str(v
ariant_score) + "," + str(variant_score - initial_score)+ ","+ str(probability)+ "\n")
	step=step+1
      gen+=1

    print '\nDONE'

    data_filename = 'premutate_rep1_bb_T_ch_T.csv'
    with open(data_filename, "w") as f:
        f.writelines(data)

#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):
 
  xi = calc_x(stab_org, beta, thresholds)
  xj = calc_x(stab_mut, beta, thresholds)

  if xj >= xi:
    return(float(1.0))
  else:
    exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))


#score considering the thresh value
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


