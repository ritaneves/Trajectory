#Beam search, pruning of the worst nodes on each level
#Gets a Lambert Leg budget and a Beam Size for the best nodes
#V1: Prunes nodes without considering its characteristics

#!/usr/bin/env python

from state_asteroids import State, MOVE_TYPE
import sys


def BeamSearch(root, beam_size):
    
    best = None
    states = [root.copy()]
    legs = 0
    flag = 0
    legs_converged = 0

    while not states == []:
        frontier = []
        for s in sorted(states, key=lambda x: x.sum):
            for move in s.moves(): 
                ns = s.copy()
		if ns.next_move == MOVE_TYPE.TOF:
		    flag = 1
		    legs += 1 #We only have a new computation when we deal with TOF
		    
                ns.move(move)

                if not ns.isterminal() or ns.tof == []:
                    frontier.append(ns)
                else:
                    if best is None or len(best.seq) < len(ns.seq):
                        best = ns.copy()

	#Sort the moves in TOF by dv
	#The frontier will only get beam_size of these
	if flag == 1:
	    frontier = sorted(frontier, key=lambda x: x.dv)[0:beam_size]

	#The new level consists on the pruned children of the last one
        states = frontier
   
    return legs, best


if __name__ == '__main__':
    
    legs, best = BeamSearch(root = State(), beam_size = 1000000)
    print legs, best
    
