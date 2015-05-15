#HOOT without greedy step, saving the entire tree
#Performs HOO when the action is continuous (T0 and TOF)
#V3: Returning ALL non-zero stats nodes or N (=5) if the former are not enough
#Expands all the planets and selects one of them at random.
#Removes nodes that reach terminal state, but not the entire path until there - may be trapped.

#TODO: Check v and rho, the number of cycles
#!/usr/bin/env python

from __future__ import division
from math import *
import sys
import time
import random
import numpy as np
import copy
import itertools

from state_rosetta_continuous import State, MOVE_TYPE, MAX_DV, MAX_EPOCH, fix_first_move, set_t_res
import state_rosetta_continuous
from tools import pretty_time


class Node:

    def __init__(self, down_lim=None, up_lim=None, state=None, last_move = None, parent = None):
        self.parent = parent
        self.children = []
        self.untried_moves = list(state.moves())
        self.state = state.copy()
        self.last_move = last_move
	self.stats = 0 # value estimate
        self.visits = 0 # visits
        self.b_value = sys.maxint
        self.down_lim = down_lim
        self.up_lim = up_lim

    #HOO Functions
    def untried(self, state):
	if state is None:
	    self.untried_moves = []
	else:
  	    self.untried_moves = list(state.moves())
	
    def Select_B(self):
	if self.children[0].b_value == self.children[1].b_value:
	    return random.choice(self.children)
	else:
            return sorted(self.children, key = lambda c: c.b_value)[-1]
           
    def ExpandBinary(self, down_lim, up_lim, s, time_type, value):
        n = Node(down_lim, up_lim, parent = self, state = s, last_move=value)
        self.children.append(n)
    
    def Update_B(self, reward, j):
	self.stats = max(self.stats, reward)
        rho = 1
        v = 1
	self.visits += 1
	#print j, np.log(j/self.visits)
	#sys.exit(0)
        U = self.stats + sqrt(abs(2*np.log((j+1)/self.visits))) + rho*v
        self.b_value = min(U, max(self.children[0].b_value, self.children[1].b_value))
        
    #UCB Functions
    def SelectUCT(self):
	return random.choice(self.children)

    def Expand(self, state, move):
        n = Node(parent=self, state=state, last_move=move)
        self.untried_moves.remove(move)
        self.children.append(n)

    def Update(self, value):
        self.visits += 1          
        self.stats = max(self.stats, value)
 
    #representation of self, if needed
    def __repr__(self):
        return str('Planets' + str(self.state.seq) + 'TOF=' + str(self.state.tof) + 'T0=' + str(self.state.t0))

    #def __str__(self, level=0):
    #    ret = "\t"*level+repr([self.last_move, self.b_value])+"\n"
    #    for child in self.children:
    #        ret += child.__str__(level+1)
    #    return ret

def MiddlePoint(down_lim, up_lim):
    return (down_lim + up_lim)/2


def HOO(rootnode, budget1, time_type):

    best = None
    finals = []
    terminals = []
    down_lim, up_lim = rootnode.state.lims()
    init = MiddlePoint(int(down_lim), int(up_lim))

    rootstate = rootnode.state.copy()
    state = rootstate #rootstate does not matter, it is being changed with state...

    #so now we have a rootnode and a rootstate. the state is the rootstate, for now. they don't matter because they don't belong to 
    #the HOO, they are just there for the appending in the end.

    #and now we create the first HOO node!

    state.move(init)
    node = Node(down_lim=down_lim, up_lim=up_lim, state=state, last_move = init)
    

    for i in range(budget1):

	state = copy.deepcopy(node.state)

    	while node.children != []:
	  #  print node.state.t0, node.children[0].state.t0, node.children[0].b_value, node.children[1].state.t0, node.children[0].b_value
            node = node.Select_B() #select with biggest b value
	   # print 'Chosen', node.state.t0, node.b_value
	   # if i == 20:
	#	sys.exit(0)
            state.UpdateLim(node.down_lim, node.up_lim)
	#print '\n'

        if node.children == []:
            state1 = state.copy()
            state2 = state.copy()

	    if time_type == 1: #T0
		t01 = MiddlePoint(state1.down_lim, node.last_move) #ta fixe
                state1.UpdateLim(state1.down_lim, node.last_move)
		state1.input_t0(t01)
                node.ExpandBinary(state1.down_lim, state1.up_lim, state1, time_type, t01) 
                
		t02 = MiddlePoint(node.last_move, state2.up_lim)
                state2.UpdateLim(node.last_move, state2.up_lim)
		state2.input_t0(t02)
                node.ExpandBinary(state2.down_lim, state2.up_lim, state2, time_type, t02) 

	    elif time_type ==2: #TOF

                tof1 = MiddlePoint(state1.down_lim, node.last_move) #ta fixe
                state1.UpdateLim(state1.down_lim, node.last_move)
		state1.input_tof(tof1)
                node.ExpandBinary(state1.down_lim, state1.up_lim, state1, time_type, tof1) 

                tof2 = MiddlePoint(node.last_move, state2.up_lim)
                state2.UpdateLim(node.last_move, state2.up_lim)
		state2.input_tof(tof2)
                node.ExpandBinary(state2.down_lim, state2.up_lim, state2, time_type, tof2) 
 
        #Simulation
        while not state.isterminal():
	    if state.moves() != []: #so we don't end up in cycles. if final, also terminal.
		state.random_move()
	    else:
		break

        value = 0
	
        if state.isfinal():
	    finals.append(node)
	    value = max(MAX_DV - state.dv, 0.)/ MAX_DV
	    #if value != 0:
	    #    print 'VALUE', value
	    node.Update_B(value, i)
	else:
	    terminals.append(node)
    
	    

	#Backpropagate
        while node.parent != None: # backpropagate from the expanded node and work back to the root 
	    if node.parent == None:
		break
            node = node.parent 

    finals = [k for k, g in itertools.groupby(sorted(finals, key = lambda c: c.last_move))]

    if len(finals) < 5:
	finals = [k for k, g in itertools.groupby(sorted(terminals, key = lambda c: c.last_move))]

    best = sorted(finals, key = lambda c: c.stats)
    best_stats = []

    for i in range(0, len(best)):
	if best[i].stats != 0:
	    best_stats = best[i:]
	    break

    if len(best_stats) > 5:
	best = best_stats #can also choose from this or IT WILL BE TOO SLOW!!! BUT WORKS
    else:
        best = best[-6:-1]

    rootnode.untried(None)
    rootnode.children = best
    for l in range(0, len(best)):
	best[l].parent = rootnode
	best[l].children = []

    return best

def HOOT(N):

    start = time.time()
    rootstate = State()
    root = Node(state=rootstate)
    n_rollouts = 0
    n_legs = 0
    best = None
    best_state = None
    node = root

    while (time.time() - start) < N:

        n_rollouts += 1
	state = copy.deepcopy(rootstate)
	n = 0
	count = 0

	#ENTER SELECTION
        while node.children != []:
	    if node.state.next_move == MOVE_TYPE.TOF:
		#print 'TOF select' 
		node = node.SelectUCT()
		state.move(node.last_move)
		node.state = state
		node.untried(state) #???
		if state.moves() == []:
		    n = 1
		    break
	    elif node.state.next_move == MOVE_TYPE.PLANET:
		#print 'Planet select' 
	    	node = node.SelectUCT()
		node.untried(state)
		state.move(node.last_move)
	    else:
		#print 'T0 select'
		node = node.SelectUCT()
		state.move(node.last_move)
	   

	#EXPANSION
	#for first node...
	if n == 0:
            if node.state.next_move == MOVE_TYPE.T0:
	        #print 'T0 expand', n_rollouts
	        nodes = HOO(node, budget1 = 1000, time_type = 1)
		node = random.choice(nodes)
		state.move(node.last_move)
	        node.untried(state)
	    elif node.state.next_move == MOVE_TYPE.TOF:
                n_legs += 1
		#print 'TOF expand'
		nodes = HOO(node, budget1 = 1000, time_type = 2)
		node = random.choice(nodes)
		state.move(node.last_move)
		node.untried(node.state)
	    elif node.state.next_move == MOVE_TYPE.PLANET:
        #ONLY if planet!!! then select again...
		#check = 0
		while len(node.untried_moves) != 0:
		    #check = 1
		    #print 'planet expand'
	            prev_state = state.copy()
                    move = node.untried_moves[0]
	            prev_state.move(move)
                    node.Expand(prev_state, move)
		#print 'planet select after expand'
		#print node, node.children
		#if check == 0:
		 #   print check
		 #   print node, state, node.untried_moves
		node = node.SelectUCT()
		state.move(node.last_move)

	# simulate

	if state.isterminal() and not state.isfinal():
	    state = node.parent.state
	    node.parent.children.remove(node)
	    node = node.parent
	else:
	    while not state.isterminal():
	        #print 'sim'
	        if state.moves() != [] and len(state.moves()) != 0: #REDUNDANCY
		    state.random_move()
	        else:
		    break
	
        # backpropagate
        value = 0
        if state.isfinal():
            value = max(MAX_DV - state.dv, 0.)/ MAX_DV
	    #if value != 0:
		#print value
            if best is None or state.dv < best.dv and state.dv != 0:
                best = state
                best_n_legs = n_legs

        while node != root:
            node.Update(value)
	    if node.parent == None:
		break
            node = node.parent 
	
    return best

if __name__=='__main__':
    
    set_t_res(32)
    fix_first_move(False)

    best = HOOT(N = 1200) #seconds!
    
    f = open('hootv3.txt', 'w')
    f.write(str('Time: ') + str(N) + str(' Sequence: ') + str(best) + str(' Function Counter ') + str(state_rosetta_continuous.FCN_COUNTER))
 
    f.close()
