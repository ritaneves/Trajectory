#!/usr/bin/env python
import sys
import time
import random
import numpy as np
import functools
import bisect
import math
import copy

from state_asteroids import State, MOVE_TYPE
from tools import pretty_time, conv_times

class Node:
    def __init__(self, parent=None, state=None, last_move=None, c_P=0.0007):
        self.Q = 0 # sum of values
        self.V = 0 # value estimate
        self.n = 0 # visits
        self.untried_moves = list(state.moves()) #HEREEEEEE ###TODO
        self.parent = parent
        self.children = []
        self.last_move = last_move
        self.state = state.copy()
        self.c_P = c_P

    def expand(self, state, move):
        n = Node(parent=self, state=state, last_move=move, c_P=self.c_P)
        self.untried_moves.remove(move)
        self.children.append(n)
        return n

    def update(self, value):
	self.n += 1
	self.V = max(self.V, value)

    def select(self, choice):
	if choice == 1: #random
	    return random.choice(self.children)
	elif choice == 2: #ucb1
            self.children.sort(key=lambda c: c.V + self.c_P * math.sqrt(np.log(self.n)/c.n))
            return self.children[-1]
	elif choice == 3: #greedy
	    self.children.sort(key=lambda c: c.V + self.c_P * self.n/c.n)
            return self.children[-1]
	elif choice == 4: #ucb1tuned
	    self.children.sort(key=ucb1tuned) 
            return self.children[-1]

    #representation of self, if needed
    def __repr__(self):
        return str('Planets' + str(self.state.seq) + 'TOF=' + str(self.state.tof) + 'T0=' + str(self.state.t0))

        
def uct(option, c_P, N):
    best = None
    best_n_legs = None
    n_rollouts = 0
    n_legs = 0
    start = time.time()

 #   f = open('test1.txt', 'w')
    rootstate = State()
    root = Node(state=rootstate, c_P=c_P)

    max_select_depth = 0
    
    while n_rollouts < N:
        n_rollouts += 1
        node = root

        # select
        select_depth = 0
        while node.untried_moves == [] and node.children != []:
	    node = node.select(option)            
            select_depth += 1

        #max_select_depth = max(max_select_depth, select_depth)

        # expand
        while node.untried_moves != []:
            move = random.choice(node.untried_moves)
            if node.state.next_move == MOVE_TYPE.TOF:
                n_legs += 1
            state = node.state.copy()
            state.move(move)
	    node = node.expand(state, move)
	    
            if node.state.next_move == MOVE_TYPE.ASTEROID and node.state.tof != []:		
	        if node.state.isterminal():
#		    print node.state.mass
		    break

	#backpropagate
           	           
        if best is None or (len(node.state.seq)-1) > best:
            best = len(node.state.seq) - 1
	    seq = node.state.seq[0:-1]

#	    f.write(str(best) + ' ' + str(node.state) + '\n')
	
        done = False
	value = (len(node.state.seq)-1)*(1/85) + (-1/85)
        
	while node is not None:
	    node.update(value)
            if node.children == [] and node.untried_moves == []:
                if node.parent is None:
                    done = True
                    break
                node.parent.children.remove(node)
            node = node.parent

        if done:
	    print 'break'
            break

    #f.close()
    return best, n_rollouts

   
if __name__=='__main__':
    
    length, n_roll = uct(option = 3, c_P=0.008, N = 1000000)


