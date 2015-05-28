#!/usr/bin/env python

import time
import random
import numpy as np
import functools
import bisect
import math

from state_asteroids import State, MOVE_TYPE
from tools import pretty_time, conv_times

MAX_DV = 10000
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

    f = open('test1.txt', 'w')
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

	#print 'expand', state.seq, state.t0, state.tof

        # backpropagate
        value = 0
        if node.state.isterminal() and node.state.dv == node.state.dv: # TODO check why dv would be NaN?
            print node.state
	           
            if best is None or len(node.state.seq) > best:
                best = len(node.state.seq) - 1
		f.write(str(best))
		
        done = False
        while node is not None:
            if node.children == [] and node.untried_moves == []:
                if node.parent is None:
                    done = True
                    break
                node.parent.children.remove(node)
            node = node.parent
        if done:
            break

	print n_rollouts
    f.close()
    return best

   
if __name__=='__main__':
    
    uct(option = 1, c_P=0.015, N = 10000000)

