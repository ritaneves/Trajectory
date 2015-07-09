#!/usr/bin/env python

import sys
import numpy as np
import time
from state_asteroids import State, MOVE_TYPE

LEGS = 0
BEST = None

def df_search(state, f):
    global BEST
    global LEGS

    if state.isterminal() and state.m_star != 0:
        if BEST is None or len(state.seq) > len(BEST.seq):
            BEST = state.copy()
	    f.write(str(BEST) + '\n')
            print '\r', BEST
            return

    moves = np.random.permutation(state.moves())

    for m in moves:
        ns = state.copy()
        if ns.next_move == MOVE_TYPE.TOF:
            LEGS += 1

        if LEGS % 100 == 0:
            print '\r', '{:,d}'.format(LEGS),
            sys.stdout.flush()
            
        ns.move(m)
        df_search(ns, f)


if __name__ == '__main__':

    f = open('sols.txt', 'w')
    start = time.time()
    BEST = None
  
    state = State()    
    df_search(state, f)

    f.close()
