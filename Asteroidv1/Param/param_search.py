#!/usr/bin/env python

import sys
import time
import random
import numpy as np
import functools
import bisect
import pickle

import multiprocessing

#from pykep_search.state_eph_grid import State, MOVE_TYPE, MAX_DV, fix_first_move, set_t_res
from state_asteroids import State, MOVE_TYPE
from tools import pretty_time

from MCTS import uct


def uct_run(c_P, N=50000):
    np.random.seed()
    return (uct(option = 2, c_P=c_P, N=N)[0], c_P)


if __name__=='__main__':
    for n_legs in [50000]:
        print '#' * 80

        from multiprocessing import Pool, Value
        from functools import partial

        N = 4000

        c_P = np.power(10, np.random.uniform(-3, 1, N))
        
        pool = Pool(7)
        res = []

        start = time.time()
        for i, r in enumerate(pool.imap_unordered(partial(uct_run, N=n_legs), c_P)):
            now = time.time()
            print '%d/%d  %.2f  %.3f  time: %s  remaining: %s' % (i+1, N, r[0], r[1],
                                                                  pretty_time(now - start),
                                                                  pretty_time((N-(i+1.))/(i+1.) * (now - start)))
            res.append(r)
            
            if (i+1) % 100 == 0:
                pickle.dump(res, open('cp_ast_%d.pkl' % n_legs, 'wb'))
