from __future__ import division
import copy
import numpy as np
import bisect
import math
import PyKEP as kep
import tools
import sys
import time

ASTEROID_TOF = [60, 90, 120, 150, 180, 210, 240, 270, 300, 330]

ASTEROID_NAMES = [523, 658, 752, 754, 786, 932, 973, 1151, 1554, 1892, 1913, 1971, 2640, 2658, 2788, 2814, 3178, 3574, 3836, 3855, 3874, 3932,
                  4066, 4439, 5032, 5402, 5758, 5765, 7713, 7879, 7889, 8255, 8426, 8487, 8499, 8511, 8868, 8978, 9224, 9622, 9634, 10533, 10751,
                  11069, 11088, 11557, 11954, 12475, 12489, 12508, 12526, 12548, 12572, 13836, 13959, 14167, 14248, 14983, 15039, 15381, 15681,
                  15773, 15776, 15948, 16000, 16005, 16040]

T0 = 8935 # launch window
MAX_MISSION_TIME = 6*365.25

T_MIN = T0
T_MAX = T0 + MAX_MISSION_TIME
T = 0.3
Isp = 3000
g = 9.81

global T_SCALE
T_SCALE = {name: np.arange(T_MIN, T_MAX, 30) for name in ASTEROID_NAMES}

#Ephemeris Calculations?
EPH = [[]]*len(ASTEROID_NAMES)

for i in range(0, len(ASTEROID_NAMES)):
    EPH[i] = []
    for time in T_SCALE[ASTEROID_NAMES[i]]:
	r1, v1 = tools.ASTEROIDS[ASTEROID_NAMES[i]].eph(time)
	EPH[i].append([r1, v1])

MOVE_TYPE = tools.enum('T0', 'ASTEROID', 'TOF')
MOVES = {
    MOVE_TYPE.T0: [0],
    MOVE_TYPE.ASTEROID: ASTEROID_NAMES,
    MOVE_TYPE.TOF: None, # defined on the fly
}

#Rosetta State    
class State():
    
    def __init__(self, seq=[658], t0=None, tof=[], dv=[], dv_damon = [], mass=2000, m_star = 0, next_move=MOVE_TYPE.T0):
        self.seq = copy.copy(seq) #Path of Planets
        self.tof = copy.copy(tof) #Path of TOF
        self.t0 = t0
        self.mass = mass
        self.dv = copy.copy(dv)
	self.dv_damon = copy.copy(dv_damon)
        self.next_move = next_move
	self.m_star = copy.copy(m_star)
        
    #Get Possible Moves
    def moves(self):
        if self.next_move is MOVE_TYPE.TOF:
            min_tof = ASTEROID_TOF[0]
            max_tof = ASTEROID_TOF[-1]
	    if len(self.tof) < 1:
		cur_t = T_SCALE[658][self.t0]
	    else: 
		cur_t = T_SCALE[self.seq[-2]][self.tof[-1]]	  

            lb = bisect.bisect(T_SCALE[self.seq[-1]], cur_t + min_tof)
            
	    if (cur_t + max_tof) > (T_SCALE[658][self.t0] + MAX_MISSION_TIME):
		ub = bisect.bisect(T_SCALE[self.seq[-1]], self.t0 + MAX_MISSION_TIME)

	    else:
	        ub = bisect.bisect(T_SCALE[self.seq[-1]], cur_t + max_tof)

	    return [i for i in range(lb-1, ub)] #return upper and lower bounds on t2

	elif self.next_move is MOVE_TYPE.ASTEROID:
	    mov = copy.copy(MOVES[self.next_move])	     
	    for ast in self.seq:
	        if ast in mov:
		    mov.remove(ast)
	    return mov		
	    
	else:	
            return MOVES[self.next_move]


    #Get to next state (planet or TOF)
    def move(self, move):
        if self.next_move == MOVE_TYPE.T0: #does not matter
            self.t0 = move
            self.next_move = MOVE_TYPE.ASTEROID
            return
        elif self.next_move == MOVE_TYPE.ASTEROID:
	    self.seq.append(move)
            self.next_move = MOVE_TYPE.TOF
            return
        elif self.next_move == MOVE_TYPE.TOF:  
	    self.tof.append(move)
	    t2 = self.tof[-1]

	    if len(self.tof) > 1:
	        t1 = self.tof[-2]
	    else:
		t1 = self.t0

	    tof = T_SCALE[self.seq[-1]][t2] - T_SCALE[self.seq[-2]][t1]

            self.next_move = MOVE_TYPE.ASTEROID

	    i = ASTEROID_NAMES.index(self.seq[-2])
	    j = ASTEROID_NAMES.index(self.seq[-1])

	    dv_lambert, dv_damon, m_star = tools.lambert_leg(self.seq[-2], self.seq[-1], i, j, t1, t2, tof)   

            self.dv.append(dv_lambert)
            self.dv_damon.append(dv_damon)

	    self.m_star = m_star

	    self.mass = self.mass * np.exp((-1*dv_lambert)/(Isp*g))

            return
        else:
            print 'unknown move type %s' % self.next_move
          
    #State is leaf, either because of constraints or because of goal reached            
    def isterminal(self): 
	if self.tof is not None and self.mass > self.m_star:
            return True

	#Extra conditions that should not be needed
	if self.seq[-1] in self.seq[0:-1]:
	    return True
	if len(self.seq)-len(self.tof) == 2:
	    if len(self.tof) > 0 and (T_SCALE[self.seq[-2]][self.tof[-1]] - T_SCALE[658][self.t0]) > MAX_MISSION_TIME:
		return True
	else: 
	    if len(self.tof) > 0 and (T_SCALE[self.seq[-1]][self.tof[-1]] - T_SCALE[658][self.t0]) > MAX_MISSION_TIME:
		return True

        return False

    def random_move(self, moves = []):
                 
        move = moves[np.random.randint(0, len(moves))]
        self.move(move)
        return move

    def copy(self):
        return State(seq=self.seq, t0=self.t0, tof=self.tof, dv=self.dv, next_move=self.next_move, m_star=self.m_star, mass=self.mass, dv_damon=self.dv_damon)

    def __key(self): #sorted
        return (self.seq, self.t0, self.tof)
	#return (self.seq[-1])

    def __eq__(s1, s2): #equality
        return s1.__key() == s2.__key()

    def __hash__(self): #identifies a particular variable to use the __eq__
        return hash(self.__key())    
        
    def __repr__(self):

	t0, tof = tools.conv_times(self.t0, self.tof, self.seq)

	if self.seq is not None:
	    s = str(self.seq) + ' '
        if self.t0 is not None:
            s += '{:8.2f} mjd2000  '.format(t0)
	if self.tof is not None:
	    s += '[' + ', '.join(['{:.2f}'.format(t) for t in tof]) + ']' + ' '
            s += str(str('M_star = ') + str(self.m_star)) + ' Mass: ' + str(self.mass)
	if len(self.dv) > 0 and len(self.dv_damon) > 0:
	    s += ' DV Lam: ' + str(self.dv[-1]) + ' DV Dam: ' + str(self.dv_damon[-1])
	return s   
