from __future__ import division
import copy
import numpy as np
import bisect
import math
import PyKEP as kep
import tools
import sys
import time

#TOF Interval definition
PLANET_TOF = {}
PLANET_TOF[('earth', 'venus')] = (30., 400.)
PLANET_TOF[('earth', 'earth')] = (340., 800.)
PLANET_TOF[('earth', 'mars')] = (100., 1000.)
PLANET_TOF[('earth', 'jupiter')] = (400., 1600.)
PLANET_TOF[('earth', '67p')] = (500., 3000.)

PLANET_TOF[('venus', 'venus')] = (224.7, 450.)
PLANET_TOF[('venus', 'mars')] = (100., 600.)
PLANET_TOF[('venus', 'jupiter')] = (400., 1600.)
PLANET_TOF[('venus', '67p')] = (500., 3000.)

PLANET_TOF[('mars', 'mars')] = (780., 1600.)
PLANET_TOF[('mars', 'jupiter')] = (400., 1600.)
PLANET_TOF[('mars', '67p')] = (500., 3000.)

PLANET_TOF[('jupiter', 'jupiter')] = (4332., 9000.)
PLANET_TOF[('jupiter', '67p')] = (500., 3000.)

#TOF the same coming and going
for (p1, p2) in PLANET_TOF.keys():
    PLANET_TOF[(p2, p1)] = PLANET_TOF[(p1, p2)]

#Comet = Final Destination
chury = kep.planet.keplerian(kep.epoch(2456879.5, 'jd'),
                   (
                       3.4630 * kep.AU, # a
                       0.64102, # e
                       7.0405 * kep.DEG2RAD, # i
                       50.147 * kep.DEG2RAD, # W
                       12.780 * kep.DEG2RAD, # w
                       303.71 * kep.DEG2RAD # M
                   ),
                   kep.MU_SUN,
                   1e-10, # mu_self
                   5.e3, # radius
                   100.e3, # save_radius
                   '67p'
)
tools.PLANETS[chury.name] = chury

#WHY?
for p in tools.PLANETS.values():
    if not p.name is 'jupiter':
        p.save_radius = 1.05
    
PLANET_NAMES = ['venus', 'earth', 'mars', 'jupiter', '67p']

T0 = (1460., 1825.) # launch window
MAX_MISSION_TIME = 4000
T_MIN = T0[0]
T_MAX = T0[-1] + MAX_MISSION_TIME

#T_res is the time sampling interval
def set_t_res(t_res):
    global T_RES
    global T_SCALE
    T_RES = t_res
    T_SCALE = {name: np.arange(T_MIN, T_MAX, tools.PLANETS[name].compute_period(kep.epoch(0))/+
kep.DAY2SEC/T_RES) for name in PLANET_NAMES}

#32 = 360/11.25
set_t_res(32)


#Ephemeris Calculations?
EPH = [[]]*len(PLANET_NAMES)

for i in range(0, len(PLANET_NAMES)):
    EPH[i] = []
    for time in T_SCALE[PLANET_NAMES[i]]:
	r1, v1 = tools.PLANETS[PLANET_NAMES[i]].eph(time)
	EPH[i].append([r1, v1])

#MAX_EPOCH = T_SCALE['67p'][-2]
#MAX_FLYBYS = 5
MAX_DV = 20000.
MAX_DV_LAUNCH = 5000.
MOVE_TYPE = tools.enum('T0', 'PLANET', 'TOF')
MOVES = {
    MOVE_TYPE.T0: [i for i in range(0, len([t for t in T_SCALE['earth'] if t >= T0[0] and t <= T0[1]]))],
    MOVE_TYPE.PLANET: PLANET_NAMES,
    MOVE_TYPE.TOF: None, # defined on the fly
}

#Fix T0
def fix_first_move(fixit):
    if fixit:
        MOVES[MOVE_TYPE.T0] = [10]
    else:
        MOVES[MOVE_TYPE.T0] = [i for i in range(0, len([t for t in T_SCALE['earth'] if t >= T0[0] and t <= T0[1]]))]
        

#Rosetta State    
class State():
    
    def __init__(self, seq=['earth'], t0=None, tof=[], vrel=None, dv=0, next_move=MOVE_TYPE.T0):
        self.seq = copy.copy(seq) #Path of Planets
        self.tof = copy.copy(tof) #Path of TOF
        self.t0 = t0
        self.vrel = copy.copy(vrel) #CHANGE? TODO
        self.dv = dv
        self.next_move = next_move
        
    #Get Possible Moves
    def moves(self):
        if self.isterminal():
            return []
        if self.next_move is MOVE_TYPE.TOF:
            min_tof, max_tof = PLANET_TOF[self.seq[-2], self.seq[-1]] #not index
	    if len(self.tof) < 1:
		cur_t = T_SCALE['earth'][self.t0]
	    else: 
		cur_t = T_SCALE[self.seq[-2]][self.tof[-1]]
	   	    	 
            lb = bisect.bisect(T_SCALE[self.seq[-1]], cur_t + min_tof)
            ub = bisect.bisect(T_SCALE[self.seq[-1]], cur_t + max_tof)

	    return [i for i in range(lb, ub)] #return upper and lower bounds on t2!! - moves that can be done on P2
	else:
            return MOVES[self.next_move]


    #Get to next state (planet or TOF)
    def move(self, move):
	# assert move in self.moves(), '%s is not in current move list' % move 
        # assert not self.isterminal(), 'current state %s is terminal' % self
        # assert not self.final(), 'current state %s is final' % self
        if self.next_move == MOVE_TYPE.T0:
            self.t0 = move
            self.next_move = MOVE_TYPE.PLANET
            return
        elif self.next_move == MOVE_TYPE.PLANET:
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
	   
            self.next_move = MOVE_TYPE.PLANET

	    i = PLANET_NAMES.index(self.seq[-2])
	    j = PLANET_NAMES.index(self.seq[-1])

	    dv, vrel_out = tools.lambert_leg(self.seq[-2], self.seq[-1], i, j, t1, t2, tof,
                                             vrel=self.vrel,
                                             rendezvous=self.isfinal())
            if len(self.tof) == 1:
                dv = max(dv - MAX_DV_LAUNCH, 0)
            self.dv += dv
            self.vrel = vrel_out
            return
        else:
            print 'unknown move type %s' % self.next_move
          
    #State is leaf, either because of constraints or because of goal reached            
    def isterminal(self):
#	if len(self.tof) == MAX_FLYBYS:
#	    return True
        if self.dv is not None and self.dv > MAX_DV:
            return True
	if len(self.seq)-len(self.tof) == 1:
            if self.t0 is not None and (len(self.tof) > 1 and (T_SCALE[self.seq[-1]][self.tof[-1]] - T_SCALE['earth'][self.t0]) > MAX_MISSION_TIME):
                return True
	elif len(self.seq)-len(self.tof) == 2:
	    if self.t0 is not None and (len(self.tof) > 1 and (T_SCALE[self.seq[-2]][self.tof[-1]] - T_SCALE['earth'][self.t0]) > MAX_MISSION_TIME):
		return True
        if self.isfinal():
            return True
        return False


    def random_move(self, continuous=False):
        # continuous random move
        if continuous and (self.next_move == MOVE_TYPE.T0 or self.next_move == MOVE_TYPE.TOF):
            move = np.random.uniform(min(MOVES[self.next_move]), max(MOVES[self.next_move]))
        else:
            moves = self.moves()
            move = moves[np.random.randint(0, len(moves))]
        self.move(move)
        return move

    #State is goal
    def isfinal(self):
        return self.seq[-1] == '67p' and len(self.seq) -1 == len(self.tof)
       
    def get_value():
        return self.dv

    def copy(self):
        return State(seq=self.seq, t0=self.t0, tof=self.tof, vrel=self.vrel, dv=self.dv, next_move=self.next_move)

    def __key(self): #sorted
        return (self.seq, self.t0, self.tof)
	#return (self.seq[-1])

    def __eq__(s1, s2): #equality
        return s1.__key() == s2.__key()

    def __hash__(self): #identifies a particular variable to use the __eq__
        return hash(self.__key())    
        
    def __repr__(self):
        s = '{:8.2f} m/s  '.format(self.dv)
	t0, tof = tools.conv_times(self.t0, self.tof, self.seq)
        if self.t0 is not None:
            s += '{:8.2f} mjd2000  '.format(t0)
            s += '-'.join([p[0] for p in self.seq]) + '  '
	    s += '[' + ', '.join(['{:.2f}'.format(t) for t in tof]) + ']'
        return s   
