import copy
import numpy as np
import bisect
import sys
import PyKEP as kep
import tools



PLANET_TOF = {}
PLANET_TOF[('earth', 'venus')] = (30., 400.)
PLANET_TOF[('earth', 'earth')] = (365.25, 800.)
PLANET_TOF[('earth', 'mars')] = (100., 500.)
PLANET_TOF[('earth', 'jupiter')] = (400., 1600.)
PLANET_TOF[('earth', 'saturn')] = (1200., 2500.)

PLANET_TOF[('venus', 'venus')] = (224.7, 450.)
PLANET_TOF[('venus', 'mars')] = (100., 600.)
PLANET_TOF[('venus', 'jupiter')] = (400., 1600.)
PLANET_TOF[('venus', 'saturn')] = (1200., 2500.)

PLANET_TOF[('mars', 'mars')] = (780., 1600.)
PLANET_TOF[('mars', 'jupiter')] = (400., 1600.)
PLANET_TOF[('mars', 'saturn')] = (1200., 2500.)

PLANET_TOF[('jupiter', 'jupiter')] = (4332., 9000.)
PLANET_TOF[('jupiter', 'saturn')] = (800., 2500.)


for (p1, p2) in PLANET_TOF.keys():
    PLANET_TOF[(p2, p1)] = PLANET_TOF[(p1, p2)]

PLANET_NAMES = ['venus', 'earth', 'mars', 'jupiter', 'saturn']

MAX_MISSION_TIME = 2200
T0 = (-900., -700.) # launch window
T_MIN = T0[0]
T_MAX = T0[-1] + MAX_MISSION_TIME

def set_t_res(t_res):
    global T_RES
    global T_SCALE
    T_RES = t_res
    T_SCALE = {name: np.arange(T_MIN, T_MAX, tools.PLANETS[name].compute_period(kep.epoch(0))/kep.DAY2SEC/T_RES) for name in PLANET_NAMES}

set_t_res(32)


print T_SCALE['earth']
sys.exit(0)

EPH = [[]]*len(PLANET_NAMES)

for i in range(0, len(PLANET_NAMES)):
    EPH[i] = []
    for time in T_SCALE[PLANET_NAMES[i]]:
	r1, v1 = tools.PLANETS[PLANET_NAMES[i]].eph(time)
	EPH[i].append([r1, v1])

#print '\n'.join('time grid - %s: %d' % (name, len(T_SCALE[name])) for name in PLANET_NAMES)

MAX_EPOCH = T_SCALE['saturn'][-2]
#MAX_FLYBYS = 5
MAX_DV = 10000

MOVE_TYPE = tools.enum('T0', 'PLANET', 'TOF')
MOVES = {
    MOVE_TYPE.T0: [i for i in range(0, len([t for t in T_SCALE['earth'] if t >= T0[0] and t <= T0[1]]))], 
    MOVE_TYPE.PLANET: PLANET_NAMES,
    MOVE_TYPE.TOF: None, # defined on the fly
}


def fix_first_move(fixit):
    if fixit:
        MOVES[MOVE_TYPE.T0] = [10]
    else:
        MOVES[MOVE_TYPE.T0] = [i for i in range(0, len([t for t in T_SCALE['earth'] if t >= T0[0] and t <= T0[1]]))]
        
    
class State():
    """State class."""
    
    def __init__(self, seq=['earth'], t0=None, tof=[], vrel=None, dv=0, next_move=MOVE_TYPE.T0):
        self.seq = copy.copy(seq)
        self.tof = copy.copy(tof)
        self.t0 = t0
        self.vrel = copy.copy(vrel)
        self.dv = dv
        self.next_move = next_move
        

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
#	    print [i for i in range(lb, ub)]
	    return [i for i in range(lb, ub)] #return upper and lower bounds on t2!! - moves that can be done on P2
	else:
            return MOVES[self.next_move]


    def move(self, move):
#        assert move in self.moves(), '%s is not in current move list' % move
#        assert not self.isterminal(), 'current state %s is terminal' % self
#        assert not self.final(), 'current state %s is final' % self
        if self.next_move == MOVE_TYPE.T0:
            self.t0 = move
            self.next_move = MOVE_TYPE.PLANET
            return
        elif self.next_move == MOVE_TYPE.PLANET:
            self.seq.append(move)
            self.next_move = MOVE_TYPE.TOF
#            if not self.isterminal():
#                self.random_move() # TODO: WOW don't use this
            return
        elif self.next_move == MOVE_TYPE.TOF:
	    self.tof.append(move)
	    t2 = self.tof[-1]
#	    print len(self.tof)
#	    print T_SCALE[self.seq[-1]][self.tof[-1]], len(self.tof)
	    if len(self.tof) > 1:
	        t1 = self.tof[-2]
	    else:
		t1 = self.t0

	    tof = T_SCALE[self.seq[-1]][t2] - T_SCALE[self.seq[-2]][t1]
	   
            self.next_move = MOVE_TYPE.PLANET

	    i = PLANET_NAMES.index(self.seq[-2])
	    j = PLANET_NAMES.index(self.seq[-1])

	    #nada errado aqui no tools...
	    dv, vrel_out = tools.lambert_leg(self.seq[-2], self.seq[-1], i, j, t1, t2, tof,
                                             vrel=self.vrel,
                                             rendezvous=self.isfinal())
            if len(self.tof) == 1:
                dv = max(dv - 5000, 0)
            self.dv += dv
            self.vrel = vrel_out
            return
        else:
            print 'unknown move type %s' % self.next_move


    def random_move(self, continuous=False):
        # continuous random move
        if continuous and (self.next_move == MOVE_TYPE.T0 or self.next_move == MOVE_TYPE.TOF):
            move = np.random.uniform(min(MOVES[self.next_move]), max(MOVES[self.next_move]))
        else:
            moves = self.moves()
            move = moves[np.random.randint(0, len(moves))]
        self.move(move)
        return move
            
            
    def isterminal(self):
        #if len(self.tof) == MAX_FLYBYS:
        #    return True
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


    def isfinal(self):
        return self.seq[-1] == 'saturn' and len(self.seq) -1 == len(self.tof)
        

    def get_value():
        return self.dv

        
    def copy(self):
        return State(seq=self.seq, t0=self.t0, tof=self.tof, vrel=self.vrel, dv=self.dv, next_move=self.next_move)

    def __key(self):
        return (self.seq, self.t0, self.tof)


    def __eq__(s1, s2):
        return s1.__key() == s2.__key()


    def __hash__(self):
        return hash(self.__key())
        
        
    def __repr__(self):
        s = '{:8.2f} m/s  '.format(self.dv)
	t0, tof = tools.conv_times(self.t0, self.tof, self.seq)

        if self.t0 is not None:
            s += '{:8.2f} mjd2000  '.format(t0)
        s += '-'.join([p[0] for p in self.seq]) + '  '
        s += '[' + ', '.join(['{:.2f}'.format(t) for t in tof]) + ']'
        return s

    
