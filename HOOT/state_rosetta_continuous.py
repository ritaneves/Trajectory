import copy
import numpy as np
import bisect

import PyKEP as kep
import tools

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
T_MIN = T0[0]
T_MAX = T0[-1] + 4000.

#Fit Grid Resolution
def set_t_res(t_res):
    global T_RES
    global T_SCALE
    T_RES = t_res
    T_SCALE = {name: np.arange(T_MIN, T_MAX, tools.PLANETS[name].compute_period(kep.epoch(0))/kep.DAY2SEC/T_RES) for name in PLANET_NAMES}

set_t_res(32)

FCN_COUNTER = 0
MAX_EPOCH = T_SCALE['67p'][-2]
MAX_FLYBYS = 5
MAX_DV = 10000.
MAX_DV_LAUNCH = 5000.
MOVE_TYPE = tools.enum('T0', 'PLANET', 'TOF')
MOVES = {
    #CHANGE HERE THE MOVE TYPE TO CONTINUOUS - FOR TOF AND T0? TODO
    MOVE_TYPE.T0: [t for t in T_SCALE['earth'] if t >= T0[0] and t <= T0[1]], 
    MOVE_TYPE.PLANET: PLANET_NAMES,
    MOVE_TYPE.TOF: None, # defined on the fly
}

#Fix T0
def fix_first_move(fixit):
    if fixit:
        MOVES[MOVE_TYPE.T0] = [min(T_SCALE['earth'], key=lambda x: abs(1550.-x))]
    else:
        MOVES[MOVE_TYPE.T0] = [t for t in T_SCALE['earth'] if t >= T0[0] and t <= T0[1]]
        

#Rosetta State    
class State():#PUT HERE DOWN AND UP LIM?? TODO
    
    def __init__(self, seq=['earth'], t0=None, tof=[], vrel=None, dv=0, next_move=MOVE_TYPE.T0, down_lim=None, up_lim=None):
        self.seq = copy.copy(seq) #Path of Planets
        self.tof = copy.copy(tof) #Path of TOF
        self.t0 = t0
        self.vrel = copy.copy(vrel)
        self.dv = dv
        self.next_move = next_move
	self.down_lim = down_lim
	self.up_lim = up_lim 

    #this is being done ok
    def lims(self):
	if self.next_move == MOVE_TYPE.TOF:
	    self.seq = [x for x in self.seq if not isinstance(x, float)]
	    self.down_lim, self.up_lim = PLANET_TOF[self.seq[-2], self.seq[-1]]	    
	elif self.next_move == MOVE_TYPE.T0:
	    self.down_lim = T0[0]
	    self.up_lim = T0[-1]
	return self.down_lim, self.up_lim

    def input_t0(self, t0):
	self.t0 = t0

    def input_tof(self, tof):
	self.tof = []
	self.tof.append(tof)

    def UpdateLim(self, down_lim, up_lim):
        self.down_lim = down_lim
        self.up_lim = up_lim

    #Get Possible Moves
    def moves(self):
        if self.isterminal():
            return []

	#RETURN NEXT MOVE
        if self.next_move is MOVE_TYPE.TOF:
	    self.seq = [x for x in self.seq if not isinstance(x, float)]
	    #PERCEBER QUE MERDA ESTA A ACONTECER AQUI TODO
            min_tof, max_tof = PLANET_TOF[self.seq[-2], self.seq[-1]]
            cur_t = self.t0 + sum(self.tof)
            lb = bisect.bisect(T_SCALE[self.seq[-1]], cur_t + min_tof)
            ub = bisect.bisect(T_SCALE[self.seq[-1]], cur_t + max_tof)
            return [t - cur_t for t in T_SCALE[self.seq[-1]][lb:ub]]
        else:
            return MOVES[self.next_move]

    #Get to next state (planet or TOF)
    def move(self, move):
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
            self.next_move = MOVE_TYPE.PLANET
            #Evaluate lambert leg
            t = self.t0 + sum(self.tof[:-1])
            tof = self.tof[-1]
	    self.seq = [x for x in self.seq if not isinstance(x, float)]
	    global FCN_COUNTER
	    FCN_COUNTER += 1
            dv, vrel_out = tools.lambert_leg(self.seq[-2], self.seq[-1], t, tof,
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
      #  if len(self.tof) == MAX_FLYBYS:
       #     return True
        if self.dv is not None and self.dv > MAX_DV:
            return True
        if self.t0 is not None and self.t0 + sum(self.tof) > MAX_EPOCH:
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
        return copy.deepcopy(self)

    def __key(self): #sorted
        return (self.seq, self.t0, self.tof)

    def __eq__(s1, s2): #equality
        return s1.__key() == s2.__key()

    def __hash__(self): #identifies a particular variable to use the __eq__
        return hash(self.__key())    
        
    def __repr__(self):
        s = '{:8.2f} m/s  '.format(self.dv)
        s += '{:7.2f} days  '.format(sum(self.tof))
        if self.t0 is not None:
            s += '{:8.2f} mjd2000  '.format(self.t0)
            s += '-'.join([p[0] for p in self.seq]) + '  '
            s += '[' + ', '.join(['{:.2f}'.format(t) for t in self.tof]) + ']'
        return s   
