import numpy as np
import PyKEP as kep
import sys
import state_asteroids

ASTEROIDS = {name: kep.planet.gtoc7(name) for name in [523, 658, 752, 754, 786, 932, 973, 1151, 1554, 1892, 1913, 1971, 2640, 2658, 2788, 2814, 3178, 3574, 3836, 3855, 3874, 3932,
                                                      4066, 4439, 5032, 5402, 5758, 5765, 7713, 7879, 7889, 8255, 8426, 8487, 8499, 8511, 8868, 8978, 9224, 9622, 9634, 10533, 10751,
                                                      11069, 11088, 11557, 11954, 12475, 12489, 12508, 12526, 12548, 12572, 13836, 13959, 14167, 14248, 14983, 15039, 15381, 15681,
                                                      15773, 15776, 15948, 16000, 16005, 16040, 658, 932, 1151, 1554, 1892, 2640, 3932, 4439, 7889, 8255, 8487, 8511, 10751, 11954,
                                                      12475, 13836, 14248, 15773, 16000]}

T_max = 0.3
g = 9.81
Isp = 3000

def enum(*sequential, **named):
    """Helper function to create enums."""
    enums = dict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['reverse_mapping'] = reverse
    return type('Enum', (), enums)


def pretty_time(secs):
    """Returns a human readable time string."""
    hs, secs = divmod(secs, 60*60)
    mins, secs = divmod(secs, 60)
    return "%d:%02d:%04.1f" % (hs, mins, secs)   


def conv_times(t0, tof, seq):
    """Returns the real tof."""
    new_tof = []
    new_t0 = state_asteroids.T_SCALE[658][t0]
    
    for i in range(0, len(tof)):
	if i == 0:
	    x = new_t0
	else:
	    x = state_asteroids.T_SCALE[seq[i]][tof[i-1]]
	y = state_asteroids.T_SCALE[seq[i+1]][tof[i]]
	new_tof.append(y-x)

    return new_t0, new_tof     


def lambert_leg(P1, P2, i, j, t1, t2, tof, vrel=None, dv_launch=0., rendezvous=False):
    """Compute a lambert leg from planet to planet.
    Arguments:
    p1 -- starting planet (str or PyKEP.planet object)
    p2 -- final planet (str or PyKEP.planet object)
    t0 -- start time of leg in MJD2000
    tof -- time of flight in days
    
    Keyword arguments:
    vrel -- caresian coordinates of the relative velocity before the flyby at p1
    dv_launch -- dv discounted at lunch (i.e. if vrel is None)
    rendezvous -- add final dv
    Returns:
    dV, vrel_out, where vrel_out is the relative velocity at the end of the leg at p2
    """

    ast1 = ASTEROIDS[P1]
    ast2 = ASTEROIDS[P2]

    r1 = state_asteroids.EPH[i][t1][0]
    v1 = state_asteroids.EPH[i][t1][1]
    r2 = state_asteroids.EPH[j][t2][0]
    v2 = state_asteroids.EPH[j][t2][1]

    lambert = kep.lambert_problem(r1, r2, tof * kep.DAY2SEC, p1.mu_central_body, False, 0)

    vrel_in = tuple(map(lambda x, y: x - y, lambert.get_v1()[0], v1))
    vrel_out = tuple(map(lambda x, y: x - y, lambert.get_v2()[0], v2))

    dv_lambert = np.linalg.norm(vrel_out) + np.linalg.norm(vrel_in)

    a, _, _, dv_damon = kep.damon(vrel_in, vrel_out, tof*kep.DAY2SEC)
    m_star = kep.max_start_mass(np.linalg.norm(a), dv_damon, T_max, Isp)
        
    return dv_lambert, dv_damon, m_star
