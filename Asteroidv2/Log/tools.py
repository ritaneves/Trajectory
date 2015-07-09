import numpy as np
import PyKEP as kep
import sys
import state_asteroids

ASTEROIDS = {name: kep.planet.gtoc7(name) for name in [523, 658, 752, 754, 786, 932, 973, 1151, 1554, 1892, 1913, 1971, 2640, 2658, 2788, 2814, 3178, 3574, 3836, 3855, 3874, 3932,
                                                      4066, 4439, 5032, 5402, 5758, 5765, 7713, 7879, 7889, 8255, 8426, 8487, 8499, 8511, 8868, 8978, 9224, 9622, 9634, 10533, 10751,
                                                      11069, 11088, 11557, 11954, 12475, 12489, 12508, 12526, 12548, 12572, 13836, 13959, 14167, 14248, 14983, 15039, 15381, 15681,
                                                      15773, 15776, 15948, 16000, 16005, 16040]}

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


def lambert_leg(P1, P2, t0, tof):
    
    ast1 = ASTEROIDS[P1]
    ast2 = ASTEROIDS[P2]

    r1, v1 = ast1.eph(kep.epoch(t0))
    r2, v2 = ast2.eph(kep.epoch(t0 + tof))

    lambert = kep.lambert_problem(r1, r2, tof * kep.DAY2SEC, ast1.mu_central_body)

    vrel_in = tuple(map(lambda x, y: -x + y, lambert.get_v1()[0], v1))
    vrel_out = tuple(map(lambda x, y: -x + y, lambert.get_v2()[0], v2))

    dv_lambert = np.linalg.norm(vrel_out) + np.linalg.norm(vrel_in)

    a, _, _, dv_damon = kep.damon(vrel_in, vrel_out, tof*kep.DAY2SEC)
    m_star = kep.max_start_mass(np.linalg.norm(a), dv_damon, T_max, Isp)
    
    return dv_lambert, dv_damon, m_star
