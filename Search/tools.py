import numpy as np
import PyKEP as kep
import sys
import state_rosetta

PLANETS = {name: kep.planet.jpl_lp(name) for name in ['venus', 'earth', 'mars', 'jupiter', 'saturn']}


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
    new_t0 = state_rosetta.T_SCALE['earth'][t0]
    
    for i in range(0, len(tof)):
	if i == 0:
	    x = new_t0
	else:
	    x = state_rosetta.T_SCALE[seq[i]][tof[i-1]]
	y = state_rosetta.T_SCALE[seq[i+1]][tof[i]]
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

    p1 = PLANETS[str(P1)]
    p2 = PLANETS[str(P2)]

    r1 = state_rosetta.EPH[i][t1][0]
    v1 = state_rosetta.EPH[i][t1][1]
    r2 = state_rosetta.EPH[j][t2][0]
    v2 = state_rosetta.EPH[j][t2][1]

    lambert = kep.lambert_problem(r1, r2, tof * kep.DAY2SEC, p1.mu_central_body, False, 0)

    vrel_in = tuple(map(lambda x, y: x - y, lambert.get_v1()[0], v1))
    vrel_out = tuple(map(lambda x, y: x - y, lambert.get_v2()[0], v2))

    if vrel is None:
        # launch
        dv = max(np.linalg.norm(vrel_in) - dv_launch, 0)
    else:
        # flyby
        #print p1.name, p2.name, np.linalg.norm(vrel_in), np.linalg.norm(vrel_out)
        dv = kep.fb_vel(vrel, vrel_in, p1)

    if rendezvous:
        dv += np.linalg.norm(vrel_out)
        
    return dv, vrel_out
