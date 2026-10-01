import sys, numpy as np
from scipy.optimize import differential_evolution
args = sys.argv[1:]
sys.argv = ['x'] + args
exec(open('adversarial_ctrl.py').read().split('best = None')[0])
res = differential_evolution(lambda th: -min_rayleigh(th), [(-0.6, 0.6)] * len(free), seed=int(args[2]),
                             popsize=25, maxiter=400, tol=1e-12, polish=True, workers=1)
print('m=%s z=%s sym=%s  DE best min normalized Delta = %.3e' % (args[0], args[1], args[4], -res.fun))
print('  theta* =', np.round(res.x, 4))
