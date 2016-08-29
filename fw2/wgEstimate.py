from scipy.optimize import minimize
import numpy as np

def argmin(f, fd, resLen, args):
    constrains = ({'type': 'eq',
             'fun': lambda x: np.sum(x) - 1.0,
             'jac': lambda x: np.ones(len(x))
             },
            {'type': 'ineq',
             'fun': lambda x: x,
             'jac': lambda x: np.identity(len(x))
             })

    init = [1.0/resLen]*resLen
    return minimize(f, init, args=args, jac=fd, options={'disp': False},
                   method='SLSQP', constraints=constrains)

class WeightPair:
    def __init__(self, id, wg):
        self.id = id
        self.weight = wg


def estimateWeight(data, guideSet):
    pass
