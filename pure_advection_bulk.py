
from scipy.linalg import solve
from matplotlib import pyplot
from typing import List


VELOCITY = 0.5


LENGTH = 0.5
N_NODES = 100
DELTA_X = LENGTH/N_NODES
DELTA_T = 0.01
TIME = 2.0
x = [0]
x += [DELTA_X/2 + i*DELTA_X for i in range(N_NODES)]
x += [LENGTH]


def get_solution(a_P:List[float],a_E:List[float],a_W:List[float],S_C_A:float,S_C_B:float,a_P0:List[float]=None):
    b0 = [0]*N_NODES
    b0[0] = S_C_A
    b0[-1] = S_C_B

    if(a_P0 is None): a_P0 = a_P
    mat_A = [([0]*N_NODES) for i in range(N_NODES)]
    mat_A[0][0], mat_A[0][1] = a_P[0], -a_E[0]
    mat_A[-1][-2], mat_A[-1][-1] = -a_W[-1], a_P[-1]
    for i in range(1,N_NODES-1): mat_A[i][i-1],mat_A[i][i],mat_A[i][i+1] = -a_W[i],a_P[i],-a_E[i]
    b = [b0[i]+a_P0[i]*rho_0[i] for i in range(len(b0))]
    return solve(mat_A, b)


plot_time = 0
PLOT_INTERVAL = 0.05

def plot_instant_solution(val_A:float, vals:List[float], val_B:float, force_plotting=False):
    global plot_time, PLOT_INTERVAL, x
    vals = list(vals)
    if(plot_time > PLOT_INTERVAL or force_plotting):
        vals.insert(0,val_A)
        vals.append(val_B)
        pyplot.plot(x,vals)
        plot_time = 0
    else: plot_time += DELTA_T


def show_balance(inflow=0,outflow=0,source=0,accumulation=0,init_amount=0):
    print('init amount:', init_amount)
    print('inflow:', inflow)
    print('outflow:', outflow)
    print('source:', source)
    print('accumulation:', accumulation)
    absolute_imbalance = outflow - inflow + accumulation
    denom = abs(init_amount)+abs(inflow)
    relative_imbalance = absolute_imbalance/denom if denom>0 else 1
    print('absolute imbalance:', absolute_imbalance)
    print('relative imbalance:', relative_imbalance)


def eps(xx): return 0 if xx<0 else 1
def pos(x): return x if x>0 else 0
def nu(x,c): return x if c>0 else 1



rho_A = 1
rho_B = 0.0
rho_0 = [0.5]*N_NODES


time = 0
inflow, outflow, accumulation = 0, 0, 0
init_amount = sum([rho_0[i] for i in range(len(rho_0))])*DELTA_X

plot_instant_solution(rho_A, rho_0, rho_B, force_plotting=True)

while(time<TIME):
    a_P0, a_P, a_E, a_W = [0]*N_NODES, [0]*N_NODES, [0]*N_NODES, [0]*N_NODES
    for i in range(N_NODES):
        a_P[i] = DELTA_X/DELTA_T + abs(VELOCITY)
        a_P0[i] = DELTA_X/DELTA_T
        a_W[i],a_E[i] = pos(VELOCITY), -pos(-VELOCITY)

    S_C_A, S_C_B = pos(VELOCITY)*rho_A, -pos(-VELOCITY)*rho_B 
    rho = get_solution(a_P,a_E,a_W,S_C_A,S_C_B,a_P0)

    if VELOCITY>0: rho_B = rho[-1]
    else: rho_A = rho[0]
    #store data for final plot
    plot_instant_solution(rho_A, rho, rho_B)

    #balance check
    inflow += VELOCITY*DELTA_T*(rho_A if VELOCITY>0 else -rho_B)
    outflow += VELOCITY*DELTA_T*(rho_B if VELOCITY>0 else -rho_A)
    accumulation += DELTA_X*sum([rho[i]-rho_0[i] for i in range(len(rho))])
    #update time and previous values
    rho_0 = rho
    time += DELTA_T


show_balance(inflow=inflow, outflow=outflow, accumulation=accumulation, init_amount=init_amount)


pyplot.show()

