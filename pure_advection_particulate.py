
from scipy.linalg import solve
from matplotlib import pyplot


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



VELOCITY = 0.5


LENGTH = 0.5
N_NODES = 50
DELTA_X = LENGTH/N_NODES
DELTA_T = 0.01
TIME = 2.0
x = [0]
x += [DELTA_X/2 + i*DELTA_X for i in range(N_NODES)]
x += [LENGTH]
plot_interval = 0.05


def eps(xx): return 0 if xx<0 else 1


rho_A = 1
m_A = 0.001
c_A = rho_A/m_A

rho_B = 0.0
m_B = 0.001
c_B = rho_B/m_B


grad_c_A = 0
grad_c_B = 0

grad_m_A = 0
grad_m_B = 0


rho_0 = [0.5]*N_NODES
m_0 = [0.001]*N_NODES
c_0 = [rho_0[i]/m_0[i] for i in range(len(m_0))]


time = 0
inflow = 0
outflow = 0
accumulation = 0
init_amount = sum([m_0[i]*c_0[i] for i in range(len(m_0))])*DELTA_X


vals = [c_A*m_A]
[vals.append(c_0[i]*m_0[i]) for i in range(len(m_0))]
vals.append(c_B*m_B if VELOCITY<0 else c_0[-1]*m_0[-1])
pyplot.plot(x,vals)
plot_time = 0


while(time<TIME):

    # =========================== particle concentration
    ac_P0 = 1/DELTA_T
    ac_P = ac_P0 + VELOCITY*(2*eps(VELOCITY)-1)/DELTA_X
    ac_E = -VELOCITY*(1-eps(VELOCITY))/DELTA_X
    ac_W = VELOCITY*eps(VELOCITY)/DELTA_X
    ac_E_A = -VELOCITY*(1-eps(VELOCITY))/DELTA_X
    ac_W_B = VELOCITY*eps(VELOCITY)/DELTA_X
    if VELOCITY>0:
        ac_P_A = ac_P
        ac_P_B = ac_P0 + VELOCITY*eps(VELOCITY)/DELTA_X
        Sc_C_A = VELOCITY*eps(VELOCITY)/DELTA_X*c_A
        Sc_C_B = -0.5*VELOCITY*(1-eps(VELOCITY))*grad_c_B
    else:
        ac_P_A = ac_P0 + VELOCITY*(eps(VELOCITY)-1)/DELTA_X
        ac_P_B = ac_P
        Sc_C_A = 0.5*VELOCITY*eps(VELOCITY)*grad_c_A
        Sc_C_B = -VELOCITY*(1-eps(VELOCITY))/DELTA_X*c_B
    bc0 = [0]*N_NODES
    bc0[0] = Sc_C_A
    bc0[-1] = Sc_C_B

    matrix_c_A = []
    first_row = [0]*N_NODES
    first_row[0], first_row[1] = ac_P_A, -ac_E_A
    last_row = [0]*N_NODES
    last_row[-2],last_row[-1] = -ac_W_B, ac_P_B
    matrix_c_A.append(first_row)
    for i in range(1,N_NODES-1):
        row = [0]*N_NODES
        row[i-1] = -ac_W
        row[i] = ac_P
        row[i+1] = -ac_E
        matrix_c_A.append(row)
    matrix_c_A.append(last_row)

    b = [bc0[i]+ac_P0*c_0[i] for i in range(len(bc0))]
    c = solve(matrix_c_A, b)
    if VELOCITY>0: c_B = c[-1]+0.5*DELTA_X*grad_c_B
    else: c_A = c[0]-0.5*DELTA_X*grad_c_A

    # =========================== particle mass
    am_P = [0]*N_NODES
    am_P0 = [0]*N_NODES
    am_E = [0]*N_NODES
    am_W = [0]*N_NODES

    cc = [c_A]
    [cc.append(c_i) for c_i in c]
    cc.append(c_B)

    C_P_MIN = 1000
    
    for i in range(N_NODES):
        c_P = max(cc[i+1], C_P_MIN)
        c_P_0 = max(c_0[i], C_P_MIN)
        am_P0[i] = c_P_0/DELTA_T
        am_P[i] = c_P/DELTA_T + VELOCITY/DELTA_X*(eps(VELOCITY)*cc[i]-(eps(VELOCITY)-1)*cc[i+2])
        am_E[i] = -VELOCITY/DELTA_X * (1-eps(VELOCITY))*cc[i+2]
        am_W[i] = VELOCITY/DELTA_X * eps(VELOCITY) * cc[i]

    am_E_A = am_E[0]
    am_W_B = am_W[-1]
    if VELOCITY>0:  
        am_P_A = am_P[0]
        am_P_B = am_P0[-1] + VELOCITY/DELTA_X*eps(VELOCITY)*cc[-1]
        Sm_C_A = VELOCITY/DELTA_X*eps(VELOCITY)*m_A*c_A
        Sm_C_B = -0.5*VELOCITY*(1-eps(VELOCITY))*grad_m_B*c_B
    else:
        am_P_A = am_P0[0] + VELOCITY/DELTA_X*(eps(VELOCITY)-1)*cc[0]
        am_P_B = am_P[-1]
        Sm_C_A = 0.5*VELOCITY*eps(VELOCITY)*grad_m_A*c_A
        Sm_C_B = -VELOCITY/DELTA_X*(1-eps(VELOCITY))*m_B*c_B

    bm0 = [0]*N_NODES
    bm0[0] = Sm_C_A
    bm0[-1] = Sm_C_B

    matrix_m_A = []
    first_row = [0]*N_NODES
    first_row[0], first_row[1] = am_P_A, -am_E_A
    last_row = [0]*N_NODES
    last_row[-2],last_row[-1] = -am_W_B, am_P_B
    matrix_m_A.append(first_row)
    for i in range(1,N_NODES-1):
        row = [0]*N_NODES
        row[i-1] = -am_W[i]
        row[i] = am_P[i]
        row[i+1] = -am_E[i]
        matrix_m_A.append(row)
    matrix_m_A.append(last_row)
 
    b = [bm0[i]+am_P0[i]*m_0[i] for i in range(len(bm0))]
    m = solve(matrix_m_A, b)
    if VELOCITY>0: m_B = m[-1]+0.5*DELTA_X*grad_m_B
    else: m_A = m[0]-0.5*DELTA_X*grad_m_A

    #store data for final plot
    if(plot_time > plot_interval):
        vals = [c_A*m_A]
        [vals.append(c[i]*m[i]) for i in range(len(m))]
        vals.append(c_B*m_B)
        pyplot.plot(x,vals)
        plot_time = 0
    else:
        plot_time += DELTA_T

    #balance check
    inflow += VELOCITY*DELTA_T*(c_A*m_A if VELOCITY>0 else -c_B*m_B)
    outflow += VELOCITY*DELTA_T*(c_B*m_B if VELOCITY>0 else -c_A*m_A)
    accumulation += DELTA_X*sum([c[i]*m[i]-c_0[i]*m_0[i] for i in range(len(m))])
    
    c_0 = c
    m_0 = m
    time += DELTA_T


show_balance(inflow=inflow, outflow=outflow, accumulation=accumulation, init_amount=init_amount)



pyplot.show()

