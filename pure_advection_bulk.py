from scipy.linalg import solve
from matplotlib import pyplot


VELOCITY = 0.5


LENGTH = 0.75
N_NODES = 100
DELTA_X = LENGTH/N_NODES
DELTA_T = 0.02
TIME = 2.0
x = [0]
x += [DELTA_X/2 + i*DELTA_X for i in range(N_NODES)]
x += [LENGTH]


def eps(xx): return 0 if xx<0 else 1


rho_A = 1
rho_B = 0

grad_A = 0
grad_B = 0


a_P0 = 1/DELTA_T


a_P = a_P0 + VELOCITY*(2*eps(VELOCITY)-1)/DELTA_X
a_E = -VELOCITY*(1-eps(VELOCITY))/DELTA_X
a_W = VELOCITY*eps(VELOCITY)/DELTA_X
a_E_A = -VELOCITY*(1-eps(VELOCITY))/DELTA_X
a_W_B = VELOCITY*eps(VELOCITY)/DELTA_X


if VELOCITY>0:
    a_P_A = a_P
    a_P_B = a_P0 + VELOCITY*eps(VELOCITY)/DELTA_X
    S_C_A = VELOCITY*eps(VELOCITY)/DELTA_X*rho_A
    S_C_B = -0.5*VELOCITY*(1-eps(VELOCITY))*grad_B
else:
    a_P_A = a_P0 + VELOCITY*(eps(VELOCITY)-1)/DELTA_X
    a_P_B = a_P
    S_C_A = 0.5*VELOCITY*eps(VELOCITY)*grad_A
    S_C_B = -VELOCITY*(1-eps(VELOCITY))/DELTA_X*rho_B


b0 = [0]*N_NODES
b0[0] = S_C_A
b0[-1] = S_C_B
rho_0 = [0.5]*N_NODES
time = 0

init_amount = sum(rho_0)*DELTA_X


absolute_imbalance = 0
inflow = 0
outflow = 0
accumulation = 0


matrix_A = []
first_row = [0]*N_NODES
first_row[0], first_row[1] = a_P_A, -a_E_A
last_row = [0]*N_NODES
last_row[-2],last_row[-1] = -a_W_B, a_P_B
matrix_A.append(first_row)
for i in range(1,N_NODES-1):
    row = [0]*N_NODES
    row[i-1] = -a_W
    row[i] = a_P
    row[i+1] = -a_E
    matrix_A.append(row)
matrix_A.append(last_row)


while(time<TIME):
    b = [b0[i]+a_P0*rho_0[i] for i in range(len(b0))]
    rho = solve(matrix_A, b)

    if VELOCITY>0: rho_B = rho[-1]+0.5*DELTA_X*grad_B
    else: rho_A = rho[0]-0.5*DELTA_X*grad_A
    time += DELTA_T
    vals = [rho_A]
    [vals.append(item) for item in rho]
    vals.append(rho_B)
    pyplot.plot(x,vals)
    inflow += VELOCITY*DELTA_T*(rho_A if VELOCITY>0 else -rho_B)
    outflow += VELOCITY*DELTA_T*(rho_B if VELOCITY>0 else -rho_A)
    accumulation += DELTA_X*sum([rho[i]-rho_0[i] for i in range(len(rho))])
    absolute_imbalance = outflow - inflow + accumulation
    relative_imbalance = absolute_imbalance/(abs(init_amount)+abs(inflow))
    rho_0 = rho


print('inflow:', inflow)
print('outflow:', outflow)
print('accumulation:', accumulation)
print('absolute imbalance:', absolute_imbalance)
print('relative imbalance:', relative_imbalance)


pyplot.show()

