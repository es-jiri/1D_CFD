from scipy.linalg import solve
from matplotlib import pyplot


DENSITY = 1e03
CONDUCTIVITY = 0
VELOCITY = 0.1


LENGTH = 0.5
N_NODES = 200
DELTA_X = LENGTH/N_NODES
DELTA_T = 0.01
TIME = 2.0
x = [0]
x += [DELTA_X/2 + i*DELTA_X for i in range(N_NODES)]
x += [LENGTH]


def eps(xx): return 0 if xx<0 else 1


phi_A = 1
phi_B = 0


a_P0 = DELTA_X/DELTA_T*DENSITY
a_P = a_P0 + DENSITY*VELOCITY*(2*eps(VELOCITY)-1) + 2*CONDUCTIVITY
a_E = CONDUCTIVITY - DENSITY*VELOCITY*(1-eps(VELOCITY))
a_W = CONDUCTIVITY + DENSITY*VELOCITY*eps(VELOCITY)
a_E_A = CONDUCTIVITY
a_W_B = CONDUCTIVITY
a_P_A = a_P + CONDUCTIVITY
a_P_B = a_P + CONDUCTIVITY
S_C_A = (2*CONDUCTIVITY+DENSITY*VELOCITY*eps(VELOCITY))*phi_A
S_C_B = (2*CONDUCTIVITY-DENSITY*VELOCITY*(1-eps(VELOCITY)))*phi_B

b0 = [0]*N_NODES
b0[0] = S_C_A
b0[-1] = S_C_B

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


phi_0 = [0]*N_NODES
time = 0


while(time<TIME):
    b = [b0[i]+a_P0*phi_0[i] for i in range(len(b0))]
    phi = solve(matrix_A, b)
    phi_0 = phi
    time += DELTA_T
    vals = [phi_A]
    [vals.append(item) for item in phi]
    vals.append(phi_B)
    pyplot.plot(x,vals)

pyplot.show()

