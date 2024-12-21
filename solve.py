from dataclasses import dataclass
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from amply import Amply
from enum import Enum
from parse import parsed

#--------------------------------------------------------------------------------------------------------------
# The data file contains (approximately):
# K => Resources (a set)
# Q => resource types (assignment done using the V parameter
# F => fronts(?), resource type needed is assigned by B param
# T => number of timeslots
# TF/TR/P/N/A => constraints on resource (availability)
# U => transit time
# C => water capacity
# S => total number of resources at each front constraints
# D/E => downloads (grabbing water)?
# W => water needed in each front
# a1/a2/a3 => weight values for the optimization function? (of we do not do lexicographical optimization)
# M => big M
# PR => priority of each front

#--------------------------------------------------------------------------------------------------------------
# Pre-compute the sizes of our sets so we can more easily index/iterate/use them

num_resources = 7 #len(aerial_resources) # len(k)
num_areas = 2 #len(areas) # len(m)
num_timeslots = 45 #len(timeslots) # len(t)
num_resource_types = 2 #max([resource.resource_type for resource in aerial_resources])

#Build matrices based on our classes (so we can easily sum in a loop)
V_qk = parsed['V']#[[1 if resource.resource_type == r_type else 0 for resource in aerial_resources] for r_type in range(num_resource_types)]
C_k = parsed['C']#[resource.water_capacity for resource in aerial_resources]
T_k = parsed['TF']#[resource.flight_duration for resource in aerial_resources]
R_k = parsed['TR']#[resource.min_rest_time for resource in aerial_resources]
N_k = parsed['N']#[resource.max_num_flights for resource in aerial_resources]
P_k = parsed['P']#[resource.crew_max_flights for resource in aerial_resources]
S_m = parsed['S']#[area.max_resources for area in areas]
I_m = [1.0 for area in range(num_areas)]  # TODO: not included in the dat files?

U_km = np.array(parsed['U']) # already indexed [k][m]
A_tk = np.array(parsed['A']) # Already [t][k]
W_tm = np.array(parsed['W']) # Already [t][m]
B_qm = np.array(parsed['B'])#np.zeros((num_resource_types, num_areas), dtype=np.int32)  # Already [q][m]

D_tkm = np.zeros((num_timeslots, num_resources, num_areas))
for t in range(num_timeslots):
    for k in range(num_resources):
        for m in range(num_areas):
            D_tkm[t, k, m] = parsed['D'][m][t][k] # ORIG: [m][t][k]


E_tkm = np.zeros((num_timeslots, num_resources, num_areas))
for t in range(num_timeslots):
    for k in range(num_resources):
        for m in range(num_areas):
            E_tkm[t, k, m] = parsed['E'][m][t][k] # ORIG: [m][t][k]

#Create the gurobi model
model = gp.Model("model")

#Decision variables
c_tkm = model.addVars(num_timeslots, num_resources, num_areas, name='c', vtype=GRB.BINARY)  # Flight towards area m

# Aux variables
s_tkm =  model.addVars(num_timeslots, num_resources, num_areas, name='s', vtype=GRB.BINARY)
e_tkm =  model.addVars(num_timeslots, num_resources, num_areas, name='e', vtype=GRB.BINARY)
x_tkm =  model.addVars(num_timeslots, num_resources, num_areas, name = 'x', vtype=GRB.BINARY)
v_tqm =  model.addVars(num_timeslots, num_resource_types, num_areas, name='v', vtype=GRB.BINARY)


TCplus_tm =  model.addVars(num_timeslots, num_areas, name='TC+', vtype=GRB.CONTINUOUS)
TCminus_tm =  model.addVars(num_timeslots, num_areas, name='TC-', vtype=GRB.CONTINUOUS)

TCplus_min = model.addVar(name='TC+min', vtype=GRB.CONTINUOUS)
TCminus_min = model.addVar(name='TC-min', vtype=GRB.CONTINUOUS)
# NOTE: TCmin = TCplus_min - TCminus_min
WO =  model.addVar(name='WO', vtype=GRB.CONTINUOUS)


# Define big M (based on the GH data)
big_M = parsed['M']
# TODO: sinle objective?
a1 = parsed['a1']
a2 = parsed['a2']
a3 = parsed['a3']



# CONSTRAINTS
# Flight initiation (3.5.1)
model.addConstrs(c_tkm.sum('*', k, '*') <= N_k[k] for k in range(num_resources))  # (6)

# (7)
for k in range(num_resources):
    for t in range(num_timeslots):
        model.addConstr(c_tkm.sum(t, k, '*') <= 1, name='(7)')

# (8)
# TODO: special attention
for k in range(num_resources):
    for t in range(num_timeslots - (T_k[k] + R_k[k] - 1)):
        # Split it up - pretty sure the first inline ranging had a bug
        lhs = 0
        for t_prime in range(T_k[k] + R_k[k]): # TODO: do we need a minus one here? We presumably want to include the -1 entry in the sum, and without it we wouldn's
            lhs += c_tkm.sum(t+t_prime, k, '*')

        model.addConstr(lhs <= 1, name='(8)')

# (9)
for k in range(num_resources):
    for m in range(num_areas):
        # print(f"x; {range(num_timeslots-(T_k[k]-2), num_timeslots)}")
        for t in range(num_timeslots-1-(T_k[k]-2), num_timeslots):  # TODO: where does the -2 come from? Also, do we need the -1?
            model.addConstr(c_tkm[t, k, m] == 0, name='(9)')


# 3.5.2: Phase of flight definitions
# (10)
for k in range(num_resources):  # NOTE: U_km msut be smaller than T_k (travel time smaller than max range, else the thing breaks)
    for m in range(num_areas):
        for t in range(num_timeslots - (T_k[k] - 1)):
            model.addConstr(e_tkm[t + int(U_km[k, m]), k, m] >= c_tkm[t, k, m], name='(10)')

# (11)
for k in range(num_resources):
    for m in range(num_areas):
        for t in range(num_timeslots - (T_k[k] - 1)):
            model.addConstr(e_tkm[t + (T_k[k] - 1) - int(U_km[k, m]), k, m] >= c_tkm[t, k, m], name='(11)')

# (12)
# TODO: special attention
for k in range(num_resources):
    for m in range(num_areas):
        for t in range(num_timeslots - (T_k[k] - 1)):
            for t_prime in range(int(U_km[k, m]+1), T_k[k] - int(U_km[k, m]) - 1): # -2 -> -1; we need to include the last 'index' of the t' range; the 0-indexing should be taken care of by our base index, iteration over num_timeslots
                model.addConstr(x_tkm[t + t_prime, k, m] >= c_tkm[t, k, m], name='(12)')

# (13)
for k in range(num_resources):
    for m in range(num_areas):
        model.addConstr(2 * c_tkm.sum('*', k, m) == e_tkm.sum('*', k, m), name='(13)')

# (14)
for k in range(num_resources):
    for m in range(num_areas):
        model.addConstr(((T_k[k] - 2 * U_km[k, m] - 2) * c_tkm.sum('*', k, m)) == x_tkm.sum('*', k, m), name='(14)')

# (15)
for t in range(num_timeslots):
    for k in range(num_resources):
        for m in range(num_areas):
            model.addConstr(x_tkm[t, k, m] + e_tkm[t, k, m] <= 1, name='(15)')

# 3.5.3: Resource activity constraints
# (16)
# TODO: special attention
for k in range(num_resources):
    for m in range(num_areas):
        for t in range(num_timeslots - (T_k[k] - 1)):
            for t_prime in range(T_k[k]): # Want to include T_k[k] - 1 in this range -> need to go one higher in the range operation (and unlike t, this is explicitly 0-idnexed in the paper)
                model.addConstr(s_tkm[t + t_prime, k, m] >= c_tkm[t, k, m], name='(16)')

# (17)
for k in range(num_resources):
    for m in range(num_areas):
        model.addConstr(T_k[k] * c_tkm.sum('*', k, m) == s_tkm.sum('*', k, m), name='(17)')

# (18)
for t in range(num_timeslots):
    for k in range(num_resources):
        model.addConstr(s_tkm.sum(t, k, '*') <= int(A_tk[t, k]), name='(18)')

# (19)
# TODO: special attention
for k in range(num_resources):
    for t in range(num_timeslots - P_k[k]):
        lhs = big_M * (1 - s_tkm.sum(t, k, '*'))

        rhs = 0
        for t_prime in range(P_k[k]-1, num_timeslots - t): # TODO: check 0-indexing
            rhs += s_tkm.sum(t+t_prime, k, '*')

        model.addConstr(lhs >= rhs, name='(19)')

# 3.5.4: Area of operation
# (20)
for m in range(num_areas):
    for t in range(num_timeslots):
        # NOTE: Split the sum into two; linear operator so should be the same

        model.addConstr((x_tkm.sum(t, '*', m) + e_tkm.sum(t, '*', m)) <= S_m[m], name='(20)')

# (21)
for m in range (num_areas):
    for t in range(num_timeslots):
        for q in range(num_resource_types):
            for k in range(num_resources):
                model.addConstr(v_tqm[t, q, m] >= V_qk[q][k] * (x_tkm[t, k, m] + e_tkm[t, k, m]), name='(21)')

# (22)
for m in range(num_areas):
    for t in range(num_timeslots):
        model.addConstr(v_tqm.sum(t, '*', m) <= 1, name='(22)')

# (23)
for m in range(num_areas):
    for q in range(num_resource_types):
        for t in range(num_timeslots):
            model.addConstr(v_tqm[t, q, m] >= int(B_qm[q, m]), name='(23)')

# 3.5.5: Target completion and water output constraints
# (24)
for t in range(num_timeslots):
    for m in range(num_areas):
        lhs = TCplus_tm[t, m] - TCminus_tm[t, m]

        rhs = 0
        for k in range(num_resources):
            rhs += C_k[k] * E_tkm[t, k, m] * e_tkm[t, k, m]
            rhs += C_k[k] * D_tkm[t, k, m] * x_tkm[t, k, m]
        rhs -= W_tm[t, m]

        model.addConstr(lhs == rhs, name='(24)')

# (25)
for t in range(num_timeslots):
    for m in range(num_areas):
        model.addConstr(TCplus_min - TCminus_min <= TCplus_tm[t, m] - TCminus_tm[t, m], name='(25)')

# (26)
WO_eqn = 0
for k in range(num_resources):
    for m in range(num_areas):
        for t in range(num_timeslots):
            WO_eqn += C_k[k] * (E_tkm[t, k, m] * e_tkm[t, k, m] + D_tkm[t, k, m] * x_tkm[t, k, m])
model.addConstr(WO == WO_eqn, name='(26)')


# Define the objective functions
f1 = 0
for t in range(num_timeslots):
    for m in range(num_areas):
        f1 += -I_m[m] * TCminus_tm[t, m]
f2 = TCplus_min - TCminus_min
f3  = WO

model.setObjective(f1, GRB.MAXIMIZE)

model.write('test.lp')
model.optimize()
print(model)

for v in model.getVars():
    print('%s %g' % (v.VarName, v.X))

print('Obj: %g' % model.ObjVal)

for t in range(num_timeslots):
    print(f"t={t}, x={x_tkm[t, 0, 1].X}, e={e_tkm[t, 0, 1].X}, s={s_tkm[t, 0, 1].X}")