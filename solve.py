from dataclasses import dataclass
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from amply import Amply
from parse import parsed

from enum import Enum
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

#------------------------------------------------INPUT DATA--------------------------------------------------
# Pre-compute the sizes of our sets so we can more easily index/iterate/use them

num_resources = 7 #K - set of aerial resources
num_areas = 2 #M - set of areas
num_timeslots = 45 #T - set of timeslots
num_resource_types = 2 #Q - set of resource types

#Build matrices based on our classes (so we can easily sum in a loop)
V_qk = parsed['V']# Vqk - Resource type: 1 if resource k is of type q, 0 otherwise
C_k = parsed['C']# Ck - Water capacity of resource k
T_k = parsed['TF']# Tk - Flight time (in time slots) of resource k
R_k = parsed['TR']# Rk - Minimum rest time (in time slots) between flights for resource k
N_k = parsed['N']# Nk - Maximum number of flights in the entire time period for resource k
P_k = parsed['P']# Maximum number of time slots the crew operating resource k can be physically present
S_m = parsed['S']# 1 if resource k is active in a flight associated with area m during time slot t
I_m = [1.0 for area in range(num_areas)]  #Importance of area m

#Ukm - Number of pure transit time slots required for resource k to travel to/from area m
U_km = np.array(parsed['U']) # Already indexed [k][m]

#Atk - Availability: 1 if resource k is available in time slot t , 0 otherwise
A_tk = np.array(parsed['A']) # Already [t][k]

#Minimum water target required at area m at time t (in liters)
W_tm = np.array(parsed['W']) # Already [t][m]

#Bqm - Area restrictions: 1 if only resources of type q can operate in area m, 0 if unrestricted
B_qm = np.array(parsed['B']) # Already [q][m]

#---------------------------------------------------------------------
#Dtkm - Average number of drops an aerial resource can make in area m in time slot t
D_tkm = np.zeros((num_timeslots, num_resources, num_areas))
for t in range(num_timeslots):
    for k in range(num_resources):
        for m in range(num_areas):
            D_tkm[t, k, m] = parsed['D'][m][t][k] # ORIG: [m][t][k]

#Etkm - Average number of drops an aerial resource k can make in area m in a transit time slot
E_tkm = np.zeros((num_timeslots, num_resources, num_areas))
for t in range(num_timeslots):
    for k in range(num_resources):
        for m in range(num_areas):
            E_tkm[t, k, m] = parsed['E'][m][t][k] # ORIG: [m][t][k]

#---------------------------------MODEL-----------------------------------
model = gp.Model("model")

#Decision variables
c_tkm = model.addVars(num_timeslots, num_resources, num_areas, name='c', vtype=GRB.BINARY)  #Flight towards area m

#Auxiliary variables -> 1 OR 0
s_tkm =  model.addVars(num_timeslots, num_resources, num_areas, name='s', vtype=GRB.BINARY) #Initiates a flight towards m in t
e_tkm =  model.addVars(num_timeslots, num_resources, num_areas, name='e', vtype=GRB.BINARY) #Active in a flight associated with m during t
x_tkm =  model.addVars(num_timeslots, num_resources, num_areas, name = 'x', vtype=GRB.BINARY) #Active if en route during t & to m
v_tqm =  model.addVars(num_timeslots, num_resource_types, num_areas, name='v', vtype=GRB.BINARY) #Active if fighting the fire during t & in m

TCplus_tm =  model.addVars(num_timeslots, num_areas, name='TC+', vtype=GRB.CONTINUOUS) #Surplus of water target completion in t & m
TCminus_tm =  model.addVars(num_timeslots, num_areas, name='TC-', vtype=GRB.CONTINUOUS) #Deficit of water target completion in t & m

TCplus_min = model.addVar(name='TC+min', vtype=GRB.CONTINUOUS) #Min. water target completion over all t & m
TCminus_min = model.addVar(name='TC-min', vtype=GRB.CONTINUOUS)#Min. water target completion over all t & m

WO =  model.addVar(name='WO', vtype=GRB.CONTINUOUS) #Water output (WO) - total water output of all m in all t

#Define big M (based on the GH data)
big_M = parsed['M']

#Weights for the single optimization function
a1 = parsed['a1']
a2 = parsed['a2']
a3 = parsed['a3']

#---------------------------------CONSTRAINTS-----------------------------------
#Flight initiation Constraints

#QUESTION THIS - SHOULD IT NOT LOOP THROUGH K?
#Equation (1.2) - Ensures no. of flights initiated  than are permitted for k
model.addConstrs(
    (c_tkm.sum('*', k, '*') <= N_k[k] for k in range(num_resources)),
    name="(1.2)" 
)

#(1.3) - Ensures a resources goes to one area during a time slot 
for k in range(num_resources):
    for t in range(num_timeslots):
        model.addConstr(c_tkm.sum(t, k, '*') <= 1, name='(1.3)')

#(1.4) - One flight within a certain time window corresponding to the respevyive flight & rest time 
# TODO: special attention
for k in range(num_resources):
    for t in range(num_timeslots - (T_k[k] + R_k[k] - 1)):
    #for t in range(num_timeslots - (T_k[k] + R_k[k])):
        # Split it up - pretty sure the first inline ranging had a bug

        lhs = 0
        #for t_prime in range(T_k[k] + R_k[k]-1):
        for t_prime in range(T_k[k] + R_k[k]):  # TODO: do we need a minus one here? We presumably want to include the -1 entry in the sum, and without it we wouldn's
            lhs += c_tkm.sum(t+t_prime, k, '*')
        model.addConstr(lhs <= 1, name='(1.4)')

#(1.5) Ensure that the resource does not initiate a flight in a posterior time slot
for k in range(num_resources):
    for m in range(num_areas):
        for t in range(num_timeslots-(T_k[k]-1), num_timeslots):  # TODO: where does the -2 come from? Also, do we need the -1?
            model.addConstr(c_tkm[t, k, m] == 0, name='(1.5)')

#---------------------------------PHASE OF FLIGHT DEFINITIONS-----------------------------------

#(1.6) - Arriving to the assigned area
for k in range(num_resources):  # NOTE: U_km must be smaller than T_k (travel time smaller than max range, else the thing breaks)
    for m in range(num_areas):
        #for t in range(num_timeslots - (T_k[k])):
        for t in range(num_timeslots - (T_k[k] - 1)):
            model.addConstr(e_tkm[t + int(U_km[k, m]), k, m] >= c_tkm[t, k, m], name='(1.6)')

#(1.7) - Leaving from assigned area
for k in range(num_resources):
    for m in range(num_areas):
        #for t in range(num_timeslots - (T_k[k])):
        for t in range(num_timeslots - (T_k[k] - 1)):
            model.addConstr(e_tkm[t + (T_k[k] - 1) - int(U_km[k, m]), k, m] >= c_tkm[t, k, m], name='(1.7)')

#(1.8) - Fighting fire in the assigned area
for k in range(num_resources):
    for m in range(num_areas):
        #for t in range(num_timeslots - (T_k[k])):
        for t in range(num_timeslots - (T_k[k] - 1)):
            #for t_prime in range(int(U_km[k, m]+1), T_k[k] - int(U_km[k, m]) - 2): # -2 -> -1; we need to include the last 'index' of the t' range; the 0-indexing should be taken care of by our base index, iteration over num_timeslots
            for t_prime in range(int(U_km[k, m]) + 1, T_k[k] - int(U_km[k, m]) -1):
                model.addConstr(x_tkm[t + t_prime, k, m] >= c_tkm[t, k, m], name='(1.8)')

#(1.9) - Ensure firefighting flight is initiated during the timeslots when the resource is arriving
for k in range(num_resources):
    for m in range(num_areas):
        model.addConstr(2 * c_tkm.sum('*', k, m) == e_tkm.sum('*', k, m), name='(1.9)')

#(1.10) - Ensure firefighting flight is initiated during the timeslots when the resource is leaving
for k in range(num_resources):
    for m in range(num_areas):
        model.addConstr((T_k[k] - 2*U_km[k, m] - 2)*c_tkm.sum('*', k, m) == x_tkm.sum('*', k, m), name='(1.10)')

#(1.11)
for t in range(num_timeslots):
    for k in range(num_resources):
        for m in range(num_areas):
            model.addConstr(x_tkm[t, k, m] + e_tkm[t, k, m] <= 1, name='(1.11)')

#---------------------------------RESOURCE CONSTRAINTS-----------------------------------
# (1.12) - Define time slots allocated for each resource to operate in a flight tied to a specific area. Arrival
# TODO: special attention
for k in range(num_resources):
    for m in range(num_areas):
        #for t in range(num_timeslots - (T_k[k])):
        for t in range(num_timeslots - (T_k[k] - 1)):
            for t_prime in range(T_k[k]): # Want to include T_k[k] - 1 in this range -> need to go one higher in the range operation (and unlike t, this is explicitly 0-indexed in the paper)
            #for t_prime in range(0, T_k[k]-1):
                model.addConstr(s_tkm[t + t_prime, k, m] >= c_tkm[t, k, m], name='(1.12)')

#(1.13) - Define time slots allocated for each resource to operate in a flight tied to a specific area. Leaving
for k in range(num_resources):
    for m in range(num_areas):
        model.addConstr(T_k[k] * c_tkm.sum('*', k, m) == s_tkm.sum('*', k, m), name='(1.13)')

#(1.14) - Resource can only be activate when it is available in a certain timeslot
for t in range(num_timeslots):
    for k in range(num_resources):
        model.addConstr(s_tkm.sum(t, k, '*') <= int(A_tk[t, k]), name='(1.14)')

#(1.15) - Limit of physical presence of the crew
# TODO: special attention
for k in range(num_resources):
    for t in range(num_timeslots - P_k[k]):
        lhs = big_M * (1 - s_tkm.sum(t, k, '*'))

        rhs = 0
        for t_prime in range(P_k[k], num_timeslots-t+1): # TODO: check 0-indexing
        #t_prime = P_k[k]
            rhs += s_tkm.sum(t+t_prime, k, '*')

        model.addConstr(lhs >= rhs, name='(1.15)')

#------------------------------ AREA OF OPERATION  -----------------------------------
#(1.16) - Allows limit for the number of resources flying in the same area at the same time
for m in range(num_areas):
    for t in range(num_timeslots):
        # NOTE: Split the sum into two; linear operator so should be the same
        model.addConstr((x_tkm.sum(t, '*', m) + e_tkm.sum(t, '*', m)) <= S_m[m], name='(1.16)')

#(1.17) - In a single area only resources of the same type are present at the same time
for m in range (num_areas):
    for t in range(num_timeslots):
        for q in range(num_resource_types):
            for k in range(num_resources):
                model.addConstr(v_tqm[t, q, m] >= V_qk[q][k]*(x_tkm[t, k, m] + e_tkm[t, k, m]), name='(1.17)')

#(1.18) - In a single area only resources of the same type are present at the same time
for m in range(num_areas):
    for t in range(num_timeslots):
        model.addConstr(v_tqm.sum(t, '*', m) <= 1, name='(1.18)')

#(1.19) - guarantees that resources are allocated to areas where they are compatible with their type of operation.
for m in range(num_areas):
    for q in range(num_resource_types):
        for t in range(num_timeslots):
            model.addConstr(v_tqm[t, q, m] >= int(B_qm[q, m]), name='(1.19)')

#---------------------------------WATER TARGET COMPLETION-----------------------------------
# (1.20) - To calculate the water target completion in each time slot and area of operation
for t in range(num_timeslots):
    for m in range(num_areas):
        lhs = TCplus_tm[t, m] - TCminus_tm[t, m]
        rhs = 0
        for k in range(num_resources):
            rhs += C_k[k] * E_tkm[t, k, m] * e_tkm[t, k, m]
            rhs += C_k[k] * D_tkm[t, k, m] * x_tkm[t, k, m]
        rhs -= W_tm[t, m]

        model.addConstr(lhs == rhs, name='(1.20)')

#(1.21) - Overall minimal water target completion in all time slots and areas of operation
for t in range(num_timeslots):
    for m in range(num_areas):
        model.addConstr((TCplus_min - TCminus_min) <= (TCplus_tm[t, m] - TCminus_tm[t, m]), name='(1.21)')

#(1.22) - Total water output
WO_eqn = 0
for k in range(num_resources):
    for m in range(num_areas):
        for t in range(num_timeslots):
            WO_eqn += C_k[k]*((E_tkm[t, k, m] * e_tkm[t, k, m]) + (D_tkm[t, k, m] * x_tkm[t, k, m]))
model.addConstr(WO == WO_eqn, name='(1.22)')

#---------------------------------OBJECTIVE-----------------------------------
# Define the objective functions
f1 = 0
for t in range(num_timeslots):
    for m in range(num_areas):
        f1 += -I_m[m] * TCminus_tm[t, m]
f2 = TCplus_min - TCminus_min
f3  = WO

# Lexicographical priority: f1 > f2 > f3
objective = [f1, f2, f3]  # List of objectives


model.setObjectiveN(-a1*objective[0], 0, priority=3)  # f1 has the highest priority
model.setObjectiveN(a2*objective[1], 1, priority=2)  # f2 has the second highest priority
model.setObjectiveN(a3*objective[2], 2, priority=1)  # f3 has the lowest priority

#model.setObjective(sum(weighted_objective), GRB.MAXIMIZE)

model.write('test.lp')
model.setParam(GRB.Param.Seed, 42)  # Set random seed for reproducibility
model.setParam(GRB.Param.TimeLimit, 60)  # Set time limit to 3600 seconds (1 hour)
model.optimize()

if model.Status == GRB.OPTIMAL:
    print("Optimal solution found within time limit.")
elif model.Status == GRB.TIME_LIMIT:
    print("Time limit reached. Best solution found:")
    print(f"Objective value: {model.ObjVal}")
else:
    print(f"Optimization stopped with status {model.Status}.")

print(f"Number of constraints: {model.numConstrs}")
print(f"Number of variables: {model.numVars}")


# for v in model.getVars():
#     print('%s %g' % (v.VarName, v.X))

# print('Obj: %g' % model.ObjVal)

# for t in range(num_timeslots):
#     print(f"t={t}, x={x_tkm[t, 0, 1].X}, e={e_tkm[t, 0, 1].X}, s={s_tkm[t, 0, 1].X}")