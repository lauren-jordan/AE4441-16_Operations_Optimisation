The data file contains (approximately):
K => Resources (a set)
Q => resource types (assignment done using the V parameter)
F => fronts(?), resource type needed is assigned by B param
T => number of timeslots
TF/TR/P/N/A => constraints on resource (availability)
U => transit time
C => water capacity
S => total number of resources at each front constraints
D/E => downloads (grabbing water)?
W => water needed in each front
a1/a2/a3 => weight values for the optimization function? (of we do not do lexicographical optimization)
M => big M
PR => priority of each front