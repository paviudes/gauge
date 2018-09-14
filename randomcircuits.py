# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 11:58:37 2018

@author: dript
"""
from sys import stdout
import numpy as np
from quantumsim.sparsedm import SparseDM
# import surface17 as surf17
import json
import time
import sys
import quantumsim.circuit as qc
#import math
#import matplotlib.pyplot as plt
# from random import *
#import pycuda

#time0 = time.time()
#error = 0.01
#cnot_length = 5
#trial_length = 10


#average_average_fidelity = 0
#average_optimal_fidelity = 0


def ConstructCircuit(qcirc, index, timing, error, s, v):
    # construct a quantum circuit and add it to the quantum circuit object
    # these are composed of the single qubit X,Y,Z rotation gates
    encoding = [["X", [0]], ["X", [1]], ["X", [2]],
    ["Y", [0]], ["Y", [1]], ["Y", [2]],
    ["S", [0]], ["S", [1]], ["S", [2]],
    ["H", [0]], ["H", [1]], ["H", [2]],
    ["CNOT", [0, 1]], ["CNOT", [0, 2]], 
    ["CNOT", [1, 0]], ["CNOT", [1, 2]],
    ["CNOT", [2, 0]], ["CNOT", [2, 1]]]
    (gate, support) = encoding[index]
    #print("Gate index")
    #print(index)
    #print("Gate")
    #print(gate)
    #print("Asking for the %s at location %s." % (gate, np.array_str(np.array(support))))
    #print("timing")
    #print(timing)
    #if (index >= 12):
    #print("CNOT on %d, %d with s = %d, v = %d" % (support[0], support[1], s, v))
    if (gate == "CNOT"):
    	#print("CNOT on %d, %d with s = %d, v = %d" % (support[0], support[1], s, v))
    	eval("qcirc.add_gate(qc.RotateY(\'%d\', angle = (1 + error) * v * np.pi/2, time = timing))" % (support[0]))
    	eval("qcirc.add_gate(qc.XX(\'%d\', \'%d\', chi = (1 + error) * (1 + error) * s * np.pi/4, time = timing + 0.1))" % (support[0], support[1]))
    	eval("qcirc.add_gate(qc.RotateX(\'%d\', angle = -s * (1 + error) * np.pi/2, time = timing + 0.2))" % (support[0]))
    	eval("qcirc.add_gate(qc.RotateX(\'%d\', angle = -s * v * (1 + error) * np.pi/2, time = timing + 0.2))" % (support[1]))
    	eval("qcirc.add_gate(qc.RotateY(\'%d\', angle = -v * (1 + error) * np.pi/2, time = timing + 0.3))" % (support[0]))
    	# print("CNOT on %d, %d with s = %d, v = %d" % (support[0], support[1], s, v))
    else:
    	angles = {"X":["np.pi/2"],"Y":["np.pi/2"],"S":["np.pi/4"], "H":["-np.pi/2", "np.pi"]}
    	rgates = {"X":["RotateX"], "Y":["RotateY"], "S":["RotateZ"], "H":["RotateY", "RotateX"]}
    	for i in range(len(rgates[gate])):
    		eval("qcirc.add_gate(qc.%s(\'%d\', angle = -(1 + error) * %s, time = timing))" % (rgates[gate][i], support[0], angles[gate][i]))
    return None


def RandomCliffGen(depth, ncnots):
	# Construct a random Clifford circuit of a given depth with a specific number of CNOTs.
	# Generate random gates at every step, stop when the number of CNOTs has reached the required amount.
	ncliff = 18
	circ = (-1) * np.ones(depth, dtype = np.int)
	count = 0
	for i in range(depth):
		circ[i] = np.random.randint(0, high = ncliff)
		# if gate is < 12, it is a single qubit gate. Else it is a CNOT.
		count = count + np.int(circ[i] >= 12)
		if (count >= ncnots):
			break
	return circ


def BellStateInitialize(qcirc):
	# Circuit for preparing the bell state.
	qcirc.add_qubit('0')
	qcirc.add_qubit('1')
	qcirc.add_qubit('2')
	qcirc.add_qubit('0p')
	qcirc.add_qubit('1p')
	qcirc.add_qubit('2p')

	qcirc.add_gate(qc.Hadamard('0',time=0))
	qcirc.add_gate(qc.Hadamard('1',time=0))
	qcirc.add_gate(qc.Hadamard('2',time=0))

	qcirc.add_gate(qc.CNOT('0p','0',time=0.5))
	qcirc.add_gate(qc.CNOT('1p','1',time=0.5))
	qcirc.add_gate(qc.CNOT('2p','2',time=0.5))
	return None


def RandomCliffordCircuit(qcirc, circ, gauge, noise):
	# Construct a random Clifford circuit using the quantumsim object
	# Preparing the Choi state of the circuit -- applying the circuit elements to one half of a bell state.
	BellStateInitialize(qcirc)
	gidx = 0
	for i in range(circ.shape[0]):
		if (circ[i] < 0):
			break
		if (circ[i] >= 12):
			# the gauges are only for the CNOT gates
			ConstructCircuit(qcirc, circ[i], i + 1, noise, (-1)**gauge[gidx], (-1)**gauge[gidx + 1])
			gidx = gidx + 2
		else:
			ConstructCircuit(qcirc, circ[i], i + 1, noise, 0, 0)
	return None

def SimulateRandomClifford(circ, gauge, noise):
	# Compute the entanglement fidelity between the noisy circuit and the perfect circuit.
	# 1. Compute the Choi state for the perfect circuit.
	# 2. Compute the Choi state for the imperfect circuit.
	# 3. Take the inner product.
	
	start = time.time()

	perfect = qc.Circuit(title='Perfect Circuit')
	RandomCliffordCircuit(perfect, circ, gauge, 0)
	pdm = SparseDM(perfect.get_qubit_names())
	perfect.apply_to(pdm)

	#print("Perfect Choi state")
	#print(sdmc.full_dm.dm.ravel())

	noisy = qc.Circuit(title='Noisy Circuit')
	RandomCliffordCircuit(noisy, circ, gauge, noise)
	ndm = SparseDM(noisy.get_qubit_names())
	noisy.apply_to(ndm)

	fidelity = np.dot(pdm.full_dm.dm.ravel(), ndm.full_dm.dm.ravel())

	# print("fidelity = %g"% (fidelity))

	end = time.time()

	# print("gauge setting\n%s\ndone in %d seconds."% (np.array_str(gauge), end - start))	

	return fidelity


def FindOptimalGauge(circ, noise):
	# Find an optimal gauge for the quantum circuit
	# For each gauge, construct the quantum circuit and compute its fidelity.
	nparams = 2 * np.count_nonzero(circ >= 12) # + np.count_nonzero(np.logical_and(circ < 12, circ > -1))
	ndegfreedom = np.power(2, nparams, dtype = np.int)
	gauges = np.zeros((ndegfreedom, nparams + 1), dtype = np.float)
	for i in range(ndegfreedom):
		gauges[i, :-1] = np.array(list(map(np.int8, np.binary_repr(i, width = nparams))), dtype = np.float)
		gauges[i, -1] = SimulateRandomClifford(circ, gauges[i, :-1].astype(np.int8), noise)
		#stdout.write("\r...%d%% done" % ((i+1)/float(ndegfreedom) * 100))
	#print("")
	optimal = gauges[np.argmax(gauges[:, -1]), :]
	print("Fidelities\n%s" % (np.array_str(gauges[:, -1])))
	return optimal


def SearchCNOTPairs(noise):
	# Construct all configurations for CNOT on two qubits.
	# 1. 15, 12
	# 2. 15, 14
	# 3. 12, 13
	# 4. 12, 17
	# 5. 12, 12
	# 6. 12, 14
	# 7. 14, 12
	# 8. 14, 14
	start = time.time()
	circs = np.array([[15, 12], [15, 14], [12, 13], [12, 17], [12, 12], [12, 14], [14, 12], [14, 14]], dtype = np.int8)
	for i in range(circs.shape[0]):
		print("Circuit element: %s"% (np.array_str(circs[i, :])))
		optimal = FindOptimalGauge(circs[i, :], noise)
		print("Optimal gauge\n%s" % (np.array_str(optimal)))
		print("Total time: %d seconds."% (time.time() - start))
	return circs

def SearchRandomCircuits(size, ncnots, noise, ncircuits):
	# Find optimal gauge values in a set of random circuits with a given number of CNOTs
	start = time.time()
	for i in range(ncircuits):
		circ = RandomCliffGen(depth, ncnots)
		print("Random circuit\n%s" % (np.array_str(circ)))
		optimal = FindOptimalGauge(circ, noise)
		print("Optimal gauge\n%s" % (np.array_str(optimal)))
		print("...%d%% done, %d seconds elapsed." % ((i+1)/float(ncircuits) * 100, time.time() - start))
	#print("")
	print("Total time: %d seconds."% (time.time() - start))
	return None

if __name__ == '__main__':
	## Find the optimal gauge for a random circuit
	size = 20
	ncnots = 4
	noise = 0.05
	ncircuits = 1000
	# SearchRandomCircuits(size, ncnots, noise, ncircuits)
	## Construct basic CNOT elements
	circs = SearchCNOTPairs(noise)


# for trial_number in range(trial_length):

#     print("Current trial number: " + str(trial_number), flush = True)

#     cnot_count = 0
#     circuit_template = np.zeros(cnot_length, dtype = np.int)

#     while(cnot_count < cnot_length):
#         current_gate = randint(0,4)
#         if(current_gate <= 3):
#             gate_placement = randint(0,2)
#         else:
#             gate_placement = randint(0,5)
#         circuit_template.append(current_gate*3 + gate_placement)
#         if(current_gate == 4):
#             cnot_count += 1

#     eps = 0

#     rand_perfect = qc.Circuit(title='3-Random Circuit')

#     rand_perfect.add_qubit('0')
#     rand_perfect.add_qubit('1')
#     rand_perfect.add_qubit('2')
#     rand_perfect.add_qubit('0p')
#     rand_perfect.add_qubit('1p')
#     rand_perfect.add_qubit('2p')

#     rand_perfect.add_gate(qc.Hadamard('0',time=0))
#     rand_perfect.add_gate(qc.Hadamard('1',time=0))
#     rand_perfect.add_gate(qc.Hadamard('2',time=0))

#     rand_perfect.add_gate(qc.CNOT('0p','0',time=0.5))
#     rand_perfect.add_gate(qc.CNOT('1p','1',time=0.5))
#     rand_perfect.add_gate(qc.CNOT('2p','2',time=0.5))

#     time = 0
#     current_param_index = 0

#     for index in circuit_template:

#         time += 1

#         if(index == 0):

#             rand_perfect.add_gate(qc.RotateX('0',angle=-(1+eps)*np.pi/2, time=time))

#         elif(index == 1):

#             rand_perfect.add_gate(qc.RotateX('1',angle=-(1+eps)*np.pi/2, time=time))

#         elif(index == 2):

#             rand_perfect.add_gate(qc.RotateX('2',angle=-(1+eps)*np.pi/2, time=time))

#         elif(index == 3):

#             rand_perfect.add_gate(qc.RotateY('0',angle=-(1+eps)*np.pi/2, time=time))

#         elif(index == 4):

#             rand_perfect.add_gate(qc.RotateY('1',angle=-(1+eps)*np.pi/2, time=time))

#         elif(index == 5):

#             rand_perfect.add_gate(qc.RotateY('2',angle=-(1+eps)*np.pi/2, time=time))

#         elif(index == 6):

#             rand_perfect.add_gate(qc.RotateZ('0',angle=-(1+eps)*np.pi/4, time=time))

#         elif(index == 7):

#             rand_perfect.add_gate(qc.RotateZ('1',angle=-(1+eps)*np.pi/4, time=time))

#         elif(index == 8):

#             rand_perfect.add_gate(qc.RotateZ('2',angle=-(1+eps)*np.pi/4, time=time))

#         elif(index == 9):

#             qb='0';ax=0;t=time;
#             rand_perfect.add_gate(qc.RotateArb(qb,ax,1-ax,0,angle=-(1+eps)*np.pi*(1+ax)/2,time=t))
#             rand_perfect.add_gate(qc.RotateArb(qb,1-ax,ax,0,angle=(1+eps)*np.pi*(2-ax)/2,time=t+0.1))

#         elif(index == 10):

#             qb='1';ax=0;t=time;
#             rand_perfect.add_gate(qc.RotateArb(qb,ax,1-ax,0,angle=-(1+eps)*np.pi*(1+ax)/2,time=t))
#             rand_perfect.add_gate(qc.RotateArb(qb,1-ax,ax,0,angle=(1+eps)*np.pi*(2-ax)/2,time=t+0.1))

#         elif(index == 11):

#             qb='2';ax=0;t=time;
#             rand_perfect.add_gate(qc.RotateArb(qb,ax,1-ax,0,angle=-(1+eps)*np.pi*(1+ax)/2,time=t))
#             rand_perfect.add_gate(qc.RotateArb(qb,1-ax,ax,0,angle=(1+eps)*np.pi*(2-ax)/2,time=t+0.1))

#         elif(index == 12):

#             con='0';tar='1';s=2*(1-0.5);v=2*(1-0.5);t=time;
#             rand_perfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#             rand_perfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#             rand_perfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#             current_param_index += 2

#         elif(index == 13):

#             con='0';tar='2';s=2*(1-0.5);v=2*(1-0.5);t=time;
#             rand_perfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#             rand_perfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#             rand_perfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#             current_param_index += 2

#         elif(index == 14):

#             con='1';tar='2';s=2*(1-0.5);v=2*(1-0.5);t=time;
#             rand_perfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#             rand_perfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#             rand_perfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#             current_param_index += 2

#         elif(index == 15):

#             con='1';tar='0';s=2*(1-0.5);v=2*(1-0.5);t=time;
#             rand_perfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#             rand_perfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#             rand_perfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#             current_param_index += 2

#         elif(index == 16):

#             con='2';tar='0';s=2*(1-0.5);v=2*(1-0.5);t=time;
#             rand_perfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#             rand_perfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#             rand_perfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#             current_param_index += 2

#         elif(index == 17):

#             con='2';tar='1';s=2*(1-0.5);v=2*(1-0.5);t=time;
#             rand_perfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#             rand_perfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#             rand_perfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#             rand_perfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#             current_param_index += 2

#     sdmc = SparseDM(rand_perfect.get_qubit_names())
#     rand_perfect.apply_to(sdmc)

#     bp=0;
#     bf=0;

#     fidelities = np.zeros(2**(2*cnot_count), dtype = np.longdouble)

#     eps = error
#     params = np.zeros(2**(2*cnot_count), dtype = np.int8)
#     for i in range(2**(2*cnot_count)):

#         # params = [float(j) for j in list(str(bin(i))[2:].zfill((2*cnot_count)))]
# 	params[i, :] = np.array(map(np.int8, np.binary_repr(i, width = (2 * cnot_count))), dtype = np.int8)

#         rand_imperfect = qc.Circuit(title='3-Random Circuit')

#         rand_imperfect.add_qubit('0')
#         rand_imperfect.add_qubit('1')
#         rand_imperfect.add_qubit('2')
#         rand_imperfect.add_qubit('0p')
#         rand_imperfect.add_qubit('1p')
#         rand_imperfect.add_qubit('2p')

#         rand_imperfect.add_gate(qc.Hadamard('0',time=0))
#         rand_imperfect.add_gate(qc.Hadamard('1',time=0))
#         rand_imperfect.add_gate(qc.Hadamard('2',time=0))

#         rand_imperfect.add_gate(qc.CNOT('0p','0',time=0.5))
#         rand_imperfect.add_gate(qc.CNOT('1p','1',time=0.5))
#         rand_imperfect.add_gate(qc.CNOT('2p','2',time=0.5))

#         time = 0;
#         current_param_index = 0;

#         for index in circuit_template:

#             time += 1;

#             if(index == 0):

#                 rand_imperfect.add_gate(qc.RotateX('0',angle=-(1+eps)*np.pi/2, time=time))

#             elif(index == 1):

#                 rand_imperfect.add_gate(qc.RotateX('1',angle=-(1+eps)*np.pi/2, time=time))

#             elif(index == 2):

#                 rand_imperfect.add_gate(qc.RotateX('2',angle=-(1+eps)*np.pi/2, time=time))

#             elif(index == 3):

#                 rand_imperfect.add_gate(qc.RotateY('0',angle=-(1+eps)*np.pi/2, time=time))

#             elif(index == 4):

#                 rand_imperfect.add_gate(qc.RotateY('1',angle=-(1+eps)*np.pi/2, time=time))

#             elif(index == 5):

#                 rand_imperfect.add_gate(qc.RotateY('2',angle=-(1+eps)*np.pi/2, time=time))

#             elif(index == 6):

#                 rand_imperfect.add_gate(qc.RotateZ('0',angle=-(1+eps)*np.pi/4, time=time))

#             elif(index == 7):

#                 rand_imperfect.add_gate(qc.RotateZ('1',angle=-(1+eps)*np.pi/4, time=time))

#             elif(index == 8):

#                 rand_imperfect.add_gate(qc.RotateZ('2',angle=-(1+eps)*np.pi/4, time=time))

#             elif(index == 9):

#                 qb='0';ax=0;t=time;
#                 rand_imperfect.add_gate(qc.RotateArb(qb,ax,1-ax,0,angle=-(1+eps)*np.pi*(1+ax)/2,time=t))
#                 rand_imperfect.add_gate(qc.RotateArb(qb,1-ax,ax,0,angle=(1+eps)*np.pi*(2-ax)/2,time=t+0.1))

#             elif(index == 10):

#                 qb='1';ax=0;t=time;
#                 rand_imperfect.add_gate(qc.RotateArb(qb,ax,1-ax,0,angle=-(1+eps)*np.pi*(1+ax)/2,time=t))
#                 rand_imperfect.add_gate(qc.RotateArb(qb,1-ax,ax,0,angle=(1+eps)*np.pi*(2-ax)/2,time=t+0.1))

#             elif(index == 11):

#                 qb='2';ax=0;t=time;
#                 rand_imperfect.add_gate(qc.RotateArb(qb,ax,1-ax,0,angle=-(1+eps)*np.pi*(1+ax)/2,time=t))
#                 rand_imperfect.add_gate(qc.RotateArb(qb,1-ax,ax,0,angle=(1+eps)*np.pi*(2-ax)/2,time=t+0.1))

#             elif(index == 12):

#                 con='0';tar='1';s=2*(params[current_param_index]-0.5);v=2*(params[current_param_index+1]-0.5);t=time;
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#                 rand_imperfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#                 rand_imperfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#                 current_param_index += 2;

#             elif(index == 13):

#                 con='0';tar='2';s=2*(params[current_param_index]-0.5);v=2*(params[current_param_index+1]-0.5);t=time;
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#                 rand_imperfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#                 rand_imperfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#                 current_param_index += 2;

#             elif(index == 14):

#                 con='1';tar='2';s=2*(params[current_param_index]-0.5);v=2*(params[current_param_index+1]-0.5);t=time;
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#                 rand_imperfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#                 rand_imperfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#                 current_param_index += 2;

#             elif(index == 15):

#                 con='1';tar='0';s=2*(params[current_param_index]-0.5);v=2*(params[current_param_index+1]-0.5);t=time;
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#                 rand_imperfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#                 rand_imperfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#                 current_param_index += 2;

#             elif(index == 16):

#                 con='2';tar='0';s=2*(params[current_param_index]-0.5);v=2*(params[current_param_index+1]-0.5);t=time;
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#                 rand_imperfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#                 rand_imperfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#                 current_param_index += 2;

#             elif(index == 17):

#                 con='2';tar='1';s=2*(params[current_param_index]-0.5);v=2*(params[current_param_index+1]-0.5);t=time;
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=v*(1+eps)*np.pi/2, time=t))
#                 rand_imperfect.add_gate(qc.XX(con, tar, chi=(1+eps)*s*(1+eps)*np.pi/4, time=t+0.1))
#                 rand_imperfect.add_gate(qc.RotateX(con,angle=-s*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateX(tar,angle=-s*v*(1+eps)*np.pi/2, time=t+0.2))
#                 rand_imperfect.add_gate(qc.RotateY(con,angle=-v*(1+eps)*np.pi/2, time=t+0.3))

#                 current_param_index += 2;

#         sdmt = SparseDM(rand_imperfect.get_qubit_names())
#         rand_imperfect.apply_to(sdmt)

#         fid = np.dot(sdmc.full_dm.dm.ravel(), sdmt.full_dm.dm.ravel())

#         fidelities[i] = fid

#         if fid>bf:
#             bp = params
#             bf = fid

#     average_fidelity = np.mean(fidelities)

#     average_average_fidelity += average_fidelity
#     average_optimal_fidelity += bf

# average_optimal_fidelity = average_optimal_fidelity/trial_length
# average_average_fidelity = average_average_fidelity/trial_length

# print()
# print("OUTPUT")
# print()
# print("Average optimal fidelity is")
# print(average_optimal_fidelity)
# print()
# print("Average average fidelity is")
# print(average_average_fidelity)
# print()
# print("Average infidelity suppression is")
# print((1 - average_average_fidelity)/(1 - average_optimal_fidelity))
# print()
# print("END OUTPUT")
# print()







