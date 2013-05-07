#!/usr/bin/env python

import subprocess

# First, test different integrators!
integrators = ["ForwardEuler","RK4","SPRK6Ruth","SPRK6SuzukiTrotter"]
systems = ["Newtonian","BarnesHut","SoftNewtonian","SoftBarnesHut"]

for i in range(4,11):
	dt = 10**(i/2.)
	dur = 1e4*dt
	for i,integrator in enumerate(integrators):
		for j,system in enumerate(systems[:2]):
			args = ["GravSim","0",str(dt),str(dur),"1",str(i),str(j)]
			with open("data/gdata-%s-%s-%s-%s-%s-%s.txt"%tuple(args[1:]),'w') as outfile:
				p = subprocess.Popen(args,stdout=outfile)
				p.wait()
				print p.returncode
				
dt = 1000
dur = 1e4*dt
for nexp in range(2,9):
	n=int(10**(i/2.))
	for i,integrator in enumerate(integrators):
		for j,system in enumerate(systems[:2]):
			args = ["GravSim",str(n),str(dt),str(dur),"1",str(i),str(j)]
			with open("data/gdata-%s-%s-%s-%s-%s-%s.txt"%tuple(args[1:]),'w') as outfile:
				p = subprocess.Popen(args,stdout=outfile)
				p.wait()
				print p.returncode