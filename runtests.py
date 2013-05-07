#!/usr/bin/env python

import subprocess

# First, test different integrators!
integrators = ["ForwardEuler","RK4","SPRK6Ruth","SPRK6SuzukiTrotter"]
systems = ["Newtonian","BarnesHut","SoftNewtonian","SoftBarnesHut"]

#for i in range(4,11):
#	dt = 10**(i/2.)
#	dur = 1e4*dt
#	for k,integrator in enumerate(integrators):
#		for j,system in enumerate(systems[:2]):
#			if k and (j%2):
#				continue
#			args = ["GravSim","0",str(dt),str(dur),"1",str(k),str(j)]
#			args_desc = ["GravSim","0",str(dt),str(dur),"1",str(k),integrator,str(j),system]
#			with open("data/gdata-%s-%s-%s-%s-%s(%s)-%s(%s).txt"%tuple(args_desc[1:]),'w') as outfile:
#				print " ".join(args_desc),"...",
#				p = subprocess.Popen(args,stdout=outfile)
#				p.wait()
#				print p.returncode
				
dt = 1000
dur = 1e4*dt
for nexp in range(2,9):
	n=int(10**(nexp/2.))
	for k,integrator in enumerate(integrators[:1]):
		for j,system in enumerate(systems):
			if j<2 or (k and (j % 2)):
				continue
			args = ["GravSim",str(n),str(dt),str(dur),"1",str(k),str(j)]
			args_desc = ["GravSim",str(n),str(dt),str(dur),"1",str(k),integrator,str(j),system]
			with open("data/gdata-%s-%s-%s-%s-%s(%s)-%s(%s).txt"%tuple(args_desc[1:]),'w') as outfile:
				print " ".join(args_desc),"...",
				p = subprocess.Popen(args,stdout=outfile)
				p.wait()
				print p.returncode