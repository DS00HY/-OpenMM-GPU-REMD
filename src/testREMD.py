from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from replicaexchange import ReplicaExchange
from mpi4py import MPI
import logging
import sys
import os
import socket

# -- Get MPI --
comm =  MPI.COMM_WORLD
rank =  comm.Get_rank()
n_proc = comm.Get_size()
#logging.basicConfig(level = logging.DEBUG)
begintime = MPI.Wtime()
#output GPU and  host

n_replicas = 4
n_iterations = 100
savelogstep = 500
savepdbstep = 500
ex_interval = 1000
input_path = "/public/home/yinghuang/test/INPUTS_Hexa/"
include_path = "/public/home/yinghuang/test/top"
_OUTPUT_PATH = "out/"

plat= Platform.getPlatformByName('CUDA')
#prop = {'DeviceIndex': str(rank%2), 'Precision': 'double'}
prop = {'Precision': 'double'}
gro = GromacsGroFile(input_path+'AFP_phy_tip3p.gro')
top = GromacsTopFile(input_path+'AFP_phy_tip3p.top', periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=include_path)
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds)
temp = 310.0 * kelvin
integrator = VerletIntegrator(0.002*picoseconds)
tstat = AndersenThermostat(310, 1/picosecond)
system.addForce(tstat)

#integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
#simulation = Simulation(top.topology, system, integrator, plat, prop)
simulation = Simulation(top.topology, system, integrator,plat)
simulation.context.setPositions(gro.positions)
simulation.context.setVelocitiesToTemperature(temp)
#simulation.minimizeEnergy()

simulation.reporters.append(PDBReporter('out/out%i.pdb' % (comm.Get_rank()), savepdbstep))
simulation.reporters.append(StateDataReporter('out/md%i.log' % (comm.Get_rank()), savelogstep, step=True, potentialEnergy=True, kineticEnergy = True, totalEnergy = True, temperature=True, speed = True))
print("now inint REMD")
remd = ReplicaExchange(this_simulation = simulation, n_replica= n_replicas, exchange_interval = ex_interval)
templist = [300.00, 303.19, 306.40, 309.63, 310.00]
replicaslist = ReplicaExchange.create_replica_paramlist(n_replicas = n_replicas, templist = templist)
remd.run(replica_parameters=replicaslist, n_iterations = n_iterations)
if rank == 0:
   time = MPI.Wtime()-begintime
   print("\n #-----  hhhh  all using time is %f  -----#\n "%(time))
print("all end!")

