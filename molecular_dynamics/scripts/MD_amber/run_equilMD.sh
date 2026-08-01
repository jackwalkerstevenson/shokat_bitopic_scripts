#!/bin/bash

source /home/khuang/amber24/amber.sh

module load sali-libraries 
module load gcc/10.2.1
module load cuda/11.5.0          # Load CUDA for GPU acceleration (pmemd.cuda)
module load mpi/openmpi-x86_64   # Load OpenMPI for CPU parallelization (pmemd.MPI)

export OMP_NUM_THREADS=20        # Allocate 20 OpenMP threads per process
export OMP_PROC_BIND=true       # Bind threads to physical cores to maximize cache locality

echo 'Remember to update the atom selections for restraints!'
endcom=400                       # Total residue index marking the boundary of your solute (Protein+DNA)
#read -p "System size (no water/ions): " endcom


# --- 01_min.in: Initial Minimization ---
# Minimizes only water and ions while keeping the complex tightly restrained
cat>01_min.in<<EOF
Initial min for solvent+ion
 &cntrl
  imin = 1, ntr = 1, restraint_wt = 10,
  maxcyc = 5000,
  ncyc = 2500,
 /
Hold complex with weak restraints
100.0
RES 1 ${endcom}
END
END
EOF

# --- 02_min.in: System-Wide Minimization ---
# Full minimization of all atoms without restrictions before heating begins
cat>02_min.in<<EOF
All min before heat 
 &cntrl
  imin   = 1,
  maxcyc = 20000,
  ncyc   = 10000,
  ntb    = 1,
  ntr    = 0,
  ntpr   = 500,
  cut    = 12.0
 /
EOF

# --- 03_heat.in: Thermalization Stage ---
# NVT Heating up to 300K over 100ps using weak classical positional restraints
cat>03_heat.in<<EOF
 Heat before MD with weak restraints for 100ps
&cntrl
  imin   = 0,
  irest  = 0,
  ntx    = 1,
  ntb    = 1,
  cut    = 12.0,
  ntr    = 1,
  ntc    = 2,
  ntf    = 2,
  tempi  = 0.0,
  temp0  = 300.0,
  ntt    = 3,
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 50000, dt = 0.002
  ntpr = 500, ntwx = 500, ntwr = 1000
 /
Hold protein+DNA with weak restraints
50.0
RES 1 ${endcom}
END
END
EOF

# --- 04_equil1.in: NPT Equilibration Step 1 ---
# Step 1 of NPT equilibration: Restraint weight set to 75.0 via group input format
cat>04_equil1.in<<EOF
Equilibrium 5ns MD with pressure
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntb = 2, 
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 10.0, 
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 500000, dt = 0.002,		
  ntpr = 5000, ntwx = 5000, ntwr = 10000
  ntr = 1,
 /
Hold protein+DNA with weak restraints
75.0
RES 1 ${endcom}
END
END
EOF

# --- 04_equil2.in: NPT Equilibration Step 2 ---
# Step 2 of NPT equilibration: Restraint weight dropped to 50.0 via group input format
cat>04_equil2.in<<EOF
Equilibrium 5ns MD with pressure
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntb = 2, 
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 10.0, 
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 500000, dt = 0.002,		
  ntpr = 5000, ntwx = 5000, ntwr = 10000
  ntr = 1,
 /
Hold protein+DNA with weak restraints
50.0
RES 1 ${endcom}
END
END
EOF

# --- 04_equil3.in: NPT Equilibration Step 3 ---
# Step 3 of NPT equilibration: Restraint weight dropped to 25.0 via group input format
cat>04_equil3.in<<EOF
Equilibrium 5ns MD with pressure
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntb = 2, 
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 10.0, 
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 500000, dt = 0.002,		
  ntpr = 5000, ntwx = 5000, ntwr = 10000
  ntr = 1,
 /
Hold protein+DNA with weak restraints
25.0
RES 1 ${endcom}
END
END
EOF

# --- 04_equil4.in: NPT Equilibration Step 4 ---
# Swaps to standard modern Amber namelist variables: restraint_wt and restraintmask on all solute atoms
cat>04_equil4.in<<EOF
Equilibrium 5ns MD with pressure
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntb = 2, 
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 10.0, 
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 500000, dt = 0.002,		
  ntpr = 5000, ntwx = 5000, ntwr = 10000
  ntr = 1,
  restraint_wt=50.0, 
  restraintmask=':1-${endcom}',
 /
EOF

# --- 04_equil5.in: NPT Equilibration Step 5 ---
# Restraint weight remains 50.0, but masks out all Hydrogens (&!@H) from positional restraints
cat>04_equil5.in<<EOF
Equilibrium 5ns MD with pressure
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntb = 2, 
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 10.0, 
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 500000, dt = 0.002,		
  ntpr = 5000, ntwx = 5000, ntwr = 10000
  restraint_wt=50.0, 
  restraintmask=':1-${endcom}&!@H'
 /
EOF

# --- 04_equil6.in: NPT Equilibration Step 6 ---
# Restraint weight dropped to 25.0, focusing exclusively on the main backbone heavy atoms (@CA,N,O,C)
cat>04_equil6.in<<EOF
Equilibrium 5ns MD with pressure
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntb = 2, 
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 10.0, 
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 500000, dt = 0.002,		
  ntpr = 5000, ntwx = 5000, ntwr = 10000
  ntr = 1,
  restraint_wt=25.0, 
  restraintmask=':1-${endcom}@CA,N,O,C'
 /
EOF

# --- 04_equil7.in: NPT Equilibration Step 7 ---
# Tightly focused restraint mask isolating only Alpha Carbons (@CA) at a weight of 25.0
cat>04_equil7.in<<EOF
Equilibrium 5ns MD with pressure
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntb = 2, 
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 10.0, 
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  nstlim = 500000, dt = 0.002,		
  ntpr = 5000, ntwx = 5000, ntwr = 10000
  ntr = 1,
  restraint_wt=25.0, 
  restraintmask=':1-${endcom}@CA'
 /
EOF

# --- 04_equil8.in: Production MD Run ---
# 5ns unrestrained production MD utilizing the Monte Carlo barostat (barostat=2) and coordinate wrapping (iwrap=1)
cat>04_equil8.in<<EOF
production for 5 ns
 &cntrl
  imin = 0, 
  irest = 1, 
  ntx = 7,
  ntxo= 2,
  barostat = 2,
  pres0 = 1.0, 
  ntp = 1,
  taup = 2.0,
  cut = 9.0, 
  ntr = 0,
  ntc = 2, 
  ntf = 2,
  tempi = 300.0, 
  temp0 = 300.0,
  ntt = 3, 
  gamma_ln = 1.0,
  ig=-1,
  iwrap=1,
  nstlim = 2500000, dt = 0.002,
  ntpr = 75000, ntwx = 2000, ntwr = 150000,
 /
EOF


# --- INITIAL CPU PROCESSING ---
# High-instability risk steps are commented out below to default directly to GPU engine instead.
# mpirun -np 20 pmemd.MPI -O -i 01_min.in -p *.prmtop -c *.inpcrd -o 01_sol_min.out -r 01_sol_min.rst -ref *.inpcrd
wait

# mpirun -np 20 pmemd.MPI  -O -i 02_min.in -p *.prmtop -c 01_sol_min.rst -o 02_all_min.out -r 02_all_min.rst -ref 01_sol_min.rst
wait

# pmemd.cuda  -O -i 03_heat.in -p *.prmtop -c 02_all_min.rst -o 03_heat.out -r 03_heat.rst  -ref 02_all_min.rst -x 03_heat.nc -inf 03_heat.mdinfo -l 03_heat.log


# --- GPU PRODUCTION AND EQUILIBRATION SEQUENCE ---
# Steps 1 through 8 run progressively, passing the resulting restart file (.rst) to the next block's input (-c)

# [Commented Out Step 1]
# pmemd.cuda -O -i 04_equil1.in -p *.prmtop -c 03_heat.rst -o 04_equil1.out -r 04_equil1.rst -x 04_equil1.nc -inf 04_equil1.mdinfo -l 04_equil1.log -ref 03_heat.rst 
wait

# NPT Equilibration 2
pmemd.cuda -O -i 04_equil2.in -p *.prmtop -c 04_equil1.rst -o 04_equil2.out -r 04_equil2.rst -x 04_equil2.nc -inf 04_equil2.mdinfo -l 04_equil2.log -ref 04_equil1.rst
wait

# NPT Equilibration 3
pmemd.cuda -O -i 04_equil3.in -p *.prmtop -c 04_equil2.rst -o 04_equil3.out -r 04_equil3.rst -x 04_equil3.nc -inf 04_equil3.mdinfo -l 04_equil3.log -ref 04_equil2.rst 
wait

# NPT Equilibration 4 (Transitioning to modern namelist masks)
pmemd.cuda -O -i 04_equil4.in -p *.prmtop -c 04_equil3.rst -o 04_equil4.out -r 04_equil4.rst -x 04_equil4.nc -inf 04_equil4.mdinfo -l 04_equil4.log -ref 04_equil3.rst
wait

# NPT Equilibration 5 (Excluding Hydrogens from restraints)
pmemd.cuda -O -i 04_equil5.in -p *.prmtop -c 04_equil4.rst -o 04_equil5.out -r 04_equil5.rst -x 04_equil5.nc -inf 04_equil5.mdinfo -l 04_equil5.log -ref 04_equil4.rst
wait

# NPT Equilibration 6 (Restraining core backbone atoms only)
pmemd.cuda -O -i 04_equil6.in -p *.prmtop -c 04_equil5.rst -o 04_equil6.out -r 04_equil6.rst -x 04_equil6.nc -inf 04_equil6.mdinfo -l 04_equil6.log -ref 04_equil5.rst
wait

# NPT Equilibration 7 (Restraining alpha carbons only)
pmemd.cuda -O -i 04_equil7.in -p *.prmtop -c 04_equil6.rst -o 04_equil7.out -r 04_equil7.rst -x 04_equil7.nc -inf 04_equil7.mdinfo -l 04_equil7.log -ref 04_equil6.rst
wait

# 5ns Production MD Simulation (Completely unrestrained, active barostat & coordinate wrapping)
pmemd.cuda -O -i 04_equil8.in -p *.prmtop -c 04_equil7.rst -o 04_equil8.out -r 04_equil8.rst -x 04_equil8.nc -inf 04_equil8.mdinfo -l 04_equil8.log -ref 04_equil7.rst
wait
