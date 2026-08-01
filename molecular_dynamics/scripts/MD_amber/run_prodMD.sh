#!/bin/bash

# Source the Amber environment configuration script to expose Amber binaries (like pmemd.cuda)
source /home/khuang/amber24/amber.sh

# Load required cluster libraries, compilers, and accelerators
module load sali-libraries 
module load gcc/10.2.1
module load cuda/11.5.0          # Load CUDA for GPU-accelerated execution
#module load mpi/openmpi-x86_64 # CPU parallelization loading (disabled for pure GPU run)

# GPU Target Directive
export CUDA_VISIBLE_DEVICES=0   # Bind the job exclusively to the first available GPU card (Device 0)

jname=linker_m1                  # Job/system prefix identifier descriptor

# Initial seeding file instruction (commented out here)
#cp 03*.rst 05_prod0.rst 

# Iteration tracking indices for coordinate parsing
start=0                         # Index marking the input coordinate restart file number
step=100                        # Index marking the output coordinate restart file number

# Generates a 500 ns long production MD input script block via a Here Document (EOF)
cat >05_prod.in<<EOF
Production 500 ns MD 
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
  nstlim = 250000000, dt = 0.002,
  ntpr = 75000, ntwx = 2000, ntwr = 150000,
 /
EOF

# Run the long production simulation step using pmemd.cuda on the targeted GPU.
# - 'nohup' keeps the process running even if the terminal logs out or drops connection.
# - '&' forks the entire execution string into a background process thread.
nohup pmemd.cuda -O \
  -i 05_prod.in \
  -p *.prmtop \
  -c 05_prod${start}.rst \
  -o 05_prod${step}.out \
  -r 05_prod${step}.rst \
  -x 05_prod${step}.nc \
  -inf 05_prod${step}.mdinfo \
  -l 05_prod${step}.log \
  -ref 05_prod${start}.rst &

wait                            # Pause execution of this script until the background process completes
rm *.mdinfo                     # Clean up temporary status info files left behind by the MD engine
