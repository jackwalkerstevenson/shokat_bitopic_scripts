#!/usr/bin/env bash

#$ -S /bin/bash       # Run the job using the Bash shell
#$ -cwd             # Execute the job in the current working directory
#$ -pe smp 50        # Request a parallel environment with 50 slots (cores)
#$ -j y              # Merge Standard Error (STDERR) into Standard Output (STDOUT)
#$ -l mem_free=1G    # Request 1 GiB of RAM per slot (50 slots = 50 GiB total)
#$ -l scratch=2G     # Request 2 GiB of local scratch space per node
#$ -l h_rt=11:59:00  # Walltime limit: Hard Runtime of ~12 hours
#$ -r n              # If the job crashes, do NOT automatically restart it

module load amber              # Load the Amber software suite (contains CPPTRAJ)
module load mpi/openmpi-x86_64 # Load OpenMPI for parallel execution

export OMP_NUM_THREADS=50      # Use 50 OpenMP threads per MPI process
export OMP_PROC_BIND=true     # Bind threads to physical cores to maximize cache efficiency

rm -f cluster* *.dat *.apf *.log
wait                           # Ensure file deletion finishes completely

# ------------------------------------------------------------------------------
# STAGE 1: Trajectory Imaging & Basic Kinetics (RMSD & Distances)
# ------------------------------------------------------------------------------
cat >ana1.ptraj<<EOF
parm ../*.prmtop               # Load the topology (structure definition) file
# Read trajectory from frame 25,000 to the end, skipping the initial equilibration
trajin 05_prod1000.nc 25000 last 1

# --- Periodic Boundary Condition (PBC) Imaging ---
# Center the system (Protein + Ligand) at the coordinate origin
center origin :1-278
# Wrap water/ions back into a familiar unit cell geometry around the centered solute
image origin center familiar

# --- Strip Solvent & Ions ---
parmstrip :WAT,Na+,Cl-         # Strips topology references to save memory
strip :WAT                     # Strips coordinates from the trajectory frames
strip :Na+
strip :Cl-

# --- Structural Alignments & Basic RMSD ---
# Align the protein backbone to the first frame to remove global translation/rotation
rms protein_fit :1-277@CA,C,N,O first
# Calculate global complex RMSD (Protein + Ligand, excluding Hydrogens) without a second fit
rms complex_rmsd :1-278@/H first nofit 
# Calculate independent ligand displacement tracking
rms ligand_rmsd :278@/H first nofit out rms.lig.dat

# --- Distance Tracking (Binding Site Kinetics) ---
# "noimage" skips PBC math here because the "image origin" command above already fixed wrapping
distance pend :278@F1,F2,F3 :118,124,147 out dist.pend.dat noimage
distance aend :278@F4,F5 :106,233,262 out dist.aend.dat noimage

run                            # Execute Stage 1
EOF

wait
mpirun -np 50 cpptraj.MPI -i ana1.ptraj
wait
rm -f ana1.ptraj

# ------------------------------------------------------------------------------
# STAGE 2: Dynamics & Thermodynamics (Ligand Quasi-Harmonic Entropy)
# ------------------------------------------------------------------------------
cat >ana1.ptraj<<EOF
parm ../*.prmtop
trajin 05_prod1000.nc 25000 last 1
center origin :1-278
image origin center familiar

# Strip everything down strictly to the ligand to calculate its configuration matrix
strip :WAT
strip :Na+
strip :Cl-
strip :1-277                   # Drop the protein completely for this step

# Calculate the mass-weighted covariance matrix of the isolated ligand
matrix mwcovar name arcov out arcov.dat 
# Diagonalize the matrix to extract vibrational frequencies and calculate classical entropy at 300K
diagmatrix arcov out qh.out thermo outthermo thermo.dat temp 300 

run                            # Execute Stage 2
EOF

wait
mpirun -np 50 cpptraj.MPI -i ana1.ptraj
wait
rm -f ana1.ptraj

# ------------------------------------------------------------------------------
# STAGE 3: Structural Properties (Fluctuations, Compactness, Hbonds)
# ------------------------------------------------------------------------------
cat >ana1.ptraj<<EOF
parm ../*.prmtop
trajin 05_prod1000.nc 25000 last 1
center origin :1-278
image origin center familiar

strip :WAT
strip :Na+
strip :Cl-

# Ensure the system is aligned to the protein backbone before calculating fluctuations
rms protein_fit :1-277@CA,C,N,O first

# --- Atomic Fluctuations (B-factors/Flexibility) ---
# Calculate per-residue fluctuation for the protein (excluding Hydrogens)
atomicfluct out fluct.pro.apf :1-277&!@H= byres
# Calculate per-atom fluctuation for the ligand (excluding Hydrogens)
atomicfluct out fluct.lig.apf :278&!@H= byatom

# --- Radius of Gyration ---
# Measure the structural compactness of the ligand over time
radgyr test :278&!(@H=) out rog.dat mass nomax

# --- Hydrogen Bonding Analysis ---
# Map intermolecular H-bonds between the protein and ligand
# Criteria: Donor-Acceptor Distance < 3.5 Å, Donor-Hydrogen-Acceptor Angle > 135°
hbond LIG_PROT :1-278 out hbond_ligprot.dat angle 135 dist 3.5 avgout hbond_avg_ligprot.dat nointramol

run                            # Execute Stage 3
EOF

wait
mpirun -np 50 cpptraj.MPI -i ana1.ptraj
wait

# Final clean up of temporary scripts and scheduler logs
rm -f ana1.ptraj analysis.sh.o*
