#!/bin/bash

# Source the AmberTools24 environment configuration script to expose Amber binaries
source /Users/kehuang/amber24/amber.sh

# Prompt the user for system-specific names
read -p "Mol. name: " mname      # Core name matching the ligand's .frcmod and .lib files
read -e -p "Input PDB: " pdb     # Filename of the input PDB (supports tab-completion via -e)
#read -p "Mol. shorthand: " mabv # [Unused] Shorthand name for directory pathing

# Define the output topology (.prmtop) and coordinate (.inpcrd) filenames
prmtop=${mname}_${pdb}.prmtop
inpcrd=${mname}_${pdb}.inpcrd

# 2. Target Salt Concentrations (Molarity)
# ------------------------------------------------------------------------------
m1=0.14 # Target concentration for Potassium ions (K+) in Molar
m2=0.01 # Target concentration for Sodium ions (Na+) in Molar
m3=0.15 # Target concentration for Chloride ions (Cl-) in Molar

# Note: Custom forcefield files copy/load commented out here
#cp ../../${mabv}/*.frcmod .
#cp ../../${mabv}/*.lib .

# 3. Step A - Run Initial tLEaP to Determine Box Volume
# ------------------------------------------------------------------------------
# Generate a temporary tLEaP script to build the box without specific salt concentrations.
# This is necessary because we need the box volume to calculate the correct number of ions.
cat>tleap.scrpt<<EOF
source leaprc.gaff2          # Load General Amber Force Field (GAFF2) for the ligand
source leaprc.water.tip3p    # Load TIP3P water model parameters
source leaprc.protein.ff19SB # Load the protein force field (ff19SB)
loadamberparams frcmod.ionsjc_tip3p # Load Joung-Cheatham ion parameters for TIP3P
loadamberparams ${mname}.mol2.frcmod # Load the ligand-specific parameter modifications
loadoff ${mname}.mol2.lib    # Load the ligand-specific library/off file
m1 = loadpdb ${pdb}          # Load the input PDB file into unit 'm1'
solvateOct m1 TIP3PBOX 15.0 iso 1.0 # Solvate in a truncated octahedron box with a 15 Å buffer
#addIonsRand m1 Na+ 0
#addIonsRand m1 Cl- 0
saveamberparm m1 test.prmtop test.inpcrd # Save a temporary topology and coordinate file
quit
EOF

wait # Ensure the script file is fully written to the disk

# Run tLEaP with the temporary script to generate the test box
tleap -f tleap.scrpt

# Clean up temporary tLEaP execution files
rm *.scrpt *.dat *.*~ leap.log

# 4. Step B - Calculate Volume & Determine Ion Counts
# ------------------------------------------------------------------------------
# Create a CPPTRAJ script to measure the exact volume of the temporary water box
cat >ptraj.vol.in<<EOF
parm test.prmtop
trajin test.inpcrd
volume test out tvol.dat     # Measures box volume over the coordinates and outputs to tvol.dat
EOF

# Run CPPTRAJ to extract the volume
cpptraj -i ptraj.vol.in

# Extract the volume value (in cubic Angstroms, Å³) from the second row, second column of tvol.dat
tvol=`awk 'NR>1{print $2}' tvol.dat` 

# Round the volume value to the nearest integer using Python
fvol=`python -c "print (round($tvol))"`

echo $fvol # Print the rounded volume to the terminal

# Convert the box volume from cubic Angstroms (Å³) to Liters (L)
# Math: 1 Å³ = 10^-24 cm³ = 10^-27 L. Done using 'bc' for high-precision decimal math.
nV=`echo "scale=50; $fvol*10^-27" | bc`

# Calculate Cl- ion count using Avogadro's number: Ions = Volume(L) * Molarity(M) * N_A
fcon1=$(echo "scale=50; $nV*$m3*6.022*10^23" | bc)
clcon=$(python3 -c "print(int(round(float('$fcon1'))))") # Round to nearest whole ion

# Calculate Na+ ion count using Avogadro's number
fcon2=$(echo "scale=50; $nV*$m2*6.022*10^23" | bc)
nacon=$(python3 -c "print(int(round(float('$fcon2'))))") # Round to nearest whole ion

# Calculate K+ ion count using Avogadro's number
fcon3=$(echo "scale=50; $nV*$m1*6.022*10^23" | bc)
kcon=$(python3 -c "print(int(round(float('$fcon3'))))") # Round to nearest whole ion

# Output the calculated number of ions to add (Note: literal single quotes prevent variable expansion here)
echo '[Cl]=$clcon [Na]=$nacon [K]=$kcon'

# 5. Step C - Build the Final Box with Calculated Ion Counts
# ------------------------------------------------------------------------------
# Create the final tLEaP script incorporating the exact number of salt ions
cat >PU1.scrpt<<EOF
source leaprc.gaff2
source leaprc.water.tip3p
source leaprc.protein.ff19SB
loadamberparams frcmod.ionsjc_tip3p
loadamberparams ${mname}.mol2.frcmod
loadoff ${mname}.mol2.lib
m1 = loadpdb ${pdb}
solvateOct m1 TIP3PBOX 15.0 iso 1.0
addIons m1 Na+ 0             # First, neutralize the system's net charge with Na+ if net-negative
addIons m1 Cl- 0             # Or neutralize the system with Cl- if net-positive (e.g., due to Zn ions)
addIonsRand m1 Na+ ${nacon}  # Randomly place the calculated molarity-based Na+ ions
addIonsRand m1 K+ ${kcon}   # Randomly place the calculated molarity-based K+ ions
addIonsRand m1 Cl- ${clcon}  # Randomly place the calculated molarity-based Cl- ions
saveamberparm m1 ${prmtop} ${inpcrd} # Save final production-ready topology and coordinate files
savepdb m1 check_$prmtop.pdb # Save a PDB copy of the finalized box for visual inspection
quit
EOF

# Run tLEaP to build the final system
tleap -f PU1.scrpt

wait # Wait for final file generation to finish completely
        
# Remove all temporary data scripts, logs, and intermediate test topologies
rm *.dat ptraj.vol.in *.*~ leap.log *.scrpt test.prmtop test.inpcrd
