# 2D_VSE_probe
Python code used for the 2D VSE probe paper

get_pdb_box_size.py: script to get the box size of each snapshot of an NPT simulation

solvent_processing.py: code to generate truncated solute-solvent clusters with a given cutoff

eliminate_pbc.py: script to get rid of the PBC of the truncated solute-solvent configurations

calculate_efield_two_bonds: code to process Q-Chem output and evaluate the electric field projections along chemical bonds (by averaging the field projections at                               two atomic positions)

calculate_efield_fromESP:   code to process Q-Chem output and evaluate the electric field projections along chemical bonds using the difference between the                                     electrostatic potentials at two atomic positions
