**Python scripts used in this project**

`get_pdb_box_size.py`: get the box size of each snapshot of an NPT simulation

`solvent_processing.py`: generate truncated solute-solvent clusters from MD trajectory with a given cutoff radius

`eliminate_pbc.py`: get rid of the PBC for the coordinates of the truncated solute-solvent clusters

`make_input_solvent_efield.py`: generate Q-Chem inputs for solvent electric field calculations using SPADE or ALMO partitioning method

`calculate_efield_two_bonds.py`: process Q-Chem outputs and evaluate the electric field projections along chemical bonds by averaging the field projections at two atomic positions

`calculate_efield_fromESP_two_bonds.py`: process Q-Chem outputs and evaluate the electric field projections along chemical bonds using the difference between the electrostatic potentials at two atomic positions

`xyzgeom.py`: supporting functions for processing XYZ geometries
