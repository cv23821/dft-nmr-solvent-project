# DFT NMR Solvent Modelling Scripts
This repisitory contains the python scripts used in a BSc Chemistry final year project investigating whether explicit solvent modelling provides meaningful improvement over implicit approaches fo rprediciting 13C NMR chemcial shifts.

# Scripts
- 	dft_nmr_explicit_15.py
    runs DFT NMR calculations using explicit solvent clusters
- 	dft_nmr_implicit.py
    runs DFT NMR calculations using implcit solvent model
-   solvent.py
    generates explicit solvent clusters around the solute structure
-   optimisation.py
    optimises molecular geometires using GFN-xTB with an implicit solvent correction

# Example inputs
-   example_SMILE.txt
    example SMILES strings used as staring structures
-   example_orca_input.inp
    example ORCA input file for DFT NMR calculation
