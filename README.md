# Iterative-protein-docking
Credit: based on 2021 Oliver Powell https://github.com/Orpowell/autodock-vina-automator
  with changes to iterate over a range of spatial coordinates instead of multiple ligands; the space vector is shifted by half box size in one direction every iteration
        the creators of autodock vina
        #################################################################
        # If you used AutoDock Vina in your work, please cite:          #
        #                                                               #
        # J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli  #
        # AutoDock Vina 1.2.0: New Docking Methods, Expanded Force      #
        # Field, and Python Bindings, J. Chem. Inf. Model. (2021)       #
        # DOI 10.1021/acs.jcim.1c00203                                  #
        #                                                               #
        # O. Trott, A. J. Olson,                                        #
        # AutoDock Vina: improving the speed and accuracy of docking    #
        # with a new scoring function, efficient optimization and       #
        # multithreading, J. Comp. Chem. (2010)                         #
        # DOI 10.1002/jcc.21334                                         #
        #                                                               #
        # Please see https://github.com/ccsb-scripps/AutoDock-Vina for  #
        # more information.                                             #
        #################################################################
        reduce is from Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.

structure for use:
Main Folder
|_iterative docking algorithm.py
|_Project Folder
  |_receptor protein.pdbqt
  |_ligand.pdbqt
  |_vina_1.2.5_win.exe

Upon running the following folders will be created within the Main Folder
conformations = each contains the autodock fitted ligand atomic spatial coordinates in pdbqt for an individual autodock run (to the corresponding space vector)
logs = each contains the results of an individual autodock run (to the corresponding space vector)
results = the aggregated extracted numerical results
In the Main Folder the parameter files for each individual space vector will be created.

The algoritm requires:
installed python (probably 3.12 or newer) with numpy and pandas packages
receptor protein template as pdbqt file
fitted ligand as pdbqt file
box size - the size of the space fragment that will be surveyed for binding (ideally cubic)
min coordinates - the lowest integer x, y, z coordinates around the receptor (i.e. the closest bottom left corner) -|_ideally the difference between min and max is an integer multiple
max coordinates - the highest integer x, y, z coordinated around the receptor (i.e. the furthest top right corner) -|  of the half box size

Getting the required files/information:
Receptor file:
1. Clean up pdb file from unnecessary molecules, heteroatoms, water, protein chains/fragments. Pymol is a good software for it https://github.com/schrodinger/pymol-open-source.
2. Run reduce to correct Asn/Gln orientations. The algorithm is part of the ADFR software package https://ccsb.scripps.edu/adfr/downloads/.
    In command prompt: "...\ADFRsuite-1.0\bin\reduce" "...\input_receptor.pdb" > "...\output_receptor.pdb"
3. Convert the pdb file into pdbqt.
    In command prompt: "...\ADFRsuite-1.0\bin\prepare_receptor" -r "...\input_receptor.pdb" -o "...\output_receptor.pdbqt"

Ligand file:
1. Get a 3D coordinate file for the molecule usully .sdf format. Can be from pubchem or other databases or de novo from openbabel or other cheminformatic software
2. Ideally optimise the structure so that it reflects a realistic conformation. (especially rigid rings and other similar structures)
3. Convert sdf to pdbqt file. Openbabel is a good software for it. https://openbabel.org/index.html
    add hydrogens make explicit
    add hydrogens appropriate for pH 7.4
    gasteiger calculate partial charges

Getting coordinates: There are two alternatives. (The space vector is always the center of the box, not the side)
- Open the receptor pdbqt file.
  Sort by the x, y, z coordinates and manually write down the highest and lowest values.
  Add quarter of the box size to the highest and subtract from the lowest.
  Adjust values so that the difference is an integer multiple of the half box size. (in theory if it isn't it wouldn't break it, but doesn't hurt)
- Open the receptor pdbqt file in AutoDockTools. https://ccsb.scripps.edu/mgltools/downloads/
  Open the "autodockFR" menu line.
  Under "Receptor" open "Set the docking box"
  Use the docking box dialoge to find the maximum and minimum in each side manually. You can right click on the scroll wheels to change values and sensitivity.
- Coordinates can be double checked with autodocktools to make sure they cover the receptor on every side.
  
