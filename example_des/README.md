ChCl-EG
======

Here are the examples of molecule files and force field database for choline chloride - ethylene glycol system. The force field files for other deep eutectic solvents can be found [here](https://github.com/kateryna-goloviznina/desff). 

Follow the [tutorial](https://github.com/kateryna-goloviznina/pol_il/tree/master) to create the input files using the polarisable force field.

1. Steps 1-2 are identical to the ones of the **Example 2**.

        fftool 1 ch.xyz 1 Cl.zmat 2 EG.zmat -b 20
        packmol <pack.inp
        fftool 1 ch.xyz 1 Cl.zmat 2 EG.zmat -b 20 -a -l
        python polarizer.py ch.xyz Cl.zmat EG.zmat -f alpha.ff -q
        python coul_tt.py -a 8 13
        python scaleLJ.py

2. Step 3. The sigma LJ parameter of all O-Cl interactions (in particular, OH and OHG atom types) should be increased from 3.37 Å to 3.70 Å. 

    The `pair_coeff` parameters of the interaction between N1 and ON atoms in the `pair-p-sc.lmp` file

        pair_coeff    7    9 lj/cut/coul/long     0.088409     3.374611  # OH Cl
        pair_coeff    9   11 lj/cut/coul/long     0.110864     3.374611  # Cl OHG

    should be replaced by

        pair_coeff    7    9 lj/cut/coul/long     0.088409     3.700000  # OH Cl
        pair_coeff    9   11 lj/cut/coul/long     0.110864     3.700000  # Cl OHG