ChCl-EG
======

Here are the examples of molecule files and force field database for choline chloride - ethylene glycol system. The force field files for other deep eutectic solvents can be found [here](https://github.com/kateryna-goloviznina/desff). 

Follow the [tutorial](https://github.com/kateryna-goloviznina/pol_il/tree/master) to create the input files using the polarisable force field.

1. Steps 1-2 are identical to the ones of the **Example 2**.

        fftool 1 ch.xyz 1 Cl.zmat 2 EG.zmat -b 20
        packmol <pack.inp
        fftool 1 ch.xyz 1 Cl.zmat 2 EG.zmat -b 20 -a -l
        polarizer -f alpha-des.ff data.lmp data-p.lmp
        coul_tt -a 8 13
        scaleLJ -f fragment-des.ff -a alpha-des.ff -s

2. Step 3. According to our recently published paper, the sigma LJ parameter of all O-Cl interactions (in particular, OH and OHG atom types) should be increased from 3.37 Å to 3.70 Å when Cl Lennard-Jones parameters from JPCB 108 (2004) 2038 are used. Hovewer, this step is not required when Cl force field parameters are taken directly from OPLS-AA force field as given in `des.ff` file.

## References

* CL&Pol for PIL, DES, electrolytes: K. Goloviznina, Z. Gong, M. Costa Gomes,
  A. A. H. Pádua, J. Chem. Theory Comput. (2021) DOI:
  [10.26434/10.1021/acs.jctc.0c01002](https://doi.org/10.1021/acs.jctc.0c01002), ChemRxiv [10.26434/chemrxiv.12999524](https://doi.org/10.26434/chemrxiv.12999524)
