pol_il
======

_Kateryna Goloviznina_ \
_[Agilio Padua](http://perso.ens-lyon.fr/agilio.padua)_

Contents
--------

* `polarizer.py`: introduces Drude induced dipoles into [LAMMPS](http://lammps.sandia.gov/) input files.

* `scaleLJ.py`: scales pair coefficients of Lennard-Jones interactions.

* `alpha.ff`: Drude induced dipole database.

* `fragment.ff`: fragment database.

* `fragment_topologies/`: structure files for typical IL fragments.

* `examples/`: examples of [C4C1im][DCA] molecule files and force field database.

Requirements
------------

* [Python](http://www.python.org/)
* [fftool](https://github.com/agiliopadua/fftool)
* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)

Obtaining
---------

Download the files or clone the repository:

git clone https://github.com/agiliopadua/pol_il.git


Tutorial
--------

These are instructions on how to build an initial configuration for a
polarisable system composed of molecules, ions or materials.

A system consisting of one [C4C1im]+ cation and one [DCA]- anion is considered as an example. All input files can be found in the `example/` folder.

1. Use `fftool` to create `data.lmp`, `in.lmp` and `pair.lmp` files. A separate `pair.lmp` file containing all i-j pair coefficients is required for future procedures and can be created using `-a` option in `fftool`. The detailed instructions on how to use `fftool` can be found [here](https://github.com/agiliopadua/fftool).

        fftool 1 c4c1im.zmat 1 dca.zmat -b 20
        packmol <pack.inp
        fftool 1 c4c1im.zmat 1 dca.zmat -b 20 -a -l

2. Add Drude induced dipoles to LAMMPS data file using `polarizer.py` script.

        python polarizer.py c4c1im.zmat dca.zmat -f alpha.ff -q -id data.lmp -od data-p.lmp

   The script requires a file containing the specification of Drude induced dipoles according to the next format and specified by `-q` option.

        # alpha.ff
        type  dm/u  dq/e  k/(kJ/molA2)  alpha/A3  thole
        CR      0.4    -1.0     4184.0   1.122   2.6    
        NA      0.4    -1.0     4184.0   1.208   2.6 
        ...
    where
    * `dm` is the mass to place on the Drude particle (taken from its core),
    * `dq` is the charge to place on the Drude particle (taken from its core),
    * `k` is the harmonic force constant of the bond between core and Drude,
    * `alpha` is the polarizability, hyrdogen aroms are not merged,
    * `thole` is a parameter of the Thole damping function.

    The harmonic constant and the charge on the Drude particle are related though

    <img src="https://latex.codecogs.com/svg.latex?\Large&space;\alpha = q_D^2/k_D" title="\Large \alpha = q_D^2/k_D" />

    Use the `-q` option to read the force constant from the input file and to calculate Drude charge from the polarizabilities or the `-k` option to read the Drude charge from `alpha.ff` and to recalculate the force constant according to the relation above.

    A Drude particle is created for each atom in the LAMMPS data file
    that corresponds to an atom type given in the Drude file.
    Since LAMMPS uses numbers for atom types in the data file, a comment
    after each line in the Masses section has to be introduced to allow
    identification of the atom types within the force field database:

          Masses
          1   14.007  # NA
          2   12.011  # CR
          ...

    This script adds new atom types, new bond types, new atoms and
    new bonds to the `data-p.lmp` file.
    It generates the commands to be included in the `LAMMPS` input script and outputs to the terminal. They are related to the topology and force field, namely `fix drude`, `pair_style` commands and examples of thermostats and barostats are proposed:

        # Commands to include in the LAMMPS input script

        # adapt the pair_style command as needed
        pair_style hybrid/overlay ... coul/long/cs 12.0 thole 2.600 12.0

        # data file with Drude oscillators added
        read_data data-p.lmp

        # pair interactions with Drude particles written to file
        # Thole damping recommended if more than 1 Drude per molecule
        include pair-drude.lmp

        # atom groups convenient for thermostats (see package documentation), etc.
        group ATOMS type 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
        group CORES type 1 2 3 4 6 9 10 12 13 14 15
        group DRUDES type 16 17 18 19 20 21 22 23 24 25 26

        # flag for each atom type: [C]ore, [D]rude, [N]on-polarizable
        fix DRUDE all drude C C C C N C N N C C N C C C C D D D D D D D D D D D

        # store velocity information of ghost atoms
        comm_modify vel yes

        # compute the temperatures of ATOMS, DC-DP pair centers of mass and DPs
        compute TATOM ATOMS temp
        compute TDRUDE all temp/drude

        # examples of termostats and barotats
        # NVT (Nose-Hoover)
        fix DTDIR all drude/transform/direct
        fix TSTAT ATOMS nvt temp 300.0 300.0 200
        fix TSTDR DRUDES nvt temp 1.0 1.0 50
        fix DTINV all drude/transform/inverse
        # NPT (Nose-Hoover)
        fix DTDIR all drude/transform/direct
        fix TSTAT ATOMS npt temp 300.0 300.0 200 iso 1.0 1.0 1000
        fix_modify TSTAT temp TATOM press thermo_press
        fix TSTDR DRUDES nvt temp 1.0 1.0 50
        fix DTINV all drude/transform/inverse

        # avoiding the flying ice cube artefact
        fix ICECUBE all momentum 1000 linear 1 1 1

        # output the temperatures of ATOMS, DC-DP pair centers of mass and DPs
        thermo_style custom step [...] c_TATOM c_TDRUDE[1] c_TDRUDE[2]

        # write Drude particles to dump file
        dump_modify ... element ... D D D D D D D D D D D

        # ATTENTION!
        #  * read_data may need 'extra/special/per/atom' keyword, LAMMPS will exit with a    message.
        #  * If using fix shake the group-ID must not include Drude particles. Use group     ATOMS, for example.
        #  * Give all I<=J pair interactions, no mixing.
        #  * Pair style coul/long/cs from CORESHELL package is used for interactions
        #    of Drude particles. Alternatively pair lj/cut/thole/long could be used,
        #    avoiding hybrid/overlay and allowing mixing. See doc pages.

    Pair i-j interactions between induced dipoles are described by `pair_coeff` in `pair-drude.lmp`, created by this script. We propose to concatenate `pair.lmp` and `pair-drude.lmp` files into `pair-p.lmp`.

3. Scale parameters of LJ interactions between particular fragments.

        python scaleLJ.py -f fragment.ff -a alpha.ff -i fragment.inp -ip pair-p.lmp -op pair-p-sc.lmp

    The script requires several input files with fragment specification, structure files of fragments in common formats (`.xyz`, `.zmat`, `.mol`, `.pdb`) and `pair-p.lmp` file.

    Format of file containing monomers and dimers specification:

        # fragment.ff
        MONOMERS
        # name       q/e       mu/D
        c2c1im       1.0      1.1558
        ..
        DIMERS
        # m1         m2       r_COM/A    k_sapt
        c2c1im       dca       2.935      0.61
        ...
    where
    * `q` is the charge of the monomer,
    * `mu` is the dipole moment of the monomer,
    * `m1` and `m2` are the monomers forming a dimer,
    * `r_COM` is the distance between the centers of mass of the monomers,
    * `k_sapt` is the scaling factor for the epsilon of LJ potential, obtained by SAPT quantum calculation (optional).

    Format of file containing fragment list with atomic indices:

        # fragment.inp
        # c4c1im dca
        c2c1im 1:8
        C4H10  9:12
        dca   13:15

    where atomic indices or/and a range of indices correspond to atomic types associating with this fragment in `data.lmp` file. In this example, NA, CR, CW, C1, HCR, C1A, HCW, H1 belong to c2c1im fragment, C2, CS, HC, CT to C4H10 fragment and N3A, CZA, NZA to dca fragment.

    The script performs modification of Lennard-Jones interaction between atoms of the fragments.
    In order to remove a double counting of the induction effects, which are included implicitly in the empirical LJ potential, the epsilon value should be scaled. By default, the scaling factor is predicted by this script on the basis of simple properties. On the other hand, it can be obtained through quantum chemistry calculation, Symmetry-Adapted Perturbation Theory (SAPT), and it can be used by the `-q` option. Scaling of the sigma value allows to adjust density of the system (if necessary), it can be enabled using `-s` option with a default value of 0.985.

    Scaled epsilon (and sigma) values for LJ interaction are printed into `pair-p-sc.lmp` file that directly can be used by LAMMPS. The scaling coefficients are outputted to the terminal.
    If they were obtained by the prediction scheme, SAPT calculated values (if availiable) are given only for the comparison.

        Epsilon LJ parameters were scaled by k_pred parameter. Changes are marked with '~'.
        Sigma LJ parameters were not scaled.
        ----------------------------------------------
         Fragment1   Fragment2    k_sapt     k_pred
            c2c1im       c4h10      0.76       0.78
            c2c1im         dca      0.61       0.68
             c4h10       c4h10      0.94       1.00
             c4h10         dca      0.69       0.72
        ----------------------------------------------
    To perform the simulation of a polarisable system, the `USER-DRUDE` package should be enabled during LAMMPS compilation.

References
----------

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/):
  L. Martinez et al. J Comp Chem 30 (2009) 2157, DOI:
  [10.1002/jcc.21224](http://dx.doi.org/10.1002/jcc.21224).
  
* [LAMMPS](http://lammps.sandia.gov/): S. Plimton, J Comp Phys
  117 (1995) 1, DOI:
  [10.1006/jcph.1995.1039](http://dx.doi.org/10.1006/jcph.1995.1039).

* CL&Pol: K. Goloviznina, J. N. Canongia Lopes, M. Costa Gomes, A. A. H. Pádua,  J. Chem. Theory Comput. 15 (2019), 5858, DOI:
  [10.1021/acs.jctc.9b00689](https://doi.org/10.1021/acs.jctc.9b00689).
