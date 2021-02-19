# pol_il

_Kateryna Goloviznina_ \
_[Agilio Padua](http://perso.ens-lyon.fr/agilio.padua)_

## Contents

* `polarizer.py`: introduces Drude induced dipoles into [LAMMPS](http://lammps.sandia.gov/) input files.

* `scaleLJ.py`: scales pair coefficients of Lennard-Jones interactions.

* `coul_tt.py`: adds Tang Toennies charge-dipole damping for densely charged atoms.

* `alpha.ff`: Drude induced dipole database.

* `fragment.ff`: fragment database.

* `fragment_topologies/`: structure files for typical IL fragments.

* `example_il/`: examples of [C4C1im][DCA] molecule files and force field database (aprotic ionic liquid).

* `example_pil/`: examples of ethylammonium nitrate (EAN) molecule files and force field database (protic ionic liquid).

* `example_des/`: examples of choline chloride - ethylene glycol (ChCl-EG) molecule files and force field database (deep eutectic solvent).

## Requirements

* [Python](http://www.python.org/)
* [fftool](https://github.com/agiliopadua/fftool)
* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)


## Obtaining

Download the files or clone the repository:

    git clone https://github.com/agiliopadua/pol_il.git


## Tutorial

These are instructions on how to build an initial configuration for a
polarisable system composed of molecules, ions or materials.

To perform the simulation of a polarisable system, the `USER-DRUDE` package should be enabled during LAMMPS compilation.

Two systems, 1-butyl-3-methylimidazolium dicyanamide ([C4C1im][DCA]) and
ethylammonium nitrate (EAN), are considered as examples of aprotic and protic ionic liquids, respectively.

### Example 1. Aprotic ionic liquid

The input files for a system consisting of one [C4C1im]+ cation and one [DCA]- anion can be found in the `example_il/` folder.


#### 1. Create input files for a non-polarizable system 

Use `fftool` to create `data.lmp`, `in.lmp` and `pair.lmp` files. A separate `pair.lmp` file containing all i-j pair coefficients is required for future procedures and can be created using the `-a` option of `fftool`. The detailed instructions on how to use `fftool` can be found [here](https://github.com/agiliopadua/fftool).

    fftool 1 c4c1im.zmat 1 dca.zmat -b 20
    packmol < pack.inp
    fftool 1 c4c1im.zmat 1 dca.zmat -b 20 -a -l


#### 2. Add Drude induced dipoles to LAMMPS data file

    python polarizer.py c4c1im.zmat dca.zmat -f alpha.ff -q -id data.lmp -od data-p.lmp

The script requires a file containing the specification of Drude induced dipoles according to the following format, and specified with the `-f` option:

    # alpha.ff
    type  dm/u  dq/e  k/(kJ/molA2)  alpha/A3  thole
    CR    0.4   -1.0     4184.0     1.122     2.6
    NA    0.4   -1.0     4184.0     1.208     2.6
    ...

where
    * `dm` is the mass to place on the Drude particle (taken from its core),
    * `dq` is the charge to place on the Drude particle (taken from its core),
    * `k` is the harmonic force constant of the bond between core and Drude,
    * `alpha` is the polarizability (hydrogen aroms are not merged),
    * `thole` is a parameter of the Thole damping function.

The harmonic constant and the charge on the Drude particle are related though the relation <img src="https://render.githubusercontent.com/render/math?math=\alpha = q_D^2/k_D">. Use the `-q` option to read the force constant from the input file and to calculate Drude charge from the polarizabilities or the `-k` option to read the Drude charge from `alpha.ff` and to recalculate the force constant according to the relation above.

A Drude particle is created for each atom in the LAMMPS data file that corresponds to an atom type given in the Drude file. Since LAMMPS uses numbers for atom types in the data file, a comment after each line in the Masses section needs to be introduced in the original data file to allow identification of the atom types within the force field database:

    Masses
    1   14.007  # NA
    2   12.011  # CR
    ...

This script then adds new atom types, new bond types, new atoms and new bonds into a new `data-p.lmp` file. It also generates commands to be included in the `LAMMPS` input script and outputs them to the terminal. They are related to the topology and the force field, namely the `fix drude` and `pair_style` commands, and include examples of thermostats:

    # Commands to include in the LAMMPS input script

    # adapt the pair_style command as needed
    pair_style hybrid/overlay ... coul/long/cs 12.0 thole 2.600 12.0

    # new data file with Drude oscillators added
    read_data data-p.lmp

    # read pair interactions involving Drude particles
    # Thole damping recommended if more than 1 Drude per molecule
    include pair-drude.lmp

    # convenient atom groups (for shake, thermostats...)
    group ATOMS type 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
    group CORES type 1 2 3 4 6 9 10 12 13 14 15
    group DRUDES type 16 17 18 19 20 21 22 23 24 25 26

    # identify each atom type: [C]ore, [D]rude, [N]on-polarizable
    fix DRUDE all drude C C C C N C N N C C N C C C C D D D D D D D D D D D

    # store velocity information of ghost atoms
    comm_modify vel yes

    variable TK equal 300.0
    variable TDRUDE equal 1.0
    variable PBAR equal 1.0

    # temperature-grouped multiple Nosé-Hoover thermostats and barostat
    fix TSTAT all tgnpt/drude temp ${TK} ${TK} 100 ${TDRUDE} 20 iso ${PBAR} ${PBAR} 1000

    # output the temperatures of molecular COM, COM of DC-DP, and DP
    thermo_style custom step [...] f_TSTAT[1] f_TSTAT[2] f_TSTAT[3]

    # write Drude particles to dump file
    dump_modify ... element ... D D D D D D D D D D D

    # ATTENTION!
    #  * read_data may need 'extra/special/per/atom' keyword, LAMMPS will exit with a message
    #  * if using fix shake the group-ID must not include Drude particles; use group ATOMS
    #  * give all I<=J pair interactions, no mixing
    #  * pair style coul/long/cs from CORESHELL package is used for interactions of DP;
    #    alternatively pair lj/cut/thole/long could be used avoiding hybrid/overlay and
    #    allowing mixing; see doc pages.

Pair i-j interactions between induced dipoles are described by `pair_coeff` in `pair-drude.lmp`. The `pair.lmp` and `pair-drude.lmp` files can be concatenated into `pair-p.lmp`, which is used for the next step.


#### 3. Scale LJ interactions between fragments

    python scaleLJ.py -f fragment.ff -a alpha.ff -i fragment.inp -ip pair-p.lmp -op pair-p-sc.lmp

The script performs modification of Lennard-Jones interaction between atoms of the fragments. To prevent double counting of the induction effects, which are included implicitly in the empirical LJ potential, the epsilon value should be scaled. By default, the scaling factor is predicted by this script on the basis of simple properties, 

<img src="https://render.githubusercontent.com/render/math?math=k_{ij} = \bigg (1 +  c_0 r_{ij}^2 \frac{ Q_i^2\alpha_j
    + Q_j^2\alpha_i}{\alpha_i\alpha_j}  +   c_1 \frac{\mu_i^2 \alpha_j +
    \mu_j^2 \alpha_i}{\alpha_i \alpha_j} \bigg )^{-1}">

It can also be obtained through quantum chemistry calculation, via Symmetry-Adapted Perturbation Theory (SAPT), which can be invoked with the `-q` option. Scaling of the sigma value allows adjustment of the density of the system (if necessary), and can be enabled using the `-s` option which has a default value of 0.985.

    python scaleLJ.py [...] -s                   - scale all fragments' sigma by 0.985
    python scaleLJ.py [...] -s 0.9               - scale all fragments' sigma by a user-defined value
    python scaleLJ.py [...] -s c2c1im c4h10      - scale the specified fragments' sigma by 0.985
    python scaleLJ.py [...] -s 0.9 c2c1im c4h10  - scale the specified fragments' sigma by a user-defined value

The script requires several input files with fragment specification, structure files of fragments in common formats (`.xyz`, `.zmat`, `.mol`, `.pdb`), and the `pair-p.lmp` file.

The file containing specification of monomers and dimers has the following format:

    # fragment.ff
    MONOMERS
    # name       q/e      mu/D
    c2c1im       1.0      1.1558
    ...
    DIMERS
    # m1         m2       r_COM/A    k_sapt
    c2c1im       dca      2.935      0.61
    ...

where
    * `q` is the charge of the monomer,
    * `mu` is the dipole moment of the monomer,
    * `m1` and `m2` are the monomers forming a dimer,
    * `r_COM` is the distance between the centers of mass of the monomers,
    * `k_sapt` is the scaling factor for the epsilon of LJ potential, obtained by SAPT quantum calculation (optional).

If equilibrium distances are missing for certain fragment dimers, these can be obtained from a geometry optimization (we used dispersion-corrected DFT, B97+D3/cc-pVTZ). Charges and dipole moments for new fragments can be calculated using quantum chemistry (we used the level MP2/cc-pVTZ). Details are given in the [CL&Pol] paper.

Format of file containing fragment list with atomic type indices:

    # fragment.inp
    # c4c1im dca
    c2c1im  1:8
    C4H10   9:12
    dca    13:15

where atomic indices or/and a range of indices correspond to atomic types associating with this fragment in the `data.lmp` file. In this example, NA, CR, CW, C1, HCR, C1A, HCW, H1 belong to the c2c1im fragment; C2, CS, HC, CT are from the C4H10 fragment and N3A, CZA, NZA are from the dca fragment. Thus, the script requires `c2c1im.zmat`, `C4H10.zmat` and `dca.zmat` structure files.

Scaled epsilon (and sigma) values for LJ interaction are printed into a `pair-p-sc.lmp` file that directly can be used by LAMMPS. The scaling coefficients are outputted to the terminal. If they are obtained by the prediction scheme, SAPT calculated values (if available) are given only for the comparison.

    Epsilon LJ parameters were scaled by k_pred parameter. Changes are marked with '~'.
    Sigma LJ parameters were not scaled.
    ------------------------------------------
    Fragment_i   Fragment_j   k_sapt    k_pred
    c2c1im       c4h10          0.76      0.78
    c2c1im       dca            0.61      0.68
    c4h10        c4h10          0.94      1.00
    c4h10        dca            0.69      0.72
    ------------------------------------------

The CL&Pol force field can be mixed with another polarisable force field, for example SWM4-NDP model of water. In this case, the scaling of epsilon should be performed only partially that can be envoked with `-p` option

    python scaleLJ.py [...] -p swm4-ndp

The scaling coefficient will depend only on the charge, dipole and molecular polarisability of this fragment

 <img src="https://render.githubusercontent.com/render/math?math=k_{ij} = \bigg (1 + c_0 r_{ij}^2 \frac{ Q_i^2}{\alpha_i}  +   c_1 \frac{\mu_i^2}{\alpha_i} \bigg )^{-1}">

that should be specified in the input files for the script.

### Example 2. Protic ionic liquid (or other strongly H-bonded systems)

The input files of a system consisting of one ethylammonium nitrate ion pair can be found in the `example_pil/` folder.

#### 1. Steps 1 to 3 are identical to Example 1

    fftool 1 N2000.zmat 1 no3.zmat -b 20
    packmol <pack.inp
    fftool 1 N2000.zmat 1 no3.zmat -b 20 -a -l
    python polarizer.py N2000.zmat no3.zmat -f alpha.ff -q
    python scaleLJ.py -q 

Here modification of the parameters of LJ interactions between N2000 and NO3 fragments is performed using the scaling factor obtained through SAPT calculation (invoked with the `-q` option).

#### 2. Add short range damping of charge-dipole Coulomb interactions

This is almos always needed between small, highly charged atoms (such as hydrogen) and induced dipoles to prevent the "polarization catastrophe".

    python coul_tt.py -d data-p.lmp -a 3

The functional form of the damping function is
    
<img src="https://render.githubusercontent.com/render/math?math=f(r) = 1 - c \cdot e^{-b r} \sum_{k=0}^4 \frac{(b r)^k}{k!}">

resulting from an adaptation to the Coulomb interaction of the damping function originally proposed by Tang Toennies for van der Waals interactions. The `b` value is set to 4.5 and the `c` value to 1.0. This function is implemented as `coul/tt` pair style in LAMMPS (version 29Oct20 or newer), the detailed description is given [here](https://lammps.sandia.gov/doc/pair_coul_tt.html). 

The script requires the `data-p.lmp` file to obtain the list of atoms and their type (polarisable or non-polarisable). 

     1   13.607  # N1 DC
     2   11.611  # C1N DC
     3    1.008  # HN
     4   11.611  # CEN DC
     5    1.008  # H1N
     6    1.008  # HCN
     7   13.601  # NO DC
     8   15.599  # ON DC
     9    0.400  # N1 DP
    10    0.400  # C1N DP
    11    0.400  # CEN DP
    12    0.400  # NO DP
    13    0.400  # ON DP

The atomic indices of small, highly charged atoms (typically, point charges without LJ sites) should be specified with the `-a` option. The short-range Coulomb interactions of those atoms with all Drude cores and Drude particles should be damped. The corresponding `pair_coeff` lines are written to the `pair-tt.lmp` file.

    pair_coeff    1    3 coul/tt 4.5 1.0
    pair_coeff    2    3 coul/tt 4.5 1.0
    pair_coeff    3    4 coul/tt 4.5 1.0
    pair_coeff    3    7 coul/tt 4.5 1.0
    pair_coeff    3    8 coul/tt 4.5 1.0
    pair_coeff    3   9* coul/tt 4.5 1.0

Here, the damped interactions are the ones of the HN atom (index 3) with Drude cores (indices 1, 2, 4, 7, 8) and Drude particles (indices 9-13).

The script prints the commands to be included by the user to the `in-p.lmp` file to declare the `coul/tt` pair style.

    To inlcude to in-p.lmp:
        pair_style hybrid/overlay ... coul/tt 4 12.0
        include pair-tt.lmp

#### 3. Modify the LJ interaction parameters of i-j pairs involved in hydrogen bonds

Hydrogen bonds (D-H...A) involving hydrogen atoms represented by 'naked' charges without Lennard-Jones sites could lead to "freezing" of a system when modelled using polarisable force field. To avoid this effect, the repulsive potential between D and A sites should be adjusted with a typical value of the sigma parameter of 3.7 to 3.8 Å.

In EAN, the hydrogen bond is formed between HN hydrogens atoms of the cation embedded into neighbouring N1 nitrogen atoms and the ON oxygen atoms of the anion. The sigma LJ parameter of the N1-ON interaction should be increased from 3.10 Å to 3.75 Å.

<img src="example_pil/ean.png" alt="ean" width="300"/>

The `pair_coeff` parameters of the interaction between N1 and ON atoms in the `pair-p-sc.lmp` file

    pair_coeff    1    8 lj/cut/coul/long     0.037789     3.101612  # N1 ON 
        
should be replaced by
        
    pair_coeff    1    8 lj/cut/coul/long     0.037789     3.750000  # N1 ON 

## References

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/):
  L. Martinez et al. J Comp Chem 30 (2009) 2157, DOI:
  [10.1002/jcc.21224](http://dx.doi.org/10.1002/jcc.21224)

* [LAMMPS](http://lammps.sandia.gov/): S. Plimton, J Comp Phys
  117 (1995) 1, DOI:
  [10.1006/jcph.1995.1039](http://dx.doi.org/10.1006/jcph.1995.1039)

* [CL&Pol](https://doi.org/10.1021/acs.jctc.9b00689): K. Goloviznina, J. N. Canongia Lopes, M. Costa Gomes, A. A. H. Pádua,  J. Chem. Theory Comput. 15 (2019), 5858, DOI:
  [10.1021/acs.jctc.9b00689](https://doi.org/10.1021/acs.jctc.9b00689), arXiv [1703.01540](https://arxiv.org/abs/1703.01540)

* CL&Pol for PIL, DES, electrolytes: K. Goloviznina, Z. Gong, M. Costa Gomes,
  A. A. H. Pádua, J. Chem. Theory Comput. (2021) DOI:
  [10.26434/10.1021/acs.jctc.0c01002](https://doi.org/10.1021/acs.jctc.0c01002), ChemRxiv [10.26434/chemrxiv.12999524](https://doi.org/10.26434/chemrxiv.12999524)
