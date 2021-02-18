#!/usr/bin/env python
# polarizer.py - add Drude oscillators to LAMMPS data file.
# Agilio Padua <agilio.padua@ens-lyon.fr>
# Alain Dequidt <alain.dequidt@uca.fr>
# Kateryna Goloviznina <kateryna.goloviznina@ens-lyon.fr>
# version 2021/02/18

import sys
import argparse
import random
import math
import shutil
from copy import deepcopy


usage = """
==============================================================================
Add Drude oscillators to LAMMPS data file.
------------------------------------------------------------------------------
Format of file containing specification of Drude oscillators (alpha.ff):
  # type  dm/u  dq/e  k/(kJ/molA2)  alpha/A3  thole
  C3H     1.0   0.0   4184.0        1.016     2.6
  ...
* dm is the mass to place on the Drude particle (taken from its core),
* dq is the charge to place on the Drude particle (taken from its core),
* k is the harmonic force constant of the bond between core and Drude,
* alpha is the polarizability, hyrdogen aroms are not merged,
* thole is a parameter of the Thole damping function.
------------------------------------------------------------------------------
Script requires structure files of ions/molecules in the system
(.zmat, .xyz, .mol, .pdb) and corresponding force field files (il.ff)
------------------------------------------------------------------------------
A Drude particle is created for each atom in the LAMMPS data file
that corresponds to an atom type given in the Drude file.
Since LAMMPS uses numbers for atom types in the data file, a comment
after each line in the Masses section has to be introduced to allow
identification of the atom types within the force field database:
  Masses
  1   12.011  # C3H
  2   12.011  # CTO
  ...
This script will add new atom types, new bond types, new atoms and
new bonds to the data file.
It will also generate some commands to be included in the LAMMPS input script,
which are related to the topology and force field, namely fix drude,
pair_style and pair_coeff commands. For information on thermostating please
read the documentation of the USER-DRUDE package.
------------------------------------------------------------------------------
This tool can also be used to revert a Drude-polarized data file to a
non-polarizable one.
==============================================================================
"""

# keywords of header and main sections (from data.py in Pizza.py)

hkeywords = ["atoms", "ellipsoids", "lines", "triangles", "bodies",
             "bonds", "angles", "dihedrals", "impropers",
             "atom types", "bond types", "angle types", "dihedral types",
             "improper types", "xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"]

skeywords = [["Masses", "atom types"],
             ["Pair Coeffs", "atom types"],
             ["Bond Coeffs", "bond types"], ["Angle Coeffs", "angle types"],
             ["Dihedral Coeffs", "dihedral types"],
             ["Improper Coeffs", "improper types"],
             ["BondBond Coeffs", "angle types"],
             ["BondAngle Coeffs", "angle types"],
             ["MiddleBondTorsion Coeffs", "dihedral types"],
             ["EndBondTorsion Coeffs", "dihedral types"],
             ["AngleTorsion Coeffs", "dihedral types"],
             ["AngleAngleTorsion Coeffs", "dihedral types"],
             ["BondBond13 Coeffs", "dihedral types"],
             ["AngleAngle Coeffs", "improper types"],
             ["Atoms", "atoms"], ["Velocities", "atoms"],
             ["Ellipsoids", "ellipsoids"],
             ["Lines", "lines"], ["Triangles", "triangles"],
             ["Bodies", "bodies"],
             ["Bonds", "bonds"],
             ["Angles", "angles"], ["Dihedrals", "dihedrals"],
             ["Impropers", "impropers"], ["Molecules", "atoms"]]


def massline(att):
    return "{0:4d} {1:8.3f}  # {2}\n".format(att['id'], att['m'], att['type'])

def bdtline(bdt):
    return "{0:4d} {1:12.6f} {2:12.6f}  {3}\n".format(bdt['id'], bdt['k'],
                                                     bdt['r0'], bdt['note'])

def atomline(at):
    return "{0:7d} {1:7d} {2:4d} {3:10.6f} {4:13.6e} {5:13.6e} {6:13.6e} "\
           " {7}\n".format(at['n'], at['mol'], at['id'], at['q'],
                           at['x'], at['y'], at['z'], at['note'])

def bondline(bd):
    return "{0:7d} {1:4d} {2:7d} {3:7d}  {4}\n".format(bd['n'], bd['id'],
                                            bd['i'], bd['j'], bd['note'])

def velline(at):
    return "{0:7d} {1:13.6e} {2:13.6e} {3:13.6e} \n".format(at['n'],
                                       at['vx'], at['vy'], at['vz'])

# --------------------------------------


class Data(object):

    def __init__(self, datafile):
        '''read LAMMPS data file (from data.py in Pizza.py)'''

        # for extract method
        self.atomtypes = []
        self.bondtypes = []
        self.atoms = []
        self.bonds = []
        self.idmap = {}

        self.nselect = 1

        f = open(datafile, "r")

        self.title = f.readline()
        self.names = {}

        headers = {}
        while 1:
            line = f.readline().strip()
            if '#' in line:
                line = line[:line.index('#')].strip()
            if len(line) == 0:
                continue
            found = 0
            for keyword in hkeywords:
                if keyword in line:
                    found = 1
                    words = line.split()
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                      keyword == "zlo zhi":
                        headers[keyword] = (float(words[0]), float(words[1]))
                    elif keyword == "xy xz yz":
                        headers[keyword] = \
                          (float(words[0]), float(words[1]), float(words[2]))
                    else:
                        headers[keyword] = int(words[0])
            if not found:
                break

        sections = {}
        while 1:
            if len(line) > 0:
                found = 0
                for pair in skeywords:
                    keyword, length = pair[0], pair[1]
                    if keyword == line:
                        found = 1
                        if length not in headers:
                            raise RuntimeError("data section {} "\
                                  "has no matching header value".format(line))
                        f.readline()
                        list_ = []
                        for _ in range(headers[length]):
                            list_.append(f.readline())
                        sections[keyword] = list_
                if not found:
                    raise RuntimeError("invalid section {} in data"\
                                       " file".format(line))
            #f.readline()
            line = f.readline()
            if not line:
                break
            if '#' in line:
                line = line[:line.index('#')]
            line = line.strip()

        f.close()
        self.headers = headers
        self.sections = sections

    def write(self, filename):
        '''write out a LAMMPS data file (from data.py in Pizza.py)'''

        with open(filename, "w") as f:
            f.write(self.title + '\n')
            for keyword in hkeywords:
                if keyword in self.headers:
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                       keyword == "zlo zhi":
                        f.write("{0:f} {1:f} {2}\n".format(
                            self.headers[keyword][0],
                            self.headers[keyword][1], keyword))
                    elif keyword == "xy xz yz":
                        f.write("{0:f} {1:f} {2:f} {3}\n".format(
                            self.headers[keyword][0],
                            self.headers[keyword][1],
                            self.headers[keyword][2], keyword))
                    else:
                        f.write("{0:d} {1}\n".format(self.headers[keyword],
                                                   keyword))
            for pair in skeywords:
                keyword = pair[0]
                if keyword in self.sections:
                    f.write("\n{}\n\n".format(keyword))
                    for line in self.sections[keyword]:
                        f.write(line)

    def extract_nonpol(self):
        """extract atom and bond info from nonpolarizable data"""

        # extract atom IDs
        missinglabels = False
        for line in self.sections['Masses']:
            tok = line.split()
            if len(tok) < 4:
                print("error: missing type for atom ID " + tok[0] +
                      " in Masses section")
                missinglabels = True
                continue
            atomtype = {}
            atomtype['id'] = int(tok[0])
            atomtype['m'] = float(tok[1])
            atomtype['type'] = tok[3]
            self.atomtypes.append(atomtype)

        if missinglabels:
            sys.exit(0)

        # extract atom registers
        for line in self.sections['Atoms']:
            tok = line.split()
            atom = {}
            atom['n'] = int(tok[0])
            atom['mol'] = int(tok[1])
            atom['id'] = int(tok[2])
            atom['q'] = float(tok[3])
            atom['x'] = float(tok[4])
            atom['y'] = float(tok[5])
            atom['z'] = float(tok[6])
            #atom['note'] = ''.join([s + ' ' for s in tok[7:]]).strip()
            atom['note'] = ' '.join(tok[7:])
            self.atoms.append(atom)
            self.idmap[atom['n']] = atom

        if 'Velocities' in self.sections:
            for line in self.sections['Velocities']:
                tok = line.split()
                atom = self.idmap[int(tok[0])]
                atom['vx'] = float(tok[1])
                atom['vy'] = float(tok[2])
                atom['vz'] = float(tok[3])

    def polarize(self, drude):
        """add Drude particles"""

        if 'Pair Coeffs' in self.sections:
            raise RuntimeError("cannot polarize a data with Pair Coeffs")

        self.extract_nonpol()

        natom = self.headers['atoms']
        nbond = self.headers['bonds']
        nattype = self.headers['atom types']
        nbdtype = self.headers['bond types']

        # create new atom types (IDs) for Drude particles and modify cores
        newattypes = []
        for att in self.atomtypes:
            att['dflag'] = 'n'
            for ddt in drude.types:
                if ddt['type'] == att['type']:
                    nattype += 1
                    newid = {}
                    newid['id'] = ddt['id'] = nattype
                    newid['m'] = ddt['dm']
                    att['m'] -= ddt['dm']
                    # label drude particles and cores
                    att['dflag'] = 'c'
                    newid['dflag'] = 'd'
                    newid['type'] = att['type'] + ' DP'
                    att['type'] += ' DC'
                    ddt['type'] += ' DC'
                    newattypes.append(newid)
                    break

        self.headers['atom types'] += len(newattypes)
        self.sections['Masses'] = []
        for att in self.atomtypes + newattypes:
            self.sections['Masses'].append(massline(att))

        # create new bond types for core-drude bonds
        newbdtypes = []
        for att in self.atomtypes:
            for ddt in drude.types:
                if ddt['type'] == att['type']:
                    nbdtype += 1
                    newbdtype = {}
                    newbdtype['id'] = ddt['bdid'] = nbdtype
                    newbdtype['k'] = ddt['k']
                    newbdtype['r0'] = 0.0
                    newbdtype['note'] = '# ' + ddt['type'] + '-DP'
                    newbdtypes.append(newbdtype)
                    break

        self.headers['bond types'] += len(newbdtypes)
        for bdt in newbdtypes:
            self.sections['Bond Coeffs'].append(bdtline(bdt))

        # create new atoms for Drude particles and new bonds with their cores
        random.seed(123)
        newatoms = []
        newbonds = []
        for atom in self.atoms:
            atom['dflag'] = ''           # [c]ore, [d]rude, [n]on-polarizable
            atom['dd'] = 0               # partner drude or core
            for att in self.atomtypes:
                if att['id'] == atom['id']:
                    break
            for ddt in drude.types:
                if ddt['type'] == att['type']:
                    natom += 1
                    newatom = deepcopy(atom)
                    newatom['n'] = natom
                    self.idmap[natom] = newatom
                    newatom['id'] = ddt['id']
                    newatom['q'] = ddt['dq']
                    newatom['note'] = atom['note']
                    if '#' not in newatom['note']:
                        newatom['note'] += ' #'
                    newatom['note'] += ' DP'
                    newatom['dflag'] = 'd'
                    newatom['dd'] = atom['n']

                    # avoid superposition of cores and Drudes
                    newatom['x'] += 0.1 * (random.random() - 0.5)
                    newatom['y'] += 0.1 * (random.random() - 0.5)
                    newatom['z'] += 0.1 * (random.random() - 0.5)
                    if 'Velocities' in self.sections:
                        newatom['vx'] = atom['vx']
                        newatom['vy'] = atom['vy']
                        newatom['vz'] = atom['vz']

                    newatoms.append(newatom)
                    atom['q'] -= ddt['dq']
                    atom['dflag'] = 'c'
                    atom['dd'] = natom
                    if '#' not in atom['note']:
                        atom['note'] += ' #'
                    atom['note'] += ' DC'

                    nbond += 1
                    newbond = {}
                    newbond['n'] = nbond
                    newbond['id'] = ddt['bdid']
                    newbond['i'] = atom['n']
                    newbond['j'] = newatom['n']
                    newbond['note'] = '# ' + ddt['type'] + '-DP'
                    newbonds.append(newbond)
                    break

        self.headers['atoms'] += len(newatoms)
        self.headers['bonds'] += len(newbonds)
        self.sections['Atoms'] = []
        for atom in self.atoms + newatoms:
            self.sections['Atoms'].append(atomline(atom))
        for bond in newbonds:
            self.sections['Bonds'].append(bondline(bond))
        if 'Velocities' in self.sections:
            self.sections['Velocities'] = []
            for atom in self.atoms + newatoms:
                self.sections['Velocities'].append(velline(atom))

        # update list of atom IDs
        for att in newattypes:
            self.atomtypes.append(att)


    def extract_pol(self, drude):
        """extract atom, drude, bonds info from polarizable data"""

        # extract atom IDs
        for line in self.sections['Masses']:
            tok = line.split()
            atomtype = {}
            atomtype['id'] = int(tok[0])
            atomtype['m'] = float(tok[1])
            if len(tok) >= 4:
                atomtype['type'] = tok[3]
                atomtype['dflag'] = 'n'
                if tok[-1] == "DC":
                    atomtype['dflag'] = 'c'
                elif tok[-1] == "DP":
                    atomtype['dflag'] = 'd'
                print(atomtype['dflag'])
            else:
                raise RuntimeError("comments in Masses section required "\
                                   "to identify cores (DC) and Drudes (DP)")
            self.atomtypes.append(atomtype)
                        
        # extract bond type data
        for line in self.sections['Bond Coeffs']:
            tok = line.split()
            bondtype = {}
            bondtype['id'] = int(tok[0])
            bondtype['k'] = float(tok[1])
            bondtype['r0'] = float(tok[2])
            bondtype['note'] = ''.join([s + ' ' for s in tok[3:]]).strip()
            self.bondtypes.append(bondtype)

        # extract atom registers
        for line in self.sections['Atoms']:
            tok = line.split()
            atom = {}
            atom['n'] = int(tok[0])
            atom['mol'] = int(tok[1])
            atom['id'] = int(tok[2])
            atom['q'] = float(tok[3])
            atom['x'] = float(tok[4])
            atom['y'] = float(tok[5])
            atom['z'] = float(tok[6])
            # atom['note'] = ''.join([s + ' ' for s in tok[7:-1]]).strip()
            if tok[-1] == 'DC':
              atom['note'] = ' '.join(tok[7:-1])
            else:
              atom['note'] = ' '.join(tok[7:])
            self.atoms.append(atom)
            self.idmap[atom['n']] = atom

        if 'Velocities' in self.sections:
            for line in self.sections['Velocities']:
                tok = line.split()
                atom = self.idmap[int(tok[0])]
                atom['vx'] = float(tok[1])
                atom['vy'] = float(tok[2])
                atom['vz'] = float(tok[3])

        # extract bond data
        for line in self.sections['Bonds']:
            tok = line.split()
            bond = {}
            bond['n'] = int(tok[0])
            bond['id'] = int(tok[1])
            bond['i'] = int(tok[2])
            bond['j'] = int(tok[3])
            bond['note'] = ''.join([s + ' ' for s in tok[4:]]).strip()
            self.bonds.append(bond)

        if 'Velocities' in self.sections:
            for line in self.sections['Velocities']:
                tok = line.split()
                atom = self.idmap[int(tok[0])]
                atom['vx'] = float(tok[1])
                atom['vy'] = float(tok[2])
                atom['vz'] = float(tok[3])


    def depolarize(self, drude):
        """remove Drude particles"""

        self.extract_pol(drude)

        atom_tp_map = {}
        bond_tp_map = {}
        atom_map = {}
        bond_map = {}
        q = {}
        atom_tp = {}
        m = {}

        for att in self.atomtypes:
            if att['dflag'] != 'd':
                atom_tp_map[att['id']] = len(atom_tp_map) + 1
            m[att['id']] = att['m']
        print(atom_tp_map)
        for atom in self.atoms:
            if atom['id'] in atom_tp_map:
                atom_map[atom['n']] = len(atom_map) + 1
            q[atom['n']] = atom['q']
            atom_tp[atom['n']] = atom['id']
        for bond in self.bonds:
            if bond['i'] in atom_map and bond['j'] in atom_map:
                bond_map[bond['n']] = len(bond_map) + 1
                if bond['id'] not in bond_tp_map:
                    bond_tp_map[bond['id']] = len(bond_tp_map) + 1
            else:
                if bond['i'] in atom_map:
                    q[bond['i']] += q[bond['j']]
                    if atom_tp[bond['j']] in m:
                        m[atom_tp[bond['i']]] += m.pop(atom_tp[bond['j']])
                else:
                    q[bond['j']] += q[bond['i']]
                    if atom_tp[bond['i']] in m:
                        m[atom_tp[bond['j']]] += m.pop(atom_tp[bond['i']])

        size = len(self.atomtypes)
        for iatom_tp in reversed(range(size)):
            att = self.atomtypes[iatom_tp]
            if att['id'] not in atom_tp_map:
                del self.atomtypes[iatom_tp]
            else:
                att['m'] = m[att['id']]
                att['id'] = atom_tp_map[att['id']]

        size = len(self.bondtypes)
        for ibond_tp in reversed(range(size)):
            bdt = self.bondtypes[ibond_tp]
            if bdt['id'] not in bond_tp_map:
                del self.bondtypes[ibond_tp]
            else:
                bdt['id'] = bond_tp_map[bdt['id']]

        size = len(self.atoms)
        for iatom in reversed(range(size)):
            atom = self.atoms[iatom]
            if atom['n'] not in atom_map:
                del self.atoms[iatom]
            else:
                atom['q'] = q[atom['n']]
                atom['n'] = atom_map[atom['n']]
                atom['id'] = atom_tp_map[atom['id']]

        size = len(self.bonds)
        for ibond in reversed(range(size)):
            bond = self.bonds[ibond]
            if bond['n'] not in bond_map:
                del self.bonds[ibond]
            else:
                bond['n'] = bond_map[bond['n']]
                bond['id'] = bond_tp_map[bond['id']]
                bond['i'] = atom_map[bond['i']]
                bond['j'] = atom_map[bond['j']]

        self.sections['Atoms'] = []
        for atom in self.atoms:
            self.sections['Atoms'].append(atomline(atom))
        self.headers['atoms'] = len(self.atoms)
        self.sections['Masses'] = []
        for att in self.atomtypes:
            self.sections['Masses'].append(massline(att))
        self.headers['atom types'] = len(self.atomtypes)
        self.sections['Bonds'] = []
        for bond in self.bonds:
            self.sections['Bonds'].append(bondline(bond))
        self.headers['bonds'] = len(self.bonds)
        self.sections['Bond Coeffs'] = []
        for bdt in self.bondtypes:
            self.sections['Bond Coeffs'].append(bdtline(bdt))
        self.headers['bond types'] = len(self.bondtypes)

        if 'Velocities' in self.sections:
            self.sections['Velocities'] = []
            for atom in self.atoms:
                self.sections['Velocities'].append(velline(atom))


    def writepairfile (self, pairfile, drude, thole, att_id):
        """print pair_style thole"""

        with open(pairfile, "w") as f:
            f.write("# interactions involving Drude particles\n")
            f.write("pair_coeff    * {0:3d}* coul/long/cs\n".format(att_id))

            f.write("# Thole damping if more than 1 Drude per molecule\n")
            # Thole parameters for I,J pairs
            ifound = False
            for atti in self.atomtypes:
                itype = atti['type'].split()[0]
                for ddt in drude.types:
                    dtype = ddt['type'].split()[0]
                    if dtype == itype:
                        alphai = ddt['alpha']
                        tholei = ddt['thole']
                        ifound = True
                        break
                jfound = False
                for attj in self.atomtypes:
                    if attj['id'] < atti['id']:
                        continue
                    jtype = attj['type'].split()[0]
                    for ddt in drude.types:
                        dtype = ddt['type'].split()[0]
                        if dtype == jtype:
                            alphaj = ddt['alpha']
                            tholej = ddt['thole']
                            jfound = True
                            break
                    if ifound and jfound:
                        alphaij = (alphai * alphaj)**0.5
                        tholeij = (tholei + tholej) / 2.0
                        if tholeij == thole:
                            f.write("pair_coeff {0:4} {1:4} thole "\
                                  "{2:7.3f}\n".format(atti['id'], attj['id'],
                                                      alphaij))
                        else:
                            f.write("pair_coeff {0:4} {1:4} thole {2:7.3f} "\
                                  "{3:7.3f}\n".format(atti['id'],attj['id'],
                                                        alphaij, tholeij))
                    jfound = False
                ifound = False
        
    def concatenatepairfile (self, inpfile, outpfile, pairfile):
        with open(outpfile, 'wb') as wfd:
            for f in [inpfile, pairfile]:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
    
    def lmpscript(self, drude, outdfile, inpfile, outpfile, thole = 2.6, cutoff = 12.0):
        """print lines for input script, including pair_style thole"""

        inpstack = "in-p.lmp"
        pairfile = "pair-drude.lmp"
        
        dfound = False
        for att in self.atomtypes:
            if att['dflag'] == 'd':
                dfound = True
                break
        if not dfound:
            print("No polarizable atoms found.")
            return

        print("Commands for the LAMMPS input script to handle "\
              "polarization written to " + inpstack)

        with open(inpstack, 'x') as f:
            f.write("# Commands to include in the LAMMPS input stack\n\n")

            f.write("# adapt the pair_style command as needed\n")
            f.write("pair_style hybrid/overlay [...] coul/long/cs {0:.1f} "\
                  "thole {1:.3f} {0:.1f}\n\n".format(cutoff, thole))

            f.write("# new data file with Drude oscillators added\n")
            f.write("read_data {0}\n\n".format(outdfile))

            f.write("# read pair interactions involving Drude particles\n")
            f.write("# Thole damping recommended if more than 1 Drude "\
                    "per molecule\n")
            f.write("include {0}\n\n".format(pairfile))

            self.writepairfile(pairfile, drude, thole, att['id'])
            self.concatenatepairfile (inpfile, outpfile, pairfile)

            f.write("# convenient atom groups (for shake, thermostats...)\n")
            gatoms = gcores = gdrudes = ""
            for att in self.atomtypes:
                if att['dflag'] != 'd':
                    gatoms += " {0}".format(att['id'])
                if att['dflag'] == 'c':
                    gcores += " {0}".format(att['id'])
                if att['dflag'] == 'd':
                    gdrudes += " {0}".format(att['id'])
            f.write("group ATOMS type" + gatoms + "\n")
            f.write("group CORES type" + gcores + "\n")
            f.write("group DRUDES type" + gdrudes + "\n\n")

            f.write("# identify each atom type: "\
                    "[C]ore, [D]rude, [N]on-polarizable\n")
            drudetypes = ""
            for att in self.atomtypes:
                drudetypes += " {0}".format(att['dflag'].upper())
            f.write("fix DRUDE all drude" + drudetypes + "\n\n")

            f.write("# store velocity information of ghost atoms\n")
            f.write("comm_modify vel yes\n\n")

            f.write("variable TK equal 300.0\n")
            f.write("variable TDRUDE equal 1.0\n")
            f.write("variable PBAR equal 1.0\n\n")
        
            f.write("# temperature-grouped multiple Nose-Hoover thermostats "\
                    " and barostat\n")
            f.write("fix TSTAT all tgnpt/drude temp ${TK} ${TK} 100 "\
                    "${TDRUDE} 20 iso ${PBAR} ${PBAR} 1000\n\n")

            f.write("# output the temperatures of molecular COM, "\
                    "COM of DC-DP, and DP\n")
            f.write("thermo_style custom step [...] "\
                    "f_TSTAT[1] f_TSTAT[2] f_TSTAT[3]\n\n")

            i = drudetypes.find("D")

            f.write("# write Drude particles to dump file\n")
            f.write("dump_modify ... element ... " + drudetypes[i:] + "\n\n")

            f.write("# ATTENTION!\n")
            f.write("#  * read_data may need 'extra/special/per/atom' keyword, "
                  "LAMMPS will exit with a message\n")
            f.write("#  * if using fix shake the group-ID must not include "
              "Drude particles; use group ATOMS\n")
            f.write("#  * give all I<=J pair interactions, no mixing\n")
            f.write("#  * pair style coul/long/cs from CORESHELL package is "\
                  "used for interactions of DP;\n")
            f.write("#    alternatively pair lj/cut/thole/long could be used "\
                  "avoiding hybrid/overlay and\n")
            f.write("#    allowing mixing; see doc pages.\n")

# --------------------------------------

kcal =  4.184                           # kJ
eV   = 96.485                           # kJ/mol
fpe0 =  0.000719756                     # (4 Pi eps0) in e^2/(kJ/mol A)
d_alpha_H = 0.323                       # default polarisability of H in A^3, PCCP 20(2018)10992

# tolerance to deduce bonds from input configuration
BondTol = 0.25                          # Angstrom

class Drude(object):
    """specification of drude oscillator types"""

    def __init__(self, sys, drudefile, polar = '', positive = False, metal = False):
        self.types = []
        self.alpha_H = 0.0

        with open(drudefile, "r") as f:
            # search for polarisability of hydrogen in force filed file
            h_found = False
            for line in f:
                line = line.strip()
                if line.startswith('#') or len(line) == 0:
                    continue
                if line.startswith('H'):
                    tok = line.split()
                    self.alpha_H = float(tok[4])
                    h_found = True
                    break

            if not h_found:
                self.alpha_H = d_alpha_H
                print('warning: polarisability of H atom is not found. Default value is {0:5.3f} (PCCP 20(2018)10992). '.format(d_alpha_H))

            # read all polarisable atoms from the beginning
            f.seek(0)
            for line in f:
                line = line.strip()
                if line.startswith('#') or len(line) == 0:
                    continue
                tok = line.split()

                #skip H
                if tok[0].startswith('H'):
                    continue

                drude = {}
                drude['type'] = tok[0]
                drude['dm'] = float(tok[1])
                dq = float(tok[2])
                k = float(tok[3])
                alpha = float(tok[4]) # polarisability without H
                drude['thole'] = float(tok[5])

                # check if there are H arrached to the atom; if yes, merge their polarisability
                alpha += self.alpha_H*sys.countH(drude['type'])
                drude['alpha'] = alpha
                
                #print(drude['type'],alpha,sys.countH(drude['type']),alpha + self.alpha_H*sys.countH(drude['type']))

                if polar == 'q':
                    dq = (fpe0 * k * alpha)**0.5
                elif polar == 'k':
                    k = dq*dq / (fpe0 * alpha)

                if positive:
                    drude['dq'] = abs(dq)
                else:
                    drude['dq'] = -abs(dq)

                if metal:
                    drude['k'] = k / (2.0 * eV)
                else:
                    drude['k'] = k / (2.0 * kcal)

                self.types.append(drude)


# --------------------------------------

class System(object):
    """simulated system"""

    def __init__(self, strfiles):
        self.mols = []
        try:
            for sfile in strfiles:
                name = sfile.split('.')[0].strip().lower()
                if len([i for i in strfiles if name in i]) > 1:
                    raise Exception('  error: repeating file name ' + name)

                m = mol(sfile)
                self.mols.append(m)
                
        except Exception as e:
            print(e)
            sys.exit(1)
    
    def countH(self, dtype):
        dr_atoms = []
        nH0 = 0
        try:
            for m in self.mols:
                dr_atoms.extend([x for x in m.atoms if x.name == dtype])
            if len(dr_atoms) > 0:
                nH0 = dr_atoms[0].nH
                if next((x for x in dr_atoms if x.nH != nH0), None) is not None:
                    raise Exception('  error: atom type %s has different number of hydrogens attached in structure files' % dtype)
            return nH0

        except Exception as e:
            print(e)
            sys.exit(1)
        
# --------------------------------------

class forcefield(object):
    """force field parameter database"""

    def __init__(self, filename):
        self.filename = filename
        self.atoms = []
        self.bonds = []

        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.strip() == '':
                        continue
                
                    if line.lower().startswith('atom'):
                        section = 'atoms'
                        continue
                    elif line.lower().startswith('bond'):
                        section = 'bonds'
                        continue
                    elif line.lower().startswith('angl'):
                        section = 'angles'
                        continue
                    elif line.lower().startswith('dihe'):
                        section = 'dihedrals'
                        continue
                    elif line.lower().startswith('impro'):
                        section = 'improper'
                        continue

                    tok = line.strip().split()
                    if section == 'atoms':
                        name = tok[0]
                        atype = tok[1]
                        self.atoms.append(atom(name))
                        self.atoms[len(self.atoms)-1].atype = tok[1]

                    elif section == 'bonds':
                        iatp = tok[0]
                        jatp = tok[1]
                        r = float(tok[3])
                        self.bonds.append(bond(iatp,jatp,r))
        
        except IOError:
            print('  error: force field file '+ filename +'not found')
            sys.exit(1)                

# --------------------------------------

class atom(object):
    """atom in a molecule or in a force field"""

    def __init__(self, name):
        self.name = name
        self.atype = ''
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.nH = 0

    def __str__(self):
        return '%s %s \t %d H' %(self.name, self.atype, self.nH)

# --------------------------------------

class bond(object):
    """covalent bond in a molecule or in a force field"""

    def __init__(self, i , j , r = 0.0):
        self.i = i
        self.j = j
        self.r = r

    def checkval(self, dist):
        delta = abs(dist - self.r)
        return delta < BondTol

    def __str__(self):
        return '%s-%s %6.4f' % (self.i.name, self.j.name, self.r)

# --------------------------------------

class zmat(object):
    """z-matrix representing a molecule, read from .zmat file"""

    def __init__(self, filename):
        self.zatom = []
        self.connect = []
        self.improper = []
        
        with open(filename, 'r') as f:

            # read molecule name
            line = f.readline()
            while line.strip().startswith('#'):
                line = f.readline()
            self.name = line.strip()

            #read z-matrix
            line = f.readline()
            while line.strip().startswith('#') or line.strip() == '':
                line = f.readline()
            
            tok = line.strip().split()
            if len(tok) > 1:   # there can be line numbers
                shift = 1
            else:
                shift = 0

            variables = False
            while line and not line.strip().lower().startswith('var'):
                tok = line.strip().split()
                if len(tok) == 0:
                    break
                name = tok[shift]
                ir = ia = id = 0
                r = a = d = 0.0
                rvar = avar = dvar = ''
                if (len(tok) - shift) > 1:
                    ir = int(tok[shift+1])
                    if tok[shift+2][0].isalpha():
                        rvar = tok[shift+2]
                        variables = True
                    else:
                        r = float(tok[shift+2])
                    if (len(tok) - shift) > 3:
                        ia = int(tok[shift+3])
                        if tok[shift+4][0].isalpha():
                            avar = tok[shift+4]
                            variables = True
                        else:
                            a = float(tok[shift+4])
                        if (len(tok) - shift) > 5:
                            id = int(tok[shift+5])
                            if tok[shift+6][0].isalpha():
                                dvar = tok[shift+6]
                                variables = True
                            else:
                                d = float(tok[shift+6])
                zatom = {'name': name,
                        'ir': ir, 'rvar': rvar, 'r': r,
                        'ia': ia, 'avar': avar, 'a': a,
                        'id': id, 'dvar': dvar, 'd': d}
                self.zatom.append(zatom)
                line = f.readline()
                
            # read variables
            if variables:
                if line.strip().lower().startswith('var') or line.strip() == '':
                    line = f.readline()
                while line:
                    tok = line.strip().split('=')
                    if len(tok) < 2:
                        break
                    key = tok[0].strip()
                    val = float(tok[1])
                    for rec in self.zatom:
                        if rec['rvar'] == key:
                            rec['r'] = val
                        if rec['avar'] == key:
                            rec['a'] = val
                        if rec['dvar'] == key:
                            rec['d'] = val
                    line = f.readline()
                        
            # read connects, improper, force field file
            self.ff = ''
            self.guessconnect = False
            while line:
                if line.strip().startswith('#') or line.strip() == '':
                    line = f.readline()
                    continue
                tok = line.strip().split()
                if tok[0] == 'reconnect':
                    self.guessconnect = True
                if tok[0] == 'connect':
                    atomi = int(tok[1])
                    atomj = int(tok[2])
                    self.connect.append([atomi, atomj])
                elif tok[0] == 'improper':
                    atomi = int(tok[1])
                    atomj = int(tok[2])
                    atomk = int(tok[3])
                    atoml = int(tok[4])
                    self.improper.append([atomi, atomj, atomk, atoml])
                else:
                    self.ff = tok[0]
                line = f.readline()

# --------------------------------------

class mol(object):
    """molecule"""

    def __init__(self, filename):
        self.atoms = []
        self.bonds = []
        self.nmol = 0
        
        try:
            with open(filename, 'r'):
                self.filename = filename
            ext = filename.split('.')[-1].strip().lower()
            if ext == 'zmat':
                self.fromzmat(filename)
            elif ext == 'mol':
                self.frommdlmol(filename)
            elif ext == 'xyz':
                self.fromxyz(filename)
            elif ext == 'pdb':
                self.frompdb(filename)
            else:
                print('  error: unsupported molecule file extension')
                sys.exit(1)
        except IOError:
            print('  error: molecule file not found')
            sys.exit(1)
    
    def fromzmat(self, filename):
        z = zmat(filename)
        self.name = z.name
        self.guessconnect = z.guessconnect
        self.ff_file = z.ff
        for zat in z.zatom:
            self.atoms.append(atom(zat['name']))
        if not self.guessconnect:
            self.atomtypes()
            for i in range(1, len(z.zatom)):
                self.bonds.append(bond( self.atoms[i],  self.atoms[z.zatom[i]['ir'] - 1], z.zatom[i]['r']))
                if self.atoms[i].name.startswith('H'):
                    self.atoms[z.zatom[i]['ir'] - 1].nH += 1
                elif self.atoms[z.zatom[i]['ir'] - 1].name.startswith('H'):
                    self.atoms[i].nH += 1
            for cn in z.connect:
                self.bonds.append(bond( self.atoms[cn[0] - 1],  self.atoms[cn[1] - 1]))
                if self.atoms[cn[0] - 1].name.startswith('H'):
                    self.atoms[cn[1] - 1].nH += 1
                elif self.atoms[cn[1] - 1].name.startswith('H'):
                    self.atoms[cn[0] - 1].nH += 1
        else:
            self.connectivity()

        return self

    def frommdlmol(self, filename):
        with open(filename, 'r') as f:
            tok = f.readline().strip().split()
            self.name = tok[0]            # molecule name
            self.guessconnect = False
            if len(tok) > 1:              # and eventually ff file
                self.ff_file = tok[1]
                if len(tok) > 2:
                    if tok[2].startswith('rec'):
                        self.guessconnect = True
            else:
                raise Exception('  error: force field file name not found in ' + filename)
            f.readline()                  # program/date info
            line = f.readline().strip()   # comment (eventually ff file)
            if line and not line.startswith('#'):
                tok = line.split()
                self.ff_file = tok[0]
                if len(tok) > 1:
                    if tok[1].startswith('rec'):
                        self.guessconnect = True
            tok = f.readline().strip().split()   # counts line
            natom = int(tok[0])
            nbond = int(tok[1])
            
            self.atoms = [None] * natom
            for i in range(natom):
                tok = f.readline().strip().split()
                self.atoms[i] = atom(tok[3])
                self.atoms[i].x = float(tok[0])
                self.atoms[i].y = float(tok[1])
                self.atoms[i].z = float(tok[2])

            if not self.guessconnect:
                self.atomtypes()
                self.bonds = [None] * nbond
                for k in range(nbond):
                    tok = f.readline().strip().split()
                    i = int(tok[0]) - 1
                    j = int(tok[1]) - 1
                    ati =  self.atoms[i]
                    atj =  self.atoms[j]
                    self.bonds[k] = bond(ati, atj)

                    if self.atoms[i].name.startswith('H'):
                        self.atoms[j].nH += 1
                    elif self.atoms[j].name.startswith('H'):
                        self.atoms[i].nH += 1
            else:
                self.connectivity()

        return self


    def frompdb(self, filename):
        with open(filename, 'r') as f:
            line = f.readline()
            while 'COMPND' not in line:
                line = f.readline()
            tok = line.strip().split()
            if len(tok) < 3:
                raise Exception('  error: force field file name not found in ' + filename)
            self.name = tok[1]
            self.ff_file = tok[2]
            line = f.readline()
            while not (line.startswith('ATOM') or line.startswith('HETATM')):
                line = f.readline()
            self.atoms = []
            i = 0
            while 'ATOM' in line or 'HETATM' in line:
                tok = line.strip().split()
                atname = tok[2]
                self.atoms.append(atom(atname))
                self.atoms[i].x = float(tok[5])
                self.atoms[i].y = float(tok[6])
                self.atoms[i].z = float(tok[7])
                i += 1
                line = f.readline()
            self.connectivity()

        return self

    def fromxyz(self, filename):
        try:
            with open(filename, 'r') as f:
                natom = int(f.readline().strip())
                self.atoms = [None] * natom
                tok = f.readline().strip().split()
                self.name = tok[0]            # molecule name
                if len(tok) > 1:              # and eventually ff file
                    self.ff_file = tok[-1]
                else:
                    raise Exception('  error: force field file name not found in ' + filename)
                for i in range(natom):
                    tok = f.readline().strip().split()
                    self.atoms[i] = atom(tok[0])
                    self.atoms[i].x = float(tok[1])
                    self.atoms[i].y = float(tok[2])
                    self.atoms[i].z = float(tok[3])
                self.connectivity()
            return self

        except Exception as e:
            print(e)
            sys.exit(1) 

    def atomtypes(self):
        self.ff = forcefield(self.ff_file)
        natom = len(self.atoms)
            
        for k in range(0, natom):
            self.atoms[k].atype = next((x for x in self.ff.atoms if x.name == self.atoms[k].name), None).atype

    def connectivity(self):    
        """determine connectivity from bond distances in force field"""
        try:
            natom = len(self.atoms)
            self.atomtypes()
            
            for i in range(0, natom-1):
                for j in range(i+1, natom):
                    r = self.dist2atoms(self.atoms[i], self.atoms[j])
                    ffbd = next((x for x in self.ff.bonds if (x.i == self.atoms[i].atype and x.j == self.atoms[j].atype) or (x.i == self.atoms[j].atype and x.j == self.atoms[i].atype)), None)                  
                    if ffbd is not None and ffbd.checkval(r):
                        self.bonds.append(bond(self.atoms[i], self.atoms[j], r))
                                               
                        if self.atoms[i].name.startswith('H'):
                            self.atoms[j].nH += 1
                        elif self.atoms[j].name.startswith('H'):
                            self.atoms[i].nH += 1

        except Exception as e:
            print(e)
            sys.exit(1)

    def dist2atoms(self, ati, atj):
        dx = atj.x - ati.x
        dy = atj.y - ati.y
        dz = atj.z - ati.z
        return math.sqrt(dx*dx + dy*dy + dz*dz)

    
# --------------------------------------

def main():
    parser = argparse.ArgumentParser(description = usage,
             formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--ffdrude', default = 'alpha.ff',
                        help = 'Drude parameter file (default: alpha.ff)')
    parser.add_argument('-t', '--thole', type = float, default = 2.6,
                        help = 'Thole damping parameter (default: 2.6)')
    parser.add_argument('-c', '--cutoff', type = float, default = 12.0,
                        help = 'distance cutoff/A (default: 12.0)')
    parser.add_argument('-q', '--qcalc', action = 'store_true',
                        help = 'Drude charges calculated from polarisability '\
                        '(default: q value from parameter file)')
    parser.add_argument('-k', '--kcalc', action = 'store_true',
                        help = 'Drude force constants calculated from '\
                        'polarisability (default: k value from parameter file)')
    parser.add_argument('-p', '--positive', action = 'store_true',
                        help = 'Drude particles have positive charge '\
                        '(default: negative charge)')
    parser.add_argument('-m', '--metal', action = 'store_true',
                        help = 'LAMMPS metal units (default: real units)')
    parser.add_argument('-d', '--depolarize', action = 'store_true',
                        help = 'remove Drude dipole polarization from '\
                        'LAMMPS data file')
    parser.add_argument('strfiles', nargs='+',
                        help = 'infile1 [infile2 ...], '\
                        '. Use extension .zmat, .mol or .xyz')
    parser.add_argument('-id', '--indfile', default = 'data.lmp',
                        help = 'input LAMMPS data file (default: data.lmp)')
    parser.add_argument('-od', '--outdfile', default = 'data-p.lmp',
                        help = 'output LAMMPS data file (default: data-p.lmp)')
    parser.add_argument('-ip', '--inpfile', default = 'pair.lmp',
                        help = 'input LAMMPS pair file (default: pair.lmp)')
    parser.add_argument('-op', '--outpfile', default = 'pair-p.lmp',
                        help = 'output LAMMPS pair file (default: pair-p.lmp)')
    
    args = parser.parse_args()

    if args.qcalc:
        polar = 'q'
    elif args.kcalc:
        polar = 'k'
    else:
        polar = ''

    data = Data(args.indfile)
    sys = System(args.strfiles)
    drude = Drude(sys, args.ffdrude, polar, args.positive, args.metal)
    if not args.depolarize:
        data.polarize(drude)
        data.lmpscript(drude, args.outdfile, args.inpfile, args.outpfile, args.thole, args.cutoff)
    else:
        data.depolarize(drude)
    data.write(args.outdfile)

if __name__ == '__main__':
    main()
