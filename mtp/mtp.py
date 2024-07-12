#!/usr/bin/env python

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import FileIOCalculator, all_changes
from ase.units import Rydberg, Bohr, GPa
from ase.calculators.singlepoint import SinglePointCalculator

def dump_cfg(frames, filename, symbol_to_type, mode='w', clean_info=True):
    """
    Transfrom ASE trajectory to MLIP's format.
    frames: list of Atoms objects
    filename: configuration file in MLIP's format
    symbol_to_type: dict, map elements to integers
    mode: mode of `open` function
    clean_info: clean `info` of an Atoms object before dumping. If you want to use energy/force/stress in Atoms.info, set it as False
    """

    with open(filename, mode) as f:
        for atoms in frames:
            ret = ''
            ret += 'BEGIN_CFG\n'
            ret += 'Size\n{}\n'.format(len(atoms))
            pbc = atoms.get_pbc()
            try:
                cell = atoms.get_cell()[:]
                ret += 'Supercell\n'
                for axis in range(3):
                    # Consider PBC. Assume pbc is [1,0,0] [1,1,0] or [1,1,1]
                    if pbc[axis]:
                        ret += '{} {} {}\n'.format(*cell[axis])
                #ret += 'Supercell\n{} {} {}\n{} {} {}\n{} {} {}\n'.format(*cell[0], *cell[1], *cell[2])
            except:
                pass
            cartes = atoms.positions
            if clean_info:
                atoms.info = dict()
            try:
                atoms.info['forces'] = atoms.get_forces()
            except:
                pass
            if 'forces' in atoms.info:
                forces = atoms.info['forces']
                ret += 'AtomData: id type cartes_x cartes_y cartes_z fx fy fz\n'
                for i, atom in enumerate(atoms):
                    ret += '{} {} {} {} {} {} {} {}\n'.format(i + 1, symbol_to_type[atom.symbol], *cartes[i], *forces[i])
            else:
                ret += 'AtomData: id type cartes_x cartes_y cartes_z\n'
                for i, atom in enumerate(atoms):
                    ret += '{} {} {} {} {}\n'.format(i + 1, symbol_to_type[atom.symbol], *cartes[i])
            try:
                atoms.info['energy'] = atoms.get_potential_energy()
            except:
                pass
            if 'energy' in atoms.info:
                ret += 'Energy\n{}\n'.format(atoms.info['energy'])
            try:
                atoms.info['stress'] = atoms.get_stress()
            except:
                pass
            if 'stress' in atoms.info:
                stress = atoms.info['stress'] * atoms.get_volume() * -1.
                ret += 'PlusStress: xx yy zz yz xz xy\n{} {} {} {} {} {}\n'.format(*stress)
            if 'identification' in atoms.info:
                ret += 'Feature identification {}\n'.format(atoms.info['identification'])
            ret += 'END_CFG\n'
            f.write(ret)


#TODO
# different cell
# pbc

def reading_state(line, state):
    # change the reading state based on current state and line
    if 'BEGIN_CFG' in line:
        state = 'begin'
    elif 'Size' in line:
        state = 'size'
    elif 'Supercell' in line:
        state = 'cell'
    elif 'AtomData' in line:
        state = 'atom'
    elif 'Energy' in line:
        state = 'en'
    elif 'PlusStress' in line:
        state = 'stress'
    elif 'END_CFG' in line:
        state = 'end'

    return state

def load_cfg(filename, type_to_symbol, save_in_info=False):
    """
    Transform MLIP's configuration file to ASE's trajectory
    filename: configuration file in MLIP's format
    type_to_symbol: dict, map integers to elements
    save_in_info: save energy/force/stress in Atoms.info
    """
    frames = []
    state = 'no'
    with open(filename) as f:
        line = 'aaaaa'
        while line:
            state = reading_state(line, state)
            #print(line)
            #print(state)
            #print(' ')

            if state == 'no':
                pass

            if state == 'begin':
                cell = np.zeros((3, 3))
                positions = []
                symbols = []
                celldim = 0
                results = dict()

            if state == 'size':
                line = f.readline()
                natoms = int(line.split()[0])

            if state == 'cell':
                #for i in range(3):
                #    line = f.readline()
                # 0D systems have no Supercell field
                if 'Supercell' not in line:
                    for j in range(3):
                        cell[celldim, j] = float(line.split()[j])
                    celldim += 1

            if state == 'atom':
                has_force = False
                if 'fx' in line:
                    has_force = True
                    forces = []

                for _ in range(natoms):
                    line = f.readline()
                    fields = line.split()
                    symbols.append(type_to_symbol[int(fields[1])])
                    positions.append(list(map(float, fields[2: 5])))
                    if has_force:
                        forces.append(list(map(float, fields[5: 8])))
                if celldim == 1:
                    pbc = [1,0,0]
                    #cell[1,1] = 30
                    #cell[2,2] = 30
                elif celldim == 2:
                    pbc = [1,1,0]
                    #cell[2,2] = 30
                elif celldim == 3:
                    pbc = [1,1,1]
                else:
                    pbc = False
                atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=pbc)
                if has_force:
                    results['forces'] = np.array(forces)
                    #atoms.info['forces'] = np.array(forces)

            if state == 'en':
                line = f.readline()
                results['energy'] = float(line.split()[0])
                #atoms.info['energy'] = float(line.split()[0])

            if state == 'stress':
                line = f.readline()
                plus_stress = np.array(list(map(float, line.split())))
                # It is possible the cell is not pbc along 3 directions
                if atoms.get_pbc().all():
                    results['stress'] = -plus_stress / atoms.get_volume()
                    #results['pstress'] = results['stress'] / GPa
                    #atoms.info['stress'] = -plus_stress / atoms.get_volume()
                    #atoms.info['pstress'] = atoms.info['stress'] / GPa

            if state == 'end':
                if save_in_info:
                    for key, val in results.items():
                         atoms.info[key] = val
                else:
                    atoms.calc = SinglePointCalculator(atoms, **results)

                frames.append(atoms)
                state = 'no'

            line = f.readline()
            #if 'identification' in line:
            #    atoms.info['identification'] = int(line.split()[2])

    return frames


class MTPCalculator(FileIOCalculator):

    # The default command to run MTP
    #command = "mlp calc-efs pot.mtp ASE_IN.CFG ASE_OUT.CFG"
    implemented_properties = ['energy', 'forces', 'stress']
    #default_parameters = dict(
    #    exe='mlp',
    #    pot='pot.mtp',
    #    symbol2type={},
    #    )




    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='mtp_ase', atoms=None, **kwargs):
        """
        Required parameters:
        symbol2type: symbol-type mapping
        exe: exe file for MTP
        pot: potential file for MTP

        """

        # command to run MTP
        command = f"{kwargs['exe']} calc-efs {kwargs['pot']} ASE_IN.CFG ASE_OUT.CFG"

        super().__init__(restart=restart, ignore_bad_restart_file=ignore_bad_restart_file,
                                  label=label, atoms=atoms, command=command, profile=None, **kwargs)

        #self.command = "{} calc-efs {} ASE_IN.CFG ASE_OUT.CFG".format(self.parameters['exe'], self.parameters['pot'])

        # type-symbol mapping
        s2t = self.parameters['symbol2type']
        t2s = dict()
        for key, val in s2t.items():
            t2s[val] = key
        self.set(type2symbol=t2s)


    def write_input(self, atoms, properties = None, system_charges = None):
        """
        """
        atoms.info = {}
        dump_cfg([atoms], 'ASE_IN.CFG', self.parameters['symbol2type'])


    def read_results(self):
        """
        """
        outAts = load_cfg('ASE_OUT.CFG', self.parameters['type2symbol'], save_in_info=True)[0]
        #self.results = {"energy" : outAts.info['energy'], "forces" : outAts.info['forces'], "stress" : outAts.info['stress']}
        self.results = {"energy" : outAts.info['energy'], "forces" : outAts.info['forces']}
        if 'stress' in outAts.info.keys():
            self.results["stress"] = outAts.info['stress']
