# -*- coding: iso-8859-1 -*-
'''Python library to manage dl_poly files'''

import logging
import gzip
import math
import itertools
import numpy as np

import molecule

class Control(object):
  '''Class representing a CONTROL file'''

  def __init__(self, filename=None):

    self.header = ''

    self.steps = None
    self.timestep = None
    self.temp = None
    self.press = None
    self.stats = None
    self.restart = False

    self.lines = []

    if filename:
      self.read_from_file(filename)

  def write_to_file(self, filename='CONTROL'):

    f = open(filename, 'w')

    f.write(self.header + '\n\n')

    if self.steps:
      f.write('STEPS %d\n' % self.steps)
    if self.timestep:
      f.write('TIMESTEP %f\n' % self.timestep)
    if self.temp:
      f.write('TEMP %f\n' % self.temp)
    if self.press:
      f.write('PRESS %f\n' % self.press)
    if self.stats:
      f.write('STATS %d\n' % self.stats)
    if self.restart:
      f.write('RESTART\n')
    f.write('\n')
    for line in self.lines:
      f.write(line + '\n')

    f.write('FINISH\n')
    f.close()

  def read_from_file(self, filename='CONTROL'):

    f = open(filename, 'r')

    self.header = f.readline().strip()

    for line in f:
      fields = line.split()
      if len(fields) < 1:
        continue
      if fields[0].upper() == 'FINISH':
        break
      elif fields[0].upper() == 'STEPS':
        self.steps = int(fields[1])
      elif fields[0].upper() == 'TIMESTEP':
        self.timestep = float(fields[1])
      elif fields[0].upper().startswith('TEMP'):
        self.temp = float(fields[1])
      elif fields[0].upper().startswith('PRES'):
        self.press = float(fields[1])
      elif fields[0].upper() == 'STATS':
        self.stats = int(fields[1])
      elif fields[0].upper() == 'RESTART':
        self.restart = True
      else:
        self.lines.append(line)

    if not self.steps or not self.timestep:
      raise ValueError('File %s does not seem to be a valid DL_POLY CONTROL file' % filename)

    f.close()

  def __str__(self):
    return 'CONTROL DL_POLY %d steps / %f timestep / %f temp / %f press' \
      % (self.steps, self.timestep, self.temp, self.press)


class Config(object):
  '''Class representing a CONFIG file'''

  def __init__(self, filename=None):

    self.atoms = []
    self.header = None
    self.levcfg = None
    self.imcon = None
    self.natms = None
    self.engcfg = None
    self.cell = [None, None, None]

    if filename:
      self.read_from_file(filename)

  def remove_mol(self, field, name):

    nb2skip = 0
    nb2rm = 0
    for mol in field.mols:
      if mol['name'] == name:
        nb2rm = mol['nummols'] * len(mol['atoms'])
        break
      else:
        nb2skip += mol['nummols'] * len(mol['atoms'])

    self.remove_atoms(nb2rm, nb2skip)

  def remove_atoms(self, nb, skip=0):

    atoms2remove = []
    nb_skip = 0
    for atom in self.atoms:
      if nb_skip < skip:
        nb_skip += 1
        continue
      if len(atoms2remove) == nb:
        break
      atoms2remove.append(atom)

    for atom in atoms2remove:
      self.atoms.remove(atom)

    self.natms = self.natms - nb

    i = 1
    for atom in self.atoms:
      atom['iatm'] = i
      i += 1

  def write_to_file(self, filename='CONFIG'):

    f = open(filename, 'w')
    f.write(self.header + '\n')

    f.write('%10d%10d%10d%20f\n' % (self.levcfg, self.imcon, self.natms, self.engcfg))

    f.write('%20f%20f%20f\n' % self.cell[0])
    f.write('%20f%20f%20f\n' % self.cell[1])
    f.write('%20f%20f%20f\n' % self.cell[2])

    for atom in self.atoms:
      f.write('%8s%10d\n' % (atom['atmnam'], atom['iatm']))
      f.write('%20f%20f%20f\n' % (atom['xxx'], atom['yyy'], atom['zzz']))
      if self.levcfg > 0:
        f.write('%20f%20f%20f\n' % (atom['vxx'], atom['vyy'], atom['vzz']))
      if self.levcfg > 1:
        f.write('%20f%20f%20f\n' % (atom['fxx'], atom['fyy'], atom['fzz']))

    f.close()

  def read_from_file(self, filename='CONFIG'):

    self.atoms = []

    f = open(filename, 'r')

    # Record 1
    self.header = f.readline().strip()

    # Record 2
    try:
      fields = f.readline().split()
      self.levcfg = int(fields[0])
      self.imcon = int(fields[1])
      self.natms = int(fields[2])
      self.engcfg = float(fields[3])
    except ValueError:
      raise ValueError('File %s does not seem to be a valid DL_POLY CONFIG file' % filename)

    if self.imcon != 0:
      # Record 3
      fields = f.readline().split()
      self.cell[0] = (float(fields[0]), float(fields[1]), float(fields[2]))
      # Record 4
      fields = f.readline().split()
      self.cell[1] = (float(fields[0]), float(fields[1]), float(fields[2]))
      # Record 5
      fields = f.readline().split()
      self.cell[2] = (float(fields[0]), float(fields[1]), float(fields[2]))

    for _ in xrange(self.natms):
      atom = {}
      # Record i
      fields = f.readline().split()
      atom['atmnam'] = fields[0]
      atom['iatm'] = int(fields[1])

      # Record ii
      fields = f.readline().split()
      atom['xxx'] = float(fields[0])
      atom['yyy'] = float(fields[1])
      atom['zzz'] = float(fields[2])

      # Record iii
      if self.levcfg > 0:
        fields = f.readline().split()
        atom['vxx'] = float(fields[0])
        atom['vyy'] = float(fields[1])
        atom['vzz'] = float(fields[2])

      # Record iv
      if self.levcfg > 1:
        fields = f.readline().split()
        atom['fxx'] = float(fields[0])
        atom['fyy'] = float(fields[1])
        atom['fzz'] = float(fields[2])

      self.atoms.append(atom)

    f.close()

  def init_from_history_conf(self, history, cnf):

    self.header = 'Config from HISTORY timestep %d' % cnf['nstep']
    self.levcfg = history.keytrj
    self.imcon = history.imcon
    self.natms = history.natms
    # Don't appear in HISTORY
    self.engcfg = 0.0
    self.cell = cnf['cell']

    self.atoms = []

    i = 0
    for atom in cnf['atoms']:
      i += 1
      self_atom = {}

      self_atom['atmnam'] = atom['atmnam']
      self_atom['iatm'] = atom['iatm']

      self_atom['xxx'] = atom['xxx']
      self_atom['yyy'] = atom['yyy']
      self_atom['zzz'] = atom['zzz']

      if history.keytrj > 0:
        self_atom['vxx'] = atom['vxx']
        self_atom['vyy'] = atom['vyy']
        self_atom['vzz'] = atom['vzz']

      if history.keytrj > 1:
        self_atom['fxx'] = atom['fxx']
        self_atom['fyy'] = atom['fyy']
        self_atom['fzz'] = atom['fzz']

      self.atoms.append(self_atom)

  def __str__(self):
    return 'CONFIG DL_POLY %s box : %s/%s/%s nbatoms : %d' \
      % (self.header, self.cell[0], self.cell[1], self.cell[2], self.natms)


class Field(object):
  '''Class representing a Field object'''

  def __init__(self, filename=None):

    self.header = None
    self.units = None

    self.mols = []

    self.vdws = []

    self.nb_3_dihedral = 0
    self.nb_4_dihedral = 0

    if filename:
      self.read_from_file(filename)

  def get_mol(self, name):
    '''Get the mol field from a name'''

    for mol in self.mols:
      if mol['name'] == name:
        return mol

    return None

  def remove_one_mol(self, name):

    mol2remove = self.get_mol(name)

    if not mol2remove:
      raise Exception('Unknown mol %s' % name)

    self.mols.remove(mol2remove)

    # Get all remaining atom names to clean vdws
    atom_names = []
    for mol in self.mols:
      for atom in mol['atoms']:
        if atom['name'] not in atom_names:
          atom_names.append(atom['name'])

    # Clean vdw to remove unused ones
    vdws2remove = []
    for vdw in self.vdws:
      if vdw['atom_type1'] not in atom_names or vdw['atom_type2'] not in atom_names:
        vdws2remove.append(vdw)

    for vdw in vdws2remove:
      self.vdws.remove(vdw)

    return mol2remove

  def write_to_file(self, filename='FIELD', is_nonpatched=False, with_comment=False):

    f = open(filename, 'w')
    f.write(self.header + '\n')
    f.write('UNITS %s\n' % self.units)
    f.write('\n')

    f.write('MOLECULAR TYPES %d\n' % len(self.mols))

    for mol in self.mols:
      f.write('%s\n' % mol['name'])
      f.write('NUMMOLS %d\n' % mol['nummols'])

      f.write('ATOMS %d\n' % len(mol['atoms']))
      for atom in mol['atoms']:
        f.write('%s %f %f 1' % (atom['name'], atom['weight'], atom['chge']))
        # Default case, don't indicate no frozen
        if not atom['frozen'] and atom['igrp'] is None:
          f.write('\n')
        else:
          if atom['frozen']:
            f.write(' 1')
          else:
            f.write(' 0')
          if atom['igrp'] is not None:
            f.write(' %d\n' % atom['igrp'])
          else:
            f.write('\n')

      if len(mol['constraints']) > 0:
        f.write('CONSTRAINTS %d\n' % len(mol['constraints']))
        for constraint in mol['constraints']:
          f.write('%d %d %f' % (constraint['atom_i1'], constraint['atom_i2'], constraint['bondlength']))
          if with_comment:
            f.write(' # %s-%s\n' % (mol['atoms'][constraint['atom_i1'] - 1]['name'], mol['atoms'][constraint['atom_i2'] - 1]['name']))
          else:
            f.write('\n')

      if len(mol['bonds']) > 0:
        f.write('BONDS %d\n' % len(mol['bonds']))
        for bond in mol['bonds']:
          f.write('%s %d %d %s' % (bond['field_type'], bond['atom_i1'], bond['atom_i2'], ' '.join(['%f' % p for p in bond['par']])))
          if with_comment:
            f.write(' # %s-%s\n' % (mol['atoms'][bond['atom_i1'] - 1]['name'], mol['atoms'][bond['atom_i2'] - 1]['name']))
          else:
            f.write('\n')

      if len(mol['rigids']) > 0:
        f.write('RIGID UNITS %d\n' % len(mol['rigids']))
        for rigid in mol['rigids']:
          f.write('%d %s' % (len(rigid), ' '.join(['%d' % i for i in rigid])))
          if with_comment:
            f.write(' # %s\n' % ('-'.join([mol['atoms'][i - 1]['name'] for i in rigid])))
          else:
            f.write('\n')

      if len(mol['angles']) > 0:
        f.write('ANGLES %d\n' % len(mol['angles']))
        for angle in mol['angles']:
          f.write('%s %d %d %d %s' % (angle['field_type'], angle['atom_i1'], angle['atom_i2'], angle['atom_i3'], ' '.join(['%f' % p for p in angle['par']])))
          if with_comment:
            f.write(' # %s-%s-%s\n' % (mol['atoms'][angle['atom_i1'] - 1]['name'], mol['atoms'][angle['atom_i2'] - 1]['name'], mol['atoms'][angle['atom_i3'] - 1]['name']))
          else:
            f.write('\n')

      if len(mol['dihedrals']) > 0:
        f.write('DIHEDRALS %d\n' % len(mol['dihedrals']))
        for dihedral in mol['dihedrals']:
          f.write('%s %d %d %d %d %s ' % (dihedral['field_type'], dihedral['atom_i1'], dihedral['atom_i2'], dihedral['atom_i3'], dihedral['atom_i4'], ' '.join(['%f' % p for p in dihedral['par']])))
          if is_nonpatched:
            for _ in xrange(3 - len(dihedral['par'])):
              f.write('%f ' % 0.0)
          else:
            for _ in xrange(4 - len(dihedral['par'])):
              f.write('%f ' % 0.0)
          # 1-4 parameters
          f.write('%f %f' % (dihedral['elec_1_4'], dihedral['vdw_1_4']))
          if with_comment:
            f.write(' # %s-%s-%s-%s\n' % (mol['atoms'][dihedral['atom_i1'] - 1]['name'], mol['atoms'][dihedral['atom_i2'] - 1]['name'], mol['atoms'][dihedral['atom_i3'] - 1]['name'], mol['atoms'][dihedral['atom_i4'] - 1]['name']))
          else:
            f.write('\n')

      f.write('FINISH\n')

    f.write('VDW %d\n' % len(self.vdws))
    for vdw in self.vdws:
      f.write('%s %s %s %s\n' % (vdw['atom_type1'], vdw['atom_type2'], vdw['field_type'], ' '.join(['%f' % p for p in vdw['par']])))
    f.write('CLOSE\n')
    f.close()

  def read_from_file(self, filename='FIELD'):

    f = open(filename, 'r')
    self.header = f.readline().strip()
    self.units = f.readline().split()[1]

    line = f.readline()
    while not (line.upper().startswith('MOLECULAR TYPES') or line.upper().startswith('MOLECULES')):
      line = f.readline()
      # If EOF
      if line == '':
        raise ValueError('File %s does not seem to be a valid DL_POLY FIELD file' % filename)

    nummoltype = int(line.split()[-1])

    # For each molecule type
    for _ in xrange(nummoltype):
      mol = {}
      namemol = f.readline().strip()
      mol['name'] = namemol
      nummols = int(f.readline().split()[1])
      mol['nummols'] = nummols

      mol['atoms'] = []
      mol['constraints'] = []
      mol['bonds'] = []
      mol['rigids'] = []
      mol['angles'] = []
      mol['dihedrals'] = []

      # ATOMS
      # Read nb
      line = f.readline()
      while not line.upper().startswith('ATOMS'):
        line = f.readline()
        # If EOF
        if line == '':
          f.close()
          return
      remain_atoms = int(line.split()[1])

      # For each atom
      while remain_atoms > 0:
        fields = f.readline().split()
        name = fields[0]
        weight = float(fields[1])
        chge = float(fields[2])
        rept = int(fields[3])
        frozen = False
        igrp = None
        try:
          if len(fields) > 4 and int(fields[4]) > 0:
            frozen = True
          if len(fields) > 5:
            igrp = int(fields[5])
        except ValueError:
          pass

        remain_atoms -= rept

        for __ in xrange(rept):
          atom = {}
          atom['name'] = name
          atom['weight'] = weight
          atom['chge'] = chge
          atom['frozen'] = frozen
          atom['igrp'] = igrp

          mol['atoms'].append(atom)

      # Skip useless stuff
      while not line.upper().startswith('CONSTRAINTS') and \
            not line.upper().startswith('BONDS') and not line.upper().startswith('RIGID') and \
            not line.upper().startswith('ANGLES') and not line.upper().startswith('DIHEDRALS') and \
            not line.upper().startswith('FINISH'):
        line = f.readline()
        # If EOF
        if line == '':
          f.close()
          return

      # CONSTRAINTS
      if line.upper().startswith('CONSTRAINTS'):
        remain_consts = int(line.split()[1])

        # For each const
        while remain_consts > 0:
          const = {}
          fields = f.readline().split()
          const['atom_i1'] = int(fields[0])
          const['atom_i2'] = int(fields[1])
          const['bondlength'] = float(fields[2])

          remain_consts -= 1

          mol['constraints'].append(const)

      # Skip useless stuff
      while not line.upper().startswith('BONDS') and not line.upper().startswith('RIGID') and \
            not line.upper().startswith('ANGLES') and not line.upper().startswith('DIHEDRALS') and \
            not line.upper().startswith('FINISH'):
        line = f.readline()
        # If EOF
        if line == '':
          f.close()
          return

      # BONDS
      if line.upper().startswith('BONDS'):
        remain_bonds = int(line.split()[1])

        # For each bond
        while remain_bonds > 0:
          bond = {}
          fields = f.readline().split()
          bond['field_type'] = fields[0]
          bond['atom_i1'] = int(fields[1])
          bond['atom_i2'] = int(fields[2])
          bond['par'] = []
          for p in fields[3:]:
            # Avoid storing comment
            try:
              p_f = float(p)
            except ValueError:
              pass
            else:
              bond['par'].append(p_f)

          remain_bonds -= 1

          mol['bonds'].append(bond)

      # Skip useless stuff
      while not line.upper().startswith('RIGID') and \
            not line.upper().startswith('ANGLES') and not line.upper().startswith('DIHEDRALS') and \
            not line.upper().startswith('FINISH'):
        line = f.readline()
        # If EOF
        if line == '':
          f.close()
          return

      # RIGIDS
      if line.upper().startswith('RIGID'):
        remain_rigids = int(line.split()[-1])

        # For each rigid
        while remain_rigids > 0:
          bond = {}
          fields = f.readline().split()
          # Skip comments
          int_fields = []
          for field in fields[1:]:
            try:
              int_fields.append(int(field))
            except ValueError:
              pass
          mol['rigids'].append([int(field) for field in int_fields])

          remain_rigids -= 1

      # Skip useless stuff
      while not line.upper().startswith('ANGLES') and not line.upper().startswith('DIHEDRALS') and \
            not line.upper().startswith('FINISH'):
        line = f.readline()
        # If EOF
        if line == '':
          f.close()
          return

      # ANGLES
      if line.upper().startswith('ANGLES'):
        remain_angles = int(line.split()[1])

        # For each angle
        while remain_angles > 0:
          angle = {}
          fields = f.readline().split()
          angle['field_type'] = fields[0]
          angle['atom_i1'] = int(fields[1])
          angle['atom_i2'] = int(fields[2])
          angle['atom_i3'] = int(fields[3])
          angle['par'] = []
          for p in fields[4:]:
            # Avoid storing comment
            try:
              p_f = float(p)
            except ValueError:
              pass
            else:
              angle['par'].append(p_f)

          remain_angles -= 1

          mol['angles'].append(angle)

      # Skip useless stuff
      while not line.upper().startswith('DIHEDRALS') and not line.upper().startswith('FINISH'):
        line = f.readline()
        # If EOF
        if line == '':
          f.close()
          return

      # DIHEDRALS
      if line.upper().startswith('DIHEDRALS'):
        remain_dihedrals = int(line.split()[1])

        # For each dihedral
        while remain_dihedrals > 0:
          dihedral = {}
          fields = f.readline().split()
          dihedral['field_type'] = fields[0]
          dihedral['atom_i1'] = int(fields[1])
          dihedral['atom_i2'] = int(fields[2])
          dihedral['atom_i3'] = int(fields[3])
          dihedral['atom_i4'] = int(fields[4])
          dihedral['par'] = []
          dihedral['elec_1_4'] = 0.0
          dihedral['vdw_1_4'] = 0.0

          # Get the nb of remain parameters (without comments)
          nb_remain = 0
          for field in fields:
            try:
              _ = float(field)
            except ValueError:
              pass
            else:
              nb_remain += 1
          nb_remain -= 4

          # 3 params + 0 1_4
          if nb_remain == 3:
            params = fields[5:8]
            self.nb_3_dihedral += 1
          # 4 params + 0 1_4
          elif nb_remain == 4:
            params = fields[5:9]
            self.nb_4_dihedral += 1
          # 3 params + 2 1_4
          elif nb_remain == 5:
            dihedral['elec_1_4'] = float(fields[8])
            dihedral['vdw_1_4'] = float(fields[9])
            params = fields[5:8]
            self.nb_3_dihedral += 1
          # 4 params + 2 1_4
          elif nb_remain == 6:
            dihedral['elec_1_4'] = float(fields[9])
            dihedral['vdw_1_4'] = float(fields[10])
            params = fields[5:9]
            self.nb_4_dihedral += 1
          else:
            params = fields[5:]

          for p in params:
            # Avoid storing comment
            try:
              p_f = float(p)
            except ValueError:
              pass
            else:
              dihedral['par'].append(p_f)

          remain_dihedrals -= 1

          mol['dihedrals'].append(dihedral)

      self.mols.append(mol)

      # skip to the end of the molecule infos
      # TODO check some missing stuff
      if line.upper().startswith('FINISH'):
        continue

      while not f.readline().upper().startswith('FINISH'):
        pass

    self.vdws = []

    line = f.readline()
    while not line.upper().startswith('VDW'):
      line = f.readline()
      # If EOF
      if line == '':
        f.close()
        return

    # VDWS
    # Read nb
    # TODO check [0]
    remain_vdws = int(line.split()[1])

    # For each vdw
    while remain_vdws > 0:
      vdw = {}
      fields = f.readline().split()
      vdw['atom_type1'] = fields[0]
      vdw['atom_type2'] = fields[1]
      vdw['field_type'] = fields[2]
      vdw['par'] = [float(p) for p in fields[3:]]

      remain_vdws -= 1

      self.vdws.append(vdw)

    f.close()

  def __str__(self):
    return 'FIELD DL_POLY %s Units %s %s' % (self.header, self.units, self.mols)


class History(object):
  '''Class representing a History file reader'''

  def __init__(self, filename='HISTORY', skip=0, xyz_only=False, field=None, raw=True, numpy=False):

    self.filename = filename

    self.is_compact = False
    self.raw = raw
    self.numpy = numpy

    self.xyz_only = xyz_only

    self.field = field

    self.i = skip + 1

    self.atom_types = {}

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Get the first line header
      self.header = self.f.next()
    except IOError:
      self.f = open(filename, 'r')
      # Get the first line header
      self.header = self.f.next()

    # Record 2
    try:
      fields = self.f.next().split()
      self.keytrj = int(fields[0])
      self.imcon = int(fields[1])
      self.natms = int(fields[2])
    except (ValueError, IndexError):
      raise ValueError('File %s does not seem to be a valid DL_POLY HISTORY file' % self.filename)

    # Check if the file is in compact mode
    fields = self.f.next().split()
    if fields[0] != 'timestep':
      self.is_compact = True
    self.f.close()

    # Reopen the file to seek to the start

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Skip the useless first line
      self.f.next()
    except IOError:
      self.f = open(filename, 'r')
      # Skip the useless first line
      self.f.next()

    self.f.next()

    # Skip the compact infos
    if self.is_compact:
      s = itertools.islice(self.f, self.natms)
      self.raw_static = list(s)
    else:
      self.raw_static = None

    # Store the size of one conf (in lines)
    self.sizeconf = 1
    if self.imcon > 0:
      self.sizeconf += 3
    if self.is_compact:
      self.sizeconf += self.natms
    else:
      self.sizeconf += self.natms * 2
    if self.keytrj > 0:
      self.sizeconf += self.natms
    if self.keytrj > 1:
      self.sizeconf += self.natms

    # Skip useless confs
    if skip:
      line2skip = skip * self.sizeconf

      for _ in xrange(line2skip):
        if self.f.next() == '':
          raise ValueError('Skipping too much conf')

      logging.info('I skipped %d lines in the HISTORY file %s for %d confs', line2skip, self.filename, skip)

    if self.field:
      # Create atom types via field
      for fmol in self.field.mols:
        for fatom in fmol['atoms']:
          key = fmol['name'] + '/' + fatom['name']
          if key in self.atom_types:
            continue
          atype = molecule.AtomType(key, fatom)
          self.atom_types[key] = atype


  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    cnf = {}

    cnf['filename'] = self.filename
    cnf['i'] = self.i
    cnf['history'] = self

    self.i += 1

    # Record i
    try:
      fields = self.f.next().split()
      cnf['nstep'] = int(fields[1])
      cnf['natms'] = int(fields[2])
      cnf['keytrj'] = int(fields[3])
      cnf['imcon'] = int(fields[4])
      cnf['tstep'] = float(fields[5])
    except IndexError:
      raise StopIteration

    if cnf['imcon'] > 0:
      cnf['cell'] = []
      # Record ii
      for _ in xrange(3):
        fields = self.f.next().split()
        cnf['cell'].append((float(fields[0]), float(fields[1]), float(fields[2])))

    # Read atoms informations
    nlines = self.natms
    if not self.is_compact:
      nlines += self.natms
    if self.keytrj > 0:
      nlines += self.natms
    if self.keytrj > 1:
      nlines += self.natms
    s = itertools.islice(self.f, nlines)
    raw_dynamic = list(s)

    cnf['numpy'] = self.numpy

    raw = (self.natms, self.keytrj, self.xyz_only, self.raw_static, raw_dynamic)
    if self.raw:
      cnf['raw'] = raw
      cnf['atoms'] = None
    else:
      cnf['atoms'] = History.raw2atoms(raw, self.numpy)
    return cnf

  @staticmethod
  def raw2atoms(raw, numpy):
    '''Transform raw lines to atoms (array numpy or dict)'''

    nbatoms, keytrj, xyz_only, static, dynamic = raw
    # Normal mode
    if static is None:
      stride = 2 + keytrj
      offset = 1
      static = dynamic[::stride]
    # Compact mode 
    else:
      stride = 1 + keytrj
      offset = 0

    if numpy:
      dtype1 = [('atmnam', '|S4'), ('iatm', int), ('weight', float), ('charge', float)]
      dtype2 = [('xxx', float), ('yyy', float), ('zzz', float)]
      dtype3 = [('vxx', float), ('vyy', float), ('vzz', float)]
      dtype4 = [('fxx', float), ('fyy', float), ('fzz', float)]
      if xyz_only:
        dtype = dtype2
      else:
        dtype = dtype1 + dtype2
        if keytrj > 0:
          dtype += dtype3
        if keytrj > 1:
          dtype += dtype4
      atoms = np.empty(nbatoms, dtype=dtype)
      if not xyz_only:
        array = np.genfromtxt(static, dtype=dtype1)
        atoms['atmnam'] = array['atmnam']
        atoms['iatm'] = array['iatm']
        atoms['weight'] = array['weight']
        atoms['charge'] = array['charge']
      array = np.genfromtxt(dynamic[offset::stride], dtype=dtype2)
      atoms['xxx'] = array['xxx']
      atoms['yyy'] = array['yyy']
      atoms['zzz'] = array['zzz']
      if keytrj > 0 and not xyz_only:
        array = np.genfromtxt(dynamic[offset+1::stride], dtype=dtype3)
        atoms['vxx'] = array['vxx']
        atoms['vyy'] = array['vyy']
        atoms['vzz'] = array['vzz']
      if keytrj > 1 and not xyz_only:
        array = np.genfromtxt(dynamic[offset+2::stride], dtype=dtype4)
        atoms['fxx'] = array['fxx']
        atoms['fyy'] = array['fyy']
        atoms['fzz'] = array['fzz']

    else:
      atoms = [{} for _ in xrange(nbatoms)]
      if not xyz_only:
        for iatom, line in enumerate(static):
          atom = atoms[iatom]
          fields = line.split()
          atom['atmnam'] = fields[0]
          atom['iatm'] = int(fields[1])
          atom['weight'] = float(fields[2])
          atom['charge'] = float(fields[3])
      for iatom, line in enumerate(dynamic[offset::stride]):
        atom = atoms[iatom]
        fields = line.split()
        atom['xxx'] = float(fields[0])
        atom['yyy'] = float(fields[1])
        atom['zzz'] = float(fields[2])
      if keytrj > 0 and not xyz_only:
        for iatom, line in enumerate(dynamic[offset+1::stride]):
          atom = atoms[iatom]
          fields = line.split()
          atom['vxx'] = float(fields[0])
          atom['vyy'] = float(fields[1])
          atom['vzz'] = float(fields[2])
      if keytrj > 1 and not xyz_only:
        for iatom, line in enumerate(dynamic[offset+2::stride]):
          atom = atoms[iatom]
          fields = line.split()
          atom['fxx'] = float(fields[0])
          atom['fyy'] = float(fields[1])
          atom['fzz'] = float(fields[2])
    return atoms


class Statis(object):
  '''Class representing a Statis file reader'''

  def __init__(self, filename='STATIS', field_filename=None, ntpatm=None, is_stress=False, is_npt=False):

    self.filename = filename
    self.ntpatm = ntpatm
    self.is_stress = is_stress
    self.is_npt = is_npt

    # Get the ntpatm values from the FIELD
    if field_filename is not None:
      field = Field(field_filename)
      self.ntpatm = len(set([pair['atom_type1'] for pair in field.vdws]))

    if (is_stress or is_npt) and not self.ntpatm:
      raise ValueError('Cannot get stress or NPT infos without ntpatm parameter')

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Read header
      self.header = self.f.readline().strip()
    except IOError:
      self.f = open(filename, 'r')
      # Read header
      self.header = self.f.readline().strip()

    line = self.f.readline()
    try:
      self.units = line.split('=')[1].strip()
    except IndexError:
      self.units = None

  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    cnf = {}

    # Record i
    try:
      fields = self.f.readline().split()
      cnf['nstep'] = int(fields[0])
      cnf['time'] = float(fields[1])
      cnf['nument'] = int(fields[2])
    except IndexError:
      raise StopIteration

    # Record ii
    fields = self.f.readline().split()
    cnf['engcns'] = float(fields[0])
    cnf['temp'] = float(fields[1])
    cnf['engcfg'] = float(fields[2])
    cnf['engsrp'] = float(fields[3])
    cnf['engcpe'] = float(fields[4])

    # Record iii
    fields = self.f.readline().split()
    cnf['engbnd'] = float(fields[0])
    cnf['engang'] = float(fields[1])
    cnf['engdih'] = float(fields[2])
    cnf['engtet'] = float(fields[3])
    cnf['enthal'] = float(fields[4])

    # Record iv
    fields = self.f.readline().split()
    cnf['tmprot'] = float(fields[0])
    cnf['vir'] = float(fields[1])
    cnf['virsrp'] = float(fields[2])
    cnf['vircpe'] = float(fields[3])
    cnf['virbnd'] = float(fields[4])

    # Record v
    fields = self.f.readline().split()
    cnf['virang'] = float(fields[0])
    cnf['vircon'] = float(fields[1])
    cnf['virtet'] = float(fields[2])
    cnf['volume'] = float(fields[3])
    cnf['tpmshl'] = float(fields[4])

    # Record vi
    fields = self.f.readline().split()
    cnf['engshl'] = float(fields[0])
    cnf['virshl'] = float(fields[1])
    cnf['alpha'] = float(fields[2])
    cnf['beta'] = float(fields[3])
    cnf['gamma'] = float(fields[4])

    # Record vii
    fields = self.f.readline().split()
    cnf['virpmf'] = float(fields[0])
    cnf['press'] = float(fields[1])

    # Remaining (non fixed) records
    # The nb of entities is counted from the second line
    remain_nument = cnf['nument'] - 27

    cur_x = 2

    # If I don't have ntpatm infos, I can't do more
    if self.ntpatm:

      # Record amsd
      cnf['amsd'] = []
      for _ in xrange(self.ntpatm):
        remain_nument -= 1
        # Line change
        if cur_x == 5:
          cur_x = 0
          fields = self.f.readline().split()
        cnf['amsd'].append(float(fields[cur_x]))
        cur_x += 1

      # Record stress
      if self.is_stress and remain_nument >= 9:
        cnf['stress'] = []
        for _ in xrange(9):
          remain_nument -= 1
          # Line change
          if cur_x == 5:
            cur_x = 0
            fields = self.f.readline().split()
          cnf['stress'].append(float(fields[cur_x]))
          cur_x += 1

      # Record cell (NPT)
      if self.is_npt and remain_nument == 9:
        cnf['cell'] = []
        cell_vect = []
        for _ in xrange(3):
          remain_nument -= 1
          # Line change
          if cur_x == 5:
            cur_x = 0
            fields = self.f.readline().split()
          cell_vect.append(float(fields[cur_x]))
          cur_x += 1
        cnf['cell'].append(tuple(cell_vect))
        cell_vect = []
        for _ in xrange(3):
          remain_nument -= 1
          # Line change
          if cur_x == 5:
            cur_x = 0
            fields = self.f.readline().split()
          cell_vect.append(float(fields[cur_x]))
          cur_x += 1
        cnf['cell'].append(tuple(cell_vect))
        cell_vect = []
        for _ in xrange(3):
          remain_nument -= 1
          # Line change
          if cur_x == 5:
            cur_x = 0
            fields = self.f.readline().split()
          cell_vect.append(float(fields[cur_x]))
          cur_x += 1
        cnf['cell'].append(tuple(cell_vect))

    # Get the remaining element on the current line
    remain_nument -= 5 - cur_x

    # Skip remaining elements
    remain_line = int(math.ceil(remain_nument / 5.))
    logging.debug('Skip %d line for %d entities', remain_line, remain_nument)
    for _ in xrange(remain_line):
      self.f.readline()

    return cnf

  def to_array(self):
    '''Convert the STATIS to a numpy array'''

    cnf = self.next()
    keys = cnf.keys()

    # All values are float except some
    dtype = [(key, float) for key in keys]
    dtype[keys.index('nstep')] = ('snetp', int)
    dtype[keys.index('nument')] = ('nument', int)
    if self.ntpatm:
      dtype[keys.index('amsd')] = ('amsd', float, self.ntpatm)
    ss = itertools.chain((cnf,), self)
    return np.array([tuple([cnf[key] for key in keys]) for cnf in ss], dtype=dtype)

class RDFDat(object):
  '''Class representing a RDFDat file reader'''

  def __init__(self, filename='RDFDAT', field_filename=None, no_intra=False):

    self.filename = filename

    # Read field if need to remove intra molecular RDF
    if field_filename and no_intra:
      field = Field(field_filename)
      atom_in_mols = {}
      for mol in field.mols:
        for atom in mol['atoms']:
          if atom['name'] not in atom_in_mols:
            atom_in_mols[atom['name']] = []
          if mol['name'] not in atom_in_mols[atom['name']]:
            atom_in_mols[atom['name']].append(mol['name'])

    f = open(filename, 'r')
    self.header = f.readline().strip()

    line = f.readline().strip()
    self.nb_rdfs = int(line.split()[0])
    self.nb_bins = int(line.split()[1])

    self.rdfs = {}

    for _ in xrange(self.nb_rdfs):
      atom_names = tuple(f.readline().split())

      # Check if atoms are only intra
      if field_filename and no_intra:
        if len(set(atom_in_mols[atom_names[0]] + atom_in_mols[atom_names[1]])) == 1:
          # Skip
          for __ in xrange(self.nb_bins):
            f.readline()
          continue

      self.rdfs[atom_names] = []

      for __ in xrange(self.nb_bins):
        line = f.readline()
        r = float(line.split()[0])
        gr = float(line.split()[1])
        self.rdfs[atom_names].append((r, gr))

    f.close()

  def __str__(self):
    return 'RDFDAT DL_POLY %s' % (self.header)
