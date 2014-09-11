# -*- coding: iso-8859-1 -*-
'''Python library to manage lammps files'''

import logging
import gzip
import itertools
import StringIO
import numpy as np
try:
  import pandas
except ImportError:
  pass

import molecule

def load_file(filename='log.lammps', run=1):
  '''Helper function to load a lammps file'''

  f = open(filename, 'r')
  line = f.readline().strip()
  f.close()

  if line.startswith('# Time-averaged data for'):
    print 'Loading LAMMPS time-averaged data from', filename
    return Ave(filename).to_array()
  elif line.startswith('# Spatial-averaged data for'):
    print 'Loading LAMMPS spatial-average data from', filename
    return Ave(filename).to_array()
  elif line.startswith('LAMMPS'):
    print 'Loading LAMMPS log file from', filename, line[6:]
    return Log(filename, run).to_array()
  else:
    raise ValueError('Unknown format for file %s' % filename)

class Data(object):
  '''Class representing a Data file'''

  def __init__(self, filename=None):

    self.header = None

    self.atoms = []
    self.bonds = []
    self.angles = []
    self.dihedrals = []
    self.impropers = []

    self.atom_types = []
    self.bond_types = []
    self.angle_types = []
    self.dihedral_types = []
    self.improper_types = []

    self.mols = []

    self.box = [(None, None), (None, None), (None, None)]
    self.tilt = (0, 0, 0)

    if filename:
      self.read_from_file(filename)
      self.filename = filename
    else:
      self.filename = None

  def read_from_file(self, filename='lammps.data'):

    nb_atoms = 0
    nb_bonds = 0
    nb_angles = 0
    nb_dihedrals = 0
    nb_impropers = 0
    nb_atom_types = 0
    nb_bond_types = 0
    nb_angle_types = 0
    nb_dihedral_types = 0
    nb_improper_types = 0

    nb_mols = 0

    f = open(filename, 'r')
    self.header = f.readline().strip()

    section = 'Header'

    for line in f:
      line = line.strip()

      if line == '':
        continue

      # Section change
      if line == 'Masses' or line == 'Pair Coeffs' or line == 'Bond Coeffs' or line == 'Angle Coeffs' or line == 'Dihedral Coeffs' or line == 'Improper Coeffs' or \
         line == 'Atoms' or line == 'Bonds' or line == 'Angles' or line == 'Dihedrals' or line == 'Impropers' or line == 'Velocities' or line == 'Extras' or line == 'EXTRA':
        section = line
        continue

      # Header
      if section == 'Header':
        if 'atoms' in line:
          nb_atoms = int(line.split()[0])
          self.atoms = [ None ] * nb_atoms
          continue
        if 'bonds' in line:
          nb_bonds = int(line.split()[0])
          continue
        if 'angles' in line:
          nb_angles = int(line.split()[0])
          continue
        if 'dihedrals' in line:
          nb_dihedrals = int(line.split()[0])
          continue
        if 'impropers' in line:
          nb_impropers = int(line.split()[0])
          continue

        if 'atom types' in line:
          nb_atom_types = int(line.split()[0])
          continue
        if 'bond types' in line:
          nb_bond_types = int(line.split()[0])
          continue
        if 'angle types' in line:
          nb_angle_types = int(line.split()[0])
          continue
        if 'dihedral types' in line:
          nb_dihedral_types = int(line.split()[0])
          continue
        if 'improper types' in line:
          nb_improper_types = int(line.split()[0])
          continue

        if 'xlo xhi' in line:
          self.box[0] = (float(line.split()[0]), float(line.split()[1]))
          continue
        if 'ylo yhi' in line:
          self.box[1] = (float(line.split()[0]), float(line.split()[1]))
          continue
        if 'zlo zhi' in line:
          self.box[2] = (float(line.split()[0]), float(line.split()[1]))
          continue
        if 'xy xz yz' in line:
          self.tilt = (float(line.split()[0]), float(line.split()[1]), float(line.split()[2]))
          continue

      # Masses
      #TODO should I append or use i ? If the lines are not sorted ???
      if section == 'Masses':
        i = int(line.split()[0])
        mass = float(line.split()[1])

        # Should I create atoms dict ?
        if len(self.atom_types) < nb_atom_types:
          atom_type = {}
          atom_type['mass'] = mass
          self.atom_types.append(atom_type)
        else:
          atom_type = self.atom_types[i - 1]
          atom_type['mass'] = mass
        continue

      # Pair Coeffs
      if section == 'Pair Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break

        # Should I create atoms dict ?
        if len(self.atom_types) < nb_atom_types:
          atom_type = {}
          atom_type['coeffs'] = coeffs
          self.atom_types.append(atom_type)
        else:
          atom_type = self.atom_types[i - 1]
          atom_type['coeffs'] = coeffs
        continue

      # Bond Coeffs
      if section == 'Bond Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        self.bond_types.append({'coeffs': coeffs})

      # Angle Coeffs
      if section == 'Angle Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        self.angle_types.append({'coeffs': coeffs})

      # Dihedral Coeffs
      if section == 'Dihedral Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        self.dihedral_types.append({'coeffs': coeffs})

      # Improper Coeffs
      if section == 'Improper Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        self.improper_types.append({'coeffs': coeffs})

      # Atoms
      #TODO assume atom_style is full
      if section == 'Atoms':
        i = int(line.split()[0])
        mol_i = int(line.split()[1])
        atom_type_i = int(line.split()[2])
        charge = float(line.split()[3])
        x = float(line.split()[4])
        y = float(line.split()[5])
        z = float(line.split()[6])

        atom = {}
        atom['i'] = i
        atom['mol_i'] = mol_i
        if mol_i > nb_mols:
          nb_mols = mol_i
        atom['atom_type_i'] = atom_type_i
        atom['atom_type'] = self.atom_types[atom_type_i - 1]
        atom['charge'] = charge
        atom['x'] = x
        atom['y'] = y
        atom['z'] = z
        # Atoms are not ordered
        self.atoms[i - 1] = atom

      # Bonds
      if section == 'Bonds':
        i = int(line.split()[0])
        bond_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])

        bond = {}
        bond['bond_type'] = self.bond_types[bond_type_i - 1]
        bond['atom1'] = self.atoms[atom1_i - 1]
        bond['atom2'] = self.atoms[atom2_i - 1]
        self.bonds.append(bond)

      # Angles
      if section == 'Angles':
        i = int(line.split()[0])
        angle_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])
        atom3_i = int(line.split()[4])

        angle = {}
        angle['angle_type'] = self.angle_types[angle_type_i - 1]
        angle['atom1'] = self.atoms[atom1_i - 1]
        angle['atom2'] = self.atoms[atom2_i - 1]
        angle['atom3'] = self.atoms[atom3_i - 1]
        self.angles.append(angle)

      # Dihedrals
      if section == 'Dihedrals':
        i = int(line.split()[0])
        dihedral_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])
        atom3_i = int(line.split()[4])
        atom4_i = int(line.split()[5])

        dihedral = {}
        dihedral['dihedral_type'] = self.dihedral_types[dihedral_type_i - 1]
        dihedral['atom1'] = self.atoms[atom1_i - 1]
        dihedral['atom2'] = self.atoms[atom2_i - 1]
        dihedral['atom3'] = self.atoms[atom3_i - 1]
        dihedral['atom4'] = self.atoms[atom4_i - 1]
        self.dihedrals.append(dihedral)

      # Impropers
      if section == 'Impropers':
        i = int(line.split()[0])
        improper_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])
        atom3_i = int(line.split()[4])
        atom4_i = int(line.split()[5])

        improper = {}
        improper['improper_type'] = self.improper_types[improper_type_i - 1]
        improper['atom1'] = self.atoms[atom1_i - 1]
        improper['atom2'] = self.atoms[atom2_i - 1]
        improper['atom3'] = self.atoms[atom3_i - 1]
        improper['atom4'] = self.atoms[atom4_i - 1]
        self.impropers.append(improper)

      # TODO ignore Velocities and Extra

    f.close()

    if None in self.atoms:
      raise ValueError('Nb of atoms is not coherent')

    # Reconstruct mols
    self.mols = []
    for i in xrange(nb_mols):
      mol = {'i': i + 1}
      self.mols.append(mol)

    for atom in self.atoms:
      mol = self.mols[atom['mol_i'] - 1]
      if 'atoms' not in mol:
        mol['atoms'] = []
      mol['atoms'].append(atom)
      atom['mol'] = mol
      del atom['mol_i']

class Log(object):
  '''Class representing a Log file reader (for thermo infos)'''

  def __init__(self, filename='log.lammps', run=1):

    self.filename = filename

    self.run = run

    self.steps = 0

    self.is_multi = False

    self.cache_lines = []

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = gzip.open(filename, 'r')
    except IOError:
      self.f = open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = open(filename, 'r')

    # Seek to the start of the thermo stuff (and skip)
    for _ in xrange(self.run):
      line = self.f.readline()
      while not line.startswith('Memory usage per processor'):
        line = self.f.readline()
        # Get and store steps (if int)
        if line.startswith('run'):
          try:
            self.steps = int(line.split()[1])
          except ValueError:
            pass
        # If EOF
        if line == '':
          if self.run != 1:
            raise ValueError('File %s does not seem to be a valid LAMMPS log file (with thermo infos) or does not contain info for run %d' % (self.filename, self.run))
          else:
            raise ValueError('File %s does not seem to be a valid LAMMPS log file (with thermo infos)' % self.filename)

    line = self.f.readline()
    if not line.startswith('-----'):
      tmp_fields = [ field.lower() for field in line.split() ]
      self.fields = []
      mult = {}
      for field in tmp_fields:
        # Manage multiple fields with same name
        if tmp_fields.count(field) > 1:
          if field not in mult:
            mult[field] = 0
          self.fields.append('%s_%d' % (field, mult[field] + 1))
          mult[field] += 1
        else:
          self.fields.append(field)
    else:
      self.next_step = int(line.split()[2])
      self.next_cpu = float(line.split()[6])
      self.is_multi = True
      self.fields = ['step', 'cpu']

      # For multi lines, fields should be get from the first record (with lines put in cache for the first read !)
      line = self.f.readline()
      # Put line in cache
      self.cache_lines.append(line)
      while not line.startswith('-----') and line != '' and not line.startswith('Memory usage per processor'):
        fields = line.split()
        # Get all the data in format: Field = Value
        # For each line, get the infos (line format: Field = Value)
        is_first = True
        for field in fields:
          if field == '=':
            continue
          # Get fields from left side
          elif is_first:
            if field not in self.fields:
              self.fields.append(field.lower())
            is_first = False
          else:
            is_first = True
        line = self.f.readline()
        # Put line in cache
        self.cache_lines.append(line)

  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    stats = {}

    if not self.is_multi:
      # One line version
      line = self.f.readline()
      values = line.split()
      if line == '' or line.startswith('Loop'):
        raise StopIteration

      # Skip garbage
      while len(values) != len(self.fields):
        line = self.f.readline()
        values = line.split()
        if line == '' or line.startswith('Loop'):
          raise StopIteration

      for (field, value) in zip(self.fields, values):
        # Special case for int values
        if field in ('step', 'elapsed'):
          stats[field] = int(value)
        else:
          stats[field] = float(value)
    else:
      # Multi lines version
      stats['step'] = self.next_step
      stats['cpu'] = self.next_cpu

      if len(self.cache_lines):
        line = self.cache_lines.pop(0)
      else:
        line = self.f.readline()
      if line == '':
        raise StopIteration
      values = line.split()
      # Get all the data in format: Field = Value
      while len(values) > 2 and values[1] == '=':
        # For each line, get the infos (line format: Field = Value)
        field = None
        for value in values:
          if value == '=':
            continue
          if field is None:
            field = value
          else:
            stats[field.lower()] = float(value)
            field = None
        if len(self.cache_lines):
          line = self.cache_lines.pop(0)
        else:
          line = self.f.readline()
        values = line.split()

      # Try to get the next record
      while not line.startswith('-----') and line != '' and not line.startswith('Memory usage per processor'):
        if len(self.cache_lines):
          line = self.cache_lines.pop(0)
        else:
          line = self.f.readline()

      # Reach next log
      if line.startswith('Memory usage per processor'):
        raise StopIteration

      # We reach the first line of the next record
      if line.startswith('-----'):
        values = line.split()
        try:
          self.next_step = int(values[2])
          self.next_cpu = float(values[6])
        except IndexError:
          raise StopIteration

    return stats

  def to_array(self):
    '''Convert the log to a numpy array'''

    nb_fields = len(self.fields)
    types = zip(self.fields, [float] * nb_fields)
    return np.array([tuple([ line[field] for field in self.fields ]) for line in self], dtype=types)


class Dump(object):
  '''Class representing a Dump file reader'''

  def __init__(self, filename='lammps.dump', skip=0, sort=True, numpy=False, data=None, raw=True):

    self.filename = filename

    self.sort = sort
    self.numpy = numpy
    self.data = data
    self.raw = raw
    self.i = skip + 1
    self.atom_types = {}

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Test the read
      line = self.f.readline().strip()
      if line != 'ITEM: TIMESTEP':
        raise IOError('Not a valid LAMMPS dump file')
      self.f.readline()
      self.f.readline()
      # And get the nb of atoms
      nbatoms = int(self.f.readline())
      self.f.close()
      # If OK, reopen it
      self.f = gzip.open(filename, 'r')
    except IOError:
      self.f = open(filename, 'r')
      # Test the read
      line = self.f.readline().strip()
      if line != 'ITEM: TIMESTEP':
        raise IOError('Not a valid LAMMPS dump file')
      self.f.readline()
      self.f.readline()
      # And get the nb of atoms
      nbatoms = int(self.f.readline())
      self.f.close()
      # If OK, reopen it
      self.f = open(filename, 'r')

    self.sizeconf = 9 + nbatoms

    # Skip useless confs
    if skip:
      line2skip = skip * self.sizeconf

      for _ in xrange(line2skip):
        if self.f.readline() == '':
          raise ValueError('Skipping too much conf')

      logging.info('I skipped %d lines in the Dump for %d confs', line2skip, skip)

    if self.data:
      # Create atom types via data
      i_at = 0
      for atom_type in self.data.atom_types:
        i_at += 1
        atype = molecule.AtomType(str(i_at), atom_type)
        self.atom_types[str(i_at)] = atype

      # Create dummy atoms from data to try to factorize into atomTypes
      for atom in self.data.atoms:
        atype = self.atom_types[str(atom['atom_type_i'])]
        atom = molecule.Atom(str(atom['i']), atom['i'], 1, atom['x'], atom['y'], atom['z'], \
                                 params={'charge': atom['charge']}, atype=atype)

      for atom_type in self.atom_types.values():
        atom_type.factorize()


  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    cnf = {}

    cnf['filename'] = self.filename
    cnf['i'] = self.i
    cnf['dump'] = self
    self.i += 1

    # Reset all atom_type
    for atom_type in self.atom_types.values():
      atom_type.reset()

    line = self.f.next().strip()
    if line == '':
      raise StopIteration
    if line != 'ITEM: TIMESTEP':
      raise ValueError('Not a valid LAMMPS dump file')
    timestep = int(self.f.next())
    cnf['timestep'] = timestep

    line = self.f.next().strip()
    if line != 'ITEM: NUMBER OF ATOMS':
      raise ValueError('Not a valid LAMMPS dump file')
    nbatoms = int(self.f.next())
    cnf['nbatoms'] = nbatoms

    # Box stuff
    line = self.f.next().strip()
    if not line.startswith('ITEM: BOX BOUNDS'):
      raise ValueError('Not a valid LAMMPS dump file')
    cnf['boundaries'] = line.split()[-3:]
    if 'xy xz yz' in line:
      is_tilt = True
    else:
      is_tilt = False

    # The format can't change, so store in raw style
    cnf['box'] = []
    line = self.f.next()
    cnf['box'].append([ float(line.split()[0]), float(line.split()[1]) ])
    if is_tilt:
      xy = float(line.split()[2])

    line = self.f.next()
    cnf['box'].append([ float(line.split()[0]), float(line.split()[1]) ])
    if is_tilt:
      xz = float(line.split()[2])

    line = self.f.next()
    cnf['box'].append([ float(line.split()[0]), float(line.split()[1]) ])
    if is_tilt:
      yz = float(line.split()[2])
    if is_tilt:
      cnf['tilt'] = [xy, xz, yz]
    else:
      cnf['tilt'] = [0.0, 0.0, 0.0]

    # Get item names
    line = self.f.next().strip()
    if not line.startswith('ITEM: ATOMS'):
      raise ValueError('Not a valid LAMMPS dump file')
    fields = line.split()[2:]

    cnf['numpy'] = self.numpy

    s = itertools.islice(self.f, cnf['nbatoms'])
    raw = (list(s), fields, self.sort)
    if self.raw:
      cnf['raw'] = raw
      cnf['atoms'] = None
    else:
      cnf['atoms'] = Dump.raw2atoms(raw, self.numpy)

    return cnf

  @staticmethod
  def raw2atoms(raw, numpy):
    '''Transform raw lines to atoms (array numpy or dict)'''

    (s, fields, sort) = raw
    if numpy:
      #atoms = np.genfromtxt(s, dtype=None, names=fields, usemask=False)
      atoms = pandas.read_csv(StringIO.StringIO(''.join(s)), ' ', header=None, index_col=False, names=fields, as_recarray=True).view(np.ndarray).copy()
      if sort:
        atoms.sort(order='id')
    else:
      atoms = []
      for line in s:
        atom = {}
        values = line.split()
        for (field, value) in zip(fields, values):
          try:
            # Try for int of float
            try:
              atom[field] = int(value)
            except ValueError:
              atom[field] = float(value)
          except ValueError:
            atom[field] = value
        atoms.append(atom)

      if sort:
        atoms = sorted(atoms, key=lambda atom: atom['id'])

    return atoms


class Ave(object):
  '''Class representing a Average file reader'''

  def __init__(self, filename):

    self.filename = filename

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = gzip.open(filename, 'r')
    except IOError:
      self.f = open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = open(filename, 'r')

    # Read the first line header
    line = self.f.readline()
    if not line.startswith('#'):
      raise ValueError('File %s does not seem to be a valid LAMMPS ave file' % self.filename)
    self.header = line.split('#', 1)[1].strip()

    # Skip the second useless lines
    line = self.f.readline()
    if not line.startswith('#'):
      raise ValueError('File %s does not seem to be a valid LAMMPS ave file' % self.filename)

    # The third ones contains field names
    line =  self.f.readline()
    if not line.startswith('#'):
      raise ValueError('File %s does not seem to be a valid LAMMPS ave file' % self.filename)

    self.fields = line.split('#')[1].split()

  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    cnf = {}

    line = self.f.readline()
    if line == '':
      raise StopIteration

    fields = line.split()
    cnf['timestep'] = int(fields[0])
    nb_line = int(fields[1])

    for field in self.fields:
      cnf[field] = []

    for _ in xrange(nb_line):
      line = self.f.readline()
      fields = line.split()

      # check the coherence between header and line
      if len(fields) != len(self.fields):
        raise ValueError('Invalid data or field names : %s' % self.fields)

      i = 0
      for field in self.fields:
        try:
          val = int(fields[i])
        except ValueError:
          val = float(fields[i])
          val = float(fields[i])
        cnf[field].append(val)
        i += 1

    return cnf

  def to_array(self):
    '''Convert the average to a numpy array'''

    nb_fields = len(self.fields)
    data = list(self)
    types = [('timestep', float)] + zip(self.fields, [float]*nb_fields, [ (len(data[0][self.fields[0]]),) ] * nb_fields)
    return np.array([ (line['timestep'],) + tuple(line[self.fields[i]] for i in xrange(nb_fields)) for line in data], dtype=types)
