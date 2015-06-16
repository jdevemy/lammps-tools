# -*- coding: iso-8859-1 -*-
'''Python library to manage mol or zmat files'''

class ZMat(object):
  '''Class representing a mol/zmat file'''

  def __init__(self, filename):

    self.filename = filename

    try:
      f = open(filename, 'r')

      # Skip comments
      line = f.readline()
      while line.startswith('#') and line != '':
        line = f.readline()

      self.name = line.strip()

      # Skip comments
      line = f.readline()
      while (line.startswith('#') or line.strip() == '') and line != '':
        line = f.readline()

      self.atoms = []

      # There can be the atom index in first column, or not !
      shift = 0
      if len(line.split()) > 1:
        shift = 1

      # First atom
      atom = {'name': line.split()[shift], 'i': 1, 'dist_atom': None, 'dist': 0.0, 'angle_atom': None, 'angle': 0.0, 'dihedr_atom': None, 'dihedr': 0.0, 'extra_connect': []}
      self.atoms.append(atom)

      i = 2
      # Read Atoms
      line = f.readline()
      while line.strip() != '':
        fields = line.split()
        name = fields[shift]
        dist_atom = self.atoms[int(fields[1 + shift]) - 1]
        try:
          dist = float(fields[2 + shift])
        except ValueError:
          dist = fields[2 + shift]

        if len(fields) > 4:
          angle_atom = self.atoms[int(fields[3 + shift]) - 1]
          angle = float(fields[4 + shift])
        else:
          angle_atom = None
          angle = 0.0

        if len(fields) > 6:
          dihedr_atom = self.atoms[int(fields[5 + shift]) - 1]
          dihedr = float(fields[6 + shift])
        else:
          dihedr_atom = None
          dihedr = 0.0

        atom = {'name': name, 'i': i, 'dist_atom': dist_atom, 'dist': dist, 'angle_atom': angle_atom, 'angle': angle, 'dihedr_atom': dihedr_atom, 'dihedr': dihedr, 'extra_connect': []}
        self.atoms.append(atom)

        i += 1

        line = f.readline()

      # Skip comments
      while (line.startswith('#') or line.strip() == '') and line != '':
        line = f.readline()

      variables = {}

      # Read variables
      while '=' in line:
        fields = line.split()
        variables[fields[0]] = float(fields[2])
        line = f.readline()

      # Reput variables values in data
      for atom in self.atoms:
        if atom['dist'] in variables:
          atom['dist'] = variables[atom['dist']]
        if isinstance(atom['dist'], basestring):
          raise Exception('Not a valid .mol/.zmat file: variable %s undefined' % atom['dist'])
        if atom['angle'] in variables:
          atom['angle'] = variables[atom['angle']]
        if isinstance(atom['angle'], basestring):
          raise Exception('Not a valid .mol/.zmat file: variable %s undefined' % atom['angle'])
        if atom['dihedr'] in variables:
          atom['dihedr'] = variables[atom['dihedr']]
        if isinstance(atom['dihedr'], basestring):
          raise Exception('Not a valid .mol/.zmat file: variable %s undefined' % atom['dihedr'])

      # Skip comments
      while (line.startswith('#') or line.strip() == '') and line != '':
        line = f.readline()

      # Read Extra Connect
      while line.startswith('connect'):
        fields = line.split()
        atom1 = self.atoms[int(fields[1]) - 1]
        atom2 = self.atoms[int(fields[2]) - 1]
        atom1['extra_connect'].append(atom2)
        atom2['extra_connect'].append(atom1)
        line = f.readline()

      # Skip comments
      while (line.startswith('#') or line.strip() == '') and line != '':
        line = f.readline()

      # Read impropers
      self.impropers = []

      while line.startswith('improper'):
        fields = line.split()
        atom1 = self.atoms[int(fields[1]) - 1]
        atom2 = self.atoms[int(fields[2]) - 1]
        atom3 = self.atoms[int(fields[3]) - 1]
        atom4 = self.atoms[int(fields[4]) - 1]
        self.impropers.append((atom1, atom2, atom3, atom4))
        line = f.readline()

      # Read ff
      self.ff = None
      if line != '' and '.ff' in line:
        self.ff = line.strip()

      self.ff_lines = []
      if self.ff is None:
        while line != '':
          self.ff_lines.append(line)
          line = f.readline()

      f.close()

    except (IndexError, ValueError, IOError), e:
      raise Exception('Not a valid .mol/.zmat file: %s' % e)

  def __str__(self):
    return 'ZMAT %s(%s)' \
      % (self.filename, self.name)

