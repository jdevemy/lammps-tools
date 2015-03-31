# -*- coding: iso-8859-1 -*-
'''Python library to manage configuration'''

import logging

import molecule
import lammps
import dlpoly

class Box(object):

  def __init__(self, cnf=None, typ=None):

    self.typ = typ

    if typ == 'dlpoly':
      # Missing things, but enough for us

      # Min Max
      self.x = (0, max(cnf['cell'][0][0], cnf['cell'][1][0], cnf['cell'][2][0]))
      self.y = (0, max(cnf['cell'][0][1], cnf['cell'][1][1], cnf['cell'][2][1]))
      self.z = (0, max(cnf['cell'][0][2], cnf['cell'][1][2], cnf['cell'][2][2]))

      # Tilt
      self.xy = self.xz = self.yz = 0.0

      # Periodic
      if cnf['imcon'] == 0:
        self.px = False
        self.py = False
        self.pz = False
      elif cnf['imcon'] == 6:
        self.px = True
        self.py = True
        self.pz = False
      else:
        self.px = True
        self.py = True
        self.pz = True

    elif typ == 'lammps_dump':
      # Min Max
      self.x = tuple(cnf['box'][0])
      self.y = tuple(cnf['box'][1])
      self.z = tuple(cnf['box'][2])

      # Tilt
      if 'tilt' in cnf:
        self.xy = cnf['tilt'][0]
        self.xz = cnf['tilt'][1]
        self.yz = cnf['tilt'][2]

      # Periodic
      if 'boundaries' in cnf:
        if 'p' in cnf['boundaries'][0]:
          self.px = True
        if 'p' in  cnf['boundaries'][1]:
          self.py = True
        if 'p' in  cnf['boundaries'][2]:
          self.pz = True

    elif typ == 'lammps_data':
      # Min Max
      self.x = tuple(cnf.box[0])
      self.y = tuple(cnf.box[1])
      self.z = tuple(cnf.box[2])

      # Tilt
      self.xy = cnf.tilt[0]
      self.xz = cnf.tilt[1]
      self.yz = cnf.tilt[2]

      # Periodic
      # Cannot detect periodic conditions

    else:
      # Default box

      # Min Max
      self.x = (-1.0, 1.0)
      self.y = (-1.0, 1.0)
      self.z = (-1.0, 1.0)

      self.xy = self.xz = self.yz = 0.0

    self.lx = abs(self.x[1] - self.x[0])
    self.ly = abs(self.y[1] - self.y[0])
    self.lz = abs(self.z[1] - self.z[0])

  def __str__(self):
    return 'Box %f/%f %f/%f %f/%f' \
      % (self.x[0], self.x[1], self.y[0], self.y[1], self.z[0], self.z[1])

  def __repr__(self):
    return '%f/%f %f/%f %f/%f' \
      % (self.x[0], self.x[1], self.y[0], self.y[1], self.z[0], self.z[1])

  def __eq__(self, other):
    return self.x == other.x and self.y == other.y and self.z == other.z \
      and self.xy == other.xy and self.xz == other.xz and self.yz == other.yz

  def __ne__(self, other):
    return not self == other

  @property
  def v1(self):
    return (self.x[1] - self.x[0], 0.0, 0.0)

  @property
  def v2(self):
    return (self.xy, self.y[1] - self.y[0], 0.0)

  @property
  def v3(self):
    return (self.xz, self.yz, self.z[1] - self.z[0])

  @property
  def len_x(self):
    return self.x[1] - self.x[0]

  @property
  def len_y(self):
    return self.y[1] - self.y[0]

  @property
  def len_z(self):
    return self.z[1] - self.z[0]

  @property
  def volume(self):
    return abs(self.x[1] - self.x[0]) * abs(self.y[1] - self.y[0]) * abs(self.z[1] - self.z[0])

class Conf(object):

  def __init__(self, cnf, typ, *args, **kwargs):

    self.box = Box(cnf, typ)
    self.mols = []
    self.timestep = 0.0
    self.dump = None
    self.typ = typ

    if typ == 'dlpoly':
      self.dump = cnf['history']
      self.timestep = cnf['nstep']
      self.filename = cnf['filename']
      self.nb_atoms = cnf['natms']
      # The constructor for a dlpoly conf, can have one more parameters
      if len(args) < 1 or args[0] is None:
        mol_names = []
      else:
        mol_names = args[0]

      if 'raw' in cnf:
        # Atoms array are generated here for parallelisability
        cnf['atoms'] = dlpoly.History.raw2atoms(cnf['raw'], cnf['numpy'])

      # Ordinary style
      if not cnf['numpy']:
        # Index to read atom in conf
        i = 0
        for fmol in self.dump.field.mols:
          # Useless molecule
          if mol_names != [] and fmol['name'] not in mol_names:
            i += fmol['nummols'] * len(fmol['atoms'])
            continue
          for i_mol in xrange(fmol['nummols']):
            # Create a molecule via an object field
            mol = molecule.Molecule(name=fmol['name'], i=i_mol + 1)
            mol.field_init(self.dump.field)
            # No atom types in this case, because infos are already in atoms: TODO ?
            for atom in mol.atoms:
              atom.x = cnf['atoms'][i]['xxx']
              atom.y = cnf['atoms'][i]['yyy']
              atom.z = cnf['atoms'][i]['zzz']
              i += 1
            self.mols.append(mol)
      # Numpy style
      else:
        self.array = cnf['atoms']

        # Index to read atom in conf
        i = 0
        for fmol in self.dump.field.mols:
          # Useless molecule
          if mol_names != [] and fmol['name'] not in mol_names:
            i += fmol['nummols'] * len(fmol['atoms'])
            continue
          for i_mol in xrange(fmol['nummols']):
            # Create a molecule via an object field
            mol = molecule.Molecule(name=fmol['name'], i=i_mol + 1)
            mol.field_init(self.dump.field)
            for j in xrange(len(mol.atoms)):
              # Replace real atoms by NPAtoms
              atom = mol.atoms[j]
              # Use atom_types to have infos from FIELD
              atype = self.dump.atom_types[mol.name + '/' + atom.name]
              npatom = molecule.NPAtom(atom.name, i + 1, j + 1, self.array, atype=atype)
              mol.atoms[j] = npatom
              i += 1
            self.mols.append(mol)

    elif typ == 'lammps_data':
      self.filename = cnf.filename
      self.nb_atoms = len(cnf.atoms)
      i = 0
      i_glob = 0
      for dmol in cnf.mols:
        i += 1
        mol = molecule.Molecule(name='M%d' % i, i=i)
        i_mol = 0
        for datom in dmol['atoms']:
          i_glob += 1
          i_mol += 1
          atom = molecule.Atom('A%d' % i_mol, i_glob, i_mol, datom['x'], datom['y'], datom['z'], \
                               {'mass': datom['atom_type']['mass'], 'charge': datom['charge']})
          # Add an extra info to help build connections
          datom['atom'] = atom

          mol.atoms.append(atom)
        self.mols.append(mol)

      for dbond in cnf.bonds:
        atom1 = dbond['atom1']['atom']
        atom2 = dbond['atom2']['atom']
        if not atom1 in atom2.neighbs:
          atom2.neighbs.append(atom1)
        if not atom2 in atom1.neighbs:
          atom1.neighbs.append(atom2)

      # Delete the useless extra info
      for datom in cnf.atoms:
        del datom['atom']

    elif typ == 'lammps_dump':
      if 'dump' in cnf:
        self.dump = cnf['dump']
      else:
        self.dump = None
      self.timestep = cnf['timestep']
      self.filename = cnf['filename']
      self.nb_atoms = cnf['nbatoms']

      if 'raw' in cnf:
        # Atoms array are generated here for parallelisability
        cnf['atoms'] = lammps.Dump.raw2atoms(cnf['raw'], cnf['numpy'])

      # Ordinary style
      if not cnf['numpy']:
        mols = []
        # Create initial Mols (too many)
        for i in xrange(cnf['nbatoms']):
          mols.append(molecule.Molecule(name='M%d' % (i + 1), i=i + 1))

        i = 0
        for catom in cnf['atoms']:
          i += 1
          i_mol = catom['mol']
          # Alone mols (i_mol == 0) are created on the fly and added after the others
          if i_mol == 0:
            mol = molecule.Molecule(name='M0', i=0)
            mols.append(mol)
          else:
            mol = mols[i_mol - 1]
          if 'element' in catom:
            name = catom['element']
          else:
            name = 'A%s' % catom['id']
          # Add params
          params = {}

          for (key, value) in catom.items():
            if key in ('x', 'y', 'z', 'id'):
              continue
            params[key] = value
          x = y = z = None
          if 'x' in catom:
            x = catom['x']
          elif 'xs' in catom:
            x = catom['xs']
          elif 'xu' in catom:
            x = catom['xu']
          elif 'xsu' in catom:
            x = catom['xsu']
          if 'y' in catom:
            y = catom['y']
          elif 'ys' in catom:
            y = catom['ys']
          elif 'yu' in catom:
            y = catom['yu']
          elif 'ysu' in catom:
            y = catom['ysu']
          if 'z' in catom:
            z = catom['z']
          elif 'zs' in catom:
            z = catom['zs']
          elif 'zu' in catom:
            z = catom['zu']
          elif 'zsu' in catom:
            z = catom['zsu']
          # If good data are available in the dump, create Atom and AtomType
          if self.dump and self.dump.data and 'type' in catom:
            atype = self.dump.atom_types[str(catom['type'])]
            atom = molecule.Atom(name, catom['id'], len(mol.atoms) + 1, x, y, z, \
                                 params=params, atype=atype)
          else:
            atom = molecule.Atom(name, catom['id'], len(mol.atoms) + 1, x, y, z, \
                                 params=params)
          mol.atoms.append(atom)

      # Numpy style
      else:
        self.array = cnf['atoms']

        # No objects needed (only array)
        if 'no_objects' in kwargs and kwargs['no_objects']:
          return

        # Create initial Mols (too many)
        mols = []
        for i in xrange(cnf['nbatoms']):
          mols.append(molecule.Molecule(name='M%d' % (i + 1), i=i + 1))

        i = 0
        for catom in cnf['atoms']:
          i += 1
          i_mol = catom['mol']
          # Alone mols (i_mol == 0) reput the right i after
          if i_mol == 0:
            mol = molecule.Molecule(name='M0', i=0)
            mols.append(mol)
          else:
            mol = mols[i_mol - 1]
          if 'element' in catom.dtype.names:
            name = catom['element']
          else:
            name = 'A%s' % catom['id']
          # If good data are available in the dump, create Atom and AtomType
          if self.dump and self.dump.data and 'type' in catom.dtype.names:
            atype = self.dump.atom_types[str(catom['type'])]
            atom = molecule.NPAtom(name, i, len(mol.atoms) + 1, self.array, atype=atype)
          else:
            atom = molecule.NPAtom(name, i, len(mol.atoms) + 1, self.array)
          mol.atoms.append(atom)

      # Keep only real mols (and skip empty mols)
      self.mols = []
      for mol in mols:
        if len(mol.atoms) > 0:
          self.mols.append(mol)

      # Should I add bonds (with data only) ?
      if 'add_bonds' in kwargs and kwargs['add_bonds'] and self.dump and self.dump.data:

        # Put all atoms in a sorted array for easy access
        sorted_atoms = [[]] * self.nb_atoms
        for mol in self.mols:
          for atom in mol.atoms:
            sorted_atoms[atom.i - 1] = atom

        # Add bonds
        for bond in self.dump.data.bonds:
          a_i1, a_i2 = bond['atom1']['i'], bond['atom2']['i']
          a1 = sorted_atoms[a_i1 - 1]
          a2 = sorted_atoms[a_i2 - 1]
          a1.neighbs.append(a2)
          a2.neighbs.append(a1)

  def get_formulas(self):
    '''Get all the formulas of the mols in this conf'''
    formulas = []
    for mol in self.mols:
      form = molecule.Formula(mol=mol)
      if form not in formulas:
        formulas.append(form)
    return formulas

  def get_clusters(self, mol_names, r):
    '''Return the list of mols in the same cluster'''

    class Cluster(object):
      clusters = []

      def __init__(self, mol):
        self.mols = [mol]
        Cluster.clusters.append(self)

      def __str__(self):
        return ' '.join([mol.name for mol in self.mols])

      def __repr__(self):
        return ' '.join([mol.name for mol in self.mols])

      def merge(self, other):
        '''Merge two cluster into one'''

        self.mols += other.mols
        for mol in other.mols:
          mol.cluster = self
        Cluster.clusters.remove(other)

    # Init stuff
    clusterable_mols = []
    for mol in self.mols:
      if str(mol.formula) not in mol_names and mol.name not in mol_names:
        continue
      clusterable_mols.append(mol)
      mol.reconstruct(self.box)
      mol.cluster = Cluster(mol)

    for mol in clusterable_mols:
      for mol2 in clusterable_mols:
        if mol == mol2 or mol.cluster == mol2.cluster:
          continue

        if mol.near(mol2, self.box, 'prox_cm', {'min_dist': 0.0, 'max_dist': r}):
          mol.cluster.merge(mol2.cluster)

    # Return group of mols
    return [cluster.mols for cluster in Cluster.clusters]

  def recenter_z(self, center_mol):
    '''Try to recenter the mols on the center mol and along z
       The recenter can take several iterations or fail'''

    # Detect if the center mol is located on the PBC
    is_inside = False
    nb_slide = 0
    while not is_inside:
      min_z = max_z = None
      for mol in self.mols:
        if mol.name != center_mol and molecule.Formula(center_mol) != mol.formula:
          continue
        mol.reconstruct(self.box)
        if max_z is None or mol.cm[2] > max_z:
          max_z = mol.cm[2]
        if min_z is None or mol.cm[2] < min_z:
          min_z = mol.cm[2]

      slide_z = 0.0
      if (max_z - min_z) > 0.9 * self.box.lz:
        nb_slide += 1
        if nb_slide > 20:
          logging.warning('Cannot center mol %s, I keep this conf unchanged...', center_mol)
          return
        slide_z = 0.1 * self.box.lz
        logging.info('I slide z from %f to move from PBC', slide_z)

        for mol in self.mols:
          for atom in mol.atoms:
            atom.z += slide_z
            # PBC MGMT
            if self.typ == 'dlpoly':
              if atom.z < -self.box.lz / 2.0:
                atom.z += self.box.lz
              elif atom.z > self.box.lz / 2.0:
                atom.z -= self.box.lz
            else:
              if atom.z < self.box.z[0]:
                atom.z += self.box.lz
              elif atom.z > self.box.z[1]:
                atom.z -= self.box.lz
          mol._cm = None
      else:
        is_inside = True

    # Recenter along z to make inside mols center
    z_cm = 0.0
    nb_z = 0
    for mol in self.mols:
      if mol.name != center_mol and molecule.Formula(center_mol) != mol.formula:
        continue
      z_cm += mol.cm[2]
      # Reset cache
      mol._cm = None
      nb_z += 1

    z_cm = z_cm / nb_z
    if self.typ == 'dlpoly':
      # Dlpoly is centered on 0
      z_shift = z_cm
    else:
      z_shift = z_cm - self.box.z[0] - (self.box.z[1] - self.box.z[0]) / 2.0
    logging.debug('Conf shifted by %f', z_shift)

    for mol in self.mols:
      for atom in mol.atoms:
        atom.z -= z_shift
        # PBC MGMT
        if self.typ == 'dlpoly':
          if atom.z < -self.box.lz / 2.0:
            atom.z += self.box.lz
          elif atom.z > self.box.lz / 2.0:
            atom.z -= self.box.lz
        else:
          if atom.z < self.box.z[0]:
            atom.z += self.box.lz
          elif atom.z > self.box.z[1]:
            atom.z -= self.box.lz
