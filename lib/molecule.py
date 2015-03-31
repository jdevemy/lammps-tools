# -*- coding: iso-8859-1 -*-
'''Python library to manage Molecules and Atoms in'''

import logging
import math
import re
import copy

import xyz
import mol
import dlpoly
import utils3d
import element

CACHE = {}

class AtomType(object):
  '''A type of atom'''

  def __init__(self, name, params={}):

    self.name = name
    self._params = params

    try:
      self.symbol = element.name2symbol(name)
    except ValueError:
      self.symbol = None

    self.fields = tuple(['name', 'symbol'] + self._params.keys())

    self.atoms = []

  def __getattr__(self, name):
    # For a strange bugfix !!!
    if name == '_params':
      raise AttributeError
    if name in self._params:
      return self._params[name]
    raise AttributeError

  def reset(self):
    '''Unlink to atoms (a new conf will arrive)'''
    self.atoms = []

  def factorize(self):
    '''Try to put in the AtomType common params of all atoms'''

    all_params = {}
    for atom in self.atoms:
      for (name, param) in atom._params.items():
        if name not in all_params:
          all_params[name] = set()
        all_params[name].update([param])
    for (name, a_p) in all_params.items():
      # Param has a uniq value, factorize it
      if len(a_p) == 1:
        self._params[name] = a_p.pop()
        self.fields = tuple(list(self.fields) + [name])
        for atom in self.atoms:
          del atom._params[name]

class BaseAtom(object):
  '''A base atom'''

  def __init__(self, name, i, i_mol, atype):

    self.name = name
    self.i = i
    self.i_mol = i_mol
    self.atype = atype
    if self.atype:
      self.atype.atoms.append(self)
    self.neighbs = []

class Atom(BaseAtom):
  '''An atom inside a molecule'''

  def __init__(self, name, i, i_mol, x, y, z, params={}, atype=None):

    super(Atom, self).__init__(name, i, i_mol, atype)
    self.x = x
    self.y = y
    self.z = z
    # Alias
    if 'q' in params:
      params['charge'] = params['q']
      del params['q']
    self._params = params

    if self.atype:
      self.fields = tuple(['name', 'i', 'i_mol', 'x', 'y', 'z'] + self._params.keys() + list(self.atype.fields))
    else:
      self.fields = tuple(['name', 'i', 'i_mol', 'x', 'y', 'z'] + self._params.keys())

    try:
      self.symbol = element.name2symbol(name)
    except ValueError:
      pass
    else:
      self.fields = tuple(list(self.fields) + ['symbol'])

  def __getattr__(self, name):
    # For a strange bugfix !!!
    if name == '_params':
      raise AttributeError
    # Alias
    if name == 'q':
      name = 'charge'
    if name in self._params:
      return self._params[name]
    elif self.atype:
      return getattr(self.atype, name)
    # Try to get mass from element
    elif name == 'mass':
      try:
        el = element.Element(self.symbol)
        return el.mass
      except ValueError:
        # Invalid symbol
        raise AttributeError
    else:
      raise AttributeError

  def __str__(self):
    if 'fr' in self.fields and 'to' in self.fields:
      return '%s(%s->%s) %d (%f %f %f) %s' \
        % (self.name, self.fr, self.to, self.i, self.x, self.y, self.z, self.neighbs)
    else:
      return '%s %d (%f %f %f) %s' \
        % (self.name, self.i, self.x, self.y, self.z, self.neighbs)

  def __repr__(self):
    return '%s %d' % (self.name, self.i)

  @property
  def coords(self):
    return (self.x, self.y, self.z)

class NPAtom(BaseAtom):
  '''An atom in numpy'''

  def __init__(self, name, i, i_mol, conf_array, atype=None):

    super(NPAtom, self).__init__(name, i, i_mol, atype)
    self.conf_array = conf_array

    if self.atype:
      self.fields = tuple(['name', 'i', 'i_mol'] + list(self.conf_array.dtype.names) + list(self.atype.fields))
    else:
      self.fields = tuple(['name', 'i', 'i_mol'] + list(self.conf_array.dtype.names))

    if 'q' in self.fields:
      self.fields += ('charge', )

    if 'element' in self.fields:
      self.fields += ('symbol', )

  def __str__(self):
    return '%s %d' \
      % (self.name, self.i)

  def __repr__(self):
    return '%s %d' % (self.name, self.i)

  def __getattr__(self, name):
    #Alias
    if name not in self.conf_array.dtype.names:
      if 'q' in self.conf_array.dtype.names and name == 'charge':
        name = 'q'
      elif 'element' in self.conf_array.dtype.names and name == 'symbol':
        name = 'element'
      elif 'xu' in self.conf_array.dtype.names and name == 'x':
        name = 'xu'
      elif 'yu' in self.conf_array.dtype.names and name == 'y':
        name = 'yu'
      elif 'zu' in self.conf_array.dtype.names and name == 'z':
        name = 'zu'

    if name in self.conf_array.dtype.names:
      return self.conf_array[self.i - 1][name]
    elif self.atype:
      return getattr(self.atype, name)
    else:
      raise AttributeError

class Formula(object):
  '''A formula for a molecule i.e. C2H6'''

  def __init__(self, form=None, mol=None):

    if form is not None:
      if 'formula_%s' % form in CACHE:
        self.fml = CACHE['formula_%s' % form].fml
        return

      self.fml = {'+': 0.}
      tokens = re.findall('[A-Z][a-z]*|[0-9]+|\+|\-|\(|\)|\[|\]', form)
      tokens.reverse()
      n = [1]
      num = 1
      # Put each token in a dict with the number as value
      for token in tokens:
        if token[0].isdigit():
          num = int(token)
        elif token[0] in '[(':
          n.pop()
        elif token[0] in ')]':
          n.append(num * n[-1])
          num = 1
        else:
          if token in self.fml:
            self.fml[token] += num * n[-1]
          else:
            self.fml[token] = num * n[-1]
          num = 1
      if '-' in self.fml:
        self.fml['+'] -= self.fml.pop('-')
      CACHE['formula_%s' % form] = self

    elif mol is not None:
      self.fml = {'+': 0.}
      for a in mol.atoms:
        if a.symbol in self.fml:
          self.fml[a.symbol] += 1
        else:
          self.fml[a.symbol] = 1
        if 'charge' in a.fields:
          self.fml['+'] += a.charge
      # Should be an integer charge
      self.fml['+'] = int(round(self.fml['+']))

  def __cmp__(self, other):
    return cmp(self.fml, other.fml)

  def __str__(self):
    s = ''
    keys = self.fml.keys()
    keys.remove('+')
    keys.sort()
    for elem in keys:
      if self.fml[elem] == 1:
        s += elem
      else:
        s += '%s%d' % (elem, self.fml[elem])
    if self.fml['+'] == -1:
      s += '-'
    elif self.fml['+'] < -1:
      s += '%d' % self.fml['+']
    elif self.fml['+'] == 1:
      s += '+'
    elif self.fml['+'] > 1:
      s += '+%d' % self.fml['+']
    return s

class Molecule(object):
  '''A molecule (with atoms inside)'''

  def xyz_init(self, mol_xyz):
    '''Constructor for a .xyz molecule'''

    self.name = mol_xyz.comment

    i = 1
    for atom in mol_xyz.atoms:
      a = Atom(atom['name'], i, i, atom['x'], atom['y'], atom['z'])
      self.atoms.append(a)
      i += 1

    # Second pass to create connections
    i = 1
    for atom in mol_xyz.atoms:
      if 'CONNECT' in atom['comment']:
        a = self.atoms[i - 1]
        i_neighbs = [int(n) for n in atom['comment'].split('CONNECT')[1].split()[0].split(',')]
        for i_n in i_neighbs:
          try:
            a.neighbs.append(self.atoms[i_n - 1])
          except KeyError:
            raise ValueError('Connection error between %s and %d', a, i_n)
      i += 1

    self.check_connections()

    # Then impropers
    for comment in mol_xyz.footer_comments:
      if comment.startswith('IMPROPER'):
        (i1, i2, i3, i4) = [int(i) for i in comment.split('IMPROPER')[1].split()[0].split(',')]
        self.impropers.append((self.atoms[i1 - 1], self.atoms[i2 - 1], \
                               self.atoms[i3 - 1], self.atoms[i4 - 1]))

    # Check if there's no alone atoms (i.e. it's not a molecule)
    for atom in self.atoms:
      if len(atom.neighbs) == 0:
        logging.info('Poor alone atom detected %s in %s, there is more than one molecule in %s (or you provide a .xyz file without CONNECT)', repr(atom), repr(self), self.filename)

  def mol_init(self, mol_mol):
    '''Constructor for a .mol molecule'''

    self.name = mol_mol.name

    i = 1
    for atom in mol_mol.atoms:
      if i == 1:
        at = Atom(atom['name'], i, i, 0.0, 0.0, 0.0)
      elif i == 2:
        dist = atom['dist']
        dist_atom = self.atoms[atom['dist_atom']['i'] - 1]
        at = Atom(atom['name'], i, i, dist, 0.0, 0.0)
        at.neighbs.append(dist_atom)
        dist_atom.neighbs.append(at)
      elif i == 3:
        dist = atom['dist']
        dist_atom = self.atoms[atom['dist_atom']['i'] - 1]
        angle = math.radians(atom['angle'])
        angle_atom = self.atoms[atom['angle_atom']['i'] - 1]

        # For this construction, the new atom is at point (x, y), atom ir
        # is at point (xr, yr) and atom ia is at point (xa, ya).  Theta
        # is the angle between the vector joining ir to ia and the
        # x-axis, a' (= theta - a) is is the angle between r and the
        # x-axis. x = xa + r cos a', y = ya + r sin a'.  From the dot
        # product of a unitary vector along x with the vector from ir to
        # ia, theta can be calculated: cos theta = (xa - xr) / sqrt((xa -
        # xr)^2 + (ya - yr)^2).  If atom ia is in third or forth quadrant
        # relative to atom ir, ya - yr < 0, then theta = 2 pi - theta.

        # delx = c[ia].x - c[ir].x;
        # dely = c[ia].y - c[ir].y;
        # theta = acos(delx / sqrt(delx * delx + dely * dely));
        # if (dely < 0.0)
        #   theta = 2.0 * M_PI - theta;
        # a = theta - a;
        # c[2].x = c[ir].x + r * cos(a);
        # c[2].y = c[ir].y + r * sin(a);
        # c[2].z = 0.0;

        delx = angle_atom.x - dist_atom.x
        dely = angle_atom.y - dist_atom.y
        theta = math.acos(delx / math.sqrt(delx * delx + dely * dely))
        if dely < 0.0:
          theta = 2.0 * math.pi - theta
        angle = theta - angle

        x = dist_atom.x + dist * math.cos(angle)
        y = dist_atom.y + dist * math.sin(angle)
        z = 0

        at = Atom(atom['name'], i, i, x, y, 0.0)
        at.neighbs.append(dist_atom)
        dist_atom.neighbs.append(at)
      else:
        dist = atom['dist']
        dist_atom = self.atoms[atom['dist_atom']['i'] - 1]
        angle = math.radians(atom['angle'])
        angle_atom = self.atoms[atom['angle_atom']['i'] - 1]
        dihedr = math.radians(atom['dihedr'])
        dihedr_atom = self.atoms[atom['dihedr_atom']['i'] - 1]

        # For this construction the new atom is at point A, atom ir
        # is at B, atom ia at C and atom id at D.  Point a is the
        # projection of A onto the plane BCD.  Point b is the
        # projection of A along the direction CB (the line defining
        # the dihedral angle between planes ABC and BCD). n = CD x CB
        # / |CD x CB| is the unit vector normal to the plane BCD. m =
        # CB x n / |CB x n| is the unit vector on the plane BCD
        # normal to the direction CB.
        #
        #                              .'A
        #                ------------.' /.-----------------
        #               /           b /  .               /
        #              /           ./    .              /
        #             /           B......a      ^      /
        #            /           /              |n    /
        #           /           /                    /
        #          /           C                    /
        #         /             \                  /
        #        /               \                /
        #       /plane BCD        D              /
        #      ----------------------------------
        #
        #             A              C------B...b
        #            /.             /        .  .
        #           / .            /    |m    . .
        #          /  .           /     V      ..
        #  C------B...b          D              a

        # BA = r;
        # vB = v_new(c[ir].x, c[ir].y, c[ir].z);
        # vC = v_new(c[ia].x, c[ia].y, c[ia].z);
        # vD = v_new(c[id].x, c[id].y, c[id].z);
        # vCB = v_sub(vC, vB);
        # vCD = v_sub(vD, vC);
        # vb = v_sub(vC, v_timesr((v_abs(vCB) - BA * cos(a)) / v_abs(vCB), vCB));
        # bA = BA * sin(a);
        # aA = bA * sin(d);
        # ba = bA * cos(d);
        # vn = v_unit(v_cross(vCD, vCB));
        # vm = v_unit(v_cross(vCB, vn));
        # va = v_add(vb, v_timesr(ba, vm));
        # vA = v_add(va, v_timesr(aA, vn));
        # c[i].x = vA.x;
        # c[i].y = vA.y;
        # c[i].z = vA.z;

        BA = dist
        vB = utils3d.Vector3d(dist_atom.x, dist_atom.y, dist_atom.z)
        vC = utils3d.Vector3d(angle_atom.x, angle_atom.y, angle_atom.z)
        vD = utils3d.Vector3d(dihedr_atom.x, dihedr_atom.y, dihedr_atom.z)
        vCB = vC - vB
        vCD = vD - vC
        vb = vC - vCB.times((vCB.length() - BA * math.cos(angle)) / vCB.length())
        bA = BA * math.sin(angle)
        aA = bA * math.sin(dihedr)
        ba = bA * math.cos(dihedr)

        # If 3 points are aligned: use an arbitrary normal vector
        try:
          vn = vCD.cross(vCB).unit()
        except ZeroDivisionError:
          vn = utils3d.Vector3d(0, 0, 1.0)
        vm = vCB.cross(vn).unit()
        va = vb + vm.times(ba)
        vA = va + vn.times(aA)

        x = vA.to_tuple()[0]
        y = vA.to_tuple()[1]
        z = vA.to_tuple()[2]

        at = Atom(atom['name'], i, i, x, y, z)
        at.neighbs.append(dist_atom)
        dist_atom.neighbs.append(at)

      self.atoms.append(at)
      i += 1

    # Manage extra connection
    for atom in mol_mol.atoms:
      at = self.atoms[atom['i'] - 1]
      for connect in atom['extra_connect']:
        neighb = self.atoms[connect['i'] - 1]
        if neighb not in at.neighbs:
          at.neighbs.append(neighb)

    # Manage impropers
    for (a1, a2, a3, a4) in mol_mol.impropers:
      self.impropers.append((self.atoms[a1['i'] - 1], self.atoms[a2['i'] - 1], \
                             self.atoms[a3['i'] - 1], self.atoms[a4['i'] - 1]))

  def field_init(self, field):
    '''Constructor for a molecule in a FIELD'''

    # In these case the coordinates stay (0, 0, 0), we only get connections

    the_mol = None
    for mol_ in field.mols:
      if self.name != mol_['name']:
        continue
      else:
        the_mol = mol_
    if not the_mol:
      raise ValueError('Unknown molecule %s in FIELD %s' % (self.name, self.filename))

    i = 1
    for atom in the_mol['atoms']:
      # Add params infos
      params = {}
      params['frozen'] = atom['frozen']
      params['igrp'] = atom['igrp']
      params['charge'] = atom['chge']
      params['mass'] = atom['weight']
      a = Atom(atom['name'], i, i, 0.0, 0.0, 0.0, params)
      self.atoms.append(a)
      i += 1

    for bond in the_mol['bonds']:
      i1 = bond['atom_i1']
      i2 = bond['atom_i2']
      a1 = self.atoms[i1 - 1]
      a2 = self.atoms[i2 - 1]
      a1.neighbs.append(a2)
      a2.neighbs.append(a1)

    for rigid in the_mol['rigids']:
      # Create connect along the rigid unit
      i1 = rigid[0]
      for i2 in rigid[1:]:
        a1 = self.atoms[i1 - 1]
        a2 = self.atoms[i2 - 1]
        a1.neighbs.append(a2)
        a2.neighbs.append(a1)
        i1 = i2

    for const in the_mol['constraints']:
      i1 = const['atom_i1']
      i2 = const['atom_i2']
      a1 = self.atoms[i1 - 1]
      a2 = self.atoms[i2 - 1]
      a1.neighbs.append(a2)
      a2.neighbs.append(a1)

    for dihedral in the_mol['dihedrals']:
      i1 = dihedral['atom_i1']
      i2 = dihedral['atom_i2']
      i3 = dihedral['atom_i3']
      i4 = dihedral['atom_i4']
      a1 = self.atoms[i1 - 1]
      a2 = self.atoms[i2 - 1]
      a3 = self.atoms[i3 - 1]
      a4 = self.atoms[i4 - 1]
      # Impropers are dihedral not along connected
      if a2 not in a1.neighbs or a3 not in a2.neighbs or a4 not in a3.neighbs:
        self.impropers.append((a1, a2, a3, a4))

  def __init__(self, filename=None, typ=None, name='', i=0):

    self.filename = filename

    self.atoms = []

    self.impropers = []

    self.name = name

    self.i = i

    self._formula = None
    self._cm = None

    if typ == 'xyz':
      mol_xyz = xyz.Xyz(self.filename)
      self.xyz_init(mol_xyz)
    elif typ == 'mol' or typ == 'zmat':
      mol_mol = mol.Mol(self.filename)
      self.mol_init(mol_mol)
    elif typ == 'FIELD':
      field = dlpoly.Field(self.filename)
      self.field_init(field)
    elif typ != None:
      raise ValueError('Unkown molecule type %s' % typ)

  def __str__(self):
    return 'Molecule %s %s : %s' \
      % (self.name, self.filename, ' '.join([str(a) for a in self.atoms]))

  def __repr__(self):
    return '%s/%s' % (self.name, self.filename)

  @property
  def formula(self):
    if self._formula is None:
      self._formula = Formula(mol=self)
    return self._formula

  @property
  def cm(self):
    '''Get the center of mass of a molecule'''

    if self._cm is None:
      x = y = z = 0.0
      mass_total = 0.0
      for atom in self.atoms:
        if 'mass' in atom.fields:
          mass = atom.mass
        else:
          # Get the mass from the element if unavailable
          try:
            mass = element.Element(atom.symbol).mass
          except ValueError:
            mass = 1.0
        x += atom.x * mass
        y += atom.y * mass
        z += atom.z * mass
        mass_total += mass
      x /= mass_total
      y /= mass_total
      z /= mass_total
      self._cm = (x, y, z)
    return self._cm

  def xyz_write(self, filename):
    '''Write the molecule into a .xyz file'''

    f = open(filename, 'w')

    f.write('%d\n'% len(self.atoms))

    f.write('%s\n'% self.name)

    for atom in self.atoms:
      if len(atom.neighbs):
        comment = ','.join([str(neighb.i) for neighb in atom.neighbs])
        f.write('%s %f %f %f # CONNECT %s\n' % (atom.name, atom.x, atom.y, atom.z, comment))
      else:
        f.write('%s %f %f %f\n' % (atom.name, atom.x, atom.y, atom.z))

    for improper in self.impropers:
      f.write('# IMPROPER %d,%d,%d,%d\n' % (improper[0].i, improper[1].i, improper[2].i, improper[3].i))

    f.close()

  def write(self, filename, typ):
    '''Write the molecule into a file'''

    if typ == 'xyz':
      self.xyz_write(filename)
    else:
      raise ValueError('Unknown molecule file type %s' % typ)

  def remove_atom(self, i):
    '''Remove one atom from the molecule'''

    atom2remove = self.atoms[i - 1]
    self.atoms.remove(atom2remove)

    for atom in self.atoms:
      if atom.i > i:
        atom.i -= 1
      if atom2remove in atom.neighbs:
        atom.neighbs.remove(atom2remove)

  def check_connections(self):
    '''Check if connections are coherent'''

    for a in self.atoms:
      for neighb in a.neighbs:
        if a not in neighb.neighbs:
          raise ValueError('Incoherent connections')

  def add_connections(self, max_r=2.0):
    '''Add connections from max radius'''

    for atom in self.atoms:
      for atom2 in self.atoms:
        if atom2 in atom.neighbs or atom is atom2:
          continue
        dx = atom2.x - atom.x
        dy = atom2.y - atom.y
        dz = atom2.z - atom.z
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        if r < max_r:
          atom.neighbs.append(atom2)
          atom2.neighbs.append(atom)

  def fusion_multi(self, max_r=0.5):
    '''Fusion atoms which are (quite) at the same place'''

    atoms2remove_i = []
    for atom in self.atoms:
      for atom2 in self.atoms:
        if atom is atom2 or atom.i in atoms2remove_i or atom2.i in atoms2remove_i or atom.name != atom2.name:
          continue
        dx = atom2.x - atom.x
        dy = atom2.y - atom.y
        dz = atom2.z - atom.z
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        if r < max_r:
          atom.x = (atom.x + atom2.x) / 2.0
          atom.y = (atom.y + atom2.y) / 2.0
          atom.z = (atom.z + atom2.z) / 2.0
          atoms2remove_i.append(atom2.i)
          for atom3 in atom2.neighbs:
            if atom3 not in atom.neighbs:
              atom.neighbs.append(atom3)

    atoms2remove_i.sort(reverse=True)
    for i in atoms2remove_i:
      self.remove_atom(i)

  def add_mol(self, other):
    nb = len(self.atoms)
    for atom in other.atoms:
      atom.i += nb
      self.atoms.append(atom)

  def duplicate(self, x=0.0, y=0.0, z=0.0):
    '''Duplicate the current molecule with moving it'''

    copy_atoms = []
    for atom in self.atoms:
      copy_atom = copy.copy(atom)
      copy_atom.neighbs = []
      copy_atom.x += x
      copy_atom.y += y
      copy_atom.z += z
      copy_atom.i += len(self.atoms)
      copy_atoms.append(copy_atom)
    self.atoms += copy_atoms

  def get_all_bonds(self):
    '''Return all the bonds of the molecules'''

    bonds = []

    for atom in self.atoms:
      for neighb in atom.neighbs:
        bond = (atom, neighb)
        bond_rev = list(bond)
        bond_rev.reverse()
        bond_rev = tuple(bond_rev)
        if bond not in bonds and bond_rev not in bonds:
          bonds.append(bond)

    return bonds

  def get_all_angles(self):
    '''Return all the angles of the molecules'''

    angles = []

    for atom in self.atoms:
      for neighb in atom.neighbs:
        for neighb2 in neighb.neighbs:
          if neighb2 == atom:
            continue
          angle = (atom, neighb, neighb2)
          angle_rev = list(angle)
          angle_rev.reverse()
          angle_rev = tuple(angle_rev)
          if angle not in angles and angle_rev not in angles:
            angles.append(angle)

    return angles

  def get_all_dihedrals(self):
    '''Return all the dihedrals of the molecules'''

    dihedrals = []

    for atom in self.atoms:
      for neighb in atom.neighbs:
        for neighb2 in neighb.neighbs:
          if neighb2 == atom:
            continue
          for neighb3 in neighb2.neighbs:
            if neighb3 == atom or neighb3 == neighb:
              continue
            dihedral = (atom, neighb, neighb2, neighb3)
            dihedral_rev = list(dihedral)
            dihedral_rev.reverse()
            dihedral_rev = tuple(dihedral_rev)
            if dihedral not in dihedrals and dihedral_rev not in dihedrals:
              dihedrals.append(dihedral)

    return dihedrals

  def get_all_impropers(self):
    '''Return all the impropers of the molecules (for coherence)'''
    return self.impropers

  def reconstruct(self, box):
    '''Reconstruct molecule with periodic conditions of box'''

    box_len_x_2 = box.len_x / 2.0
    box_len_y_2 = box.len_y / 2.0
    box_len_z_2 = box.len_z / 2.0

    with_neighbs = False
    for atom in self.atoms:
      atom.is_reco = False
      atom.x_shift = 0.0
      atom.y_shift = 0.0
      atom.z_shift = 0.0
      if atom.neighbs != []:
        with_neighbs = True

    for atom in self.atoms:
      atom.is_reco = True

      # Should I use neighbs, or do I take all atoms ?
      if with_neighbs:
        neighbs = atom.neighbs
      else:
        neighbs = self.atoms

      for neighb in neighbs:
        if neighb.is_reco:
          continue
        neighb.is_reco = True

        # Reco the neighbour from the root atom
        d_x = neighb.x - atom.x
        d_y = neighb.y - atom.y
        d_z = neighb.z - atom.z

        neighb.x_shift += atom.x_shift
        neighb.y_shift += atom.y_shift
        neighb.z_shift += atom.z_shift

        if abs(d_x) > box_len_x_2:
          neighb.x_shift -= math.copysign(box.len_x, d_x)
        if abs(d_y) > box_len_y_2:
          neighb.y_shift -= math.copysign(box.len_y, d_y)
        if abs(d_z) > box_len_z_2:
          neighb.z_shift -= math.copysign(box.len_z, d_z)

    for atom in self.atoms:
      atom.x += atom.x_shift
      atom.y += atom.y_shift
      atom.z += atom.z_shift

  def rotate(self, origin, angle_x=0.0, angle_y=0.0, angle_z=0.0):
    '''Rotate the whole molecule with an angle x, y, and z'''

    for atom in self.atoms:
      v = utils3d.Vector3d(atom.x - origin[0], atom.y - origin[1], atom.z - origin[2])
      v.rotate(angle_x, angle_y, angle_z)
      atom.x = origin[0] + v.x
      atom.y = origin[1] + v.y
      atom.z = origin[2] + v.z

  def rotate_vect(self, origin, v_rot, angle):
    '''Rotate the whole molecule along a vector'''

    for atom in self.atoms:
      v = utils3d.Vector3d(atom.x - origin[0], atom.y - origin[1], atom.z - origin[2])
      v.rotate_vect(v_rot, angle)
      atom.x = origin[0] + v.x
      atom.y = origin[1] + v.y
      atom.z = origin[2] + v.z

  def near(self, other, box, typ, params):
    '''Check if an other molecule is near or not and return the connexion type'''

    if typ == 'hbond':
      # Self is acceptor
      # Other is donor
      acceptors = params['acceptors'][self.name]['ranks']
      donors = params['donors'][other.name]['ranks']
      for acceptor in acceptors:
        a = self.atoms[acceptor - 1]
        for donor in donors:
          d = other.atoms[donor[0] - 1]
          h = other.atoms[donor[1] - 1]

          # Length criterion
          v_ha = utils3d.Vector3d(a.x - h.x, a.y - h.y, a.z - h.z)

          # Optimization
          if (v_ha.x > params['max_dist'] and v_ha.x < box.lx/2.0) or \
             (v_ha.y > params['max_dist'] and v_ha.y < box.ly/2.0) or \
             (v_ha.z > params['max_dist'] and v_ha.z < box.lz/2.0):
            continue

          if v_ha.x > box.lx/2.0:
            v_ha.x -= box.lx
          if v_ha.y > box.ly/2.0:
            v_ha.y -= box.ly
          if v_ha.z > box.lz/2.0:
            v_ha.z -= box.lz
          if v_ha.x < -box.lx/2.0:
            v_ha.x += box.lx
          if v_ha.y < -box.ly/2.0:
            v_ha.y += box.ly
          if v_ha.z < -box.lz/2.0:
            v_ha.z += box.lz

          # Optimization
          if v_ha.x > params['max_dist'] or v_ha.y > params['max_dist'] or v_ha.z > params['max_dist']:
            continue

          length = v_ha.length()
          if length < params['min_dist'] or length > params['max_dist']:
            continue

          v_hd = utils3d.Vector3d(d.x - h.x, d.y - h.y, d.z - h.z)

          if v_hd.x > box.lx/2.0:
            v_hd.x -= box.lx
          if v_hd.y > box.ly/2.0:
            v_hd.y -= box.ly
          if v_hd.z > box.lz/2.0:
            v_hd.z -= box.lz
          if v_hd.x < -box.lx/2.0:
            v_hd.x += box.lx
          if v_hd.y < -box.ly/2.0:
            v_hd.y += box.ly
          if v_hd.z < -box.lz/2.0:
            v_hd.z += box.lz

          # Angle criterion
          cross = v_ha.unit().dot(v_hd.unit())
          angle = math.degrees(math.acos(cross))
          if angle < params['min_angle'] or angle > params['max_angle']:
            continue
          return (d.name, h.name, a.name)

      # No hbond found
      return False
    elif typ == 'prox':
      for s in self.atoms:
        for o in other.atoms:
          v = utils3d.Vector3d(s.x - o.x, s.y - o.y, s.z - o.z)
          if v.x > box.lx/2.0:
            v.x -= box.lx
          if v.y > box.ly/2.0:
            v.y -= box.ly
          if v.z > box.lz/2.0:
            v.z -= box.lz
          if v.x < -box.lx/2.0:
            v.x += box.lx
          if v.y < -box.ly/2.0:
            v.y += box.ly
          if v.z < -box.lz/2.0:
            v.z += box.lz

          length = v.length()
          if length > params['min_dist'] and length < params['max_dist']:
            return True
      # No proximity detected
      return False
    elif typ == 'prox_cm':
      v = utils3d.Vector3d(self.cm[0] - other.cm[0], self.cm[1] - other.cm[1], self.cm[2] - other.cm[2])
      if v.x > box.lx/2.0:
        v.x -= box.lx
      if v.y > box.ly/2.0:
        v.y -= box.ly
      if v.z > box.lz/2.0:
        v.z -= box.lz
      if v.x < -box.lx/2.0:
        v.x += box.lx
      if v.y < -box.ly/2.0:
        v.y += box.ly
      if v.z < -box.lz/2.0:
        v.z += box.lz

      length = v.length()
      if length > params['min_dist'] and length < params['max_dist']:
        return True
      else:
        return False
    else:
      return False

  def fusion(self, other, ff):
    '''Fusion 2 molecules into one, with dummy atoms for the second one'''

    def continue_fusion(ff, mol1, mol2, i1, i2, already_fusions1, already_fusions2):
      '''Recursive function to search for fusion point between two mols with a fusion point'''
      atom1 = mol1.atoms[i1 - 1]
      atom2 = mol2.atoms[i2 - 1]
      logging.debug('Already fusion %s %s', already_fusions1, already_fusions2)
      l = [(i1, i2)]
      # Search for neigbors affinity
      for n1 in atom1.neighbs:
        # Avoid previous fusion point
        if n1.i_mol in already_fusions1:
          continue
        for n2 in atom2.neighbs:
          # Avoid previous fusion point
          if n2.i_mol in already_fusions2:
            continue
          # Found fusion point
          # TODO check if the two fusion points are located "near" !!!
          if n1.symbol == n2.symbol and ff.atoms[n1.name].mass == ff.atoms[n2.name].mass:
            logging.debug('Found fusion on %d %d', n1.i_mol, n2.i_mol)
            already_fusions1 += [n1.i_mol]
            already_fusions2 += [n2.i_mol]
            l += continue_fusion(ff, mol1, mol2, n1.i_mol, n2.i_mol, already_fusions1, already_fusions2)
            break
      return l

    atoms_symb_mass_1 = set([(atom.symbol, ff.atoms[atom.name].mass) for atom in self.atoms])
    atoms_symb_mass_2 = set([(atom.symbol, ff.atoms[atom.name].mass) for atom in other.atoms])

    common_atoms = atoms_symb_mass_1.intersection(atoms_symb_mass_2)

    if common_atoms == set():
      raise ValueError('Unable to fusion molecules: no common atoms')

    best_fusion = None

    for (symbol, mass) in common_atoms:

      # Search for fusion point in mol1
      mol1_fusion_points = []
      for atom in self.atoms:
        if atom.symbol == symbol and ff.atoms[atom.name].mass == mass:
          mol1_fusion_points.append(atom.i_mol)

      # Search for fusion point in mol2
      mol2_fusion_points = []
      for atom in other.atoms:
        if atom.symbol == symbol and ff.atoms[atom.name].mass == mass:
          mol2_fusion_points.append(atom.i_mol)

      # Try the fusion and search for the best one
      for point1 in mol1_fusion_points:
        for point2 in mol2_fusion_points:
          logging.debug('Try fusion point %d %d', point1, point2)
          l = continue_fusion(ff, self, other, point1, point2, [point1], [point2])
          logging.debug('Result: %s', l)
          if best_fusion is None or len(l) > len(best_fusion[1]):
            best_fusion = ((point1, point2), l)

    fusion_point = best_fusion[0]
    dx = self.atoms[fusion_point[1] - 1].x - self.atoms[fusion_point[0] - 1].x
    dy = self.atoms[fusion_point[1] - 1].y - self.atoms[fusion_point[0] - 1].y
    dz = self.atoms[fusion_point[1] - 1].z - self.atoms[fusion_point[0] - 1].z

    # Create the fusion mol
    fusion_mol = Molecule(filename='%s-%s' % (self.filename, other.filename), name='%s-%s' % (self.name, other.name))

    fusioned1 = [f[0] for f in best_fusion[1]]

    fusioned2 = [f[1] for f in best_fusion[1]]

    i = 0
    # Create common atoms (with only need atom1 ref)
    for (i1, i2) in best_fusion[1]:
      i += 1
      atom1 = self.atoms[i1 - 1]
      atom2 = other.atoms[i2 - 1]
      atomf = Atom(atom1.name, i, 1, atom1.x, atom1.y, atom1.z, params={'fr': atom1.name, 'to': atom2.name, 'i1': i1, 'i2': i2})
      fusion_mol.atoms.append(atomf)

    i = 0
    # Create common atoms neighbouring (with only atom1 ref)
    for (i1, i2) in best_fusion[1]:
      i += 1
      atom1 = self.atoms[i1 - 1]
      atomf = fusion_mol.atoms[i - 1]
      for neighb1 in atom1.neighbs:
        if neighb1.i_mol in fusioned1:
          atomf.neighbs.append(fusion_mol.atoms[fusioned1.index(neighb1.i_mol)])

    # Add atoms not in common (included dummy)
    while len(fusioned1) < len(self.atoms) or len(fusioned2) < len(other.atoms):
      for atomf in fusion_mol.atoms:
        if atomf.i1 is not None:
          atom1 = self.atoms[atomf.i1 - 1]
          for neighb1 in atom1.neighbs:
            if neighb1.i_mol not in fusioned1:
              atom = Atom(neighb1.name, len(fusion_mol.atoms) + 1, 1, neighb1.x, neighb1.y, neighb1.z, params={'fr': neighb1.name, 'to': None, 'i1': neighb1.i_mol, 'i2': None})
              fusion_mol.atoms.append(atom)
              atom.neighbs.append(atomf)
              atomf.neighbs.append(atom)
              fusioned1.append(neighb1.i_mol)

        if atomf.i2 is not None:
          atom2 = other.atoms[atomf.i2 - 1]
          for neighb2 in atom2.neighbs:
            if neighb2.i_mol not in fusioned2:
              atom = Atom(neighb2.name, len(fusion_mol.atoms) + 1, 1, neighb2.x - dx, neighb2.y - dy, neighb2.z - dz, params={'fr': None, 'to': neighb2.name, 'i1': None, 'i2': neighb2.i_mol,})
              fusion_mol.atoms.append(atom)
              atom.neighbs.append(atomf)
              atomf.neighbs.append(atom)
              fusioned2.append(neighb2.i_mol)

    # Remove useless infos and Change names for transformation atoms
    for atom in fusion_mol.atoms:
      del atom._params['i1']
      del atom._params['i2']
      atom.fields = list(atom.fields)
      atom.fields.remove('i1')
      atom.fields.remove('i2')
      atom.fields = tuple(atom.fields)

    return fusion_mol
