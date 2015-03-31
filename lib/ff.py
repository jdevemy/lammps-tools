# -*- coding: iso-8859-1 -*-
'''Python library to manage force fields'''

import logging

class AtomFF(object):
  '''A force field for an atom'''

  def __init__(self, atom_type, mass, charge, polarizability, thole, field_type, par, comment=''):

    # Force field is lj or q (lj special case)
    self.atom_type = atom_type
    self.mass = mass
    self.charge = charge
    self.polarizability = polarizability
    self.thole = thole
    self.field_type = field_type
    self.par = par
    self.comment = comment

    logging.debug('Add %s', self)

  def __str__(self):
    return 'Atom FF %s %f %f %f %f %s %s %s' \
      % (self.atom_type, self.mass, self.charge, self.polarizability, self.thole, self.field_type, self.par, self.comment)

  def __eq__(self, other):
    return self.atom_type == other.atom_type and self.mass == other.mass and self.charge == other.charge \
      and self.polarizability == other.polarizability and self.thole == other.thole and self.field_type == other.field_type

  def __ne__(self, other):
    return not(self.__eq__(other))

class BondFF(object):
  '''A force field for a bond'''

  def __init__(self, atom_type1, atom_type2, field_type, par, comment=''):

    # Force field is harm or cstr (constraint)
    self.atom_type1 = atom_type1
    self.atom_type2 = atom_type2
    self.field_type = field_type
    self.par = par
    self.comment = comment

    logging.debug('Add %s', self)

  def __str__(self):
    return 'Bond FF %s %s %s %s %s' \
      % (self.atom_type1, self.atom_type2, self.field_type, self.par, self.comment)

class AngleFF(object):
  '''A force field for an Angle'''

  def __init__(self, atom_type1, atom_type2, atom_type3, field_type, par, comment=''):

    self.atom_type1 = atom_type1
    self.atom_type2 = atom_type2
    self.atom_type3 = atom_type3
    self.field_type = field_type
    self.par = par
    self.comment = comment

    logging.debug('Add %s', self)

  def __str__(self):
    return 'Angle FF %s %s %s %s %s %s' \
      % (self.atom_type1, self.atom_type2, self.atom_type3, self.field_type, self.par, self.comment)

class DihedralFF(object):
  '''A force field for a Dihedral'''

  def __init__(self, atom_type1, atom_type2, atom_type3, atom_type4, field_type, par, elec_1_4, vdw_1_4, comment=''):

    self.atom_type1 = atom_type1
    self.atom_type2 = atom_type2
    self.atom_type3 = atom_type3
    self.atom_type4 = atom_type4
    self.field_type = field_type
    self.par = par
    self.elec_1_4 = elec_1_4
    self.vdw_1_4 = vdw_1_4
    self.comment = comment

    logging.debug('Add %s', self)

  def __str__(self):
    return 'Dihedral FF %s %s %s %s %s %s %s' \
      % (self.atom_type1, self.atom_type2, self.atom_type3, self.atom_type4, self.field_type, self.par, self.comment)

class ImproperFF(object):
  '''A force field for an Improper'''

  def __init__(self, atom_type1, atom_type2, atom_type3, atom_type4, field_type, par, elec_1_4, vdw_1_4, comment=''):

    self.atom_type1 = atom_type1
    self.atom_type2 = atom_type2
    self.atom_type3 = atom_type3
    self.atom_type4 = atom_type4
    self.field_type = field_type
    self.par = par
    self.elec_1_4 = elec_1_4
    self.vdw_1_4 = vdw_1_4
    self.comment = comment

    logging.debug('Add %s', self)

  def __str__(self):
    return 'Improper FF %s %s %s %s %s %s %s' \
      % (self.atom_type1, self.atom_type2, self.atom_type3, self.atom_type4, self.field_type, self.par, self.comment)


class ForceField(object):
  '''Class representing a force field'''

  def nff_init(self):
    '''Constructor for a .nff force fields'''

    f = open(self.filename, 'r')

    in_description = False

    for line in f:
      # Special case for multi lines Description
      if in_description:
        if line.startswith('/FF_DESCRIPTION'):
          in_description = False
          continue
        self.description += line
        continue

      if line.startswith('FF_NAME'):
        self.name = line.split('"')[1]
      elif line.startswith('FF_DESCRIPTION'):
        in_description = True
        self.description = ''
        continue

      if line.startswith('ENERGY_UNIT'):
        self.energy_unit = line.split('=')[1].strip().lower()

      if line.startswith('ATOM'):
        par = {}
        atom_type = line.split()[1]
        mass = float(line.split('MASS=')[1].split()[0])
        #TODO understand X
        try:
          charge = float(line.split('CHARGE=')[1].split()[0])
        except (ValueError, IndexError):
          charge = 0.0

        par['sigma'] = float(line.split('SIG=')[1].split()[0].replace(',', '.'))
        par['epsilon'] = float(line.split('EPS=')[1].split()[0].replace(',', '.'))

        try:
          comment = line.split('>')[1].strip()
        except IndexError:
          comment = ''

        self.atoms[atom_type] = AtomFF(atom_type, mass, charge, 0.0, 0.0, 'lj', par, comment)

      if line.startswith('BOND'):
        par = {}
        atom_type1 = line.split()[1]
        atom_type2 = line.split()[2]
        field_type = line.split('TYPE=')[1].split()[0].lower()
        par1 = float(line.split('PAR1=')[1].split()[0].replace(',', '.'))
        try:
          par2 = float(line.split('PAR2=')[1].split()[0].replace(',', '.'))
        except IndexError:
          par2 = None

        try:
          comment = line.split('>')[1].strip()
        except IndexError:
          comment = ''

        if field_type == 'harm':
          par['k'] = par1
          par['r'] = par2
        elif field_type == 'cstr':
          par['r'] = par1

        self.bonds[(atom_type1, atom_type2)] = BondFF(atom_type1, atom_type2, field_type, par, comment)
        self.bonds[(atom_type2, atom_type1)] = BondFF(atom_type2, atom_type1, field_type, par, comment)

      if line.startswith('ANGLE'):
        par = {}
        atom_type1 = line.split()[1]
        atom_type2 = line.split()[2]
        atom_type3 = line.split()[3]
        field_type = line.split('TYPE=')[1].split()[0].lower()
        par1 = float(line.split('PAR1=')[1].split()[0].replace(',', '.'))
        par2 = float(line.split('PAR2=')[1].split()[0].replace(',', '.'))

        try:
          comment = line.split('>')[1].strip()
        except IndexError:
          comment = ''
        if field_type == 'harm':
          par['k'] = par1
          par['theta'] = par2

        self.angles[(atom_type1, atom_type2, atom_type3)] = AngleFF(atom_type1, atom_type2, atom_type3, field_type, par, comment)
        self.angles[(atom_type3, atom_type2, atom_type1)] = AngleFF(atom_type3, atom_type2, atom_type1, field_type, par, comment)

      if line.startswith('DIHEDRAL'):
        par = {}
        atom_type1 = line.split()[1]
        atom_type2 = line.split()[2]
        atom_type3 = line.split()[3]
        atom_type4 = line.split()[4]
        field_type = line.split('TYPE=')[1].split()[0].lower()
        par1 = float(line.split('PAR1=')[1].split()[0].replace(',', '.'))
        par2 = float(line.split('PAR2=')[1].split()[0].replace(',', '.'))
        par3 = float(line.split('PAR3=')[1].split()[0].replace(',', '.'))
        # TODO manage VDW=X -> bond "-126"
        vdw_1_4 = float(line.split('VDW=')[1].split()[0].replace(',', '.'))
        elec_1_4 = float(line.split('ELEC=')[1].split()[0].replace(',', '.'))

        if field_type == 'cos':
          par['A'] = par1
          par['delta'] = par2
          par['m'] = par3

        #TODO manage multiple params for a same dihedral in .nff

        try:
          comment = line.split('>')[1].strip()
        except IndexError:
          comment = ''

        self.dihedrals[(atom_type1, atom_type2, atom_type3, atom_type4)] = DihedralFF(atom_type1, atom_type2, atom_type3, atom_type4, field_type, par, elec_1_4, vdw_1_4, comment)
        self.dihedrals[(atom_type4, atom_type3, atom_type2, atom_type1)] = DihedralFF(atom_type4, atom_type3, atom_type2, atom_type1, field_type, par, elec_1_4, vdw_1_4, comment)

    f.close()

  def ff_init(self):
    '''Constructor for a .ff force fields'''

    f = open(self.filename, 'r')

    # First line is the name of the ff
    self.name = f.readline().split('#')[1].strip()

    self.description = ''

    current_section = None

    for line in f:
      if line == '' or not line.startswith('#'):
        break
      self.description += line.split('#')[1].strip() + '\n'

    for line in f:

      if line.startswith('#') or line.strip() == '':
        continue

      # Section change
      if line.startswith('ATOMS'):
        current_section = 'ATOMS'
        continue
      elif line.startswith('TRANSLATION'):
        current_section = 'TRANSLATION'
        continue
      elif line.startswith('BONDS'):
        current_section = 'BONDS'
        continue
      elif line.startswith('ANGLES'):
        current_section = 'ANGLES'
        continue
      elif line.startswith('DIHEDRALS'):
        current_section = 'DIHEDRALS'
        continue
      elif line.startswith('IMPROPER'):
        current_section = 'IMPROPERS'
        continue

      if current_section == 'ATOMS':
        par = {}
        atom_type = line.split()[0]
        # Test if translations are Here (New version)
        try:
          mass = float(line.split()[1])
        except ValueError:
          # Translation
          atom_type2 = line.split()[1]
          self.translations[atom_type] = atom_type2

          shift = 1
          mass = float(line.split()[1 + shift])
        else:
          shift = 0

        charge = float(line.split()[2 + shift])
        field_type = line.split()[3 + shift].lower()
        # For field type q, no params
        if field_type == 'lj':
          par['sigma'] = float(line.split()[4 + shift])
          par['epsilon'] = float(line.split()[5 + shift])
          try:
            polarizability = float(line.split()[6 + shift])
          except (ValueError, IndexError):
            polarizability = 0.0
          try:
            thole = float(line.split()[7 + shift])
          except (ValueError, IndexError):
            thole = 0.0
        elif field_type == 'mart':
          par['sigma'] = float(line.split()[4 + shift])
          par['epsilon'] = float(line.split()[5 + shift])
          try:
            polarizability = float(line.split()[6 + shift])
          except (ValueError, IndexError):
            polarizability = 0.0
          try:
            thole = float(line.split()[7 + shift])
          except (ValueError, IndexError):
            thole = 0.0
        elif field_type == 'q':
          par['sigma'] = 0.0
          par['epsilon'] = 0.0
          field_type = 'lj'
          try:
            polarizability = float(line.split()[6 + shift])
          except (ValueError, IndexError):
            polarizability = 0.0
          try:
            thole = float(line.split()[7 + shift])
          except (ValueError, IndexError):
            thole = 0.0
        elif field_type == 'nm':
          par['Eo'] = float(line.split()[4 + shift])
          par['n'] = float(line.split()[5 + shift])
          par['m'] = float(line.split()[6 + shift])
          par['r0'] = float(line.split()[7 + shift])
          try:
            polarizability = float(line.split()[8 + shift])
          except (ValueError, IndexError):
            polarizability = 0.0
          try:
            thole = float(line.split()[9 + shift])
          except (ValueError, IndexError):
            thole = 0.0
        elif field_type == 'buck':
          par['A'] = float(line.split()[4 + shift])
          par['rho'] = float(line.split()[5 + shift])
          par['C'] = float(line.split()[6 + shift])
          try:
            polarizability = float(line.split()[7 + shift])
          except (ValueError, IndexError):
            polarizability = 0.0
          try:
            thole = float(line.split()[8 + shift])
          except (ValueError, IndexError):
            thole = 0.0

        self.atoms[atom_type] = AtomFF(atom_type, mass, charge, polarizability, thole, field_type, par)

      elif current_section == 'TRANSLATION':
        atom_type1 = line.split()[0]
        atom_type2 = line.split()[1]
        self.translations[atom_type2] = atom_type1

      elif current_section == 'BONDS':
        par = {}
        atom_type1 = line.split()[0]
        atom_type2 = line.split()[1]
        field_type = line.split()[2].lower()
        par1 = float(line.split()[3])
        par2 = float(line.split()[4])

        # Special case indicating a constraint with a - before k (!!!) or cons
        if par2 < 0:
          field_type = 'cstr'
          par['r'] = par1
          par['k'] = -par2
        elif field_type == 'cons':
          field_type = 'cstr'
          par['r'] = par1
          par['k'] = par2
        elif field_type == 'harm':
          par['r'] = par1
          par['k'] = par2
        else:
          raise ValueError('Unknown bond type %s' % field_type)

        self.bonds[(atom_type1, atom_type2)] = BondFF(atom_type1, atom_type2, field_type, par)
        self.bonds[(atom_type2, atom_type1)] = BondFF(atom_type2, atom_type1, field_type, par)

      elif current_section == 'ANGLES':
        par = {}
        atom_type1 = line.split()[0]
        atom_type2 = line.split()[1]
        atom_type3 = line.split()[2]
        field_type = line.split()[3].lower()
        par1 = float(line.split()[4])
        par2 = float(line.split()[5])

        if field_type == 'harm':
          par['theta'] = par1
          par['k'] = par2
        elif field_type == 'hcos':
          par['theta'] = par1
          par['k'] = par2
        elif field_type == 'cons':
          field_type = 'cstr'
          par['theta'] = par1
          par['k'] = par2
        else:
          raise ValueError('Unknown angle type %s' % field_type)

        self.angles[(atom_type1, atom_type2, atom_type3)] = AngleFF(atom_type1, atom_type2, atom_type3, field_type, par)
        self.angles[(atom_type3, atom_type2, atom_type1)] = AngleFF(atom_type3, atom_type2, atom_type1, field_type, par)

      elif current_section == 'DIHEDRALS':
        par = {}
        atom_type1 = line.split()[0]
        atom_type2 = line.split()[1]
        atom_type3 = line.split()[2]
        atom_type4 = line.split()[3]
        field_type = line.split()[4].lower()
        par1 = float(line.split()[5])
        par2 = float(line.split()[6])
        try:
          par3 = float(line.split()[7])
        except (IndexError, ValueError):
          par3 = None
        try:
          par4 = float(line.split()[8])
        except (IndexError, ValueError):
          par4 = None

        # Mix all in OPLS
        if field_type == 'cos3' or field_type == 'cos4' or field_type == 'opls':
          if field_type == 'cos3':
            par4 = 0.0
          par['A1'] = par1
          par['A2'] = par2
          par['A3'] = par3
          par['A4'] = par4
          field_type = 'opls'
        elif field_type == 'charmm':
          par['K'] = par1
          par['n'] = par2
          par['d'] = par3
        else:
          raise ValueError('Unknown dihedral type %s' % field_type)

        # Default 1-4 is 0.5
        elec_1_4 = 0.5
        vdw_1_4 = 0.5

        self.dihedrals[(atom_type1, atom_type2, atom_type3, atom_type4)] = DihedralFF(atom_type1, atom_type2, atom_type3, atom_type4, field_type, par, elec_1_4, vdw_1_4)
        self.dihedrals[(atom_type4, atom_type3, atom_type2, atom_type1)] = DihedralFF(atom_type4, atom_type3, atom_type2, atom_type1, field_type, par, elec_1_4, vdw_1_4)

      elif current_section == 'IMPROPERS':
        par = {}
        atom_type1 = line.split()[0]
        atom_type2 = line.split()[1]
        atom_type3 = line.split()[2]
        atom_type4 = line.split()[3]
        field_type = line.split()[4].lower()
        par1 = float(line.split()[5])
        par2 = float(line.split()[6])
        try:
          par3 = float(line.split()[7])
        except (IndexError, ValueError):
          par3 = None
        try:
          par4 = float(line.split()[8])
        except (IndexError, ValueError):
          par4 = None

        # Mix all in OPLS
        if field_type == 'cos3' or field_type == 'cos4' or field_type == 'opls':
          if field_type == 'cos3':
            par4 = 0.0
          par['A1'] = par1
          par['A2'] = par2
          par['A3'] = par3
          par['A4'] = par4
          field_type = 'opls'
        else:
          raise ValueError('Unknown improper type %s' % field_type)

        # Default 1-4 is 0.5
        elec_1_4 = 0.5
        vdw_1_4 = 0.5

        self.impropers[(atom_type1, atom_type2, atom_type3, atom_type4)] = ImproperFF(atom_type1, atom_type2, atom_type3, atom_type4, field_type, par, elec_1_4, vdw_1_4)
        self.impropers[(atom_type4, atom_type3, atom_type2, atom_type1)] = ImproperFF(atom_type4, atom_type3, atom_type2, atom_type1, field_type, par, elec_1_4, vdw_1_4)

    f.close()

  def __init__(self, filename):

    self.filename = filename

    self.atoms = {}

    self.bonds = {}

    self.angles = {}

    self.dihedrals = {}

    self.impropers = {}

    self.translations = {}

    self.description = ''

    self.name = ''

    self.energy_unit = 'kj/mol'

    if filename.endswith('.nff'):
      self.nff_init()
    elif filename.endswith('.ff'):
      self.ff_init()
    else:
      raise ValueError('Unknown force field file type %s' % filename)

  def __str__(self):
    return 'Force Field %s' \
      % (self.filename)

  def translate(self, atom_name):
    '''Translate the atom with translations from this ff'''

    if atom_name in self.translations:
      return self.translations[atom_name]
    else:
      return atom_name
