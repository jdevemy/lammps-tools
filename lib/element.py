# -*- coding: iso-8859-1 -*-
'''Python library to manage chemistry things'''

import csv

ELEMENTS_FILE = '/usr/local/data/elements.csv'

def name2symbol(name):
  '''Convert an atom name to the corresponding symbol'''

  if len(name) == 0 or name is None:
    raise ValueError('Empty name')

  if len(name) == 1:
    return name
  elif len(name) == 2:
    if name[1].islower() and name != 'Hw' and name != 'Ow':
      return name
    else:
      return name[0]
  else:
    if name[1].islower():
      if name[2].islower():
        return name[0] + name[1] + name[2]
      else:
        return name[0] + name[1]
    else:
      return name[0]

class Element(object):
  '''An element'''

  elements = None

  def __init__(self, name):

    # Read the data file if needed
    if not Element.elements:
      Element.elements = {}
      reader = csv.reader(open(ELEMENTS_FILE, 'r'), delimiter=',', quotechar='"')
      for row in reader:
        if row[0].startswith('#'):
          continue
        symbol = row[0]
        atomic_nb = int(row[1])
        mass = float(row[2])
        try:
          b_coeff = complex(row[3])
        except ValueError:
          b_coeff = None

        Element.elements[symbol] = {'atomic_nb': atomic_nb, 'mass': mass, 'b_coeff': b_coeff}

    if name in Element.elements:
      symbol = name
    else:
      symbol = name2symbol(name)
      if symbol not in Element.elements:
        raise ValueError('Not a correct symbol %s' % name)

    self.symbol = symbol
    self.atomic_nb = Element.elements[symbol]['atomic_nb']
    self.mass = Element.elements[symbol]['mass']
    self.b_coeff = Element.elements[symbol]['b_coeff']

  def __str__(self):
    return '%s %d %f' \
      % (self.symbol, self.atomic_nb, self.mass)

  def __repr__(self):
    return '%s' % (self.symbol)

