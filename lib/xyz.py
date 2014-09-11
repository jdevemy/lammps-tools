# -*- coding: iso-8859-1 -*-
'''Python library to manage xyz files'''

import logging

class Xyz(object):
  '''Class representing a xyz file'''

  def __init__(self, filename):

    self.filename = filename

    try:
      f = open(filename, 'r')

      self.nb_atoms = int(f.readline())

      self.comment = f.readline().strip()

      self.footer_comments = []

      self.atoms = []

      i = 0
      for line in f:
        i += 1
        name = line.split()[0]
        x = float(line.split()[1])
        y = float(line.split()[2])
        z = float(line.split()[3])

        try:
          comment = line.split('#')[1].strip()
        except IndexError:
          comment = ''

        atom = {}
        atom['name'] = name
        atom['x'] = x
        atom['y'] = y
        atom['z'] = z
        atom['comment'] = comment
        self.atoms.append(atom)

        if i == self.nb_atoms:
          break

      # Footer comments
      for line in f:
        if line.startswith('#'):
          self.footer_comments.append(line.split('#', 1)[1].strip())

      f.close()

    except (IndexError, ValueError, IOError), e:
      raise Exception('Not a valid .xyz file: %s' % e)

  def __str__(self):
    return 'XYZ %s(%s)' \
      % (self.filename, self.comment)

