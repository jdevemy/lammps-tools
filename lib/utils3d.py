# -*- coding: iso-8859-1 -*-
'''Python library with utils for 3d stuff'''

import math

class Vector3d(object):
  '''A vector in 3d'''

  def __init__(self, *coords):
    if len(coords) == 3:
      self.x = coords[0]
      self.y = coords[1]
      self.z = coords[2]
    else:
      self.x = coords[0][0]
      self.y = coords[0][1]
      self.z = coords[0][2]

  def __str__(self):
    return '%f %f %f' % (self.x, self.y, self.z)

  def __repr__(self):
    return '%f %f %f' % (self.x, self.y, self.z)

  def to_tuple(self):
    return (self.x, self.y, self.z)

  def __add__(self, other):
    x = self.x + other.x
    y = self.y + other.y
    z = self.z + other.z
    return Vector3d(x, y, z)

  def __sub__(self, other):
    x = self.x - other.x
    y = self.y - other.y
    z = self.z - other.z
    return Vector3d(x, y, z)

  def cross(self, other):
    '''Cross product'''
    x = self.y * other.z - self.z * other.y
    y = self.z * other.x - self.x * other.z
    z = self.x * other.y - self.y * other.x
    return Vector3d(x, y, z)

  def dot(self, other):
    '''Dot product'''
    return self.x * other.x + self.y * other.y + self.z * other.z

  def times(self, number):
    '''Product with a number'''
    x = self.x * number
    y = self.y * number
    z = self.z * number
    return Vector3d(x, y, z)

  def unit(self):
    '''Unit a vector'''
    return self.times(1/self.length())

  def length(self):
    return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

  def change_axis(self, v1, v2, v3):
    # Axis change (v1, v2, v3 are the new axis)
    new_x = self.dot(v1)
    new_y = self.dot(v2)
    new_z = self.dot(v3)
    self.x = new_x
    self.y = new_y
    self.z = new_z
    return self

  def rotate(self, angle_x, angle_y, angle_z):

    angle_x = math.radians(angle_x)
    angle_y = math.radians(angle_y)
    angle_z = math.radians(angle_z)

    if angle_x != 0.0:
      new_y = self.y * math.cos(angle_x) - self.z * math.sin(angle_x)
      new_z = self.y * math.sin(angle_x) + self.z * math.cos(angle_x)
      self.y = new_y
      self.z = new_z
    if angle_y != 0.0:
      new_x = self.x * math.cos(angle_y) + self.z * math.sin(angle_y)
      new_z = - self.x * math.sin(angle_y) + self.z * math.cos(angle_y)
      self.x = new_x
      self.z = new_z
    if angle_z != 0.0:
      new_x = self.x * math.cos(angle_z) - self.y * math.sin(angle_z)
      new_y = self.x * math.sin(angle_z) + self.y * math.cos(angle_z)
      self.x = new_x
      self.y = new_y
    return self

  def rotate_vect(self, v, angle):

    angle = math.radians(angle)

    u = v.unit()

    c = math.cos(angle)
    s = math.sin(angle)
    ux2 = u.x * u.x
    uy2 = u.y * u.y
    uz2 = u.z * u.z

    x = (ux2 + (1 - ux2) * c) * self.x + (u.x * u.y * (1 - c) - u.z * s) * self.y + (u.x * u.z * (1 - c) + u.y * s) * self.z
    y = (u.x * u.y * (1 - c) + u.z * s) * self.x + (uy2 + (1 - uy2) * c) * self.y + (u.y * u.z * (1 - c) - u.x * s) * self.z
    z = (u.x * u.z * (1 - c) - u.y * s) * self.x + (u.y * u.z * (1 - c) + u.x * s) * self.y + (uz2 + (1 - uz2) * c) * self.z

    self.x = x
    self.y = y
    self.z = z
    return self

def get_dihedral_angle(coord_a, coord_b, coord_c, coord_d):
  '''Compute the dihedral angle (in deg) with the (4) coordinates'''

  vA = Vector3d(coord_a)
  vB = Vector3d(coord_b)
  vC = Vector3d(coord_c)
  vD = Vector3d(coord_d)
  vAB = vB - vA
  vBC = vC - vB
  vCD = vD - vC
  vABC = vAB.cross(vBC)
  vBCD = vBC.cross(vCD)
  phi = math.degrees(math.acos(vABC.dot(vBCD) / (vABC.length() * vBCD.length())))

  if vABC.cross(vBCD).dot(vBC) < 0.0:
    phi = 360.0 - phi

  return phi
