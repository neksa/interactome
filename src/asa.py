####
From https://github.com/boscoh/pdbremix/blob/master/pdbremix/asa.py
###

try:
  import numpy
  from v3numpy import *
except:
  from v3list import *


import math
import random

import numpy as np

def is_similar_mag(a, b, small=1E-5):
  return abs(abs(a)-abs(b)) <= small


def vector(*args):
  n_arg = len(args)
  if n_arg == 0:
    return np.zeros(3, dtype=np.float)
  if n_arg == 1:
    data = args[0]
    if len(data) == 3:
      return np.array(data, copy=True)
    raise TypeError('vector() with 1 argument must have 3 elements')
  if n_arg == 3:
    return np.array(args, dtype=np.float, copy=True)
  else:
    raise TypeError('vector() takes 0, 1 or 3 arguments')


def set_vector(*args):
  "Changes values of a vector in place"
  vector = args[0]
  if len(args) == 2:
    vector[:] = args[1]
  elif len(args) == 4:
    vector[:] = args[1:4]


def crd(vector):
  "Returns values of vector as a sequence of floats"
  return vector


def is_similar_matrix(m1, m2, small=1E-5):
  iter1 = np.ndenumerate(m1)
  iter2 = np.ndenumerate(m2)
  for (i1, val1), (i2, val2) in zip(iter1, iter2):
    if not is_similar_mag(val1, val2, small):
      return False
  return True


is_similar_vector = is_similar_matrix


def mag(vector):
  return np.sqrt(np.dot(vector, vector))


def scale(vector, s):
  return  s*vector


def norm(vector):
  return vector/mag(vector)


radians = np.radians

degrees = np.degrees

cross = np.cross

dot = np.dot


def identity():
  m = np.zeros((4, 3))
  m[:3,:3] = np.eye(3)
  return m


def transform(matrix, vector):
  return np.dot(matrix[:3,:3], vector) + matrix[3,:]  


def left_inverse(matrix):
  inverse = identity()
  r = matrix[:3,:3].transpose()
  inverse[:3,:3] = r
  inverse[3,:] = -np.dot(r, matrix[3,:])
  return inverse


# from http://stackoverflow.com/a/6802723
# uses the right hand screw rule for theta
def rotation(axis, theta):
  m = identity()
  a = np.cos(theta/2)
  b, c, d = norm(axis) * np.sin(theta/2)
  m[0] = [a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)    ]
  m[1] = [2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b)    ]
  m[2] = [2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-b*b-c*c]
  return m


def translation(displacement):
  m = identity()
  m[3,:] = displacement
  return m


def combine(m1, m2):
  m3 = identity()
  m3[:3,:3] = np.dot(m1[:3,:3], m2[:3,:3])
  m3[3,:] = np.dot(m1[:3,:3], m2[3,:]) + m1[3,:]
  return m3




# common to v3list and v3numpy


def parallel(v, axis):
  l = mag(axis)
  if is_similar_mag(l, 0):
    return v
  else:
    return scale(axis, dot(v, axis)/l/l) 


def perpendicular(v, axis):
  return v - parallel(v, axis)


def normalize_angle(angle):
  while abs(angle) > math.pi:
    if angle > math.pi:
      angle -= math.pi*2
    if angle < -math.pi:
      angle += 2*math.pi
  if is_similar_mag(abs(angle + math.pi), 0):
    angle = math.pi
  return angle


def vec_angle(a, b):
  a_len = mag(a)
  b_len = mag(b)
  if is_similar_mag(a_len, 0) or is_similar_mag(b_len, 0):
    return 0.0
  c = dot(a, b) / a_len / b_len
  if c >=  1.0:
    return 0.0
  elif c <= -1.0:
    return math.pi
  else:
    return math.acos(c)  


def vec_dihedral(a, axis, c):
  ap = perpendicular(a, axis)
  cp = perpendicular(c, axis)
  angle = vec_angle(ap, cp)
  if dot(cross(ap, cp), axis) > 0:
    angle = -angle
  return angle


def dihedral(p1, p2, p3, p4):
  return vec_dihedral(p1-p2, p2-p3, p4-p3)


def distance(p1, p2):
  return mag(p1 - p2)


def rotation_at_center(axis, theta, center):
  t = translation(-center)
  r = rotation(axis, theta)
  t_inv = translation(center)
  return combine(t_inv, combine(r, t))


def get_center(crds):
  center = vector()
  for crd in crds:
    center += crd
  return scale(center, 1.0/float(len(crds)))


def get_width(crds):
  center = get_center(crds)
  max_diff = 0
  for crd in crds:
    diff = v3.distance(crd, center)
    if diff > max_diff:
      max_diff = diff
  return 2*max_diff


def random_mag(): 
  return random.uniform(0, 90)


def random_real(): 
  return random.uniform(-90, 90)


def random_vector():
  return vector(random_real(), random_real(), random_real())


def random_rotation():
  return rotation(random_vector(), radians(random_real()))


def random_matrix():
  return combine(random_rotation(), translation(random_vector()))


############


class Atom:
  def __init__(
      self, pos=None,
      atom_type="", res_num=None):
    self.is_hetatm = False
    self.pos = v3.vector() if pos is None else pos
    self.vel = v3.vector()
    self.mass = 0.0
    self.type = ""
    self.element = ""
    self.chain_id = " "
    self.res_type = ""
    self.res_num = ""
    self.res_insert = ""
    self.bfactor = 0.0
    self.occupancy = 0.0
    self.num = 0
    self.alt_conform = " "

  def copy(self):
    return copy.deepcopy(self)

  def type_str(self):
    atom_type = self.type
    if len(atom_type) == 1:
      atom_type = " %s  " % atom_type
    elif len(atom_type) == 2:
      atom_type = " %s " % atom_type
    elif len(atom_type) == 3:
      if atom_type[0].isdigit():
        atom_type = "%s " % atom_type
      else:
        atom_type = " %s" % atom_type
    return atom_type    

  def pdb_str(self):
    if self.is_hetatm:
      field = "HETATM"
    else:
      field = "ATOM  "
    x, y, z = self.pos
    s = "%6s%5s %4s %-4s%1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f" \
            % (field, 
               str(self.num)[-5:], 
               self.type_str(),
               self.res_type, 
               self.chain_id,
               str(self.res_num)[-4:], 
               self.res_insert,
               x, y, z,
               self.occupancy, 
               self.bfactor)
    return s
               
  def __str__(self):
    x, y, z = self.pos
    return "%s%s-%s (% .1f % .1f % .1f)" \
            %  (self.res_type, self.res_num, 
                self.type, x, y, z)

  def transform(self, matrix):
    new_pos = v3.transform(matrix, self.pos)
    v3.set_vector(self.pos, new_pos)


def AtomFromPdbLine(line):
  """Returns an Atom object from an atom line in a pdb file."""
  atom = Atom()
  if line.startswith('HETATM'):
    atom.is_hetatm = True
  else:
    atom.is_hetatm = False
  atom.num = int(line[6:11])
  atom.type = line[12:16].strip(" ")
  element = ''
  for c in line[12:15]:
    if not c.isdigit() and c != " ":
      element += c
  if element[:2] in two_char_elements:
    atom.element = element[:2]
  else:
    atom.element = element[0]
  atom.alt_conform = line[16]
  atom.res_type = line[17:21].strip()
  atom.chain_id = line[21]
  atom.res_num = int(line[22:26])
  atom.res_insert = line[26]
  if atom.res_insert == " ":
    atom.res_insert = ""
  x = float(line[30:38])
  y = float(line[38:46])
  z = float(line[46:54])
  v3.set_vector(atom.pos, x, y, z)
  try:
    atom.occupancy = float(line[54:60])
  except:
    atom.occupancy = 100.0
  try:
    atom.bfactor = float(line[60:66])
  except:
    atom.bfactor = 0.0
  return atom
  
  
def cmp_atom(a1, a2):
  if a1.num < a2.num:
    return -1
  else:
    return 0


radii = { 
 'H': 1.20,
 'N': 1.55,
 'NA': 2.27,
 'CU': 1.40,
 'CL': 1.75,
 'C': 1.70,
 'O': 1.52,
 'I': 1.98,
 'P': 1.80,
 'B': 1.85,
 'BR': 1.85,
 'S': 1.80,
 'SE': 1.90,
 'F': 1.47,
 'FE': 1.80,
 'K':  2.75,
 'MN': 1.73,
 'MG': 1.73,
 'ZN': 1.39,
 'HG': 1.8,
 'XE': 1.8,
 'AU': 1.8,
 'LI': 1.8,
 '.': 1.8
}
two_char_elements = [e for e in radii.keys() if len(e) == 2]


def add_radii(atoms):
  for atom in atoms:
    if atom.element in radii:
      atom.radius = radii[atom.element]
    else:
      atom.radius = radii['.']


def get_center(atoms):
  center = v3.vector()
  for atom in atoms:
    center += atom.pos
  return v3.scale(center, 1.0/float(len(atoms)))


def get_width(atoms, center):
  max_diff = 0
  for atom in atoms:
    diff = v3.distance(atom.pos, center)
    if diff > max_diff:
      max_diff = diff
  return 2*max_diff


class AtomList:

  def __init__(self, pdb=""):
    self.id = ''
    self._atoms = []
    if pdb:
      self.read_pdb(pdb)

  def copy(self):
    return copy.deepcopy(self)

  def n_atom(self):
    return len(self._atoms)

  def atoms(self):
    return self._atoms

  def atom(self, i):
    return _atoms[i]
    
  def clear(self):
    for atom in self._atoms:
      del atom
    del self._atoms[:]

  def transform(self, matrix):
    for atom in self._atoms:
      atom.transform(matrix)

  def insert_atom(self, atom):
    self._atoms.append(atom)
    
  def erase_atom(self, atom_type):
    for atom in self._atoms:
      if atom.type == atom_type:
        self._atoms.remove(atom)
        del atom
        return

  def read_pdb(self, fname):
    self.clear()
    for line in open(fname, 'r').readlines():
      if line.startswith(("ATOM", "HETATM")):
        atom = AtomFromPdbLine(line);
        if len(self._atoms) == 1:
          self.id = atom.chain_id
        self.insert_atom(atom)
      if line.startswith("ENDMDL"):
        return

  def write_pdb(self, pdb):
    with open(pdb, 'w') as f:
      for atom in sorted(self._atoms, cmp=cmp_atom):
        f.write(atom.pdb_str() + '\n')

  def set_id(self, new_id):
    self.id = new_id
    for a in self.atoms():
      a.chain_id = new_id



class Residue:

  def __init__(self, in_type, in_chain_id, in_num, in_insert=''):
    self.type = in_type
    self.chain_id = in_chain_id
    self.num = in_num
    self.insert = in_insert
    self._atom_dict = {}
 
  def name(self):
    tag = ""
    if self.chain_id != " " and self.chain_id != "":
      tag += self.chain_id + ":"
    tag += str(self.num)
    if self.insert:
      tag += self.insert
    return tag  

  def __str__(self):
    atom_name_list = [a.type for a in self.atoms()]
    atom_name = " ".join(atom_name_list)
    return "%s-%s { %s }" % (self.type, self.num, atom_name)

  def copy(self):
    return copy.deepcopy(self)
  
  def n_atom(self):
    return len(self._atom_dict)
    
  def atom(self, atom_type):
    return self._atom_dict[atom_type]
    
  def has_atom(self, atom_type):
    return atom_type in self._atom_dict.keys()
    
  def atoms(self):
    return self._atom_dict.values()
  
  def atom_name(self, atom_type):
    return self.type + self.num + ":" + atom_type

  def insert_atom(self, atom):
    self._atom_dict[atom.type] = atom
    atom.chain_id = self.chain_id
    atom.res_num = self.num
    atom.res_type = self.type
  
  def erase_atom(self, atom_type):
    del self._atom_dict[atom_type]
    
  def set_num(self, i, insert=""):
    self.num = i
    self.insert = insert
    for atom in self.atoms():
      atom.res_num = self.num
      atom.res_insert = insert
     
  def inc_num(self):
    self.set_num(self.num+1, self.insert)

  def dec_num(self):
    self.set_num(self.num-1, self.insert)
    
  def dec_insert(self):
    l = self.insert;
    if l == "A" or l == "a":
      self.insert = ''
    else:
      i = string.ascii_letters.find(l)
      self.insert = string.ascii_letters[i-1]

  def transform(self, matrix):
     for atom in self.atoms():
       atom.transform(matrix)

  def set_chain_id(self, chain_id):
    self.chain_id = chain_id
    for a in self.atoms():
      a.chain_id = chain_id

  def set_type(self, res_type):
    self.type = res_type
    for a in self.atoms():
      a.res_type = res_type


class Polymer(AtomList):

  def __init__(self, fname=""):
    AtomList.__init__(self)
    self._residues = []
    if fname:
      self.read_pdb(fname)

  def copy(self):
    return copy.deepcopy(self)

  def residue(self, i):
    return self._residues[i]
    
  def residues(self):
    return self._residues

  def insert_atom(self, i, atom):
    self._atoms.append(atom)
    self.residue(i).insert_atom(atom)
    
  def erase_atom(self, i, atom_type):
    atom = self.residue(i).atom(atom_type)
    self._atoms.remove(atom)
    self.residue(i).erase_atom(atom_type)
    del atom
    
  def clear(self):
    del self._residues[:]
    AtomList.clear(self)
    
  def n_residue(self):
    return len(self._residues)
    
  def insert_residue(self, i, res):
    is_insertion = False
    if i < self.n_residue()-1:
      save_res_num = self.residue(i).num
      if self.residue(i+1).num == save_res_num:
        is_insertion = True

    if self.n_residue() == 0:
      res.set_num(res.num, res.insert)
    elif i < self.n_residue():
      res.set_num(self.residue(i).num, self.residue(i).insert)
    else:
      res.set_num(self.residue(i-1).num, "")
      res.inc_num()

    self._residues.insert(i, res)
    for atom in res.atoms():
      self.insert_atom(i, atom)

    for j in range(i+1, self.n_residue()):
      self.residue(j).inc_num()

    if is_insertion:
      while self.residue(i+1).insert:
        for j in range(i+1, self.n_residue()):
          if self.residue(j).res_num == save_res_num:
            self.residue(k).dec_insert()
    
  def append_residue(self, res):
    self._residues.append(res)
    for atom in res.atoms():
      self.insert_atom(self.n_residue()-1, atom)

  def erase_residue(self, i):  
    save_res_num = self.residue(i).num

    for atom in self.residue(i).atoms():
      self._atoms.remove(atom)
      del atom
    self._residues.pop(i)  
    
    if i < self.n_residue():
      if self.residue(i).num == save_res_num:
        # erasing residue in an insertion
        for j in range(i, self.n_residue()):
          if self.residue(j).num == erase_res_num_int:
            self.residue(j).dec_insert()
      else:
        for j in range(i, self.n_residue()):
          self.residue(j).dec_num()
    
  def get_i_residue(self, tag):
    # clean up tag
    tag = tag.strip()
    if tag[0] == ":":
      tag = tag[1:]
    if not tag[0].isdigit() and tag[1].isdigit():
      tag = tag[0] + ":" + tag[1:]
    for i, residue in enumerate(self.residues()):
      if tag.lower() == residue.tag().lower():
        return i
    raise Exception("Can't find residue " + tag)
  
  def extract_polymer(self, i, j):
    extract = Polymer()
    for res in self.residues()[i:j]:
      extract.append_residue(res.copy())
    return extract
 
  def chain_ids(self):
    chain_id = [r.chain_id for r in self.residues()]
    return list(set(chain_id))

  def extract_chain(self, chain_id):
    extract = Polymer()
    for res in self.residues():
      if res.chain_id == chain_id:
        extract.append_residue(res.copy())
    return extract
 
  def insert_polymer(self, i, insert):
    for res in reversed(insert.residues()):
      self.insert_residue(i, res.copy())
    
  def load_residue_bfactors(self, res_bfactors):
    for i, r in enumerate(self.residues()):
      for atom in r.atoms():
        if i >= len(res_bfactors):
          return
        else:
          atom.bfactor = res_bfactors[i]

  def __str__(self):
    res_name_list = [str(res) for res in self._residues]
    return "\n".join(res_name_list)
 
  def read_pdb(self, fname):
    self.clear()
    res_num = -1
    res_insert = " "
    for line in open(fname, 'r').readlines():
      if line.startswith("ATOM") or line.startswith("HETATM"):
        atom = AtomFromPdbLine(line);
        if (res_num != atom.res_num) or (res_insert != atom.res_insert):
          residue = Residue(atom.res_type, atom.chain_id,
                            atom.res_num, atom.res_insert)
          self.append_residue(residue)
          res_num = atom.res_num
          res_insert = atom.res_insert
        self.insert_atom(-1, atom)
      if line.startswith("ENDMDL"):
        return

  def write_pdb(self, pdb):
    f = open(pdb, 'w')
    n_atom = 0
    for res in self.residues():
      res_atoms = res.atoms()
      res_atoms.sort(cmp_atom)
      for atom in res_atoms:
        f.write(atom.pdb_str() + '\n')
    f.close()



############


import math
import v3
import pdbatoms


def generate_sphere_points(n):
  """
  Returns list of 3d coordinates of points on a sphere using the
  Golden Section Spiral algorithm.
  """
  points = []
  inc = math.pi * (3 - math.sqrt(5))
  offset = 2 / float(n)
  for k in range(int(n)):
    y = k * offset - 1 + (offset / 2)
    r = math.sqrt(1 - y*y)
    phi = k * inc
    points.append(v3.vector(math.cos(phi)*r, y, math.sin(phi)*r))
  return points


def find_neighbor_indices(atoms, probe, k):
  """
  Returns list of indices of atoms within probe distance to atom k. 
  """
  neighbor_indices = []
  atom_k = atoms[k]
  radius = atom_k.radius + probe + probe
  indices = range(k)
  indices.extend(range(k+1, len(atoms)))
  for i in indices:
    atom_i = atoms[i]
    dist = v3.distance(atom_k.pos, atom_i.pos)
    if dist < radius + atom_i.radius:
      neighbor_indices.append(i)
  return neighbor_indices


def calculate_asa(atoms, probe, n_sphere_point=960):
  """
  Returns list of accessible surface areas of the atoms,
  using the probe and atom radius to define the surface.
  """
  sphere_points = generate_sphere_points(n_sphere_point)
 
  const = 4.0 * math.pi / len(sphere_points)
  areas = []
  for i, atom_i in enumerate(atoms):
    
    neighbor_indices = find_neighbor_indices(atoms, probe, i)
    n_neighbor = len(neighbor_indices)
    j_closest_neighbor = 0
    radius = probe + atom_i.radius
    
    n_accessible_point = 0
    for point in sphere_points:
      is_accessible = True
      test_point = v3.scale(point, radius) + atom_i.pos
      cycled_indices = range(j_closest_neighbor, n_neighbor)
      cycled_indices.extend(range(j_closest_neighbor))
      
      for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff = v3.distance(atom_j.pos, test_point)
        if diff*diff < r*r:
          j_closest_neighbor = j
          is_accessible = False
          break
      if is_accessible:
        n_accessible_point += 1
    
    area = const*n_accessible_point*radius*radius 
    areas.append(area)
  return areas

##################################
  
#!/usr/bin/env python

__doc__ = """
Calculates the Accessible Surface Area of atoms in a PDB file. 
Algorithm adapted from Rose lab's chasa.py implementation of
Shrake & Rupley (19743) JMB:351. (c) 2007 Bosco Ho.

Usage: asa.py [-n <n_sphere>] <in_pdb> [<out_pdb>]

  <out_pdb>   PDB file in which the atomic ASA values are written 
              to the b-factor column.
  -n <n_sphere>  number of points used in generating the spherical
              dot-density for the calculation (default=960). The 
              more points, the more accurate (but slower) the 
              calculation. 
"""

from pdbremix.asa import calculate_asa
from pdbremix import pdbatoms
from pdbremix.docopt import docopt

if __name__ == "__main__":
  arg = docopt(__doc__)

  n_sphere = 960
  if arg['-n']:
    n_sphere = int(arg['-n'])
  print "Points on sphere: ", n_sphere

  mol = pdbatoms.AtomList(arg['<in_pdb>'])
  atoms = mol.atoms()
  pdbatoms.add_radii(atoms)

  asas = calculate_asa(atoms, 1.4, n_sphere)
  print "%.1f angstrom squared." % sum(asas)

  if arg['<out_pdb>']:
    for asa, atom in zip(asas, atoms):
      atom.bfactor = asa
    mol.write_pdb(arg['<out_pdb>'])




