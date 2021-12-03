import numpy as np

def parse_xyz_file(filename):
   AtomList = []
   CoordList = []
   fr = open(filename, 'r')
   for line in fr.readlines():
      l = re.search('^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$', line)
      if l!=None:
         AtomList.append(l.group(1))
         coord = float(l.group(2)),float(l.group(3)),float(l.group(4))
         CoordList.append(coord)
   fr.close()
   Coords = np.array(CoordList)
   return AtomList, Coords

def compute_distance(Coords, idx1, idx2):
   vec_12 = Coords[idx2-1] - Coords[idx1-1]
   return np.linalg.norm(vec_12)

def compute_angle(Coords, idx1, idx2, idx3):
   vec_21 = Coords[idx1-1] - Coords[idx2-1]
   vec_21 /= np.linalg.norm(vec_21)
   vec_23 = Coords[idx3-1] - Coords[idx2-1]
   vec_23 /= np.linalg.norm(vec_23)
   theta = np.arccos(np.dot(vec_21, vec_23)) * 180.0/np.pi
   return theta
