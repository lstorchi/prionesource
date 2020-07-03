#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Read part of a matrix from a binary file making at most one copy of array
# in memory.
#
# The code array.frombuffer(file.read()) creates two copies in memory.
# The numpy.fromfile() routine can't read into an existing array.
#
# Shailesh Kumar Panday, PhD
def read_array(path, byte_offset, ijk_origin, ijk_size, ijk_step,
               full_size, type, byte_swap):

    if (tuple(ijk_origin) == (0,0,0) and
        tuple(ijk_size) == tuple(full_size)):
        m = read_full_array(path, byte_offset, full_size,
                            type, byte_swap)
        return m


# -----------------------------------------------------------------------------
# Read an array from a binary file making at most one copy of array in memory.
#
def read_full_array(path, byte_offset, size, type, byte_swap,
                    block_size = 2**20):

    a = allocate_array(size, type)
    file = open(path, 'rb')
    file.seek(byte_offset)
    file.readinto(a)
    file.close()

    if byte_swap:
        a.byteswap(True)

    return a
  
# -----------------------------------------------------------------------------
#
from numpy import float32
def allocate_array(size, value_type = float32, step = None, 
                   reverse_indices = True, zero_fill = False):

    if step is None:
        msize = size
    else:
        msize = [1+(sz-1)/st for sz,st in zip(size, step)]

    shape = list(msize)
    if reverse_indices:
        shape.reverse()

    if zero_fill:
        from numpy import zeros as alloc
    else:
        from numpy import empty as alloc

    try:
        m = alloc(shape, value_type)
    except ValueError:
        # numpy 1.0.3 sometimes gives ValueError, sometimes MemoryError
        report_memory_error(msize, value_type)
    except MemoryError:
        report_memory_error(msize, value_type)
  
    return m
  
# -----------------------------------------------------------------------------
#
def report_memory_error(size, value_type):
    from numpy import dtype, product, float
    vtype = dtype(value_type)
    tsize = vtype.itemsize
    bytes = product(size, dtype=float)*float(tsize)
    mbytes = bytes / 2**20
    sz = ','.join(['%d' % s for s in size])
    e = ('Could not allocate %.0f Mbyte array of size %s and type %s.\n'
         % (mbytes, sz, vtype.name))
    
    raise (Exception, e)
# -----------------------------------------------------------------------------
# Read DelPhi or GRASP unformatted phi file.  This was derived from the
# Chimera DelphiViewer extension file reading code.
#

# -----------------------------------------------------------------------------
#
class DelPhiGridMap:

  def __init__(self, path):

    self.path = path

    file = open(path, 'rb')

    file.seek(0,2)                              # goto end of file
    self.file_size = file.tell()
    file.seek(0,0)                              # goto beginning of file

    if self.file_size == 0:
      raise (SyntaxError, 'Empty file')

    swap = self.need_to_swap_bytes(file)
    uplbl = self.read_record(file, swap)
    morelabels = self.read_record(file, swap)
    self.data_offset = file.tell()
    self.skip_record(file, swap)
    botlbl = self.read_record(file, swap)
    params = self.read_record(file, swap)
    file.close()

    from numpy import float32, int32, float64
    if len(params) == 16:	# GRASP Phi file
      size = 65
      self.value_type = float32
      sc = params
    elif len(params) == 20:	# DelPhi Phi file
      size = string_values(params[16:20], int32, swap)[0]
      self.value_type = float32
      sc = params[:16]
    elif len(params) == 36:	# 2008 Mac DelPhi Phi file
      size = string_values(params[32:36], int32, swap)[0]
      self.value_type = float64
      sc = params[:32]
    else:
      raise (SyntaxError, ('Parameter record size %d must be 16 or 20 or 36'
                          % len(params)))
    pval = string_values(sc, self.value_type, swap)
    self.scale = pval[0]
    self.xyz_center = pval[1:4]

    step = 1.0/self.scale
    half_size = step * ((size - 1) / 2)
    xyz_origin = [c - half_size for c in self.xyz_center]
    #print("while  parsing: xyz_origin", xyz_origin)
    self.swap = swap
    self.size = (size, size, size)
    self.xyz_step = (step, step, step)
    self.xyz_origin = xyz_origin
    
  # ---------------------------------------------------------------------------
  #
  def need_to_swap_bytes(self, file):

    from numpy import frombuffer, int32
    v = frombuffer(file.read(4), int32)[0]
    file.seek(0,0)
    return (v < 0 or v >= 65536)
    
  # ---------------------------------------------------------------------------
  #
  def read_record(self, file, swap, skip = False):

    from numpy import int32
    size = string_values(file.read(4), int32, swap)[0]
    if size < 0:
      raise (SyntaxError, 'Negative record size %d' % size)
    if size > self.file_size:
      raise (SyntaxError, ('Record size %d > file size %d'
                          % (size, self.file_size)))

    if skip:
      file.seek(size, 1)
      string = ''
    else:
      string = file.read(size)

    from numpy import int32
    esize = string_values(file.read(4), int32, swap)[0]
    if esize != size:
      raise (SyntaxError, ('Record size at end of record %d' % esize + 
                          ' != size at head of record %d' % size))
      
    return string
    
  # ---------------------------------------------------------------------------
  #
  def skip_record(self, file, swap):

    self.read_record(file, swap, skip = True)
    
  # ---------------------------------------------------------------------------
  #
  def matrix(self, ijk_origin, ijk_size, ijk_step):
    data = read_array(self.path, self.data_offset + 4,
                      ijk_origin, ijk_size, ijk_step, self.size,
                      self.value_type, self.swap)
    return data

  # ---------------------------------------------------------------------------
  #
  def get_data(self):
    matrix = self.matrix((0, 0, 0), self.size, self.xyz_step)
    return matrix

  # ---------------------------------------------------------------------------
  #
  def delphimap_to_text_file(self, filename):
    """
    This method writes out the associated .phi file to a text file where
    first 4 lines are header 
    5th line onwards there are nX * xY * nZ data lines.
    data format: index x y z phi
    """
    matrix = self.get_data()
    fout = open(filename, 'w')
    fout.write("scale = {:8.3f} per angstrom\n".format(self.scale))
    fout.write("number of grids  (nX, nY, nZ) = {:8d} {:8d} {:8d}\n".format(*self.size))
    fout.write("origin of grid   (X0, Y0, Z0) = {:8.3f} {:8.3f} {:8.3f}\n".format(*self.xyz_origin))
    fout.write("coordinate steps (dX, dY, dZ) = {:8.3f} {:8.3f} {:8.3f}\n".format(*self.xyz_step))
    idx = 0
    for i in range(matrix.shape[0]):
      x = self.xyz_origin[0] + i * self.xyz_step[0]
      for j in range(matrix.shape[1]):
        y = self.xyz_origin[1] + j * self.xyz_step[1]
        for k in range(matrix.shape[2]):
          z = self.xyz_origin[2] + k * self.xyz_step[2]
          fout.write("%9d %8.3f %8.3f %8.3f %10.5f\n"%(idx, x, y, z, matrix[i,j,k]))
          idx += 1
      fout.flush()
    fout.close()
    
    def get_ijk_coord_and_phi(self, matrix, i, j, k):
      x = self.xyz_origin[0] + i * self.xyz_step[0]
      y = self.xyz_origin[1] + j * self.xyz_step[1]
      y = self.xyz_origin[2] + k * self.xyz_step[1]
      phi = matrix[i, j, k]
      return (x , y, z, phi)

  def delphimap_to_kont_file(self, filename):
    """
    This method writes out the associated .phi file to a ASCII kont file 
    """
    matrix = self.get_data()
    fout = open(filename, 'w')

    idx = 1
    for k in range(matrix.shape[2]):
      z = self.xyz_origin[2] + k * self.xyz_step[2]
      for i in range(matrix.shape[0]):
        x = self.xyz_origin[0] + i * self.xyz_step[0]
        for j in range(matrix.shape[1]):
          y = self.xyz_origin[1] + j * self.xyz_step[1]

          fout.write("%9d %8.3f %8.3f %8.3f\n"%(idx, x, y, z))
          idx += 1
    
    fout.write("Fake \n")
    
    for k in range(matrix.shape[2]):
      for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
          fout.write("%8.3f\n"%(matrix[i, j, k]))

    fout.close()
    
# -----------------------------------------------------------------------------
#
def string_values(string, type, swap):
  from numpy import frombuffer
  values = frombuffer(string, type)
  if swap:
    values = values.byteswap()
  return values


'''
def DelPhiMap_scipy(filename, txtoutfile):
  from scipy.io import FortranFile

  f = FortranFile(filename, 'r' )
    
  uplbl = f.read_record('a20')
  nxtolbl = f.read_record('a70')
  epmap = f.read_reals(dtype='float32').reshape((65,65,65))
  botlbl = f.read_record('a16')
  scalemin = f.read_reals(dtype='float32')
    
  scale = scalemin[0]
  oldmid = scalemin[1:4]
  fout = open(txtoutfile, 'w')
  idx = 1
  for iz in range(65):
    IZ = iz + 1
    for ix in range(65):
      IX = ix + 1
      for iy in range(65):
        IY = iy + 1
    
        x = (IX - 33)/scale + oldmid[0]
        y = (IY - 33)/scale + oldmid[1]
        z = (IZ - 33)/scale + oldmid[2]
    
        fout.write("%6d %8.3f %8.3f %8.3f %10.5f\n"%(idx, x, y, z, epmap[ix, iy, iz]))
        idx += 1
  fout.close()
'''    
  
if __name__ == "__main__":
  """
  creating a DelPhiGridMap object from a given .phi potential file
  >>> delphi_data = DelPhiGridMap("/path/to/test.phi")
  
  obtaining the gridmatrix data
  >>> matrix = delphi_data.get_data()
  
  Writing phimap to a text format
  >>> delphi_data.delphimap_to_text_file("/path/to/test-1.txt")
  """
  delphi_data = DelPhiGridMap("wt_90-231_aligned_2.phi")
        
  matrix = delphi_data.get_data()
    
  delphi_data.delphimap_to_kont_file("delphi.kont")
    
  #DelPhiMap_scipy("/media/shailesh/DATA/Delphi_ALL/phimap_read_example/test.phi", "/media/shailesh/DATA/Delphi_ALL/phimap_read_example/test-2.txt")
