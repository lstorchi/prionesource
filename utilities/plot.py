import numpy
import argparse
import matplotlib.pyplot as plt

if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  parser.add_argument("-f","--file", help="input file", \
    required=True, default="", type=str)
  args = parser.parse_args()

  X1 = 0
  Y1 = 1
  S1 = 2
  X2 = 21
  Y2 = 22
  S2 = 23


  X3 = 7
  Y3 = 8
  S3 = 9
  X4 = 14
  Y4 = 15
  S4 = 16

  fp = open(args.file, "r")
  fp.readline()

  x1 = []
  y1 = []
  s1 = []

  x2 = []
  y2 = []
  s2 = []

  x3 = []
  y3 = []
  s3 = []

  x4 = []
  y4 = []
  s4 = []
  for l in fp:
      sl = l.split(",")
      if sl[X1] != "":
        x1.append(numpy.float64(sl[X1]))
      if sl[X2] != "":
        x2.append(numpy.float64(sl[X2]))

      if sl[Y1] != "":
        y1.append(numpy.float64(sl[Y1]))
      if sl[Y2] != "":
        y2.append(numpy.float64(sl[Y2]))
      
      if sl[S1] != "":
        s1.append(numpy.float64(sl[S1]))
      if sl[S2] != "":
        s2.append(numpy.float64(sl[S2]))

      if sl[X3] != "":
        x3.append(numpy.float64(sl[X3]))
      if sl[X4] != "":
        x4.append(numpy.float64(sl[X4]))

      if sl[Y3] != "":
        y3.append(numpy.float64(sl[Y3]))
      if sl[Y4] != "":
        y4.append(numpy.float64(sl[Y4]))
      
      if sl[S3] != "":
        s3.append(numpy.float64(sl[S3]))
      if sl[S4] != "":
        s4.append(numpy.float64(sl[S4]))
 
 
  fp.close()

  plt.errorbar(x1, y1, s1,  linestyle='None', \
      marker='^', label="APBS standard")
  plt.plot(x1, y1, linestyle='--')

  plt.errorbar(x2, y2, s2,  linestyle='None', \
      marker='^', label="Coulomb damped dielectric")
  plt.plot(x2, y2, linestyle='--')
  

  plt.errorbar(x3, y3, s3,  linestyle='None', \
      marker='^', label="APBS flat")
  plt.plot(x3, y3, linestyle='--')

  plt.errorbar(x4, y4, s4,  linestyle='None', \
      marker='^', label="Coulomb standard")
  plt.plot(x4, y4, linestyle='--')

  plt.legend(loc="lower left")
  plt.show()
