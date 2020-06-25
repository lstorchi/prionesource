import sys
import re

filename = ""

if len(sys.argv) < 2:
    print "usage: ", sys.argv[0] , " filename "
    exit(1)

xmin = [] 
ymin = [] 
zmin = [] 
xmax = [] 
ymax = [] 
zmax = []

for i in range(1,len(sys.argv)):
    filename = sys.argv[i]

    fp = open(filename, "r")
    
    xc = 0.0 
    yc = 0.0 
    zc = 0.0 
    xbm = 0.0 
    ybm = 0.0 
    zbm = 0.0
    xbx = 0.0 
    ybx = 0.0 
    zbx = 0.0
    
    for l in fp:
    
        if l.find("GRID_CENTER") != -1:
            line = l.replace("GRID_CENTER", "")
    
            p = re.compile(r'\s+')
            line = p.sub(' ', line)
            line = line.lstrip()
            line = line.rstrip()
                                   
            plist = line.split(",")
    
            xc = float(plist[0])
            yc = float(plist[1])
            zc = float(plist[2])
    
        if l.find("INNERBOX") != -1:
            line = l.replace("INNERBOX", "")
    
            p = re.compile(r'\s+')
            line = p.sub(' ', line)
            line = line.lstrip()
            line = line.rstrip()
                                   
            plist = line.split(",")
    
            xbm = float(plist[0])
            ybm = float(plist[1])
            zbm = float(plist[2])

        if l.find("OUTERBOX") != -1:
            line = l.replace("OUTERBOX", "")
    
            p = re.compile(r'\s+')
            line = p.sub(' ', line)
            line = line.lstrip()
            line = line.rstrip()
                                   
            plist = line.split(",")
    
            xbx = float(plist[0])
            ybx = float(plist[1])
            zbx = float(plist[2])

    xb = xbm + ((xbx - xbm)/2.0)
    yb = ybm + ((ybx - ybm)/2.0)
    zb = zbm + ((zbx - zbm)/2.0)

    print str(xc-(xb/2.0))+";"+str(xc+(xb/2.0))+";"+ \
        str(yc-(yb/2.0))+";"+str(yc+(yb/2.0))+";"+ \
        str(zc-(zb/2.0))+";"+str(zc+(zb/2.0))

    xmin.append(xc-(xb/2.0))
    xmax.append(xc+(xb/2.0))
    ymin.append(yc-(yb/2.0))
    ymax.append(yc+(yb/2.0))
    zmin.append(zc-(zb/2.0))
    zmax.append(zc+(zb/2.0))
    
    fp.close()


print str(min(xmin))+";"+str(max(xmax))+";"+ \
        str(min(ymin))+";"+str(max(ymax))+";"+ \
        str(min(zmin))+";"+str(max(zmax))


