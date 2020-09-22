Samples:

$ python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-fit.txt -v 
$ python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-fit.txt -v  --ddielectric

$ python3 mol2todx.py -s 2.0 -f lista-fit.txt 
$ python3 mol2todx.py -s 2.0 -f lista-fit.txt -c

$ python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-e200k-110NaCl.txt -v 
$ python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-e200k-110NaCl.txt -v  --ddielectric

$ python3 mol2todx.py -s 2.0 -f lista-e200k-110NaCl.txt 
$ python3 mol2todx.py -s 2.0 -f lista-e200k-110NaCl.txt -c

This is printing also the list to be used in the carbo'

$ for name in *.dx 
  do 
    echo $name  
    python3 normalizedx.py -f $name
   done

