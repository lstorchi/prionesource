Samples:

python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-fit.txt -v 
python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-fit.txt -v --ddielectric --interpolate lista-fit_coulomb_mean.dx

python3 mol2todx.py -s 2.0 -f lista-fit.txt --interpolate lista-fit_coulomb_mean.dx
python3 mol2todx.py -s 2.0 -f lista-fit.txt -c --interpolate lista-fit_coulomb_mean.dx

python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-e200k-110NaCl.txt -v --interpolate lista-fit_coulomb_mean.dx
python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista-e200k-110NaCl.txt -v --ddielectric --ddielectric --interpolate lista-fit_coulomb_mean.dx

python3 mol2todx.py -s 2.0 -f lista-e200k-110NaCl.txt --interpolate lista-fit_coulomb_mean.dx
python3 mol2todx.py -s 2.0 -f lista-e200k-110NaCl.txt -c --interpolate lista-fit_coulomb_mean.dx

This is printing also the list to be used in the carbo'

for name in *.dx 
do 
  echo $name  
  python3 normalizedx.py -f $name
done

echo "lista-fit-dx-vs-lista-e200k-110NaCl-dx.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx.txt -f2 lista-e200k-110NaCl-dx.txt -v --show > lista-fit-dx-vs-lista-e200k-110NaCl-dx.out

echo "lista-fit-dx-vs-lista-fit-dx.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx.txt -f2 lista-fit-dx.txt -v --show > lista-fit-dx-vs-lista-fit-dx.out

echo "lista-e200k-110NaCl-dx-vs-lista-e200k-110NaCl-dx.png"
python3 dxcarbosimilarity.py -f1 lista-e200k-110NaCl-dx.txt -f2 lista-e200k-110NaCl-dx.txt -v --show > lista-e200k-110NaCl-dx-vs-lista-e200k-110NaCl-dx.out


echo "lista-fit-dx-ddielectric-vs-lista-e200k-110NaCl-dx-ddielectric.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx-ddielectric.txt -f2 lista-e200k-110NaCl-dx-ddielectric.txt -v --show > lista-fit-dx-ddielectric-vs-lista-e200k-110NaCl-dx-ddielectric.out

echo "lista-fit-dx-ddielectric-vs-ista-fit-dx-ddielectric.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx-ddielectric.txt -f2 lista-fit-dx-ddielectric.txt -v --show > lista-fit-dx-ddielectric-vs-lista-fit-dx-ddielectric.out

echo "lista-e200k-110NaCl-dx-ddielectric-vs-lista-e200k-110NaCl-dx-ddielectric.png"
python3 dxcarbosimilarity.py -f1 lista-e200k-110NaCl-dx-ddielectric.txt -f2 lista-e200k-110NaCl-dx-ddielectric.txt -v --show > lista-e200k-110NaCl-dx-ddielectric-vs-lista-e200k-110NaCl-dx-ddielectric.out


echo "lista-fit-dx-apbs-vs-lista-e200k-110NaCl-dx-apbs.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx-apbs.txt -f2 lista-e200k-110NaCl-dx-apbs.txt -v --show > lista-fit-dx-apbs-vs-lista-e200k-110NaCl-dx-apbs.out

echo "lista-fit-dx-apbs-vs-ista-fit-dx-apbs.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx-apbs.txt -f2 lista-fit-dx-apbs.txt -v --show > lista-fit-dx-apbs-vs-lista-fit-dx-apbs.out

echo "lista-e200k-110NaCl-dx-apbs-vs-lista-e200k-110NaCl-dx-apbs.png"
python3 dxcarbosimilarity.py -f1 lista-e200k-110NaCl-dx-apbs.txt -f2 lista-e200k-110NaCl-dx-apbs.txt -v --show > lista-e200k-110NaCl-dx-apbs-vs-lista-e200k-110NaCl-dx-apbs.out


echo "lista-fit-dx-apbs-flat-vs-lista-e200k-110NaCl-dx-apbs-flat.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx-apbs-flat.txt -f2 lista-e200k-110NaCl-dx-apbs-flat.txt -v --show > lista-fit-dx-apbs-flat-vs-lista-e200k-110NaCl-dx-apbs-flat.out

echo "lista-fit-dx-apbs-flat-vs-ista-fit-dx-apbs-flat.png"
python3 dxcarbosimilarity.py -f1 lista-fit-dx-apbs-flat.txt -f2 lista-fit-dx-apbs-flat.txt -v --show > lista-fit-dx-apbs-flat-vs-lista-fit-dx-apbs-flat.out

echo "lista-e200k-110NaCl-dx-apbs-flat-vs-lista-e200k-110NaCl-dx-apbs-flat.png"
python3 dxcarbosimilarity.py -f1 lista-e200k-110NaCl-dx-apbs-flat.txt -f2 lista-e200k-110NaCl-dx-apbs-flat.txt -v --show > lista-e200k-110NaCl-dx-apbs-flat-vs-lista-e200k-110NaCl-dx-apbs-flat.out

