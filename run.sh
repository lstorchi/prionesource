python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista_wtype.txt -v 2> lista_wtype_dx_coulomb.txt
python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista_wtype.txt -v --ddielectric --interpolate lista_wtype_coulomb_mean.dx 2> lista_wtype_dx_coulomb_ddielectric.txt

python3 mol2todx.py -s 2.0 -f lista_wtype.txt --interpolate lista_wtype_coulomb_mean.dx 2> lista_wtype_dx_apbs.txt
python3 mol2todx.py -s 2.0 -f lista_wtype.txt -c --interpolate lista_wtype_coulomb_mean.dx 2> lista_wtype_dx_apbs_flat.txt

python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista_e200k.txt -v --interpolate lista_wtype_coulomb_mean.dx 2> lista_e200k_dx_coulomb.txt
python3 mol2tocoulombdx.py -s 2.0 -d 30.0 -f lista_e200k.txt -v --ddielectric --ddielectric --interpolate lista_wtype_coulomb_mean.dx 2> lista_e200k_dx_coulomb_ddielectric.txt

python3 mol2todx.py -s 2.0 -f lista_e200k.txt --interpolate lista_wtype_coulomb_mean.dx 2> lista_e200k_dx_apbs.txt
python3 mol2todx.py -s 2.0 -f lista_e200k.txt -c --interpolate lista_wtype_coulomb_mean.dx 2> lista_e200k_dx_apbs_flat.txt

for name in *.dx 
do 
  echo $name  
  python3 normalizedx.py -f $name
done


for method in coulomb coulomb_ddielectric apbs  apbs_flat
do
  for mol1 in wtype e200k
  do 
    for mol2 in wtype e200k
    do
      echo $mol1 " vs " $mol2  " using " $method
    
      echo "$mol1"_dx_"$method"_vs_"$mol2"_dx_"$method".png
      python3 dxcarbosimilarity.py -f1 lista_"$mol1"_dx_"$method".txt -f2 lista_"$mol2"_dx_"$method".txt -v > "$mol1"_dx_"$method"_vs_"$mol2"_dx_"$method".out
      mv finalcarbo.png "$mol1"_dx_"$method"_vs_"$mol2"_dx_"$method".png

    done
  done
done

