> model1_vs_model1_DRY.txt
> model1_vs_model1_OH2.txt

for name in $(ls /pub3/redo/docking-pp_e200k-110NaCl/e200k_110NaCl_model1_A-B/output_files/*.pdb)
do 
  filename="${name%.*}"
  python3 centrodsranker.py -f $name -d 20.0 -m 100 -c 30 -s 1.0 -e -0.8 2> /dev/null 1> "$filename"_DRY.out
  python3 centrodsranker.py -f $name -d 20.0 -c 40 -s 1.0 -e -6.0 -p OH2 2> /dev/null 1> "$filename"_OH2.out
  VALUE=$(grep Pair "$filename"_DRY.out)
  echo $name $VALUE >> model1_vs_model1_DRY.txt
  VALUE=$(grep Pair "$filename"_OH2.out)
  echo $name $VALUE >> model1_vs_model1_OH2.txt
done 

