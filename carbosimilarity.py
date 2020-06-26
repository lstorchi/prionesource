import sys
import argparse

sys.path.append("./common")
import carbo

if __name__ == "__main__":
    STEPVAL = 1.4
    DELTAVAL = 20.0
    coulombconst = 1.0
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f1","--file1", help="input the first mol2 file and weight filename  \n"+ \
      "   filenames are comma separated \"file1.mol,weight1.txt", required=True, default="", type=str)
    parser.add_argument("-f2","--file2", help="input the second mol2 file and weight filename  \n"+ \
      "   filenames are comma separated \"file2.mol,weight2.txt", required=True, default="", type=str)
    parser.add_argument("-s","--stepvalue", help="input step value defaul="+str(STEPVAL), \
      required=False, default=STEPVAL, type=float)
    parser.add_argument("-d","--deltaval", help="input delta value defaul="+str(DELTAVAL), \
      required=False, default=DELTAVAL, type=float)

    args = parser.parse_args()

    if len(args.file1.split(",")) != 2:
      print("Error in --file1 option two comma separated filenames expected")
      exit(1)

    if len(args.file2.split(",")) != 2:
      print("Error in --file2 option two comma separated filenames expected")
      exit(1)

    filename1 = args.file1.split(",")[0]
    weightsname1 = args.file1.split(",")[1]

    filename2 = args.file2.split(",")[0]
    weightsname2 = args.file2.split(",")[1] 
      
    carbo.carbo_similarity (filename1, weightsname1, filename2, weightsname2, \
        STEPVAL, DELTAVAL, coulombconst)

