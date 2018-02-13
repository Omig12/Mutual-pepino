##########################################################################################
#                                                                                        #
# Copyright 2018 Megaprobe-Lab                                                           #
#                                                                                        #
# This is software created by the megaprobe lab under the GPL3 license.                  #
#                                                                                        #
# This program parses any ole GFA file you may have with read information for two        #
# organisms.                                                                             #
#                                                                                        #
# To run the program, utilize the command:                                               #
# $ "Python3 Second_Program.py GFA -f <filename>"                               #
# The GFA file utilized must be located on the same directory as the parser or full path #
# must be specified.                                                                     #
#                                                                                        #
##########################################################################################

################
# Load Modules # 
################

import re
import os
import sys
import argparse

def main():
    regex = r"^L\t(\d+:\d+):\(A:(\d+),B:(\d+)\)\t[+-]\t(\d+:\d+):\(A:(\d+),B:(\d+)\)"
    prog = re.compile(regex, flags = re.M)
    with open(args.file,'r') as file:
        for line in file.readlines():
            if line[0] != "L":
                continue
            else:
                # print(line)
                g = prog.match(line)
                A = [int(g.group(2)),int(g.group(5))]
                B = [int(g.group(3)),int(g.group(6))]
                coef = DiffExp(A,B)
                if  coef <= args.threshold:
                    with open(args.output ,'a+') as out:
                        out.write(line + "Differential Coefficient: " + str(coef)+ "\n\n" )           
                else:
                    continue
                    # if os.path.isfile(args.output):
                    #     print(line[:-1])
                    #     # print({g.group(1):{"A":g.group(2), "B": g.group(3)},g.group(4):{"A":g.group(5), "B": g.group(6)} })
                    #     print(coef)
                    #     print("Did not write to file")
                    # else:
                        # with open(args.output ,'a+') as out:
                        #     out.write(line + "Differential Coefficient: " + str(coef)+ "\n\n" )           

# A Delta B / Over k-mer total 
def DiffExp(kA, kB):
    return (abs((kA[0] + kA[1]) - (kB[0] + kB[1]))/(kA[0] + kA[1] + kB[0] + kB[1]))


parser = argparse.ArgumentParser(description="Differential Expression Coefficient")
parser.add_argument('-f', '--file', required=True, help="The GFA input file with Differential Expression info, ouput from dbg.py")
parser.add_argument('-t','--threshold', default = 1.0, type= float, required=False, help="The threshold of differential expression, only show links with similarity less than float value given. (default: 1.0 'Identical')")
parser.add_argument('-o','--output', default="diff_exp.txt",required=False,help="Output GFA file name. (default: 'diff_exp.txt')")

args = parser.parse_args()


if __name__ == "__main__":
    main()