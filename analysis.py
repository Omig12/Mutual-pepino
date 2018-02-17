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
# $ "python3 analysis.py -f <filename> [-t <threshold>] [-o <outputfilename>]"           #
# The GFA file utilized must be located on the same directory as the parser or full path #
# must be specified.                                                                     #
#                                                                                        #
##########################################################################################

################
# Load Modules # 
################

import argparse
import re
from sys import exit
# import os

# # Dominant node:
# def Dominant(kA, kB):
# 	'''Calculates differential expression inside a node on the gfa graph that are linked.'''
# 	return ("A by " if kA[0] > kB[0] else "B by " )+str(abs(kA[0]-kB[0])/(kA[0] + kB[0])), ("A by " if kA[1] > kB[1] else "B by " ) + str(abs(kA[1]-kB[1])/(kA[1] + kB[1]))
	

# # Dominant node:
# def Dominant(kA, kB):
# 	'''Calculates differential expression inside a node on the gfa graph that are linked.'''
# 	d = "A"
# 	if (kA[0] > kB[0]):
# 	   	d = "A: "
# 	elif (kA[0] == kB[0]):
# 		d = "None: "
# 	else:
# 		d = "B: "
# 	return	d + str(abs(kA[0]-kB[0])/(kA[0] + kB[0]))



# |A1 Delta B1| + |A2 Delta B2| / Over k-mer total
def DiffExp(kA, kB):
	'''Calculates differential expression between two node on the gfa graph that are linked.'''
	return ( (abs(kA[0]-kB[0]) + abs(kA[1]-kB[1]))/ (kA[0] + kA[1] + kB[0] + kB[1]) )

# |A1 Delta B1| + |A2 Delta B2| / Over k-mer total
# def DiffExp(kA, kB):
# 	'''Calculates differential expression between two node on the gfa graph that are linked.'''
# 	return ( abs(abs(kA[0]-kB[0]) - abs(kA[1]-kB[1]))/ (kA[0] + kA[1] + kB[0] + kB[1]) )

# |Delta A| + |Delta B| / Over k-mer total
# def DiffExp(kA, kB):
# 	'''Calculates differential expression between two node on the gfa graph that are linked.'''
# 	return ( (abs(kA[0]-kA[1]) + abs(kB[0]-kB[1])) / (kA[0] + kA[1] + kB[0] + kB[1]) )

# |A Delta B| / Over k-mer total 
# def DiffExp(kA, kB):
# 	'''Calculates differential expression between two node on the gfa graph that are linked.'''
# 	return (abs((kA[0] + kA[1]) - (kB[0] + kB[1]))/(kA[0] + kA[1] + kB[0] + kB[1]))

def main():
	'''Main Program function runs parser, differential expression and output'''
	# Regular expresion capture groups descriptions: g0 == all captures [g1:g6] 
	#   		 L	   0:1:(A:1,B:5)				+	  0:60:(A:1,B:501)				+	9M
	###############   g1  ####### g2 #### g3 ###########   g4  ####### g5 #### g6 ###
	# regex = r"^L\t(\d+:\d+):\(A:(\d+),B:(\d+)\)\t[+-]\t(\d+:\d+):\(A:(\d+),B:(\d+)\)"
	################   g1   ######### g2 ###### g3 #############   g4   ######### g5 ###### g6 ###############
	regex = r"^L\s+(\d+\:\d+)\:\(A\:(\d+)\,B\:(\d+)\)\s+[+-]\s+(\d+\:\d+)\:\(A\:(\d+)\,B\:(\d+)\)\s+\+\s+\d\w"
	prog = re.compile(regex, flags=re.M)
	with open(args.file,'r') as file:
		for line in file.readlines():
			if line[0] != "L":
				continue
			else:
				g = prog.match(line)
				# print(g.group(0))
				A = [int(g.group(2)), int(g.group(5))]
				B = [int(g.group(3)), int(g.group(6))]
				coef = DiffExp(A,B)
				print(line + "Differential Coefficient: " + str(coef)+ "\n\n")
				if  coef >= args.threshold:
					with open(args.output,"a+") as out:
						out.write(line + "Differential Coefficient: " + str(coef)+ "\n\n" )           
				else:
					continue
				# if args.dominant:
				# 	d1, d2 = Dominant(A,B)
				# 	print("Dominant in %s: %s \nDominant in %s: %s\n\n"%(g.group(1), d1, g.group(4), d2))
				# else:
				# 	continue
	# if os.path.isfile(args.output):
	#     # print({g.group(1):{"A":g.group(2), "B": g.group(3)},g.group(4):{"A":g.group(5), "B": g.group(6)} })
	#     print("Did not write to file")
	# else:
	# with open(args.output ,'a+') as out:
	#     out.write(line + "Differential Coefficient: " + str(coef)+ "\n\n" )           

parser = argparse.ArgumentParser(description="Differential Expression Coefficient Analysis Script")
parser.add_argument('-f', '--file', required=True, help="The GFA input file with Differential Expression info, ouput from dbg.py.")
parser.add_argument('-t','--threshold', default = 0.0, type= float, required=False, help="The threshold of differential expression, only show links with more differentially expressed than the float value t given, where 0 =< t <= 1. (default: 0.0 'Identical')")
# parser.add_argument('-d','--dominant', default = 0, type= bool, required=False, help="Calculate differential expression for individual kamers or nodes.")
parser.add_argument('-o','--output', default="diff_exp.txt",required=False,help="Output GFA file name. (default: 'diff_exp.txt')")
args = parser.parse_args()

if __name__ == "__main__":
	if args.file[-3:] != "gfa":
		if input("File does not appear to be a GFA, Do you wish to continue (y/n)? ") == "y":
			main()
		else:
			exit(0)
	else:
		main()
