#!env python

import argparse
import warnings

def interface_standard():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-if', '--input-file-list', dest='ifname', help='tree file list',required=True, nargs='+')
    parser.add_argument('-of', '--output-file-name', dest='ofname', help='output file name',required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = interface_standard()
    writer = open(args.ofname,'w')
    for fname in args.ifname :
        with open(fname) as handler :
            writer.write("".join([ line.strip() for line in handler ]) + "\n")
    writer.close()
            
    
