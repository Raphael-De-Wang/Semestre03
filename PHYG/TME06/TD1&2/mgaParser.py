import sys
import os
import getopt
import commands

def main(argv=None):
    try:
        if argv is None:
            argv = sys.argv;
            if len(argv) <= 1:
                print "No Parameters provided";
        try:
            help_message = 'MetaGeneAnnotator Parser Script  \n usage: python pp.py -i [MetaGeneAnnotator output] -o [fasta output] ';
            opts, args = getopt.getopt(argv[1:], "h:i:o:g:", ["help","input=","output=","genome="]);
        except getopt.error:
            print 'MetaGeneAnnotator Parser Script  \n usage: python pp.py -i [MetaGeneAnnotator output] -o [fasta output] ';
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message);
            if option in ("-i", "--input"):
                inputFile  = value;
            if option in ("-o", "--output"):
                outputFile  = value
            if option in ("-g", "--genome"):
                genomeFile = value
        try:
            inputFile #checks if input exists otherwise raises exception
        except:
            print "Input file not specified";
        try:
            outputFile #checks if input exists otherwise raises exception
        except:
            print "Output file not specified";
        try:
            genomeFile #checks if input exists otherwise raises exception
        except:
            print "Genome file not specified";
    except:
        print "\tfor help use --help or -h";
        return 2;
    
	apply_everything(genomeFile,inputFile,outputFile)
	
    return 1;

if __name__ == "__main__":
    sys.exit(main())
