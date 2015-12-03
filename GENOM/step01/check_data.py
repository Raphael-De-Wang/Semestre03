#!env python

import argparse

def interface_standard():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-pheno', dest='pheno', help='input phenotype fiel name',required=True)
    parser.add_argument('-ped',   dest='ped'  , help='input ped file name',required=True)
    parser.add_argument('-of',    dest='ofname', help='output file name',required=True)
    return parser.parse_args()


def read_pheno(phenofile):
	handler = open(phenofile)
	pheno_id = []
	for line in handler:
		pheno_id.append(line.split('\t')[0])
        del pheno[0]
	return pheno_id

    
def read_ped(pedfile,pheno_id):
	handler = open(pedfile)
	ped_data = []
	missing_id = []
	for line in handler:
		current_id = line.split()[0]
		if current_id in pheno_id:
			ped_data.append(line[:-1])
		else:
			missing_id.append(current_id)
	return (ped_data,missing_id)

def write_ped(ofname,ped):    
    fo = open(ofname,"w")
    for i in ped:
	fo.write(i)
	fo.write('\n')
        fo.close()

    
if __name__ == "__main__":
    args = interface_standard()
    pheno = read_pheno(args.pheno)
    (ped,miss) = read_ped(args.ped,pheno)
    write_ped(args.ofname,ped)

