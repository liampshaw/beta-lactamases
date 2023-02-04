import argparse
import re


def get_options():
    parser = argparse.ArgumentParser(description='Convert prokka output to abricate-style output (for downstream plotting)')
    parser.add_argument('--prokka', help='Input prokka gff', type=str)
    return parser.parse_args()

def main():
	args = get_options()
	print("#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE")
	with open(args.prokka, 'r') as f:
		for line in f.readlines():
			if line.startswith("##FASTA"):
				break
			elif line.startswith("##"):
				pass
			else:
				line = line.strip('\n').split('\t')
				output_string = '\t'.join(["NA", line[0], line[3], line[4], line[6], re.sub("ID=", "", re.sub(";.*", "", re.sub(".*;gene=", "", line[8]))), \
										"NA", "NA", "NA", "NA", "NA", "NA", "NA",\
										re.sub(";.*", "", re.sub(".*;product=", "", line[8])),"NA"])
				print(output_string)





if __name__=="__main__":
	main()
