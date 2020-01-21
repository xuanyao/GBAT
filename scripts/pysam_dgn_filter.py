import pysam
import sys
import argparse
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--bfile',default=None, type=str,help='input bam file')
	parser.add_argument('--outdir',default=None, type=str,help='output directory')
	args = parser.parse_args()
	defaults = vars(parser.parse_args(''))
	opts = vars(args)

	bfile=args.bfile
	outfolder=args.outdir
	f = pysam.AlignmentFile(bfile,'rb')
	outfile=outfolder+"/"+bfile.split('/')[-1]
	out = pysam.AlignmentFile(outfile, "wb", template=f)


	for read in f:
		if(read.cigarstring):
			if(read.mapq==255):
				out.write(read)

	out.close()
	f.close()
