import os
import glob
import argparse
from itertools import islice

# Function 1: Set up input arguments
def get_args():
	parser = argparse.ArgumentParser(description="Will take in a set of files, the .log and .tre files from Jenny's introgression PAUP runs and parse them into a usable format.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-l', '--logfilesdir', type=str, help='Path to a set of .log files (i.e. test1.log, test2.log, test3.log, etc.). Easiest to name this as a subdiretory within the working directory.', required=True)
	parser.add_argument('-t', '--treefilesdir', type=str, help='Path to a set of .tre files (i.e. test1.tre, test2.tre, test3.tre, etc.). Easiest to name this as a subdiretory within the working directory.', required = True)
	parser.add_argument('-o', '--outfile', type=str, help='Name of the output file with parsed data.', required = True)
	parser.add_argument('-w', '--workpath', type=str, help='Directory path where work is to be done. Must have the trailing "/".', required = True)

	args = parser.parse_args()
	LOGDIR = args.logfilesdir
	TREEDIR = args.treefilesdir
	OUTFILE = args.outfile
	WORKPATH = args.workpath

	return LOGDIR, TREEDIR, OUTFILE, WORKPATH
	
	####For testing with ipython only
#	LOGDIR = 'logs'
#	TREEDIR = 'trees'
#	OUTFILE = 'output.txt'
#	WORKPATH = '/lustre/scratch/daray/introgression_data_files/'
	
#Function 2: Extract data from log files
def LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR):
	os.chdir(WORKPATH + '/' + LOGDIR)
	LOGFILENAME = os.path.basename(LOGFILE)
	ITER = LOGFILENAME.split('.')[0]
	ITER = ITER.replace('test', '')
	#Get informative characters value for this log file
	with open(LOGFILE) as SEARCHLOG:
		for LINE in SEARCHLOG:
			LINE = LINE.rstrip()
#			print(LINE)
			if 'Number of included characters = ' in LINE:
				INFORMATIVE = LINE.split()[5] 
	#Get tree lengths
			if 'Length  ' in LINE:
				LENGTH1 = LINE.split()[1]
				LENGTH2 = LINE.split()[2]
				LENGTH3 = LINE.split()[3]

	#Get information from tree file
	os.chdir(WORKPATH + '/' + TREEDIR)
	TREEFILENAME = 'test' + ITER + '.tre'
#	print('Working on tree file: ' + TREEFILENAME)
	with open (TREEFILENAME) as SEARCHTREE:
	#Get taxon IDs
		for LINE in SEARCHTREE:
			if 'Translate' in LINE:
				TAXON1, TAXON2, TAXON3, TAXON4 = islice(SEARCHTREE, 4)
				TAXON1 = TAXON1.split()[1]
				TAXON1 = TAXON1.replace(',', '')
				TAXON2 = TAXON2.split()[1]
				TAXON2 = TAXON2.replace(',', '')
				TAXON3 = TAXON3.split()[1]
				TAXON3 = TAXON3.replace(',', '')
				TAXON4 = TAXON4.split()[1]
				TAXON4 = TAXON4.replace(',', '')
				QUARTET = TAXON1 + ',' + TAXON2 + ',' + TAXON3 + ',' + TAXON4
	return [INFORMATIVE, LENGTH1, LENGTH2, LENGTH3, TAXON1, TAXON2, TAXON3, TAXON4, QUARTET]

def GETTOPO(WORKPATH, LOGDIR, LOGFILE, TREEDIR, TAXON1, TAXON2, TAXON3, TAXON4):
	os.chdir(WORKPATH + '/' + LOGDIR)
	LOGFILENAME = os.path.basename(LOGFILE)
	ITER = LOGFILENAME.split('.')[0]
	ITER = ITER.replace('test', '')
#	print('Working on tree file: ' + ITER)
	os.chdir(WORKPATH + '/' + TREEDIR)
	TREEFILE = 'test' + ITER + '.tre'
	with open(TREEFILE) as SEARCHTREE:
		for LINE in SEARCHTREE:
			if "tree 'PAUP_1' = [&U] " in LINE:
				TOPO1 = LINE.split()[4]
				TOPO1 = TOPO1.replace(';', '')
				TOPO1 = TOPO1.replace('1', TAXON1)
				TOPO1 = TOPO1.replace('2', TAXON2)
				TOPO1 = TOPO1.replace('3', TAXON3)
				TOPO1 = TOPO1.replace('4', TAXON4)
#				print('Tree 1: ' + TOPO1)
			if "tree 'PAUP_2' = [&U] " in LINE:
				TOPO2 = LINE.split()[4]
				TOPO2 = TOPO2.replace(';', '')
				TOPO2 = TOPO2.replace('1', TAXON1)
				TOPO2 = TOPO2.replace('2', TAXON2)
				TOPO2 = TOPO2.replace('3', TAXON3)
				TOPO2 = TOPO2.replace('4', TAXON4)
#				print('Tree 2: ' + TOPO2)
			if "tree 'PAUP_3' = [&U] " in LINE:
				TOPO3 = LINE.split()[4]
				TOPO3 = TOPO3.replace(';', '')
				TOPO3 = TOPO3.replace('1', TAXON1)
				TOPO3 = TOPO3.replace('2', TAXON2)
				TOPO3 = TOPO3.replace('3', TAXON3)
				TOPO3 = TOPO3.replace('4', TAXON4)
#				print('Tree 3: ' + TOPO3)
	return [TOPO1, TOPO2, TOPO3]
                
def GETASTRAL(WORKPATH, LOGDIR, LOGFILE, TREEDIR, TAXON1, TAXON2, TAXON3, TAXON4):
	os.chdir(WORKPATH + '/' + LOGDIR)
	LOGFILENAME = os.path.basename(LOGFILE)
	ITER = LOGFILENAME.split('.')[0]
	ITER = ITER.replace('test', '')
#	print('Working on tree file: ' + ITER)
	os.chdir(WORKPATH + '/' + TREEDIR)
	TREEFILE = 'test' + ITER + '_accepted.tre'
	with open(TREEFILE) as SEARCHTREE:
		for LINE in SEARCHTREE:
			if "tree SIMPLE = [&U] " in LINE:
				ASTRAL = LINE.split()[4]
				ASTRAL = ASTRAL.replace(';', '')
				ASTRAL = ASTRAL.replace('1', TAXON1)
				ASTRAL = ASTRAL.replace('2', TAXON2)
				ASTRAL = ASTRAL.replace('3', TAXON3)
				ASTRAL = ASTRAL.replace('4', TAXON4)
	return ASTRAL

def GETSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, TOPO1, TOPO2, TOPO3):
#	print('TOPO1 is ' + TOPO1)
#	print('TOPO2 is ' + TOPO2)
#	print('TOPO3 is ' + TOPO3)
	SPLITS = [TAXON1 + ',' + TAXON2, TAXON1 + ',' + TAXON3, TAXON1 + ',' + TAXON4, TAXON2 + ',' + TAXON3, TAXON3 + ',' + TAXON4, TAXON2 + ',' + TAXON1, TAXON3 + ',' + TAXON1, TAXON4 + ',' + TAXON1, TAXON3 + ',' + TAXON2, TAXON4 + ',' + TAXON3, TAXON2 + ',' + TAXON4, TAXON4 + ',' + TAXON2]
	ALTSPLITS = [TAXON3 + ',' + TAXON4, TAXON2 + ',' + TAXON4, TAXON2 + ',' + TAXON3, TAXON1 + ',' + TAXON4, TAXON1 + ',' + TAXON2, TAXON3 + ',' + TAXON4, TAXON2 + ',' + TAXON4, TAXON3 + ',' + TAXON2, TAXON1 + ',' + TAXON4, TAXON1 + ',' + TAXON2, TAXON3 + ',' + TAXON1, TAXON3 + ',' + TAXON1]
	for X in range(0, 12):	
		TESTSPLIT = SPLITS[X]
#		print(TESTSPLIT)
		TESTALTSPLIT = ALTSPLITS[X] 
#		print(TESTALTSPLIT)
		if TESTSPLIT in TOPO1:
#			print(TESTSPLIT + ' is in TOPO1')
			TAX1TOPO1SPLIT1 = TESTSPLIT.split(',')[0]
			TAX2TOPO1SPLIT1 = TESTSPLIT.split(',')[1]
			TAX1TOPO1SPLIT2 = TESTALTSPLIT.split(',')[0]
			TAX2TOPO1SPLIT2 = TESTALTSPLIT.split(',')[1]
			TOPO1SPLIT1 = TAX1TOPO1SPLIT1 + ' ' + TAX2TOPO1SPLIT1
			TOPO1SPLIT2 = TAX1TOPO1SPLIT2 + ' ' + TAX2TOPO1SPLIT2
#			print(TOPO1SPLIT1)
#			print(TOPO1SPLIT2)
#			return [TOPO1SPLIT1, TOPO1SPLIT2]
	for X in range(0, 12):	
		TESTSPLIT = SPLITS[X]
#		print(TESTSPLIT)
		TESTALTSPLIT = ALTSPLITS[X] 
#		print(TESTALTSPLIT)
		if TESTSPLIT in TOPO2:
#			print(TESTSPLIT + ' is in TOPO2')
			TAX1TOPO2SPLIT1 = TESTSPLIT.split(',')[0]
			TAX2TOPO2SPLIT1 = TESTSPLIT.split(',')[1]
			TAX1TOPO2SPLIT2 = TESTALTSPLIT.split(',')[0]
			TAX2TOPO2SPLIT2 = TESTALTSPLIT.split(',')[1]
			TOPO2SPLIT1 = TAX1TOPO2SPLIT1 + ' ' + TAX2TOPO2SPLIT1
			TOPO2SPLIT2 = TAX1TOPO2SPLIT2 + ' ' + TAX2TOPO2SPLIT2
#			return [TOPO2SPLIT1, TOPO2SPLIT2]
	for X in range(0, 12):	
		TESTSPLIT = SPLITS[X]
#		print(TESTSPLIT)
		TESTALTSPLIT = ALTSPLITS[X] 
#		print(TESTALTSPLIT)
		if TESTSPLIT in TOPO3:
#			print(TESTSPLIT + ' is in TOPO3')
			TAX1TOPO3SPLIT1 = TESTSPLIT.split(',')[0]
			TAX2TOPO3SPLIT1 = TESTSPLIT.split(',')[1]
			TAX1TOPO3SPLIT2 = TESTALTSPLIT.split(',')[0]
			TAX2TOPO3SPLIT2 = TESTALTSPLIT.split(',')[1]
			TOPO3SPLIT1 = TAX1TOPO3SPLIT1 + ' ' + TAX2TOPO3SPLIT1
			TOPO3SPLIT2 = TAX1TOPO3SPLIT2 + ' ' + TAX2TOPO3SPLIT2
#			return [TOPO3SPLIT1, TOPO3SPLIT2]
	return [TOPO1SPLIT1, TOPO1SPLIT2, TOPO2SPLIT1, TOPO2SPLIT2, TOPO3SPLIT1, TOPO3SPLIT2]

def ASTRALSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, ASTRALTREE):
#	print('ASTRALTREE is ' + ASTRALTREE)
	SPLITS = [TAXON1 + ',' + TAXON2, TAXON1 + ',' + TAXON3, TAXON1 + ',' + TAXON4, TAXON2 + ',' + TAXON3, TAXON3 + ',' + TAXON4, TAXON2 + ',' + TAXON1, TAXON3 + ',' + TAXON1, TAXON4 + ',' + TAXON1, TAXON3 + ',' + TAXON2, TAXON4 + ',' + TAXON3, TAXON2 + ',' + TAXON4, TAXON4 + ',' + TAXON2]
	ALTSPLITS = [TAXON3 + ',' + TAXON4, TAXON2 + ',' + TAXON4, TAXON2 + ',' + TAXON3, TAXON1 + ',' + TAXON4, TAXON1 + ',' + TAXON2, TAXON3 + ',' + TAXON4, TAXON2 + ',' + TAXON4, TAXON3 + ',' + TAXON2, TAXON1 + ',' + TAXON4, TAXON1 + ',' + TAXON2, TAXON3 + ',' + TAXON1, TAXON3 + ',' + TAXON1]
	for X in range(0, 12):	
		TESTSPLIT = SPLITS[X]
#		print(TESTSPLIT)
		TESTALTSPLIT = ALTSPLITS[X] 
#		print(TESTALTSPLIT)
		if TESTSPLIT in ASTRALTREE:
#			print(TESTSPLIT + ' is in ASTRALTREE')
			TAX1ASTRALSPLIT1 = TESTSPLIT.split(',')[0]
			TAX2ASTRALSPLIT1 = TESTSPLIT.split(',')[1]
			TAX1ASTRALSPLIT2 = TESTALTSPLIT.split(',')[0]
			TAX2ASTRALSPLIT2 = TESTALTSPLIT.split(',')[1]
			ASTRALSPLIT1 = TAX1ASTRALSPLIT1 + ' ' + TAX2ASTRALSPLIT1
			ASTRALSPLIT2 = TAX1ASTRALSPLIT2 + ' ' + TAX2ASTRALSPLIT2
	return [ASTRALSPLIT1, ASTRALSPLIT2]

		
def main():	
##Get input arguments
	LOGDIR, TREEDIR, OUTFILE, WORKPATH = get_args()
	print('LOGDIR = ' + LOGDIR)
	print('TREEDIR = ' + TREEDIR)
	print('OUTFILE = ' + OUTFILE)
	print('WORKPATH = ' + WORKPATH)

#Go to working appropriate directory and get the right log file to work with
	OUT = open(OUTFILE, 'w+')
	OUT.write('Test run \t Quartet \t Split \t Informative characters \t tree length \tASTRAL\n')
	os.chdir(WORKPATH + '/' + LOGDIR)
	for LOGFILE in glob.glob('*.log'):
		print('Working on ' + LOGFILE)
		INFORMATIVE = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[0]
		QUARTET = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[8]
		LENGTH1 = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[1]
		LENGTH2 = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[2]
		LENGTH3 = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[3]
		TAXON1 = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[4]
		TAXON2 = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[5]
		TAXON3 = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[6]
		TAXON4 = LOGDATA(WORKPATH, LOGDIR, LOGFILE, TREEDIR)[7]
		print('Number of informative characters: ' + INFORMATIVE)
		print('This quartet: ' + QUARTET)
		print('Lentgh of tree 1: ' + LENGTH1)
		print('Length of tree 2: ' + LENGTH2)
		print('Length of tree 3: ' + LENGTH3)
        
		TOPO1 = GETTOPO(WORKPATH, LOGDIR, LOGFILE, TREEDIR, TAXON1, TAXON2, TAXON3, TAXON4)[0]
		TOPO2 = GETTOPO(WORKPATH, LOGDIR, LOGFILE, TREEDIR, TAXON1, TAXON2, TAXON3, TAXON4)[1]
		TOPO3 = GETTOPO(WORKPATH, LOGDIR, LOGFILE, TREEDIR, TAXON1, TAXON2, TAXON3, TAXON4)[2]
		print('TOPO1: ' + TOPO1)
		print('TOPO2: ' + TOPO2)
		print('TOPO3: ' + TOPO3)
        
		TOPO1SPLIT1 = GETSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, TOPO1, TOPO2, TOPO3)[0]
		TOPO1SPLIT2 = GETSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, TOPO1, TOPO2, TOPO3)[1]
		TOPO2SPLIT1 = GETSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, TOPO1, TOPO2, TOPO3)[2]
		TOPO2SPLIT2 = GETSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, TOPO1, TOPO2, TOPO3)[3]
		TOPO3SPLIT1 = GETSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, TOPO1, TOPO2, TOPO3)[4]
		TOPO3SPLIT2 = GETSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, TOPO1, TOPO2, TOPO3)[5]
		print('TOPO1SPLIT1: ' + TOPO1SPLIT1)
		print('TOPO1SPLIT2: ' + TOPO1SPLIT2)
		print('TOPO2SPLIT1: ' + TOPO2SPLIT1)
		print('TOPO2SPLIT2: ' + TOPO2SPLIT2)
		print('TOPO3SPLIT1: ' + TOPO3SPLIT1)
		print('TOPO3SPLIT2: ' + TOPO3SPLIT2)
        
		ASTRALTREE = GETASTRAL(WORKPATH, LOGDIR, LOGFILE, TREEDIR, TAXON1, TAXON2, TAXON3, TAXON4)
		ASTRALSPLIT1 = ASTRALSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, ASTRALTREE)[0]
		ASTRALSPLIT2 = ASTRALSPLITS(TAXON1, TAXON2, TAXON3, TAXON4, ASTRALTREE)[1]
		print('ASTRALTREE: ' + ASTRALTREE)
		print('ASTRALSPLIT1: ' + ASTRALSPLIT1)
		print('ASTRALSPLIT2: ' + ASTRALSPLIT2 + '\n')
        
		os.chdir(WORKPATH)       
		if ((TOPO1SPLIT1 == ASTRALSPLIT1) and (TOPO1SPLIT2 == ASTRALSPLIT2)):
			OUT.write(LOGFILE + '\t' + QUARTET + '\t' + TOPO1SPLIT1 + '|' + TOPO1SPLIT2 + '\t' + INFORMATIVE + '\t' + LENGTH1 + '\tASTRAL\n')
		else:
			OUT.write(LOGFILE + '\t' + QUARTET + '\t' + TOPO1SPLIT1 + '|' + TOPO1SPLIT2 + '\t' + INFORMATIVE + '\t' + LENGTH1 + '\n')        
		if ((TOPO2SPLIT1 == ASTRALSPLIT1) and (TOPO2SPLIT2 == ASTRALSPLIT2)):
			OUT.write(LOGFILE + '\t' + QUARTET + '\t' + TOPO2SPLIT1 + '|' + TOPO2SPLIT2 + '\t' + INFORMATIVE + '\t' + LENGTH2 + '\tASTRAL\n')
		else:
			OUT.write(LOGFILE + '\t' + QUARTET + '\t' + TOPO2SPLIT1 + '|' + TOPO2SPLIT2 + '\t' + INFORMATIVE + '\t' + LENGTH2 + '\n')        
		if ((TOPO3SPLIT1 == ASTRALSPLIT1) and (TOPO3SPLIT2 == ASTRALSPLIT2)):
			OUT.write(LOGFILE + '\t' + QUARTET + '\t' + TOPO3SPLIT1 + '|' + TOPO3SPLIT2 + '\t' + INFORMATIVE + '\t' + LENGTH3 + '\tASTRAL\n')
		else:
			OUT.write(LOGFILE + '\t' + QUARTET + '\t' + TOPO3SPLIT1 + '|' + TOPO3SPLIT2 + '\t' + INFORMATIVE + '\t' + LENGTH3 + '\n')        
	OUT.close()
	
if __name__ =="__main__":main()	
