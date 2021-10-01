# start
# round 1-3 run genome assembly and map genomes to a reference genome
import glob
import os
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS of each species",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/selected_species/round1/allfasta/',
                      metavar='input/')
required.add_argument("-seq",
                      help="file extension of sequencing files",
                      type=str, default='_final.scaffolds.fasta',
                      metavar='.fa or .fasta')
required.add_argument("-ref",
                      help="kegg reference",
                      type=str, default='/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/AHR.hmm',
                      metavar='AHR.hmm')

# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='/scratch/users/anniz44/scripts/1MG/AHR/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='/scratch/users/anniz44/genomes/AHR/',
                      metavar='AHR/')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running hmmsearch (default 1)",
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-hmm',
                          help="Optional: complete path to hmmsearch if not in PATH",
                          metavar="/usr/local/bin/hmmsearch",
                          action='store', default='hmmsearch', type=str)
optional.add_argument('-pro',
                      help="Optional: complete path to prodigal if not in PATH, None for no prodigal (default)",
                      metavar="/usr/local/bin/prodigal",
                      action='store', default='prodigal', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
input_script = args.s
input_script_hmm = input_script + '/AHR_search'
output_dir = args.o + '/AHR_search'
database = args.ref
try:
    os.mkdir(args.o)
except IOError:
    pass

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(input_script_hmm)
except IOError:
    pass


################################################## Function ########################################################
def run_hmm(input_fasta):
    # run kegg
    cutoff = 0.01
    fastaname = os.path.split(input_fasta)[-1].split(args.seq)[0]
    input_faa = input_fasta.replace(args.seq,'.faa')
    cmds = ''
    try:
        f1 = open(input_faa,'r')
    except IOError:
        cmds +=(args.pro + ' -q -i %s -a %s\n' % (input_fasta, input_faa))
    cmds += ('%s --tblout %s/%s.ahr.txt --cpu 40 -E %s %s %s\n') % (
        args.hmm,output_dir,fastaname, cutoff, database, input_faa)
    return cmds


################################################## Main ########################################################
# generate code
allfasta = glob.glob('%s/*/*%s'%(args.i,args.seq))+glob.glob('%s/*/fasta/*%s'%(args.i,args.seq))
print(len(allfasta))
i = 0
allscript = []
for fasta in allfasta:
    newscript = os.path.join(input_script_hmm, str(i%100) + '.sh')# totally 100 subscripts
    cmds = run_hmm(fasta)
    if newscript not in allscript:
        f1 = open(newscript, 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
        f1.close()
        allscript.append(newscript)
    else:
        f1 = open(newscript, 'a')
        f1.write('%s' % (cmds))
        f1.close()
    i += 1

# sum all codes
f1 = open(os.path.join(input_script, 'allAHRsearch.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in allscript:
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allAHRsearch.sh'%(input_script))
################################################### END ########################################################
