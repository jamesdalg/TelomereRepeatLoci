# Author: Lina Sieverling

#!/usr/bin/python

# Usage: python /home/sieverli/Code/telomere_insertion_analysis/snakemake_telomere_insertions/src/add_mate_mapq.py \
#                -i <*_tumor_discordant_reads.tsv> \
#                -b <*_merged.mdup.bam> \
#                -o <*_discordant_reads_filtered_with_mapq.tsv>


# Description: - parses tables with telomere insertion reads
#              - skips all reads where mates are mapped to decoy sequences
#              - retrieves mapping quality of mates from original BAM file and adds it to output table
#              - if mate is not found, mapping quality is empty


import os
import sys, getopt
import re
import numpy

# ----------------------------------------------------------------
# read command line args 
# ----------------------------------------------------------------
myopts, args = getopt.getopt(sys.argv[1:],"i:b:o:")
 
for opt, arg in myopts:
  if opt == '-i':
    telomere_insertion_table_file = arg
  elif opt == '-b':
    alignment_bam_file = arg
  elif opt == '-o':
    outfile_path = arg
  else:
    print("Usage: %s -i input_table -b bam_file -o output_file" % sys.argv[0])


# list of chromosomes accepted for output
chromosome_list = [ str(i) for i in range(1,22+1) ] + ["X","Y"]
chromosome_list_with_chr = ["chr" + s for s in chromosome_list]
chromosome_list = chromosome_list_with_chr

#########################################################################################################################################

# ----------------------------------------------------------------
# read in table containing read names from telomere insertions
# ----------------------------------------------------------------

telomere_insertion_table = numpy.genfromtxt(telomere_insertion_table_file, skip_header=1, delimiter="\t", dtype=None, comments=None)


# ----------------------------------------------------------------
# get mapping quality and strand of mate from original BAM file
# ----------------------------------------------------------------

output = '\t'.join(["read_name", "mate_chr", "mate_position", "mate_mapq", "mate_strand"])
for read in telomere_insertion_table:

  read_name = read[0]
  chromosome = read[1]
  position = str(read[2])

  # skip mates mapped to decoy sequences

  chromosome=chromosome.decode('ascii')
  if chromosome not in chromosome_list:
    # if "chrUn" not in str(chromosome):
    #   if "random" not in str(chromosome):
    #     if "chrM" not in str(chromosome):
    #       import ipdb
    #       ipdb.set_trace()
    #       print(type(chromosome));
    #       print(type(str(chromosome)));
    #       print(type(chromosome_list));
    #       print("lost in decoy" + str(chromosome));
    #       print(str(chromosome_list));
    continue
  read_name = read_name.decode('ascii')
  #import ipdb
  #ipdb.set_trace()
  #print(type(chromosome));
  # get original read of mate (skip secondary and supplementary alignments)
  extract_mq_command = "samtools view -F 2304 " + alignment_bam_file + " " + chromosome + ":" + position + "-" + position + "| grep  \"" + read_name + "\" #   | cut -f 5";
  #print(extract_mq_command);

  original_read = os.popen(extract_mq_command).read().rstrip()
#  import ipdb
#  ipdb.set_trace()
 # print(original_read);
 # print()
  original_read_list = original_read.split("\t")
  #print(original_read_list);

  try:
    mapq = original_read_list[4]
  except:
    mapq = ''

  try:
    flag = int(original_read_list[1])
    if flag & 0x10:
      strand = "-"
    else:
      strand = "+"
  except:
    strand = ''

  #print(mapq)
  read_list = [read_name, chromosome, position, mapq, strand]

  output += '\n' + '\t'.join(read_list) 




# ----------------------------------------------------------------
# write output
# ----------------------------------------------------------------

outfile = open( outfile_path, "w")  
outfile.write(output)       
outfile.close



