#!/usr/bin/env python

'''

*** Input options

-i, --input=file            The alignment file to be used for to generate the count table. The file
                                may be in the sam (zipped or not)format. The extension of
                                the file should contain either the '.sam' indication. Bam files are not
                                supported at the moment due to soome instability in the BAM reader 
                                regarding certain aligner formats.

-a, --annotation=file        The annotation file that will be used to generate the counts. The file should
                                be in the gff format (see https://www.sanger.ac.uk/resources/software/gff/spec.html
                                 for details).

-g, --genome_details=file    A tabulated file containing the names and sizes of each chromosome. !!! The
                                the chromosome names should be identical to the ones used to generate
                                the alignment file !!! The file should look like this example (no header):
                                    chromI    1304563
                                    chromII    6536634
                                    ...



*** Output options

-v, --verbose                When specified, the option will result in the printing of informations on the
                                alignments and the filtering steps (number of ambiguous alignments...etc).
                                Default: OFF

-o, --output_file=file        The output file where the count table should be stored. If the specified file
                                does not already exist it will be created automatically. Otherwise, it
                                will be overwritten



*** Filtering options

-q, --quality_threshold=INT    The quality threshold below which alignments will be discarded. The alignment 
                                quality index typically ranges between 1 and 41. Default: 30.

-l, --lowest_size=INT        The lowest possible size accepted for DNA fragments. Any pair of reads with an
                                insert size below that value will be discarded. Default: 130.

-L, --longest_size=INT        The longest possible size accepted for DNA fragments. Any pair of reads with an
                                insert size above this value will be discarded. Default: 180.

-d, --duplicate_filter=INT    The number of estimated PCR duplicates for accepted for a given genomic location.
                                Default: 1.



*** Count table options

-f, --feature_type=STRING    The feature types to be used when generating the count table. The feature type
                                 will be matched 3rd column of the GFF file. Default: 'Transcript'

-w, --downstream_window=INT    The window size used to obtain the counts downstream of the TSS. Default: 1000bp.

-u, --upstream_window=INT    The window size used to obtain the counts upstream of the TSS. Default: 1500bp.


'''



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#    Import relevant libraries and packages
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

from __future__ import print_function
import argparse, HTSeq, numpy, sys



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#    Get all the information on the pair of reads - and the genomic features
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def get_feature_name(feature):
    """Function to obtain a usable ID for the genomic features considered in the gff file
    """
    name = feature.attr['Name']
    return name


def alignmentQC_test(bundle,th_QC):
    """Function to check the alignment quality of the bundle of reads considered. The threshold may be 
    specified in the command line using the -qc option. The default value is 30 (based on illumina 1.8 system, ranging from 0 to 41)
    """
    bool = False
    for alignment in bundle:
        if alignment.aQual < th_QC:
            bool = True
    return bool


def is_divergent(strand1,strand2,start1,start2,stop1,stop2):
    """Function that checks if the pair of reads considered are mapped in a divergent fashion (inward, facing)
    to their reference genome
    """
    bool = True
    if strand1 == '+' and strand2 == '-':
        if (start1 <= start2) and (stop1 <= stop2):
            bool = False
    elif strand1 == '-' and strand2 == '+':
        if (start2 <= start1) and (stop2 <= stop1):
            bool = False
    return bool


def get_spanning_size(strand1,strand2,start1,start2,stop1,stop2):
    """Function that computes the insert size between paired end reads
    """
    size = ''
    fragment_start = ''
    fragment_stop = ''
    if strand1 == '+' and strand2 == '-':
        size = stop2 - start1
        fragment_start = start1
        fragment_stop = stop2
    elif strand1 == '-' and strand2 == '+':
        size = stop1 - start2
        fragment_start = start2
        fragment_stop = stop1
    return fragment_start, fragment_stop, size


def pair_information(read1, read2):
    """Function that returns all the important information concerning the pairs of reads.These include:
        - chromosomes, strands, interval boundaries
        - the presence of mapping aberrations (reads mapped on different chromosomes or on the same strand)
        - the assessment of the convergence of the pair (when applicable)
        - the insert size together with the start and stop coordinate of the asscoiated DNA fragment (when applicable)
    """
    ## initialize variables
    aberrations = ''
    divergent = ''
    insertSize = ''
    fragment_start = ''
    fragment_stop = ''

    ##    Retrieve the preleminary information on the pair
    chrom1 = read1.iv.chrom
    chrom2 = read2.iv.chrom

    strand1 = read1.iv.strand
    strand2 = read2.iv.strand

    start1 = read1.iv.start
    start2 = read2.iv.start

    stop1 = read1.iv.end
    stop2 = read2.iv.end

    ## check for an alignment of the pair without aberrations
    aberrations = not((chrom1 == chrom2) and (strand1 != strand2))

    ## check if the pair is divergent
    divergent = is_divergent(strand1,strand2,start1,start2,stop1,stop2)
    if not(divergent):## convergent pair
        fragment_start, fragment_stop, insertSize = get_spanning_size(strand1,strand2,start1,start2,stop1,stop2)
    return aberrations, divergent, insertSize, fragment_start, fragment_stop, chrom1


def make_PCR_index(chrom,start,stop):
    """Function thats generate a string unique for unique molecule and later used to estimate PCR duplicates
    The string is organised as followed:
        [chromosome]_[start coordinate]_[stop coordinate]
    """
    return '_'.join([chrom,str(start),str(stop)])





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#    Set up the input and output flows
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def load_genomic_array(chromsomeSizes):
    """ Creates a genomic array with the appropriate size to store the count information
    """
    chromosomeSizes = open(chromsomeSizes,'r')
    dic_counts_by_chrom = {}# Initialize variables

    for line in chromosomeSizes:
        name=line[:-1].split('\t')[0]
        size=int(line[:-1].split('\t')[1])
        dic_counts_by_chrom[name]=numpy.zeros(int(size),'i')
        print('{chrom}\t{S}'.format(chrom = name, S = size))
    chromosomeSizes.close()
    return dic_counts_by_chrom


def set_up_IO(fileIN,fileOUT,gff,downstream,upstream):
    '''Function that will open all the file required for the alignment processing
    '''
    ## Open alignment
    alignIN = HTSeq.SAM_Reader(fileIN)
    alignIN = HTSeq.bundle_multiple_alignments(alignIN)

    ## Open GFF file
    annotation = HTSeq.GFF_Reader(gff,end_included = True)

    ## Open output file - write the header
    countTable = open(fileOUT,'w')
    coordinates = '\t'.join(i for i in map(str,range(-upstream,downstream)))
    countTable.write('name\t{coord}\n'.format(coord = coordinates))
    return alignIN, annotation, countTable





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#    Print the processing reports
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def print_report_alignment(C, C_multiple, C_alignmentQC, C_aberration, C_size, C_PCR, C_used):
    '''Function that prints out a report when specified
    '''
    print('{count_all} pairs of reads processed'.format(count_all = C),file=sys.stderr)
    print('{count_multiple} pairs with multiple alignments'.format(count_multiple = C_multiple),file=sys.stderr)
    print('{count_alignmentQC} pairs with low quality alignments'.format(count_alignmentQC = C_alignmentQC),file=sys.stderr)
    print('{count_aberration} pairs with aberrations'.format(count_aberration = C_aberration),file=sys.stderr)
    print('{count_size} pairs with an insert size out of the specifiec boundaries'.format(count_size = C_size),file=sys.stderr)
    print('{count_PCR} estimated PCR duplicates'.format(count_PCR = C_PCR),file=sys.stderr)
    print('{count_retained} pairs used to generate the final count table'.format(count_retained = C_used),file=sys.stderr)


def print_report_countTable(processedFeatures, retainedFeature):
    '''Function that will print a report on the count table writing process
    '''
    print('{proc} annotated features processed in the annotation file provided'.format(proc = processedFeatures),file=sys.stderr)
    print('{ret} annotated features retained to write the count table'.format(ret = retainedFeature),file=sys.stderr)





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#    Processing of the alignment and writing to the count table
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def get_profile(counts,feature,downstream,upstream):
    '''Computes the counts around the annotated features
    '''
    strand = feature.iv.strand
    chrom = feature.iv.chrom
    profile = ''
    if strand == '+':
        START = feature.iv.start - upstream
        STOP = feature.iv.start + downstream
        try:
            profile = counts[chrom][START:STOP]
        except KeyError:
            print('!!! Chromosome names in the annotation file do not match the names provided in the Genome file',file=sys.stderr)
            print('Skipping feature {name}'.format(name = feature.name),file=sys.stderr)
    elif strand == '-':
        START = feature.iv.end - downstream
        STOP = feature.iv.end + upstream
        try:
            profile = numpy.flipud(counts[chrom][START:STOP])
        except KeyError:
            print('!!! Chromosome names in the annotation file do not match the names provided in the Genome file',file=sys.stderr)
            print('Skipping feature {name}'.format(name = feature.name),file=sys.stderr)
    return profile


def write_count_table(counts,annotation,featureType,downstream,upstream,countTable):
    '''Function that will write the final count table to a file, using the count information and the GFF file provided.
    '''
    ## Go through the GFF file
    processedFeatures = 0
    retainedFeature = 0
    for feature in annotation:
        processedFeatures += 1

        if not(feature.type in featureType):## Check that the feature type is correct
            continue

        retainedFeature += 1

        ## Obtain the occupancy profile
        profile = get_profile(counts,feature,downstream,upstream)
        profileSTR = map(str,profile)
        profileSTR_for_file = '\t'.join(profileSTR)

        ## Write the profile to a file
        if profile != '':
            countTable.write('{name}\t{counts}\n'.format(name = get_feature_name(feature), counts = profileSTR_for_file))
    return processedFeatures, retainedFeature


def parse_alignment(counts,alignIN,th_QC,th_low,th_high, th_PCR):
    '''Function that parses the main alignment file. Reads with ambiguous or low quality alignments, divergent orientation,
    or insert size not matching the indicated criteria will be filtered out. The center of the nucleosomal particle associated
    with the pair of reads will be evaluated using the center of the insert size
    '''

    ## Counter of aligned reads
    C = 0
    C_multiple = 0
    C_alignmentQC = 0
    C_aberration = 0
    C_size = 0
    C_PCR = 0
    C_used = 0
    PCR_duplicates = {}

    ## Go through the alignment file and update the count information
    for reads in alignIN:
        C += 1
        if C % 500000 == 0:
            print('{N} reads processed'.format(N = C),file=sys.stderr)

        if len(reads) != 2:# Multiple alignments
            C_multiple += 1
            continue

        ## Retreive the read informations
        read1 = reads[0]
        read2 = reads[1]

        if (read1 is None) or (read2 is None):#Check if the reads are aligned
            C_aberration += 1
            continue

        if (not read1.aligned) or (not read2.aligned):
            C_aberration += 1
            continue

        if alignmentQC_test(reads,th_QC):#Check the quality of the alignments
            C_alignmentQC += 1
#            C_aberration += 1
            continue

        ## Obtain information on the pair
        aberrations, divergent, insertSize, fragment_start, fragment_stop, chrom = pair_information(read1, read2)

        if aberrations or divergent:#Chek the proper orientation and alignment
            C_aberration += 1
            continue

        if (insertSize < th_low) or (insertSize > th_high):#Filter pairs with unadequate size
            C_size += 1
            continue

        ## Filter the PCR duplicates
        fragmentID = make_PCR_index(chrom, fragment_start, fragment_stop)
        if (fragmentID in PCR_duplicates):
            if PCR_duplicates[fragmentID] == th_PCR:
                C_PCR += 1
                continue
            else:
                PCR_duplicates[fragmentID] += 1
        else:
            PCR_duplicates[fragmentID] = 1

        ## Get the nucleosome center coordinate and update the counts
        center = fragment_start + ((fragment_stop - fragment_start)/2)
        counts[chrom][center] += 1
        C_used += 1
    return counts, C, C_multiple, C_alignmentQC, C_aberration, C_size, C_PCR, C_used


def get_ChIP_count_table(fileIN,fileOUT,gff,chromosomeSizes,downstream,upstream,th_QC,th_low,th_high,th_PCR,report,featureType):
    """Main function to generate the count table. Each pair of reads is examined independantly,
    and after all filtering steps, the center of the genomic interval covered is used as an approximation of the
    center of the nucleosome.
    """
    ## Open the files for input and output
    alignIN, annotation, countTable = set_up_IO(fileIN,fileOUT,gff,downstream,upstream)

    ## Load the genome in a numpy array to store the counts
    counts = load_genomic_array(chromosomeSizes)
    print('Genome file loaded',file=sys.stderr)

    ## Go through the alignment
    counts, C, C_multiple, C_alignmentQC, C_aberration, C_size, C_PCR, C_used = parse_alignment(counts,alignIN,th_QC,th_low,th_high,th_PCR)
    print('Alignment file processed',file=sys.stderr)

    ## Print a count report
    if report:
        print_report_alignment(C, C_multiple, C_alignmentQC, C_aberration, C_size, C_PCR, C_used)

    ## Go through the annotation and write the count table
    processedFeatures, retainedFeature = write_count_table(counts,annotation,featureType,downstream,upstream,countTable)
    print('Counts written to file {output}'.format(output = fileOUT),file=sys.stderr)

    ## Print a count table report
    if report:
        print_report_countTable(processedFeatures, retainedFeature)

    alignIN.close()
#    annotation.close()
    countTable.close()
    return ''





#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#    Main
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def main():## Command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_file', type=str, required=True, metavar='SAM/BAM',
    help="""The alignment file to be used for to generate the count table. The file
                                may be in the sam (zipped or not)format. The extension of
                                the file should contain either the '.sam' indication. Bam files are not
                                supported at the moment due to soome instability in the BAM reader 
                                regarding certain aligner formats.""")

    parser.add_argument('-a', '--annotation', dest='annot',type=str, required=True, metavar='GFF', 
    help="""The annotation file that will be used to generate the counts. The file should
                                be in the gff format (see https://www.sanger.ac.uk/resources/software/gff/spec.html
                                 for details).""")

    parser.add_argument('-g', '--genome_details', dest='genome',type=str, required=True, 
    metavar='Chromosome Sizes File', 
    help="""A tabulated file containing the names and sizes of each chromosome. !!! The
                                the chromosome names should be identical to the ones used to generate
                                the alignment file !!! The file should look like this example (no header):
                                    chromI    1304563
                                    chromII    6536634
                                    ...""")

    parser.add_argument('-o', '--output_file', dest = 'out', type=str, required=True, metavar='Count Table', 
    help="""The output file where the count table should be stored. If the specified file
                                does not already exist it will be created automatically. Otherwise, it
                                will be overwritten""")


    parser.add_argument('-v', '--verbose', dest = 'verbose', action='store_true', 
    help="""When specified, the option will result in the printing of informations on the
                                alignments and the filtering steps (number of ambiguous alignments...etc).
                                Default: OFF""")

    parser.add_argument('-q', '--quality_threshold', dest = 'th_qc', type=int, default=30, 
    help=""" The quality threshold below which alignments will be discarded. The alignment 
                                quality index typically ranges between 1 and 41. Default: 30.""")
                                
    parser.add_argument('-l', '--lowest_size', dest = 'th_low', type=int, default=130, 
    help="""The lowest possible size accepted for DNA fragments. Any pair of reads with an
                                insert size below that value will be discarded. Default: 130.""")

    parser.add_argument('-L', '--longest_size', dest = 'th_high', type=int, default=180, 
    help="""The longest possible size accepted for DNA fragments. Any pair of reads with an
                                insert size above this value will be discarded. Default: 180.""")


    parser.add_argument('-d', '--duplicate_filter', dest = 'th_PCR', type=int, default=1, 
    help="""The number of estimated PCR duplicates for accepted for a given genomic location.
                                Default: 1.""")
    
    parser.add_argument('-f', '--feature_type', dest = 'featureType', action='store', default='ORF', 
    help="""The feature types to be used when generating the count table. 
            The feature type will be matched 3rd column of the GFF file. 
            Default: 'Transcript'""")
    
    parser.add_argument('-w', '--downstream_window', dest = 'downstream', type=int, default=1500, 
    help="""The window size used to obtain the counts downstream of the TSS. Default: 1000bp""")
    
    parser.add_argument('-u', '--upstream_window', dest = 'upstream', type=int, default=1000, 
    help="""The window size used to obtain the counts upstream of the TSS. Default: 1500bp""")





    args = parser.parse_args()
    print(args,file=sys.stderr)

    ## Call the function to generate the count table
    get_ChIP_count_table(
        args.input_file,
        args.out,
        args.annot,
        args.genome,
        args.downstream,
        args.upstream,
        args.th_qc,
        args.th_low,
        args.th_high,
        args.th_PCR,
        args.verbose,
        args.featureType)


if __name__ == '__main__':
    main()

