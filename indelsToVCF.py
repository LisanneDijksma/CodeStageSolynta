import subprocess
import re
import sys


def openbed(firstrun, **configALL):
    """

    This function opens the BED file, retrieved from filelist object.
    Also checks whether the bed file is not empty.

    :param firstrun: First run of script (True or False)
    :param configALL: Dictionary containing all variables
    :return: the bed file

    """
    if firstrun:
        file = "%(outprfx)s.indels.bed" % configALL
        bed = open(file).readlines()
        return bed
    if not firstrun:
        file = "%(pathData)s/%(outprfx)s.bam.indels.bed" % configALL
        bed = open(file).readlines()
        return bed


def changeVCFname(**configALL):
    """
    This function will create VCF with indel name if the indels.bed is empty

    :param configALL: Dictionary containing all variables

    """
    convertName = "%(pathBCF)s view %(outprfx)s.goodheader.calls.vcf.gz -Oz -o %(outprfx)s.indels.vcf.gz \
    ; %(pathTabix)s -f %(outprfx)s.indels.vcf.gz" % configALL

    executer = subprocess.Popen(convertName, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                cwd="%(pathData)s" % configALL)
    executer.communicate()


def convertlines(bed, **configALL):
    """

    FIXME: check why sometimes an offset is detected in indel length
    This function will identify whether the INDELS are hetero- or homozygous based on threshold of > 0.8 en < 0.8.
    If there is an INDEL with score lower than 0.4 a warning is given

    :param bed: bed file
    :param configALL: Dictionary containing all variables
    :rtype: list
    :return: all lines containing indels that need to be added to VCF file

    """
    newline_list = []

    if len(bed) >= 1:
        for line in bed:
            oldline = line.split()
            indel_start_pos = int(oldline[1]) + 1
            indel_end_pos = int(oldline[2])
            alt_base = oldline[6]
            alt_base = list(re.sub("^[ACGT]*[ACGT]*", "", alt_base))

            if alt_base:
                alt_base = "".join(alt_base)
            else:
                alt_base = "-"


            if not re.search("^[ACGT].[0-9]*[ACGT]$", oldline[6]):
                #Means that INDEL is heterozygous
                if float(oldline[4]) < 0.8 and float(oldline[4]) >= 0.4:

                    if oldline[6].startswith("+"):
                        refseq = "%(pathSAM)s faidx %(reffasta)s %(chr)s:" % configALL + str(indel_start_pos) + "-" \
                                 + str(indel_start_pos) + " | sed 1d | tr -d '\n' "
                        seq = subprocess.Popen(refseq, stdout=subprocess.PIPE, shell=True, cwd="%(pathData)s"
                                                                                               % configALL)
                        base = seq.stdout.read().decode("utf-8").rstrip()

                        newline = """ "{0[0]}\t" """.format(oldline) + str(indel_start_pos) + """ "\t.\t" """ + base \
                                  + """ "\t" """ + base + alt_base \
                                  + """ "\t225\t.\tINDEL;DP={0[5]};AC=1;AN=2\tGT:PL\t0/1:.\n" """.format(oldline)
                        newline_list.extend([newline])

                    elif oldline[6].startswith("-"):

                        refseq = "%(pathSAM)s faidx %(pathRef)s %(chr)s:" % configALL + str(indel_start_pos) + "-" \
                                 + str(indel_end_pos) + " | sed 1d | tr -d '\n'"
                        seq = subprocess.Popen(refseq, stdout=subprocess.PIPE, shell=True, cwd="%(pathData)s"
                                                                                               % configALL)
                        base = seq.stdout.read().decode("utf-8").rstrip()

                        if re.search(r"-1[ACGT]$", oldline[6]):
                            alt_base = "-"
                            newline = """ "{0[0]}\t" """.format(oldline) + str(indel_start_pos) + """ "\t.\t" """ + \
                                      base + """ "\t" """ + alt_base \
                                      + """ "\t225\t.\tINDEL;DP={0[5]};AC=1;AN=2\tGT:PL\t0/1:.\n" """.format(oldline)
                            newline_list.extend([newline])
                        else:
                            alt_base = list(re.sub("[^A-Z]", "", alt_base))
                            alt_base = alt_base[0]
                            newline = """ "{0[0]}\t" """.format(oldline) + str(indel_start_pos) + """ "\t.\t" """ + \
                                      base + """ "\t" """ + alt_base \
                                      + """ "\t225\t.\tINDEL;DP={0[5]};AC=1;AN=2\tGT:PL\t0/1:.\n" """.format(oldline)
                            newline_list.extend([newline])

                #Means that INDEL is homozygous
                elif float(oldline[4]) >= 0.8 or float(oldline[4]) <= 1.0:

                    if oldline[6].startswith("+"):

                        refseq = "%(pathSAM)s faidx %(reffasta)s %(chr)s:" % configALL + str(indel_start_pos) + "-" \
                                 + str(indel_start_pos) + " | sed 1d | tr -d '\n' "

                        seq = subprocess.Popen(refseq, stdout=subprocess.PIPE, shell=True, cwd="%(pathData)s"
                                                                                               % configALL)
                        base = seq.stdout.read().decode("utf-8").rstrip()

                        newline = """ "{0[0]}\t" """.format(oldline) + str(indel_start_pos) + """ "\t.\t" """ + base \
                                  + """ "\t" """ + base + alt_base \
                                  + """ "\t225\t.\tINDEL;DP={0[5]};AC=2;AN=2\tGT:PL\t1/1:.\n" """.format(oldline)
                        newline_list.extend([newline])

                    elif oldline[6].startswith("-"):
                        refseq = "%(pathSAM)s faidx %(reffasta)s %(chr)s:" % configALL + str(indel_start_pos) + "-" \
                                 + str(indel_end_pos) + " | sed 1d | tr -d '\n'"
                        seq = subprocess.Popen(refseq, stdout=subprocess.PIPE, shell=True, cwd="%(pathData)s"
                                                                                               % configALL)
                        base = seq.stdout.read().decode("utf-8").rstrip()

                        if re.search(r"-1[ACGT]$", oldline[6]):
                            alt_base = "-"
                            newline = """ "{0[0]}\t" """.format(oldline) + str(indel_start_pos) + """ "\t.\t" """ \
                                      + base + """ "\t" """ + alt_base \
                                      + """ "\t225\t.\tINDEL;DP={0[5]};AC=2;AN=2\tGT:PL\t1/1:.\n" """.format(oldline)
                            newline_list.extend([newline])
                        else:
                            alt_base = list(re.sub("[^A-Z]", "", alt_base))
                            alt_base = alt_base[0]
                            newline = """ "{0[0]}\t" """.format(oldline) + str(indel_start_pos) + """ "\t.\t" """ +\
                                      base + """ "\t" """ + alt_base \
                                      + """ "\t225\t.\tINDEL;DP={0[5]};AC=2;AN=2\tGT:PL\t1/1:.\n" """.format(oldline)
                            newline_list.extend([newline])

                #Means that there is a phasing error/possible polyploid INDEL
                else:
                    log("Below zygosity threshold in %(file)s at position " % configALL, indel_start_pos, indel_end_pos)

                    raise NotImplementedError

    return newline_list


def addlinestovcf(newline_list, **configALL):
    """

    This function executes the commandlines which will add all the newlines to VCF file

    :param newline_list: list containing all the lines that need to be added to VCF
    :param configALL: Dictionary containing all variables

    """

    add_newlines = "".join(newline_list).replace(" ", "")

    createnewvcf = "%(pathBCF)s view %(outprfx)s.goodheader.calls.vcf.gz -o %(outprfx)s.addlines.vcf " \
                   "; cat %(outprfx)s.addlines.vcf " \
                   "| echo " % configALL + add_newlines \
                   + " >> %(outprfx)s.addlines.vcf " \
                     "; %(pathBCF)s sort %(outprfx)s.addlines.vcf | " \
                     "%(pathBCF)s view - -Oz -o %(outprfx)s.indels.vcf.gz " \
                     "; %(pathTabix)s -f %(outprfx)s.indels.vcf.gz" % configALL
    executer = subprocess.Popen(createnewvcf, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                cwd="%(pathData)s" % configALL)
    executer.communicate()


def main(firstrun, reffasta, **configALL):
    """

    This script will process indels from a bed file and add them to a VCF

    :param firstrun: script running for first time or not (True or False)
    :param reffasta: fasta from reference sequence
    :param configALL: Dictionary containing all variables
    """
    bed = openbed(firstrun, **configALL)
    #Checks if bed file is not empty, else will change the name for consistency
    if not bed:
        log("No indels detected in %(file)s in %(region)s" % configALL)
        changeVCFname(**configALL)
    else:
        configALL["reffasta"] = reffasta
        newline_list = convertlines(bed, **configALL)
        addlinestovcf(newline_list, **configALL)

def log(*args):
    """
    Writes to stdout
    :param args: thing to write to stdout

    """
    sys.stdout.write(str(args) + "\n")


