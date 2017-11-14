RTEventsCounter script usage and requirements
Copyright 2017 Peter Y. Wang
File name: RTEventsCounter.py
Version: 1.3
-------------------------------------------------
GPL statement:
This file is part of RTEventsCounter.

RTEventsCounter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RTEventsCounter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RTEventsCounter.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------
Requirement:
Python 2.7
    with modules: sys, os, re, traceback, time
-------------------------------------------------
Usage:

Non-random primers mode:
python RTEventsCounter.py samFilePath.sam refFilePath.fa primerFilePath.csv

Random primers mode:
python RTEventsCounter.py samFilePath.sam refFilePath.fa
-------------------------------------------------
Description:

RTEventsCounter was based partially on the ShapeMapper software code by Steven
Busan (2014) from the Weeks Lab in UNC, available from
http://www.chem.unc.edu/rna/software.html.
RTEventsCounter was designed primarily for RNA chemical probing experiments.
Specifically, RTEventsCounter was made for identifying and tallying reverse
transcription events from aligned next-generation sequencing data, and
organizing them into a detailed table. It also includes respective probability
values for each nucleotide.
RTEventsCounter supports mutations (total as well as broken down into each
outcome base), terminations (stops), deletions, and insertions. It can be used
with targeted probing with specific known primers, or data from random primers;
it is compatible with both paired-end and unpaired reads. It also includes a
limited range of available filters.
Outputs from separate runs can be concatenated into a single output CSV. It can
also log each run info for reference.
Example files have been included for reference.
-------------------------------------------------
Files needed:

conf.py
    Configuration file; must be within the same directory and named "conf.py".
    The formatting and default values are available in the file distributed
    alongside the main script file. All values must be present.

refFile.fa
    A FASTA file containing the reference sequences of targets being aligned to
    is required.
    Each target name should be on a single line and preceded with a ">", and the
    target names must match exactly, case-sensitive, the target names in the SAM
    and primer files. Whitespaces around (not within) the target names are
    allowed.
    The sequences should immediately follow on a new line after their respective
    target names. Line breaks and whitespaces within the sequences are allowed.
    The sequences should be in A/T/G/C; lower cases and U are automatically
    fixed by the script, but other non-whitespace characters may cause errors.
    The script will attempt to call and count ALL and ONLY the targets specified
    in the reference file. Use the appropriate FASTA to specify the targets to
    analyze in a data set.

samFile.sam
    A SAM file of the alignments of reads is required. The SAM file must contain
    at least all the required fields as outlined in SAM specifications.
    The Phred quality scores should be in Phred+33 encoding.
    The script can handle alignments of paired reads, non-paired reads, or a mix
    of both.
    If aligning pairs, the --fr option (not -rf or --ff) must be used in Bowtie
    2.

primerFile.csv
    A CSV file containing the primers information is required. The primer file
    can contain either (in sequences mode) the sequences of the parts in primer
    oligos that pair with the RNA, or (in coordinates mode) the coordinates of
    such regions (inclusive, 1-based, 5'->3' in the reference).
    The primer file must contain at least all the targets being aligned to, but
    can contain information of other targets as well. Therefore, it is possible
    to have a master primer file for different experiments and select the
    alignment targets in each read using respective FASTA reference files.
    The formatting must be exactly as follows:
    On each line, in sequences mode:
        target,primer_name,sequence
    On each line, in coordinates mode:
        target,primer_name,start,end
    Empty lines, and whitespaces not within each field, except from the
    sequences, are allowed, for formatting and visual appearance purposes.
    Note that, if using sequences mode, the script would first search for the
    coordinates using the sequences, and would take the first match only even if
    the sequence occurs multiple times.
-------------------------------------------------
Outputs:

output.csv
    The main output table of the counts. File name is changeable using the
    configuration file conf.py.
    If appending mode is off, the existing file with the same name will be
    overwritten. If appending is on, the script will append the output to the
    end without new headers, but only if it detects a non-empty file (it checks
    the first line to see if it contains any non-whitespace characters).

counterLog.txt
    If enabled in the configuration file, the script will print all stdout and
    stderr lines to a log file instead, with a timestamp using the operating
    system time. It includes all error notices, target and primer names, summary
    counts of reads, and the average reverse transcription lengths excluding
    primer regions.
