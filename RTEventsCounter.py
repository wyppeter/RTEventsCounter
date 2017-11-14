#!/usr/bin/env python
#
# Python script for parsing through sequence alignments to map and count RT mutations, indels, and stop events
# Copyright 2017 Peter Y. Wang
# Simon Lab, Yale University, New Haven, CT
# File name: RTEventsCounter.py
# Version: 1.3
# Python version: 2.7
#
# Parts of this code are based on the ShapeMapper software code by Steven Busan (2014) from the Weeks Lab in UNC
# Available from http://www.chem.unc.edu/rna/software.html
# --------------------------------
# GPL statement:
# This file is part of RTEventsCounter.
#
# RTEventsCounter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RTEventsCounter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RTEventsCounter.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import re
import traceback
import time

# Get configurations
try:
    import conf
except:
    sys.exit("Error importing configuration file conf.py: %s\n" % traceback.format_exc())

# Check file paths arguments inputs
if not conf.randomPrimers and len(sys.argv) != 4:
    sys.stdout.write("Usage: python RTEventsCounter.py samFilePath.sam refFilePath.fa primerFilePath.csv\n")
    sys.exit(1)
elif conf.randomPrimers and len(sys.argv) != 3:
    sys.stdout.write("Usage: python RTEventsCounter.py samFilePath.sam refFilePath.fa\n")
    sys.exit(1)

# Prepare logging
try:
    if conf.saveLog:
        logFile = sys.stdout = sys.stderr = open("counterLog.txt", "a")  # Redirect outputs
        logFile.write("\n>> %s\n" % time.asctime())  # Timestamp
except:
    sys.exit("Error preparing log: %s\n" % traceback.format_exc())

try:
    # Constants and functions
    NUCS = ["A", "T", "G", "C"]
    fileOut = None

    def safeDivide(num, den):
        # Divide numbers as floats; NA for zero division
        try:
            if num == "NA" or den == "NA":
                return "NA"
            else:
                return float(num)/float(den)
        except ZeroDivisionError:
            return "NA"
        except:
            sys.stderr.write("Unexpected error at division: %s\n" % sys.exc_info()[:2])
            return "NA"

    def checkSeq(lineinput):
        # Check and fix FASTA sequence lines while parsing
        line = "".join(lineinput.split())  # Remove all whitespace
        uFound = unknownFound = False
        for n in line:
            if n in ["a", "t", "g", "c"]:  # Lower case
                line = line.replace(n, n.upper())
            elif n in ["U", "u"]:  # U to T
                uFound = True
                line = line.replace(n, "T")
            elif n not in NUCS:  # Non-nucleotide characters
                unknownFound = True
                line = line.replace(n, "N")
        if uFound:
            sys.stderr.write("FASTA notice: U found and substituted as T\n")
        if unknownFound:
            sys.stderr.write("FASTA warning: Non-nucleotide character found!\n")
        return line

    def parseFasta(fastaFile):
        # Takes a FASTA file and returns a dict of sequences associated with their names
        seqs = {}
        lines = fastaFile.readlines()
        seqName = ""
        for line in lines:
            if line[0] == ">":
                # Sequence name
                seqName = line.strip()[1:]
                if seqName in seqs:
                    sys.stderr.write("FASTA error: Duplicate seq name (%s) found in FASTA!\n" % seqName)
                else:
                    seqs[seqName] = ""
            elif seqName != "":
                # Sequence, strip whitespace
                seqs[seqName] += checkSeq("".join(line.split()))
            else:
                sys.stderr.write("FASTA error: Skipped FASTA line, no seq name associated.\n")
        return seqs

    def primerSearch(primerDict, refDict):
        # Search for primers and return a dict of primer coordinates (inclusive, 1-based)
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        coordinatesDict = {}
        for primerSet in primerDict:
            if primerSet in refDict:  # Only use primers for targets given in FASTA
                coordinatesDict[primerSet] = {}
                for primerName in primerDict[primerSet]:
                    primerSeq = primerDict[primerSet][primerName][::-1]  # Reverse sequence since primer is rev comp
                    refSeq = "".join([comp.get(n, "N") for n in refDict[primerSet]])  # Complement sequence since primer is rev comp
                    match = re.search(primerSeq, refSeq)
                    if match:
                        # Found primer!
                        coordinatesDict[primerSet][primerName] = [match.start()+1, match.end()]
                        match = None
                    else:
                        sys.stderr.write("Primer search error: Primer seq (%s, %s) not found in sequence.\n" % (primerSet, primerName))
        return coordinatesDict

    def parsePrimers(primersLines, mode):
        # Takes a .csv table primer coordinates or sequences and returns them in a dict
        primers = {}
        lines = primersLines.readlines()
        for line in lines:
            if line.strip() != "":  # Allow empty lines for formatting in the csv
                lineSplit = line.split(",")
                if (mode == "seqs" and len(lineSplit) == 3) or (mode == "coords" and len(lineSplit) == 4):
                    setName = lineSplit[0].strip()
                    if setName not in primers:
                        # Make a new primer dict key
                        primers[setName] = {}
                    if mode == "seqs":
                        primers[setName][lineSplit[1].strip()] = checkSeq(lineSplit[2].strip())
                    else:  # coord mode
                        primers[setName][lineSplit[1].strip()] = [int(lineSplit[2].strip()), int(lineSplit[3].strip())]
                else:
                    sys.stderr.write("Primer parse error: Skipped primer line, invalid format.\n\tError line: %s\n" % line)
        return primers

    def parseCigarString(cigarString, read, quals, refSeq, startIndex):
        # Parse CIGAR string and generate raw events and quals
        # Resolve CIGAR strings into lists of numbers and event types
        cigarRegRE = re.compile(r"[0-9]+[M|I|D|S]")  # Note: Only CIGAR strings with M, I, D, S are allowed
        numberRE = re.compile(r"(^[0-9]+)")
        splitCigarString = cigarRegRE.findall(cigarString.strip())
        cigarList = [numberRE.split(reg)[1:] for reg in splitCigarString]  # Break CIGAR into each event

        # Prepare output dicts
        events = {}
        alignedQuals = {}

        # Prepare variables
        lastQual = "#"  # For assigning qual to deleted region
        refIndex = startIndex  # Index for along the ref seq; begin at leftmost matching position; 0-based
        readQualIndex = 0  # Index to select read and qual nucs by event region chunks

        # Go through every CIGAR region
        try:
            for i in xrange(len(cigarList)):
                region = cigarList[i]
                regionLength = int(region[0])
                regionType = region[1]

                # Match/mismatch
                if regionType == "M":
                    matchRegion = read[readQualIndex: readQualIndex + regionLength]
                    qualsRegion = quals[readQualIndex: readQualIndex + regionLength]
                    readQualIndex += regionLength
                    for mIndex in xrange(regionLength):
                        nuc = matchRegion[mIndex]
                        qual = qualsRegion[mIndex]
                        lastQual = qual
                        alignedQuals[refIndex] = qual
                        if nuc != refSeq[refIndex]:
                            events[refIndex] = nuc
                        else:
                            events[refIndex] = "|"
                        refIndex += 1
                # Insertion
                elif regionType == "I":
                    if (conf.countInserts and
                            events[refIndex-1] == "|" and refSeq[refIndex] == read[readQualIndex + regionLength] and
                            read[readQualIndex] != refSeq[refIndex] and read[readQualIndex + regionLength-1] != refSeq[refIndex-1]):
                        # Record on 5 prime neighboring nuc if unambiguous, and not flanked by mismatches
                        events[refIndex-1] = "^"
                    readQualIndex += regionLength
                # Deletion
                elif regionType == "D":
                    for delIndex in xrange(refIndex, refIndex+regionLength):
                        events[delIndex] = "-"
                        alignedQuals[delIndex] = lastQual
                    refIndex += regionLength
                # Soft clipping
                elif regionType == "S":
                    try:
                        readQualIndex += regionLength
                    except IndexError:
                        pass
                    clipRange = None
                    if i == len(cigarList)-1:  # Rightmost end of read/CIGAR
                        clipRange = xrange(refIndex, refIndex+regionLength)
                    elif i == 0:  # Leftmost end of read/CIGAR
                        clipRange = xrange(refIndex, refIndex-regionLength-1, -1)
                    if clipRange is not None:
                        for clipIndex in clipRange:
                            if clipIndex >= 0 and clipIndex < len(refSeq):
                                events[clipIndex] = "s"
                                alignedQuals[clipIndex] = "#"
        except (IndexError, KeyError):
            sys.stderr.write("CIGAR parsing index/key error!\n")
            return None, None

        return events, alignedQuals

    def cullEvents(events, refSeq):
        # Cull events
        eventChars = ["s", "-", "A", "T", "G", "C", "N", "^"]
        culledEvents = dict(events)

        # Scan left to right
        for i in xrange(min(events.keys()), max(events.keys())+1):
            if i not in events.keys():
                # Absent gaps in events
                culledEvents[i] = "~"
            elif events[i] in NUCS:
                # Mismatch
                try:
                    if events[i+1] in eventChars:
                        # With adjacent downstream mutation events
                        # Take only the most 3 prime event when consecutive events
                        culledEvents[i] = "|"
                except KeyError:
                    pass
            elif events[i] == "-":
                # Deletion
                if events[i+1] not in eventChars:
                    # With no adjacent downstream mutation events (|/~)
                    # Handle ambiguous deletions
                    k = 1
                    while events[i-k] == "-":
                        k += 1
                    ambiguous = False
                    # Slide deletion upstream and downstream
                    # Note: 5 prime start is i-k+1; 3 prime end is i; maximum offset is k
                    notDelSeq = refSeq[:i-k+1] + refSeq[i+1:]
                    for offset in xrange(-k, k+1):
                        notSubSeq = ""
                        if offset != 0:
                            try:
                                notSubSeq = refSeq[:i-k+1+offset] + refSeq[i+offset+1:]
                            except IndexError:
                                pass
                        if notSubSeq == notDelSeq:
                            ambiguous = True
                    if not ambiguous:
                        culledEvents[i] = "-"
                    else:
                        # Remove ambiguous deletion
                        culledEvents[i] = "|"
                        try:
                            offset = -1
                            while events[i+offset] == "-":
                                culledEvents[i+offset] = "|"
                                offset -= 1
                        except KeyError:
                            pass
                elif events[i+1] != "-":
                    # At the edge with soft clip (s), or with downstream mutation (A/T/G/C/N), or insert (^)
                    # Takes only the most 3 prime as event
                    culledEvents[i] = "~"
                else:
                    # Long deletion (-)
                    culledEvents[i] = "~"
                    if events[i+2] != "-":
                        # start of long deletion at i+1; record as stop if opted
                        if conf.longDelsAs == "stops":
                            culledEvents[i+1] = "!"
                        elif conf.longDelsAs == "none":
                            culledEvents[i+1] = "~"
            elif events[i] not in ["~", "|", "s", "^"]:
                # Unwanted events; replace with match character
                culledEvents[i] = "|"

        # Record stop events
        b = min(culledEvents.keys())
        if b != 0:
            # Skip over soft clips
            if culledEvents[b] == "s":
                while culledEvents[b] == "s":
                    b += 1
            # Must not be 3 prime of events
            if culledEvents[b] not in eventChars:
                culledEvents[b-1] = "!"

        return culledEvents

    def parseSamLine(line):
        # Takes a SAM file line of alignment and returns a list of the components
        splitLine = line.split()
        # Get relevant fields
        try:
            refName = splitLine[2]
            if refName not in refSeqs.keys():
                return ["skip", None, None, None, None, None]
            clusterName = splitLine[0]
            startIndex = int(splitLine[3])-1  # 0-based
            mapQual = int(splitLine[4])
            cigarString = splitLine[5]
            rawRead = splitLine[9]
            rawQual = splitLine[10]
        except IndexError:
            sys.stderr.write("SAM line parsing index error.\n")
            return ["skip", None, None, None, None, None]

        # Consolidate grabbed values
        if cigarString != "*":
            events, qualities = parseCigarString(cigarString, rawRead, rawQual, refSeqs[refName], startIndex)
            if events is not None:
                return [clusterName, refName, events, qualities, mapQual, rawRead]
            else:
                return ["skip", None, None, None, None, None]
        else:
            return ["skip", None, None, None, None, None]

    def combineEvents(R1, R2, Q1, Q2):
        # Combine two events/quals dicts for F/R seqs
        combinedEvents = dict(R1)
        combinedQuals = dict(Q1)
        # Put seq2 in first, and then slot seq1 in together
        for i in R2.keys():
            combinedEvents[i] = R2[i]
            combinedQuals[i] = Q2[i]
            if i in R1.keys():
                combinedQuals[i] = chr(int((ord(Q1[i]) + ord(Q2[i]))/2.0))
                # Handle positions that do not agree between mate pairs
                if R1[i] != R2[i]:
                    try:
                        # Disregard soft clips
                        if R1[i] == "s":
                            combinedEvents[i] = R2[i]
                        elif R2[i] == "s":
                            combinedEvents[i] = R1[i]
                        # Take higher qual events
                        elif ord(Q2[i])-33 > ord(Q1[i])-33:
                            combinedEvents[i] = R2[i]
                        else:
                            combinedEvents[i] = R1[i]
                    except IndexError:
                        sys.stdout.write("len(R1) = %i, len(Q1) = %i, len(R2) = %i, len(Q2) = %i\n" % (len(R1), len(Q1), len(R2), len(Q2)))
                        raise
        return combinedEvents, combinedQuals

    # Start main code

    # Note down arguments if logging
    if conf.saveLog:
        sys.stdout.write("Arguments: %s\n" % " ".join(sys.argv))

    # Note down key configs
    sys.stdout.write("Configurations: ")
    # Pairing
    if conf.alignPaired:
        if conf.ignoreUnpaired:
            sys.stdout.write("paired (exc. unpaired)")
        else:
            sys.stdout.write("paired (inc. unpaired)")
    else:
        sys.stdout.write("unpaired")
    # Primer Type
    if conf.randomPrimers:
        sys.stdout.write(", random primers\n\n")
    else:
        if conf.includeMisprimed:
            sys.stdout.write(", non-random primers (inc. misprimed)")
        else:
            sys.stdout.write(", non-random primers (exc. misprimed)")
        # Limited range
        if conf.eventsRange != [0, 0]:
            sys.stdout.write(", limited range to %s\n\n" % str(conf.eventsRange))
        else:
            sys.stdout.write(", unlimited range\n\n")

    # Get files and input data
    # Get SAM and prepare iterator
    with open(sys.argv[1], "rU") as samFile:
        samFilelines = samFile.readlines()
    samFileiter = iter(samFilelines)
    # Save sample name; remove path and filename extension
    sampleName = os.path.splitext(os.path.split(sys.argv[1])[1])[0]

    # Get reference seqs
    with open(sys.argv[2], "rU") as refFile:
        refSeqs = parseFasta(refFile)

    # Get primer seqs if not in random primer mode
    if not conf.randomPrimers:
        try:
            with open(sys.argv[3], "rU") as primerFile:
                # Check primer import mode and import accordingly
                # Take only primer entries of targets included in the FASTA ref
                if conf.primerInputModeSeqs:
                    # Use sequences; need to search for coordinates first
                    primerAll = primerSearch(parsePrimers(primerFile, "seqs"), refSeqs)
                    primerSeqs = {k: primerAll[k] for k in refSeqs.keys()}
                else:
                    # Use coordinates
                    primerAll = parsePrimers(primerFile, "coords")
                    primerSeqs = {k: primerAll[k] for k in refSeqs.keys()}
        except:
            sys.stderr.write("Primer list import failed, file at '%s':\n" % sys.argv[3])
            raise

    sys.stdout.write("Sample: %s\nReference sequence name(s):\n" % sampleName)

    # Set up primers
    primerKeyList = {}
    if conf.randomPrimers:
        for seqName in refSeqs:
            primerKeyList[seqName] = ["NA"]
            sys.stdout.write("- %s\t(%i nt)\n" % (seqName, len(refSeqs[seqName])))
    else:
        if conf.includeMisprimed:
            baseKeyList = ["NA"]
        else:
            baseKeyList = []
        for seqName in refSeqs:
            # For primers management
            primerKeyList[seqName] = sorted(primerSeqs[seqName].keys()) + baseKeyList
            # For primers printing
            printPrimerList = ",".join([(" %s (%s)" % (key, primerSeqs[seqName][key][0])) for key in sorted(primerSeqs[seqName].keys())])
            sys.stdout.write("- %s\t(%i nt)\n\t[Primers: %s]\n" % (seqName, len(refSeqs[seqName]), printPrimerList))

    # Build data structure and summary
    if conf.eventsRange != [0, 0]:  # Determine limited range mode
        # Truncating events; the lists storing events will store the only relevant range of indices, and so will not be in sync with the ref seq indices
        # Index shift determined later for each primer set
        eventsLimited = True  # Boolean to track if shifting will be needed
        if not conf.randomPrimers:  # Random primer mode will not limit range
            eventsBrackets = {}
            rangeLen = conf.eventsRange[1] - conf.eventsRange[0]
    else:
        eventsLimited = False
    tempFullRange = False  # For flow control below for misprimed reads
    mismatch = {}  # Mismatch event: mismatch[seqName][primerName][fromNuc][toNuc][i]
    mutationCount = {}  # Total mut count across mismatch types: mutationCount[seqName][primerName][i]
    deletion = {}  # Deletion events: deletion[seqName][primerName][fromNuc][i]
    stops = {}  # Stop events: stops[seqName][primerName][i]
    insertion = {}  # Insertion events: insertion[seqName][primerName][i]
    depth = {}  # Coverage depth: depth[seqName][primerName][i]
    RT = {}  # RT read through: RT[seqName][primerName][i]
    RTLength = {}  # Recording length of each RT: RTLength[seqName][primerName][...]
    qualsList = {}  # Average alignment quals per nuc: qualsList[seqName][primerName][i][...]
    sumSamLines = 0
    sumSamLowQual = 0
    sumSamDiscard = 0
    sumUnpairedReads = 0
    sumNoPrimer = 0
    sumNoExtension = 0
    sumCulledEventsPerGroup = {}
    for seqName in refSeqs:
        if eventsLimited and not conf.randomPrimers:
            eventsBrackets[seqName] = {}  # Prepare events range limit data
        else:
            rangeLen = len(refSeqs[seqName])  # Event lists should simply follow ref seq indices
        sumCulledEventsPerGroup[seqName] = {}
        mismatch[seqName] = {}
        mutationCount[seqName] = {}
        deletion[seqName] = {}
        stops[seqName] = {}
        insertion[seqName] = {}
        depth[seqName] = {}
        RT[seqName] = {}
        RTLength[seqName] = {}
        qualsList[seqName] = {}
        for primerName in primerKeyList[seqName]:
            if eventsLimited and primerName != "NA":
                # Work out the respective allowed nuc ranges for each primer if eventsRange defined
                eventsBrackets[seqName][primerName] = [primerSeqs[seqName][primerName][0] - conf.eventsRange[1]-1,
                                                       primerSeqs[seqName][primerName][0] - conf.eventsRange[0]-1]
            elif not conf.randomPrimers and primerName == "NA":
                tempFullRange = True  # Full range only for misprimed bucket
                holdRange = rangeLen
                rangeLen = len(refSeqs[seqName])
            sumCulledEventsPerGroup[seqName][primerName] = 0
            mismatch[seqName][primerName] = {}
            mutationCount[seqName][primerName] = [0]*rangeLen
            deletion[seqName][primerName] = {}
            stops[seqName][primerName] = [0]*rangeLen
            insertion[seqName][primerName] = [0]*rangeLen
            depth[seqName][primerName] = [0]*rangeLen
            RT[seqName][primerName] = [0]*rangeLen
            RTLength[seqName][primerName] = []
            qualsList[seqName][primerName] = [[] for n in xrange(rangeLen)]
            for fromNuc in NUCS:
                deletion[seqName][primerName][fromNuc] = [0]*rangeLen
                mismatch[seqName][primerName][fromNuc] = {}
                for toNuc in [mutNuc for mutNuc in NUCS if mutNuc != fromNuc]:  # Only dict for mismatches
                    mismatch[seqName][primerName][fromNuc][toNuc] = [0]*rangeLen
            if tempFullRange:
                rangeLen = holdRange  # Restore

    # Start initial parsing
    sys.stdout.write("\nParsing and analyzing events...\n")
    # Go through sequences and produce counts
    covered = ["|", "-", "A", "T", "G", "C", "^"]  # for depth
    readThrough = ["|", "-", "A", "T", "G", "C", "^", "~"]  # for RT
    # For recording previous data for pair combination
    prevClusterName = ""
    prevTargetName = ""
    prevLine = ""
    prevEvents = []
    prevQualities = []
    prevRawRead = ""
    # Control booleans
    pairFound = False
    updateCounts = False

    while True:
        # Get SAM line data and generate events and quals
        # Loop goes on until it hits the end of the SAM (StopIteration of iterator of SAM lines)
        try:
            line = samFileiter.next()
        except StopIteration:
            break
        if line[0] == "@":
            # Skip headers
            continue
        parsedLine = parseSamLine(line)
        sumSamLines += 1
        clusterName = parsedLine[0]
        targetName = parsedLine[1]
        events = parsedLine[2]
        qualities = parsedLine[3]
        mappingQual = parsedLine[4]
        rawRead = parsedLine[5]

        isCombinedRead = False  # Track combined pairs

        if clusterName == "skip":
            # Loop through each line of data
            sumSamDiscard += 1
        elif mappingQual < conf.minMapQual:
            # Exclude reads whose mapping quals are too low
            sumSamLowQual += 1
        else:
            # First, handle events and combine pair reads
            if conf.alignPaired:
                if not pairFound:
                    if clusterName == prevClusterName and targetName == prevTargetName:
                        # Found matching pair; combine reads
                        pairFound = True
                        updateCounts = True
                        isCombinedRead = True
                        eventsToWrite, qualsToWrite = combineEvents(prevEvents, events, prevQualities, qualities)
                    elif prevClusterName != "":
                        # Non-matching read
                        # Update with prev info if counting unpaired, but not in the 1st line case (no 0th line)
                        sumUnpairedReads += 1
                        if not conf.ignoreUnpaired:
                            eventsToWrite = dict(prevEvents)
                            qualsToWrite = dict(prevQualities)
                            pairFound = False
                            updateCounts = True
                else:
                    # Pair found
                    pairFound = False
                    updateCounts = False
            elif prevClusterName != "":  # Non-pairing mode
                updateCounts = True
                eventsToWrite = dict(prevEvents)
                qualsToWrite = dict(prevQualities)

            # Second, cull events dicts and update counts
            if updateCounts:
                updateCounts = False

                # Cull events
                refSeq = refSeqs[prevTargetName]
                culledEventsToWrite = cullEvents(eventsToWrite, refSeq)

                # Check and mark primers
                primer = None
                eventStart = max(culledEventsToWrite.keys())  # This includes both soft-clips and primers
                eventEnd = min(culledEventsToWrite.keys())
                if conf.randomPrimers:
                    primer = "NA"
                    # Mask off a designated absolute number of bases from the end for primers
                    for p in xrange(max(eventStart+1 - conf.randPrimerMaskLen, eventEnd), eventStart+1):
                        culledEventsToWrite[p] = "<"  # Mask events within primer region
                    sumCulledEventsPerGroup[prevTargetName]["NA"] += 1
                else:
                    # Use defined or found primer ranges to mask primer regions
                    for primerName in primerSeqs[prevTargetName]:
                        primerStart = primerSeqs[prevTargetName][primerName][1]-1  # 0-based; the supposed start of the candidate primer being checked
                        if eventStart in xrange(primerStart - conf.primerStartTolerance, primerStart + conf.primerStartTolerance+1):
                            # Appropriate primer assignment found
                            for p in xrange(primerSeqs[prevTargetName][primerName][0]-1, eventStart+1):  # 0-based
                                culledEventsToWrite[p] = "<"  # Mask events within primer region
                            primer = primerName  # Assign primer
                            sumCulledEventsPerGroup[prevTargetName][primerName] += 1
                            break
                    if primer is None:
                        sumNoPrimer += 1
                        # No primers assigned
                        if conf.includeMisprimed:
                            # Assign NA as primer group
                            primer = "NA"
                            sumCulledEventsPerGroup[prevTargetName]["NA"] += 1
                        else:
                            # Skip this read
                            continue

                # Generate RT excluding soft-clip and primer regions
                RTEvents = {ind: eve for ind, eve in culledEventsToWrite.items() if eve not in ["<", "s"]}
                RTLengthThis = len(RTEvents.keys())
                if RTLengthThis:  # Include this check to avoid error on empty sequence
                    RTend = min(RTEvents.keys())  # End of reverse transcription
                    RTstart = max(RTEvents.keys())  # Start of reverse transcription
                    # Take care of mid-sequence stops in long del stops mode
                    if RTEvents[RTend] == "!":
                        RTLengthThis -= 1
                else:
                    # Empty string (primer only)
                    sumNoExtension += 1
                    continue
                RTLength[prevTargetName][primer] += [RTLengthThis]

                # Cut off events beyond range for misprimed reads
                if eventsLimited and primer == "NA":
                    RTEvents = {inRange: RTEvents[inRange]
                                for inRange in RTEvents.keys() if inRange > RTstart - conf.eventsRange[1] and inRange <= RTstart - conf.eventsRange[0]}

                # Recording events into data structure
                for i in RTEvents.keys():
                    if (not eventsLimited or primer == "NA" or
                            (primer != "NA" and i >= eventsBrackets[prevTargetName][primer][0] and i < eventsBrackets[prevTargetName][primer][1])):
                        # Records if this event is within range unless we are not limiting range, or if it is in a random primer or misprimed read

                        c = RTEvents[i]  # The event char

                        # Check offset in index
                        if eventsLimited and primer != "NA":
                            shiftedInd = i - eventsBrackets[prevTargetName][primer][0]
                        else:
                            shiftedInd = i

                        # Coverage and read-through
                        if c in covered:  # covered cases
                            depth[prevTargetName][primer][shiftedInd] += 1
                            # Get lists of quals
                            qualsList[prevTargetName][primer][shiftedInd] += [ord(qualsToWrite[i])-33]
                        if c in readThrough:  # read-through frequency
                            RT[prevTargetName][primer][shiftedInd] += 1

                        # Specific events
                        if c in NUCS:  # mismatch
                            mutationCount[prevTargetName][primer][shiftedInd] += 1
                            mismatch[prevTargetName][primer][refSeq[i]][c][shiftedInd] += 1
                        elif c == "-":  # deletion
                            deletion[prevTargetName][primer][refSeq[i]][shiftedInd] += 1
                        elif c == "!":  # stop
                            stops[prevTargetName][primer][shiftedInd] += 1
                        elif c == "^":  # 5 prime of insert
                            insertion[prevTargetName][primer][shiftedInd] += 1

            # Record data for pairing
            prevClusterName = clusterName
            prevTargetName = targetName
            prevLine = line
            prevEvents = events
            prevRawRead = rawRead
            prevQualities = qualities

    # Build output
    # Prepare output
    headers = ",".join([
        "sample",
        "target",
        "primer",
        "nt",
        "avgQual",
        "ref",
        "toA",
        "toT",
        "toG",
        "toC",
        "stop",
        "mut",
        "in",
        "del",
        "depth",
        "RT",
        "iPtoA",
        "iPtoT",
        "iPtoG",
        "iPtoC",
        "iPstop",
        "iPmut",
        "iPin",
        "iPdel\n"
    ])
    # Determine behavior if output file already exists
    if os.path.isfile(conf.outputFileName):
        fileOut = open(conf.outputFileName, "r+U")  # Read + write, universal
        if conf.appendIfOutputFileExists and fileOut.readline().strip():
            # Append and non-empty
            sys.stderr.write("Output file already exists; appending current output to file.\n")
            fileOut.seek(0, 2)  # Go to the end to append
        else:
            if conf.appendIfOutputFileExists:
                # Append but empty
                sys.stderr.write("Output file already exists but is empty; overwriting file.\n")
            else:
                # Overwrite
                sys.stderr.write("Output file already exists; overwriting file.\n")
            fileOut.seek(0)  # Go to the beginning
            fileOut.truncate()  # Wipe file to prevent incomplete overwrite
            fileOut.write(headers)
    else:
        # Create new file if it doesn't already exist
        fileOut = open(conf.outputFileName, "w")
        fileOut.write(headers)

    # Prepare printing nt number shift if 1-based; placed here to avoid looped lookup of conf
    if conf.ntBase1:
        baseShift = 1
    else:
        baseShift = 0

    # Loop through rows
    for seqName in refSeqs:
        seqLen = len(refSeqs[seqName])
        for primerName in primerKeyList[seqName]:
            # Set primer names to print
            if primerName == "NA":
                printPrimerName = "NA"  # Misprimed or random primers
            elif conf.startSitePrimerName:
                printPrimerName = str(primerSeqs[seqName][primerName][0])  # Primer 5 primer nt number as names
            else:
                printPrimerName = primerName  # Defined primer names

            # Count and print over appropriate range
            if conf.randomPrimers or primerName == "NA" or not eventsLimited:
                printRange = xrange(seqLen)
            else:
                printRange = xrange(rangeLen)

            # Get nuc number shift due to truncated events
            if eventsLimited and primerName != "NA":
                eventsShift = eventsBrackets[seqName][primerName][0]
            else:
                eventsShift = 0

            # Count and print
            for i in printRange:
                # Out of range nucleotides (due to range limit shifts) omitted
                if i + eventsShift < 0 or i + eventsShift >= seqLen:
                    continue

                # Build total counts and depth
                stopCount = stops[seqName][primerName][i]
                inCount = insertion[seqName][primerName][i]
                delCount = sum([deletion[seqName][primerName][fromNuc][i] for fromNuc in NUCS])
                mutCount = mutationCount[seqName][primerName][i]
                nucDepth = depth[seqName][primerName][i]
                nucRT = RT[seqName][primerName][i]

                # Build by nucleotide mut counts
                toCounts = {"A": 0, "T": 0, "C": 0, "G": 0}
                for fromNuc in mismatch[seqName][primerName]:
                    for toN in mismatch[seqName][primerName][fromNuc]:
                        toCounts[toN] += mismatch[seqName][primerName][fromNuc][toN][i]
                if conf.sameNucAsNA:
                    toCounts[refSeqs[seqName][i + eventsShift]] = "NA"
                else:
                    toCounts[refSeqs[seqName][i + eventsShift]] = nucDepth - mutCount - delCount

                # Build Pstop
                Pstop = safeDivide(stopCount, nucRT + stopCount)

                # Get average quals per nuc
                qualListNuc = qualsList[seqName][primerName][i]
                avgQual = (float(sum(qualListNuc))/len(qualListNuc) if len(qualListNuc) > 0 else 0.0)

                # Produce csv output
                fileOut.write(",".join([
                    sampleName,
                    seqName,
                    printPrimerName,
                    str(i + baseShift + eventsShift),
                    "%.1f" % avgQual,
                    refSeqs[seqName][i + eventsShift],
                    str(toCounts["A"]),
                    str(toCounts["T"]),
                    str(toCounts["G"]),
                    str(toCounts["C"]),
                    str(stopCount),
                    str(mutCount),
                    str(inCount),
                    str(delCount),
                    str(nucDepth),
                    str(nucRT),
                    str(safeDivide(toCounts["A"], nucDepth)),
                    str(safeDivide(toCounts["T"], nucDepth)),
                    str(safeDivide(toCounts["G"], nucDepth)),
                    str(safeDivide(toCounts["C"], nucDepth)),
                    str(Pstop),
                    str(safeDivide(mutCount, nucDepth)),
                    str(safeDivide(inCount, nucDepth)),
                    str(safeDivide(delCount, nucDepth)) + "\n"
                ]))
    sys.stdout.write("Complete! Output at %s.\n\n" % conf.outputFileName)

    # Print summary info
    sys.stdout.write(("Summary of run: \n" +
                      "\tReads parsed: %i\n" +
                      "\tReads skipped due to low mapping qual: %i\n" +
                      "\tReads discarded: %i\n" +
                      "\tReads of primer only: %i\n" +
                      "\tUnpaired reads: %i\n") % (sumSamLines, sumSamLowQual, sumSamDiscard, sumNoExtension, sumUnpairedReads))
    if not conf.randomPrimers:
        if conf.includeMisprimed:
            sys.stdout.write("\tEvents analyzed without primer assignment: %i\n" % sumNoPrimer)
        else:
            sys.stdout.write("\tEvents skipped due to no primer assignment: %i\n" % sumNoPrimer)
    sys.stdout.write("Events analyzed counts:\n")
    if conf.randomPrimers:
        for seqName in refSeqs:
            sys.stdout.write("\t%s: %i\n" % (seqName, sumCulledEventsPerGroup[seqName]["NA"]))
    else:
        for seqName in refSeqs:
            for primerName in primerKeyList[seqName]:
                sys.stdout.write("\t%s[%s]: %i\n" % (seqName, (primerName if primerName != "NA" and not conf.randomPrimers else "No primer"), sumCulledEventsPerGroup[seqName][primerName]))
    sys.stdout.write("Average RT lengths:\n")  # Does not include no extension reads
    if conf.randomPrimers:
        for seqName in refSeqs:
            RTLengthList = RTLength[seqName]["NA"]
            avgRTLength = (float(sum(RTLengthList))/len(RTLengthList) if len(RTLengthList) > 0 else 0.0)
            sys.stdout.write("\t%s: %.2f\n" % (seqName, avgRTLength))
    else:
        for seqName in refSeqs:
            for primerName in primerKeyList[seqName]:
                RTLengthList = RTLength[seqName][primerName]
                avgRTLength = (float(sum(RTLengthList))/len(RTLengthList) if len(RTLengthList) > 0 else 0.0)
                sys.stdout.write("\t%s[%s]: %.2f\n" % (seqName, (primerName if primerName != "NA" and not conf.randomPrimers else "No primer"), avgRTLength))
except SystemExit:
    # sys.exit() calls
    pass
except:
    # Print out error
    sys.stderr.write(traceback.format_exc())
    sys.exit(1)
finally:
    if conf.saveLog:
        # Close log file
        logFile.close()
    if fileOut:
        # Close output if opened
        fileOut.close()
    sys.exit(0)
