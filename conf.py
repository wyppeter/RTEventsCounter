# Configuration file for RTEventsCounter.py
# Copyright 2017 Peter Y. Wang
# Simon Lab, Yale University, New Haven, CT
# File name: conf.py
# Version: 1.3
# Python version: 2.7
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


# >> Files controls

# Set output file name (.csv)
outputFileName = "output.csv"

# Overwrite or append when output file already exists
appendIfOutputFileExists = True

# Save log to "counterLog.txt"
saveLog = True


# >> Output format controls

# Nucleotide numbering base in output (0-based or 1-based)
ntBase1 = True

# Output same-nucleotide toN as actual count or "NA"
sameNucAsNA = True


# >> Special event controls

# Long deletions (>1nts)
# "dels"/"stops"/"none"
longDelsAs = "dels"

# Include insertions counts
countInserts = True


# >> Data selection and Paired-reads

# Paired reads behavior
alignPaired = True
ignoreUnpaired = True  # Useful only when alignPaired is True

# Choose to only use events within a certain range (For non-random primer mode only)
# 0 = first nucleotide of start of aligned RT sequence (3 prime side in ref seq direction), excluding primer nucs
# [inclusive-start, exclusive-end]; in reverse to event strings (3 prime to 5 prime)
# e.g. [3, 9] => ...~~~~876543~~~<<<<<~~~~...
# Both must be non-negative integers
# Should be shorter than the total length of the ref seq
# [0, 0] to toggle off
eventsRange = [0, 1000]

# Mapping qual threshold to include a SAM line
minMapQual = 0


# >> Primers settings

# Set primer input mode (For non-random primer mode only)
# True = sequences / False = coordinates
# In sequences mode, coordinates are first auto-generated through sequence searching
primerInputModeSeqs = True

# Toggle whether to use 5 prime nt number as primer name automatically (For non-random primer mode only)
startSitePrimerName = True

# Include misprimed reads as primer "NA" (For non-random primer mode only)
includeMisprimed = True

# Start site tolerance +/- for primer designation (For non-random primer mode only)
# e.g. 3 = ...NNNNYYYSYYYNNNN...
primerStartTolerance = 5


# >> Random primer settings

# Primer mode
# Toggle random primer mode
# Note that arguments to run are different in random primer mode
randomPrimers = False

# Absolute number of bases to remove from the end for primers (For random primer mode only)
randPrimerMaskLen = 8
