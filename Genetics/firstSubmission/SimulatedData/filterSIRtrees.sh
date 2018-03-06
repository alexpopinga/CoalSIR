#!/bin/bash
#
# Filter nexus and newick tree files produced by SIR_1000sims.xml to
# produce new files containing only those trees with >= 100 leaves.
# Note that phylostat (from github.com/tgvaughan/PhyloPaint) is used
# for the leaf counting.

if ( ! [ -e SIR_1000sims.newick ] ) || ( ! [ -e SIR_1000sims.nexus ] ) || ( ! [ -e SIR_1000sims.json ] ); then
    echo SIR_1000sims.xml output files could not be found.
    exit 1
fi

# Count leaves
echo -n "Counting leaves..."
phylostat -n SIR_1000sims.newick nleaves > nleaves.txt
echo "done."

# Filter newick
echo -n "Filtering newick..."
paste nleaves.txt SIR_1000sims.newick \
    | awk '{if ($1>=100) { $1=""; print $0 }}' \
    | sed 's/^ *//' > SIR_1000sims_filtered.newick
echo "done."

# Filter nexus
echo -n "Filtering nexus..."
echo -e "#nexus\nbegin trees;" > SIR_1000sims_filtered.nexus
grep TREE SIR_1000sims.nexus | paste nleaves.txt - \
    | awk '{if ($1>=100) { $1=""; print $0 }}' \
    | sed 's/^ *//' >> SIR_1000sims_filtered.nexus
echo "end;" >> SIR_1000sims_filtered.nexus
echo "done."

# Extract filtered trajectory numbers
echo -n "Extracting filtered trajectory numbers..."
grep TREE SIR_1000sims_filtered.nexus | cut -d= -f1 | cut -d_ -f2 | sed 's/ *$//' > trajNumbersFiltered.txt
echo "done."

# Write file containing useful statistics
echo -n "Writing statistics file..."
ntrees=`wc -l < SIR_1000sims_filtered.newick`
phylostat SIR_1000sims_filtered.newick nleaves height origin >statsFiltered.txt
( echo trajNum; cat trajNumbersFiltered.txt ) | paste -d' ' - statsFiltered.txt > temp.txt
( echo treeNum; for ((i=1; i<=$ntrees; i++)); do echo $i; done ) \
    | paste -d' ' - temp.txt > statsFiltered.txt
rm temp.txt
echo "done."