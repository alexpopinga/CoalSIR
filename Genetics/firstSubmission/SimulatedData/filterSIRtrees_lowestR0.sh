#!/bin/bash
#
# Filter nexus and newick tree files produced by SIR_1000sims.xml to
# produce new files containing only those trees with >= 100 leaves.
# Note that phylostat (from github.com/tgvaughan/PhyloPaint) is used
# for the leaf counting.

if ( ! [ -e SIR_tree_lowestR0.newick ] ) || ( ! [ -e SIR_tree_lowestR0.nexus ] ) || ( ! [ -e SIR_tree_lowestR0.json ] ); then
    echo SIR_tree_lowestR0.xml output files could not be found.
    exit 1
fi

# Count leaves
echo -n "Counting leaves..."
phylostat -n SIR_tree_lowestR0.newick nleaves > nleaves_lowestR0.txt
echo "done."

# Filter newick
echo -n "Filtering newick..."
paste nleaves_lowestR0.txt SIR_tree_lowestR0.newick \
    | awk '{if ($1>=100) { $1=""; print $0 }}' \
    | sed 's/^ *//' > SIR_tree_lowestR0_filtered.newick
echo "done."

# Filter nexus
echo -n "Filtering nexus..."
echo -e "#nexus\nbegin trees;" > SIR_tree_lowestR0_filtered.nexus
grep TREE SIR_tree_lowestR0.nexus | paste nleaves_lowestR0.txt - \
    | awk '{if ($1>=100) { $1=""; print $0 }}' \
    | sed 's/^ *//' >> SIR_tree_lowestR0_filtered.nexus
echo "end;" >> SIR_tree_lowestR0_filtered.nexus
echo "done."

# Extract filtered trajectory numbers
echo -n "Extracting filtered trajectory numbers..."
grep TREE SIR_tree_lowestR0_filtered.nexus | cut -d= -f1 | cut -d_ -f2 | sed 's/ *$//' > trajNumbersFiltered_lowestR0.txt
echo "done."

# Write file containing useful statistics
echo -n "Writing statistics file..."
ntrees=`wc -l < SIR_tree_lowestR0_filtered.newick`
phylostat SIR_tree_lowestR0_filtered.newick nleaves height origin >statsFiltered_lowestR0.txt
( echo trajNum; cat trajNumbersFiltered_lowestR0.txt ) | paste -d' ' - statsFiltered_lowestR0.txt > temp.txt
( echo treeNum; for ((i=1; i<=$ntrees; i++)); do echo $i; done ) \
    | paste -d' ' - temp.txt > statsFiltered_lowestR0.txt
rm temp.txt
echo "done."