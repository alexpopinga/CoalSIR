#!/bin/bash
#
# Filter nexus and newick tree files produced by SIR_1000sims.xml to
# produce new files containing only those trees with >= 100 leaves.
# Note that phylostat (from github.com/tgvaughan/PhyloPaint) is used
# for the leaf counting.

if ( ! [ -e SIR_tree_R0equals1.newick ] ) || ( ! [ -e SIR_tree_R0equals1.nexus ] ) || ( ! [ -e SIR_tree_R0equals1.json ] ); then
    echo SIR_tree_R0equals1.xml output files could not be found.
    exit 1
fi

# Count leaves
echo -n "Counting leaves..."
phylostat -n SIR_tree_R0equals1.newick nleaves > nleaves_R0equals1.txt
echo "done."

# Filter newick
echo -n "Filtering newick..."
paste nleaves_R0equals1.txt SIR_tree_R0equals1.newick \
    | awk '{if ($1>=100) { $1=""; print $0 }}' \
    | sed 's/^ *//' > SIR_tree_R0equals1_filtered.newick
echo "done."

# Filter nexus
echo -n "Filtering nexus..."
echo -e "#nexus\nbegin trees;" > SIR_tree_R0equals1_filtered.nexus
grep TREE SIR_tree_R0equals1.nexus | paste nleaves_R0equals1.txt - \
    | awk '{if ($1>=100) { $1=""; print $0 }}' \
    | sed 's/^ *//' >> SIR_tree_R0equals1_filtered.nexus
echo "end;" >> SIR_tree_R0equals1_filtered.nexus
echo "done."

# Extract filtered trajectory numbers
echo -n "Extracting filtered trajectory numbers..."
grep TREE SIR_tree_R0equals1_filtered.nexus | cut -d= -f1 | cut -d_ -f2 | sed 's/ *$//' > trajNumbersFiltered_R0equals1.txt
echo "done."

# Write file containing useful statistics
echo -n "Writing statistics file..."
ntrees=`wc -l < SIR_tree_R0equals1_filtered.newick`
phylostat SIR_tree_R0equals1_filtered.newick nleaves height origin >statsFiltered_R0equals1.txt
( echo trajNum; cat trajNumbersFiltered_R0equals1.txt ) | paste -d' ' - statsFiltered_R0equals1.txt > temp.txt
( echo treeNum; for ((i=1; i<=$ntrees; i++)); do echo $i; done ) \
    | paste -d' ' - temp.txt > statsFiltered_R0equals1.txt
rm temp.txt
echo "done."