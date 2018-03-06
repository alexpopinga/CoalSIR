#!/usr/bin/env python

from sys import argv, exit, stdout
import os, os.path
from argparse import ArgumentParser, FileType
from math import log10
import json

def getLeafAndOriginTimes(j, trajNum):
    """Extract leaf times from sample times in trajectory."""

    traj = j['trajectories'][trajNum]
    n = len(traj['t'])

    isSample = [(traj['Rh'][i]-traj['Rh'][i-1])>0.0 for i in range(1,n)]
    timesForward = [traj['t'][i] for i in range(1,n) if isSample[i-1]]

    origin = max(timesForward)
    timesBackward = [origin-t for t in timesForward]

    return (timesBackward,origin)


def generateXML(j, trajNum, trajFileName, treefname, beta, minus, out=stdout):
    """Assembles BEAST 2 XML for simulating a single
    stochastic coalescent tree corresponding to the
    chosen trajectory."""

    leafTimes, origin = getLeafAndOriginTimes(j, trajNum)

    # Write header
    out.write("""
<beast version ='2.0' namespace ='
\tbeast.evolution.tree
\t:beast.evolution.alignment
\t:beast.phylodynamics.epidemiology
\t:beast.core.parameter
\t:master.utilities'>\n\n""")

    # Write taxon set
    out.write("\t<taxonset id='taxa' spec='TaxonSet'>\n")
    for i in range(len(leafTimes)):
        out.write("\t\t<taxon id='t{}' spec='Taxon'/>\n".format(i+1))
    out.write('\t</taxonset>\n\n')

    # Write tip ages
    out.write("\t<trait id='tipDates' spec='TraitSet' traitname='date-backward'\n\t\tvalue='")
    for i,t in enumerate(leafTimes):
        if (i>0):
            out.write(",\n\t\t")
        out.write("t{}={}".format(i+1, t))

    out.write("'>\n\t\t<taxa idref='taxa'/>\n\t</trait>\n\n")

    # Define population size expression
    if minus:
        expression='(I-1)/(2*{}*S)'.format(beta)
    else:
        expression='I/(2*{}*S)'.format(beta)

    # Specify population function 
    out.write(
        """\t<populationFunction spec='PopulationFunctionFromJSON' id='popFun'
\t\tfileName='../{}'
\t\tpopSizeExpression='{}'
\t\torigin='{}'
\t\ttrajNum='{}' />\n\n""".format(trajFileName, expression, origin, trajNum))

    # Specify simulation
    out.write("""\t<run spec='CoalescentSimulatorBasic'
\t\treplicates='1'
\t\toutputFileName='{}'
\t\ttrait='@tipDates'
\t\tmaxHeight='{}'
\t\tpopulationFunction="@popFun" />\n""".format(treefname, origin));

    # Close beast tag
    out.write("</beast>\n")


#### MAIN ####

parser = ArgumentParser("genSCTreeXMLs",
                        description="""Generates XMLs necessary
for simulating stochastic coalescent trees.""",
                        epilog="Origin times are obtained from trajectory file.")

parser.add_argument("jsonfile", metavar="trajectory_file.json", type=FileType('r'),
                    help="MASTER-generated trajectory file to use.")

parser.add_argument("trajnumfile", metavar="traj_numbers.txt", type=FileType('r'),
                    help="File containing trajectory numbers to use.")

parser.add_argument("dirname", type=str,
                    help="Name of new dir for generated XML files.")

parser.add_argument("--beta","-b", type=float, default=0.00075,
                    help="Infection rate. (Default 0.00075.)")

parser.add_argument("--gamma","-g", type=float, default=0.3,
                    help="Removal rate. (Default 0.3.)")

parser.add_argument("--minus","-m", action="store_true",
                    help="Use (I-1)/(2*beta*S) instead of I/(2*beta*S).")

args = parser.parse_args(argv[1:])

# Extract trajectory numbers
trajNums = [int(x.strip()) for x in args.trajnumfile.readlines()]
print "{} trajectory numbers found.".format(len(trajNums))

# Load JSON file
j = json.load(args.jsonfile)
print "JSON file '{}' loaded.".format(args.jsonfile.name)


# Create output directory
if os.path.exists(args.dirname):
    print "Cannot create directory '{}': file already exists.".format(args.dirname)
    exit(1)
os.mkdir(args.dirname)

# Generate XML files
i=1
ndigits = int(log10(len(trajNums)))+1
outfname = "{}/treeSim_{:0" + str(ndigits) + "d}.xml"
treefname = "simulated_{:0" + str(ndigits) + "d}.tree"
for trajNum in trajNums:
    with open(outfname.format(args.dirname, i),"w") as outFile:
        stdout.write("\rGenerating file {} of {}...".format( i, len(trajNums)))
        generateXML(j, trajNum, args.jsonfile.name, treefname.format(i), args.beta, args.minus, out=outFile)
    i += 1
print "\nDone."
