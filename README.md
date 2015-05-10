smilefinder
===========


SmileFinder is a python script built to detect signatures of positive selection in diploide genomes.

Before deciding to use SmileFinder, we recommend to read the note: 

SmileFinder: a resampling-based Platform to evaluate signatures of selection from genome-wide sets of matching allele frequency data between populations, Guiblet W.M. et al..


Requirements:

- Python2.x
- Numpy
- Scipy


Originally, we used SmileFinder with data from the Human Genome Diversity Project. A first script, name count.py, gets data from the HGDP_FinalReport_Forward.txt to calculate the values of Heterozygosity and Fst at each locus.

To run count.py, you will need:

- HGDP_FinalReport_Forward.txt
- hgdp-popfile.txt
- HGDP_Map.txt

HGDP_FinalReport_Forward.txt and HGDP_Map.txt are available on http://www.hagsc.org/hgdp/files.html.

hgdp-popfile.txt needs to be customized according to your experiment. In this file you list the populations and individuals used in your calculations. We recommand at least 15 individuals (30 genomes) for each population.

Using count.py:

python count.py count.py HGDP_FinalReport_Forward.txt hgdp-popfile.txt HGDP_Map.txt  > output.csv

In the output, you will find in this order: Name of the marker, Chromosome, Position, Major Allele in the first Population, Frequency of the Major Allele, Number of Individuals in the first population, Expected Heterozygosity in the first population, Observed Heterozygosity in the first population,  Major Allele in the second Population, Frequency of the Major Allele, Number of Individuals in the second population, Expected Heterozygosity in the second population, Observed Heterozygosity in the second population, Fst between the two populations


This output is used as input for smilefinder.py. Be careful, the markers must be sorted in order according to their chromosome and position. If you wish to create your own input without using count.py, the only fields required are:

- Name of the marker
- Chromosome
- Position
- Observed Heterozygosity in both populations
- Fst

The other fields are not required and can be let empty.

When running smilefinder.py, you will be asked to provide:


- The name of the input file.
- Chose a name for the report file
- How many sliding windows you wish to have (incremented by 2 positions)
- Size of the shortest window (odd values)
- How much resampling (size of the radomized distribution - up to 100,000,000)
- Sensitivity of the report ( < 1 and > 0)


SmileFinder will provide two outputs. First, the SmileFinderCompleteTable: this is the most important table, containing all the percentiles calculated at each marker position. The second output is the report file. It filters the CompleteTable according to the Sensitivity, returning the markers with the oddest values. We recommand the users to filter the CompleteTable according to their needs (example: graphing the percentiles values in a region of interest for vizualizing the behavior of the genetic diversity).



We hope SmileFinder will be able to provide researches powerful insights about natural selection. Any questions about the program can be addressed to wilfried.guiblet@upr.edu 
