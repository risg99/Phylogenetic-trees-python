# Importing necessary libraries
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

# Read the sequences and align
align = AlignIO.read('msa.phy','phylip')

# Calculate the distance matrix
calculator = DistanceCalculator('identity')
distMatrix = calculator.get_distance(align)

# Creating a DistanceTreeConstructor object
constructor = DistanceTreeConstructor()

####################################################
					#  UPGMA
####################################################

# Construct the phlyogenetic tree using UPGMA algorithm
UPGMATree = constructor.upgma(distMatrix)

# Draw the phlyogenetic tree
Phylo.draw(UPGMATree)

# Printing the phlyogenetic tree using terminal
Phylo.draw_ascii(UPGMATree)

####################################################
					#  NJ
####################################################

# Construct the phlyogenetic tree using NJ algorithm
NJTree = constructor.nj(distMatrix)

# Draw the phlyogenetic tree
Phylo.draw(NJTree)

# Printing the phlyogenetic tree using terminal
Phylo.draw_ascii(NJTree)


####################################################
			#  Terminal Output
####################################################
'''
                                    _______________________________ sumatran
  _________________________________|
 |                                 |_______________________________ orangutan
_|
 |            _____________________________________________________ gorilla
 |           |
 |___________|                                  ___________________ bonobo
             |        _________________________|
             |_______|                         |___________________ chimpanzee
                     |
                     |_____________________________________________ homo_sapie

                                         __________________________ sumatran
  ______________________________________|
 |                                      |________________________ orangutan
 |
_|____________________________________________ gorilla
 |
 |     _____________________________________ homo_sapie
 |____|
      |                     ______________ bonobo
      |____________________|
                           |________________ chimpanzee

[Finished in 30.0s]
'''