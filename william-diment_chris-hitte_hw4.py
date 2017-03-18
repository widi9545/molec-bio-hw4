##########################################################################################
#	File: william-diment_chris-hitte_hw4.py
#	
#	Purpose: The purpose of the program is to output a phylogenetic tree in Newick format, using a distance matrix that we calculate
#	from a given set of sequences. We then use a UPGMA clustering algorithm to build the tree, and parse a Newick output from that tree.
#
#	Developers: Chris Hitte, William Diment
#	CSCI 4314, Spring 2017
#	Homework 4
#
##########################################################################################
#
#	Sample command line to run the program
#
#	python william-diment_chris-hitte_hw4.py -f dataset1.fasta -t dataset1output -s scoring1.txt -g -1
#
#	Usage: -F FASTA File -T Newick Output File -S Scoring Matrix -G gap penalty
#	-F File, -f File, use this specify the fasta sequence file that you want to use
#	-S ScoringMatrix, -s, use this to specify the scoring matrix that you wish to use
#	-T NewickOutputFile, -t, specify the file in which you want to output the Newick tree
#	-G Gap, -g, specify the gap penalty that you wish to apply
#
#########################################################################################
#
#	References: Okeson Jan 2016
#               alexokeson_hw1.py
#               Formatting of the header comment and function comments
#
#
#########################################################################################


import argparse
import numpy

#########################################################################################
#
#	All global variables and datas structures are initialized here
#
#########################################################################################
FASTA_FILE = ""
SCORING_FILE = ""
GAP_PENALTY = -1 
NEWICK_TREE_FILE = "" 

SEQUENCE_LIST = []
SEQUENCE_INDEX = {}

DISTANCE_MATRIX = None

SCORING_MATRIX = numpy.zeros(shape=(4, 4)) #We initialize our scoring matrix
SCORING_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3} #We use the scoring index to easily index into the scoring matrix

TREE_NODE_INDEX = {} #This dictionary will contain all the indicies for each node and the clusters associated with that node.
TREE_ROOT = None #This will represent the root node of the tree - the common ancestor of all the sequences

###################################################################################################
#
#	We define our necessary classes here. 
#	We store our sequences in the Sequence class, and our global alignment results
#	in the appropriately named classes.
#	We define two classes for use in the tree - the TreeNode class for each individual node
#	The TreeCluster class is used to keep track of clusters of nodes
#
###################################################################################################

class Sequence:
    def __init__(self, id, name):
        self.id = id
        self.name = name.replace(' ', '_').replace('\t', '_')
        self.seq_data = ""


class GlobalAlignmentResult:
    def __init__(self, seq_1, seq_2, distance, length, align_1, align_2):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        self.distance = distance
        self.length = length
        self.align_1 = align_1
        self.align_2 = align_2


class TreeCluster:
    def __init__(self, score):
        self.score = score
        self.left = None
        self.right = None


class TreeNode:
    def __init__(self, label, score):
        self.label = label
        self.score = score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-F', '-f', '--fasta', type=file, help="fasta filename", required=True)
    parser.add_argument('-T', '-t', '--tree', type=str, help="destination for the generated Newick tree", required=True)
    parser.add_argument('-S', '-s', '--score', type=file, help="scoring matrix filename", required=True)
    parser.add_argument('-G', '-g', '--gap', type=int, help="gap penalty", required=True)

    args = parser.parse_args()

    global FASTA_FILE
    FASTA_FILE = args.fasta

    global SCORING_FILE
    SCORING_FILE = args.score

    global GAP_PENALTY
    GAP_PENALTY = args.gap

    global NEWICK_TREE_FILE
    NEWICK_TREE_FILE = args.tree

    #
    parse_score_matrix()

    # We parse the given FASTA file
    parse_fasta_file()

    # We calculate our global alignments
    compute_global_distances()

    # We print the distance matrix!
    print_distance_matrix()

    # We use the UPGMA clustering algorithm here
    compute_upgma()

    #We output the newick tree for later use in visualization
    print_newick_tree()

    print "exiting..."

####################################################################################
#
#	We parse a given scoring matrix for use in computing our global alignments
#
####################################################################################
def parse_score_matrix():
    print "parsing scoring matrix..."

    lines = SCORING_FILE.readlines()

    i = 0
    j = 0
    for line in lines:
        scores = line.split()

        for score in scores[1:]:
            SCORING_MATRIX[i][j] = int(score)
            j += 1

        i += 1
        j = 0

    print "done"
    print


######################################################################################
#
#	We parse a given FASTA file for our use.
#	After the FASTA file has been parsed, we initialize a global matrix
#	that we will use in calculating the distances from our global alignment algorithm
#
######################################################################################
def parse_fasta_file():
    print "parsing fasta file..."

    read_sequences()

    global DISTANCE_MATRIX
    DISTANCE_MATRIX = numpy.zeros((len(SEQUENCE_LIST), len(SEQUENCE_LIST)))

    print "done"
    print

#########################################################################################
#
#	We read and store each individual sequence from the FASTA file that we want to use.
#	We update SEQUENCE_LIST with each new sequence, where we store the given ID for the sequence
#	and the bases associated with that sequence ID
#	
##########################################################################################
def read_sequences():
    id = 0

    line = FASTA_FILE.readline()
    while line != '' and line[0] == '>':
        new_sequence = Sequence(id, line[1:-1])
        SEQUENCE_LIST.append(new_sequence)

        seq_data = ""
        line = FASTA_FILE.readline()
        while line != '' and line[0] != '>':
            seq_data += line[:-1]
            line = FASTA_FILE.readline()

        #We store the data from the sequence in the 'data' field of the sequence data structure
        new_sequence.seq_data = seq_data

        id += 1

########################################################################################
#
#	We compute our global distance matrix here. 
#	We use two functions - compute_global_distance_matrix to initialize our distance matrix according to the
#	input sequences. We then use the input sequences and the initialized distance matrix to calculate
#	the global alignment
#
#######################################################################################

def compute_global_distances():
    print "calculating global distance matrix..."

    i = 0
    j = 0
    while i < len(SEQUENCE_LIST):
        seq_1 = SEQUENCE_LIST[i]

        while j < len(SEQUENCE_LIST):
            if j != i:
                seq_2 = SEQUENCE_LIST[j]

                matrix = compute_global_distance_matrix(seq_1, seq_2)

                calculate_global_distance(seq_1, seq_2, matrix)

            j += 1

        i += 1
        j = i

    print "done"
    print
########################################################################################
#
#
#	Algorithm to initialize a distance matrix. 
#	Input: Two sequences that we wish to compare
#	
#	We initialize an empty m x n matrix according to the length of the two sequences If there is a match or mismatch in the sequences, 
#	we index into our scoring matrix and input the appropriate match/mismatch score. if there is a gap, we input the gap penalty that was specified on the command line.
#	We then choose the max from these options to maximize the distance between the two sequences.
#
#			
########################################################################################
def compute_global_distance_matrix(seq_1, seq_2):
    seq_data_1 = ' ' + seq_1.seq_data
    seq_data_2 = ' ' + seq_2.seq_data

    matrix = numpy.zeros((len(seq_data_2), len(seq_data_1)))

    i = 1
    j = 1
    while i < len(seq_data_2):
        while j < len(seq_data_1):
            diagonal = matrix[i - 1][j - 1] + SCORING_MATRIX[SCORING_INDEX[seq_data_1[j]]][SCORING_INDEX[seq_data_2[i]]] #We index into the scoring matrix and update the match/mismatch accordingly

            top = matrix[i - 1][j] + GAP_PENALTY	#If there is a gap, we input the gap penalty as specified on the command line.

            left = matrix[i][j - 1] + GAP_PENALTY

            matrix[i][j] = max(top, left, diagonal, 0)

            j += 1

        i += 1
        j = 1

    return matrix

#####################################################################################
#
#
#	Algorithm to calculate a global alignment.
#	Input: Two sequences we wish to align, a given distance matrix.
#	
#	We traceback through the matrix that we initialized in order to calculate our global alignment
#	If the diagonal in the matrix is greater than or equal to the value of the top or left sides, we choose the diagonal
#	as this represents a match or mismatch between the two sequences that we are comparing.
#	if the left side or top side is greater than the diagonal, this represents a gap in the two sequences that we must fill. 
#	If the left side is greater than the top side, we choose the left side, otherwise we choose the topside.
#	We prefer to use the diagonal as much as we can, as we want to minimize gaps in the alignment between the two sequences
#	
#
#####################################################################################
def calculate_global_distance(seq_1, seq_2, matrix):
    seq_data_1 = seq_1.seq_data
    seq_data_2 = seq_2.seq_data

    align_1 = ""
    align_2 = ""

    i = len(seq_data_2)
    j = len(seq_data_1)
    while i > 0 and j > 0:
        u = matrix[i - 1][j]
        l = matrix[i][j - 1]
        d = matrix[i - 1][j - 1]

        if d >= l and d >= u: #Previously this set of variables was flipped, we change it here to maximize the distance matrix
            align_1 = seq_data_1[j-1] + align_1
            align_2 = seq_data_2[i-1] + align_2
            i -= 1
            j -= 1

        elif l >= u:
            align_1 = seq_data_1[j-1] + align_1
            align_2 = '-' + align_2
            j -= 1

        else:
            align_1 = '-' + align_1
            align_2 = seq_data_2[i-1] + align_2
            i -= 1

    while j >= 0:
        align_1 = seq_data_1[j] + align_1
        j -= 1

    while i >= 0:
        align_2 = seq_data_2[i] + align_2
        i -= 1
   
    if len(align_1) > len(align_2):  #if the length of the alignments are mismatched, we input the appropriate gaps here
        align_2 = ((len(align_1) - len(align_2)) * '-') + align_2
    elif len(align_2) > len(align_1):
        align_1 = ((len(align_2) - len(align_1)) * '-') + align_1

    distance = 0
    length = len(align_1) 
    
    i = 0
    while i < length:
        if align_1[i] != align_2[i]: #if the alignments do not match this represents the distance between the two sequences
            distance += 1 #as a result we increase the distance here

        i += 1

        #We input the distance associated with the two sequences into the global distance matrix.
    DISTANCE_MATRIX[seq_1.id][seq_2.id] = round(float(distance) / length, 5)
    DISTANCE_MATRIX[seq_2.id][seq_1.id] = round(float(distance) / length, 5)
#############################################################################
#
#	We iterate through each row of the distance matrix in order to print the distance matrix
#
#############################################################################

def print_distance_matrix():
    print "Distance Matrix for all Sequences:"
    print

    line = ""
    for row in DISTANCE_MATRIX:
        for column in row:
            line += "{0:.5f} ".format(round(column, 5))

        print line
        line = ""

    print
###########################################################################
#
#	Algorithm to implement UPGMA clustering.
#	Input: None
#	
#	For each sequence that we have compared, we append it to an initial cluster list. 
#	We then sequentially find clusters i,j to append to the tree, checking to see if it is initialized,
#	and updating the distance if it is a root node. We then initialize an instance of the TreeCluster class, which has an effect
#	of representing the root node (ancestor) of two sequences, and append the appropriate clusters to that root node. We repeat this process
#	until we have processed every cluster.
#	
############################################################################

def compute_upgma():
    print "computing UPGMA for global distance matrix..."

    node = None

    cluster_list = []
    for seq in SEQUENCE_LIST:
        SEQUENCE_INDEX[seq.name] = seq.id
        cluster_list.append([seq.name])

    while len(cluster_list) > 1:
        next_c_i, next_c_j, score = find_next_cluster(cluster_list) #We find the next clusters and the distance associated with them

        if str(cluster_list[next_c_i]) in TREE_NODE_INDEX:  
            node_1 = TREE_NODE_INDEX[str(cluster_list[next_c_i])] #If the node is in our treeIndex, we take the values from the tree index
        else:
            node_1 = TreeNode(str(cluster_list[next_c_i][0]), score / 2) #else we initialize the node

        if str(cluster_list[next_c_j]) in TREE_NODE_INDEX:
            node_2 = TREE_NODE_INDEX[str(cluster_list[next_c_j])]
        else:
            node_2 = TreeNode(str(cluster_list[next_c_j][0]), score / 2)

        if isinstance(node_1, TreeCluster): #We set the score of each node here
            node_1.score = score / 2 - node_1.score

        if isinstance(node_2, TreeCluster):
            node_2.score = score / 2 - node_2.score

        node = TreeCluster(score / 2) #We set the root node/its distance and assign it left and right children
        node.left = node_1
        node.right = node_2

        TREE_NODE_INDEX[str(cluster_list[next_c_i] + cluster_list[next_c_j])] = node #We update the treeNode index with the root node and the clusters associated with the root node

        new_cluster_list_item = cluster_list[next_c_i] + cluster_list[next_c_j]

        cluster_list[next_c_i] = new_cluster_list_item
        del cluster_list[next_c_j] #We remove the old cluster from consideration

    global TREE_ROOT
    TREE_ROOT = node #This represents the farthest distance in the tree - IE the common ancestor of each sequence we have compared.

    print "done"
    print

####################################################################################
#
#	Algorithm to find the next cluster for our hierarchy.
#	Input: a given cluster list
#	
#	We compute the distance for the 'compute_cluster_distance' algorithm and compare it to the other distances that 
#	already exist in the tree. If the cluster is not initialized; we set it as the minimum value. 
#	If the tree is initialized, we compare the distance to the minimum value. If it's smaller than the minimum value,
#	this represents the next cluster, and so we return the next two clusters along with the distance score.
#	Otherwise, we continue on until we find the next cluster
#	
#
#####################################################################################

def find_next_cluster(cluster_list):
    initialized = False
    min_score = -1

    i = 0
    next_c_i = i

    j = 1
    next_c_j = j

    while i < len(cluster_list):
        while j < len(cluster_list):
        	
            score = compute_cluster_distance(cluster_list[i], cluster_list[j])

            if not initialized: #we check to see if its initialized - if it isn't, we initialize it and set the minimum score
                min_score = score
                initialized = True
            elif score < min_score: #We compare the distance score - if its smaller then the min score we choose these clusters as our next two clusters
                next_c_i = i
                next_c_j = j
                min_score = score #we update the minimum distance score with the new, lower score

            j += 1

        i += 1
        j = i + 1

    return next_c_i, next_c_j, min_score #we return the next two clusters and the distance value associated with them

######################################################################
#
#	Algorithm to compute the cluster distance.
#	Input: Two given clusters or collection of clusters
#	We compare the distance of the two given clusters or collection of clusters to the distance matrix
#	and calculate the score by updating their distance with the value from the given distance matrix, which we
#	then divide by the length (magnitude) of the clusters or cluster list that we are comparing.
#	We then return the distance (score) for use in choosing the next cluster to fit within our tree.
#
######################################################################

def compute_cluster_distance(c_i, c_j):
    score = 0.0

    for x in c_i:
        for y in c_j:
            score += DISTANCE_MATRIX[SEQUENCE_INDEX[x]][SEQUENCE_INDEX[y]] #We update the distance with the values from the distance matrix

    score /= (len(c_i) * len(c_j)) #This is actually where we divide by the magnitude of the clusters we are comparing

    return score

#########################################################################
#
#	We write the Newick representation of the tree to the file that we specified on the command line
#	
##########################################################################

def print_newick_tree():
    print "Newick Tree:"
    print

    tree = calculate_newick_tree(TREE_ROOT)

    out = file(NEWICK_TREE_FILE, 'w')
    out.write(tree)

    print tree
    print
#########################################################
#
#	We calculate the Newick representation of the hierarchical tree that we have built
#
########################################################

def calculate_newick_tree(node):
    if isinstance(node, TreeNode):
        return "{0}:{1}".format(node.label, "{0:.5f}".format(round(node.score, 5)))

    elif node != TREE_ROOT:
        return "({0}, {1}):{2}".format(calculate_newick_tree(node.left),
                                       calculate_newick_tree(node.right),
                                       "{0:.5f}".format(round(node.score, 5)))
    else:
        return "({0}, {1});".format(calculate_newick_tree(node.left), calculate_newick_tree(node.right))


if __name__ == "__main__":
    main()
