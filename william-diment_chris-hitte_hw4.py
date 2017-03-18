import argparse
import numpy

FASTA_FILE = ""
SCORING_FILE = ""
GAP_PENALTY = -1
NEWICK_TREE_FILE = ""

SEQUENCE_LIST = []
SEQUENCE_INDEX = {}

DISTANCE_MATRIX = None

SCORING_MATRIX = numpy.zeros(shape=(4, 4))
SCORING_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

TREE_NODE_INDEX = {}
TREE_ROOT = None


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

    # we parse the file
    parse_fasta_file()

    # we calculate the alignments and build the distance matrix
    compute_global_distances()

    # we print the distance matrix!
    print_distance_matrix()

    #
    compute_upgma()

    #
    print_newick_tree()

    print "exiting..."

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


def parse_fasta_file():
    print "parsing fasta file..."

    read_sequences()

    global DISTANCE_MATRIX
    DISTANCE_MATRIX = numpy.zeros((len(SEQUENCE_LIST), len(SEQUENCE_LIST)))

    print "done"
    print


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

        new_sequence.seq_data = seq_data

        id += 1


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


def compute_global_distance_matrix(seq_1, seq_2):
    seq_data_1 = ' ' + seq_1.seq_data
    seq_data_2 = ' ' + seq_2.seq_data

    matrix = numpy.zeros((len(seq_data_2), len(seq_data_1)))

    i = 1
    j = 1
    while i < len(seq_data_2):
        while j < len(seq_data_1):
            diagonal = matrix[i - 1][j - 1] + SCORING_MATRIX[SCORING_INDEX[seq_data_1[j]]][SCORING_INDEX[seq_data_2[i]]]

            top = matrix[i - 1][j] + GAP_PENALTY

            left = matrix[i][j - 1] + GAP_PENALTY

            matrix[i][j] = max(top, left, diagonal, 0)

            j += 1

        i += 1
        j = 1

    return matrix


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

        if d >= l and d >= u:
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

    if len(align_1) > len(align_2):
        align_2 = ((len(align_1) - len(align_2)) * '-') + align_2
    elif len(align_2) > len(align_1):
        align_1 = ((len(align_2) - len(align_1)) * '-') + align_1

    distance = 0
    length = len(align_1)

    i = 0
    while i < length:
        if align_1[i] != align_2[i]:
            distance += 1

        i += 1

    DISTANCE_MATRIX[seq_1.id][seq_2.id] = round(float(distance) / length, 5)
    DISTANCE_MATRIX[seq_2.id][seq_1.id] = round(float(distance) / length, 5)


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


def compute_upgma():
    print "computing UPGMA for global distance matrix..."

    node = None

    cluster_list = []
    for seq in SEQUENCE_LIST:
        SEQUENCE_INDEX[seq.name] = seq.id
        cluster_list.append([seq.name])

    while len(cluster_list) > 1:
        next_c_i, next_c_j, score = find_next_cluster(cluster_list)

        if str(cluster_list[next_c_i]) in TREE_NODE_INDEX:
            node_1 = TREE_NODE_INDEX[str(cluster_list[next_c_i])]
        else:
            node_1 = TreeNode(str(cluster_list[next_c_i][0]), score / 2)

        if str(cluster_list[next_c_j]) in TREE_NODE_INDEX:
            node_2 = TREE_NODE_INDEX[str(cluster_list[next_c_j])]
        else:
            node_2 = TreeNode(str(cluster_list[next_c_j][0]), score / 2)

        if isinstance(node_1, TreeCluster):
            node_1.score = score / 2 - node_1.score

        if isinstance(node_2, TreeCluster):
            node_2.score = score / 2 - node_2.score

        node = TreeCluster(score / 2)
        node.left = node_1
        node.right = node_2

        TREE_NODE_INDEX[str(cluster_list[next_c_i] + cluster_list[next_c_j])] = node

        new_cluster_list_item = cluster_list[next_c_i] + cluster_list[next_c_j]

        cluster_list[next_c_i] = new_cluster_list_item
        del cluster_list[next_c_j]

    global TREE_ROOT
    TREE_ROOT = node

    print "done"
    print


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

            if not initialized:
                min_score = score
                initialized = True
            elif score < min_score:
                next_c_i = i
                next_c_j = j
                min_score = score

            j += 1

        i += 1
        j = i + 1

    return next_c_i, next_c_j, min_score


def compute_cluster_distance(c_i, c_j):
    score = 0.0

    for x in c_i:
        for y in c_j:
            score += DISTANCE_MATRIX[SEQUENCE_INDEX[x]][SEQUENCE_INDEX[y]]

    score /= (len(c_i) * len(c_j))

    return score


def print_newick_tree():
    print "Newick Tree:"
    print

    tree = calculate_newick_tree(TREE_ROOT)

    out = file(NEWICK_TREE_FILE, 'w')
    out.write(tree)

    print tree
    print


def calculate_newick_tree(node):
    if isinstance(node, TreeNode):
        return "{0}:{1}".format(node.label, node.score)

    elif node != TREE_ROOT:
        return "({0}, {1}):{2}".format(calculate_newick_tree(node.left),
                                       calculate_newick_tree(node.right),
                                       node.score)
    else:
        return "({0}, {1});".format(calculate_newick_tree(node.left), calculate_newick_tree(node.right))


if __name__ == "__main__":
    main()
