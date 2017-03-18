import argparse
import numpy


FASTA_FILE = ""
SCORING_FILE = ""
GAP_PENALTY = -1
NEWICK_TREE_FILE = ""

SCORING_MATRIX = numpy.zeros((4, 4))
SEQUENCE_LIST = []
DISTANCE_MATRIX = None

# use scoringdict['character'] to look up value associated with that character for use in scoring matrix
scoringDict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class Sequence:
    def __init__(self, id, name):
        self.id = id
        self.name = name
        self.seq_data = ""


class GlobalAlignmentResult:
    def __init__(self, seq_1, seq_2, distance, length, align_1, align_2):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        self.distance = distance
        self.length = length
        self.align_1 = align_1
        self.align_2 = align_2


class TreeNodes:
    def __init__(self):
        self.dist = 0
        self.label = 0
        self.parent = 0
        self.right = None
        self.left = None


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
    parse_file()

    global DISTANCE_MATRIX
    DISTANCE_MATRIX = numpy.zeros((len(SEQUENCE_LIST), len(SEQUENCE_LIST)))

    # we calculate the alignments and build the distance matrix
    compute_global_distances()

    # we print the distance matrix!
    print_distance_matrix()

    #
    #UPGMA()


def parse_score_matrix():
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


def parse_file():
    # we read the sequences from the file and initialize a global distance matrix that we will use
    # the matrix is initialized to the size of the sequence lists, as we will need a N x M matrix to account for
    # each of the sequences
    id = 0

    # we read the lines from the file
    line = FASTA_FILE.readline()
    while line != '' and line[0] == '>':
        # we get the title of the sequence here and append it to a global sequence list
        new_sequence = Sequence(id, line[1:-1])
        SEQUENCE_LIST.append(new_sequence)

        # we look for every character between > and the end of the line so that we get every base in there
        # this is then passed to the class that we have declared above for ease of access
        seq_data = ""
        line = FASTA_FILE.readline()
        while line != '' and line[0] != '>':
            seq_data += line[:-1]
            line = FASTA_FILE.readline()

        new_sequence.seq_data = seq_data
        # we ID each sequence to associate it with a given number
        id += 1


def compute_global_distances():
    i = 0
    j = 0

    while i < len(SEQUENCE_LIST):
        seq_1 = SEQUENCE_LIST[i]

        while j < len(SEQUENCE_LIST):
            if j != i:
                seq_2 = SEQUENCE_LIST[j]

                # we return the matrix calculated below, having associated the sequence from the sequence list with the
                # ID in the class list to get the sequence
                matrix = compute_global_distance_matrix(seq_1, seq_2)

                # we then calculate the alignment from the matrix above
                calculate_global_distance(seq_1, seq_2, matrix)

            j += 1

        i += 1
        j = i


def compute_global_distance_matrix(seq_1, seq_2):
    # here we begin building our scoring matrix to calculate the alignment, initialize the data from the sequences
    # that we have in the sequence class
    seq_data_1 = ' ' + seq_1.seq_data
    seq_data_2 = ' ' + seq_2.seq_data

    # we initialize a matrix to the size of the sequences
    matrix = numpy.zeros((len(seq_data_2), len(seq_data_1)))

    # for the top row and column, we initialize their values according to the standard scoring system for glbobal
	# alignment
    i = 0
    while i < len(seq_data_1):
        matrix[0][i] = i
        i += 1

    i = 0
    while i < len(seq_data_2):
        matrix[i][0] = i
        i += 1

    # we begin to iterate through the matrix to calculate our alignment. as is standard, we try to look for diagonal
	# matches if possible
    i = 1
    j = 1
    while i < len(seq_data_2):
        while j < len(seq_data_1):

            # calc diagonal
            # if characters are aligned, we say its a match and give it a score of 0 so we can minimize it
            # if characters arent aligned we say its a mismatch so that we increase the distance
            diagonal = matrix[i - 1][j - 1]
            if seq_data_1[j] != seq_data_2[i]:
                x = scoringDict[seq_data_1[j]]
                y = scoringDict[seq_data_2[i]]
                diagonal += SCORING_MATRIX[x, y]

            # if its a match we will have a score of the diagonal plus 0 - again so we minimize it
            # if its a mismatch, it will be score of the diagonal plus 1 - increasing the distance
            # calc top
            # same concept for the top, only we specifically use indel here
            x = scoringDict[seq_data_1[j]]
            y = scoringDict[seq_data_2[i]]
            mismatch = SCORING_MATRIX[x, y]
            top = matrix[i - 1][j] + mismatch

            # calc left
            # same concept for the left, again only using indel here
            x = scoringDict[seq_data_1[j]]
            y = scoringDict[seq_data_2[i]]
            mismatch = SCORING_MATRIX[x, y]
            left = matrix[i][j - 1] + mismatch
            # we calculate the minmimum to minmize the distance and then put that value in for the matrix
            matrix[i][j] = max(top, left, diagonal)

            j += 1

        i += 1
        j = 1

    return matrix


def calculate_global_distance(seq_1, seq_2, matrix):
    # again, we associate the sequences with the sequences in the Sequence class
    seq_data_1 = seq_1.seq_data
    seq_data_2 = seq_2.seq_data
    # we initialize our alignment strings
    align_1 = ""
    align_2 = ""
    # we get the length of the sequences here to properly iterate through the matrix. we will use a backtrace through
	# the matrix
    i = len(seq_data_2) - 1
    j = len(seq_data_1) - 1
    while i > 0 and j > 0:
        # we take the values from three directions, above, left, diagonal. We are looking for the minimum here, and will
		#  choose from these values below
        u = matrix[i - 1][j]
        l = matrix[i][j - 1]
        d = matrix[i - 1][j - 1]
        # if the diagonal is at least less than or equal to both the left and top values, we choose the diagonal for our
		#  value as it represents a match between the two sequences
        if d >= l and d >= u:
            align_1 = seq_data_1[j] + align_1
            align_2 = seq_data_2[i] + align_2
            i -= 1
            j -= 1
        # if the left side is greater than the diagonal but less then the top, we choose the left side for our value,
		# and it represents a gap in the second sequence
        elif l >= u:
            align_1 = seq_data_1[j] + align_1
            align_2 = '-' + align_2
            j -= 1
        # if the top is smaller then the left and the right, we choose the top and it also represents a gap in the first
		#  sequence
        else:
            align_1 = '-' + align_1
            align_2 = seq_data_2[i] + align_2
            i -= 1
            # if we reach the end of the first sequence we append the end of the sequence to the alignment here
    while j >= 0:
        align_1 = seq_data_1[j] + align_1
        j -= 1
    # if we reach the end of the second sequence we append the end of the sequence to the alignment here
    while i >= 0:
        align_2 = seq_data_2[i] + align_2
        i -= 1

    # we compare the lengths of the alignments here - if one is longer than the other, we correct it so that it will
	# align properly
    if len(align_1) > len(align_2):
        align_2 = ((len(align_1) - len(align_2)) * '-') + align_2
    elif len(align_2) > len(align_1):
        align_1 = ((len(align_2) - len(align_1)) * '-') + align_1

    distance = 0
    length = len(align_1)

    # if the two sequences do not align, there is either a mismatch or an indel, so we increase the 'distance' by the
	# gap provided on the command prompt
    i = 0
    while i < length:
        if align_1[i] != align_2[i]:
            distance += GAP_PENALTY

        i += 1
    # we build the distance matrix here. the formula is the mismatches/length, portrayed here as distance/length
    # e compare each sequence to every other sequence, using only the top half of the matrix
    # we do not need to use the full matrix, as the top half/bottom half of the matrix are simply mirrored across the
	# diagonal
    DISTANCE_MATRIX[seq_1.id][seq_2.id] = round(float(distance) / length, 5)
    DISTANCE_MATRIX[seq_2.id][seq_1.id] = round(float(distance) / length, 5)


def print_distance_matrix():
    print "Distance Matrix for all Sequences:"
    print

    line = ""
    # we build up each line in the distance matrix here, and format it according to the specifications that it be within
	# a 0-1 decimal range
    for row in DISTANCE_MATRIX:
        for column in row:
            line += "{0:.5f} ".format(round(column, 5))
        # we print the line, and then clear it for the nextl ine
        print line
        line = ""


def updateDistances(distMatrix, i, j):
    # idea: index into distance matrix
    # then, take the row/column that we saved down below
    # and compare the sequence distances associated with that row and column
    # and compare them to the distances with all the other sequences
    # if we wanted to update the distance for c and we have sequences a,b we want (D1(a,c) + d1(b,c))/2
    # then, we will need to delete those rows/colums from the matrix
    # deleting them is the easy part, gotta finangle with the updating distance part
    return 1


def UPGMA():
    print


if __name__ == "__main__":
    main()
