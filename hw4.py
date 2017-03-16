import argparse
import numpy
sequence_list = []
distance_matrix = None
scoringMatrix = numpy.zeros(shape=(4,4))
#use scoringdict['character'] to look up value associated with that character for use in scoring matrix
scoringDict = {'A': 0, 'C': 1, 'G':2, 'T':3 }



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

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-F', '-f', '--file', type=file, help="fasta filename", required=True)
	parser.add_argument('-T', '-t', '--width', type=str, help="output file", default="output.fasta")
	parser.add_argument('-S', '-s', '--score', type=file, help="scoring matrix", default=None)
	parser.add_argument('-G', '-g', '--g', type=int, help="gap penalty", default=None)
	args = parser.parse_args()

	global WIDTH_CONST
	WIDTH_CONST = args.width
    #we parse the file
	parse_file(args.file)
	global g
	g = args.g
	scoringMatrix = parseScoreMatrix(args.score)

    #we calculate the alignments and build the distance matrix
	compute_global_distances()
    #we print the distance matrix!
	print_distance_matrix()


def parse_file(file):
    read_sequences(file)
    #we read the sequences from the file and initialize a global distance matrix that we will use 
    #the matrix is initialized to the size of the sequence lists, as we will need a N x M matrix to account for
    #each of the sequences
    global distance_matrix
    distance_matrix = numpy.zeros((len(sequence_list), len(sequence_list)))

def parseScoreMatrix(scoreFile):
	counter = 0
	file = scoreFile
	matrixSeq = []
	fo = file.read()
	for x in range(0, len(fo)):
		if fo[x].isdigit() and fo[x-1] != '-':
			y = int(fo[x])
			matrixSeq.append(y)
		if fo[x].isdigit() and fo[x-1] == '-':
			y = int(fo[x])
			y = y * -1
			matrixSeq.append(y)
	for x in range(0, 4):
		for y in range(0, 4):
			scoringMatrix[x, y] = matrixSeq[counter]
			counter = counter + 1
	return scoringMatrix


def read_sequences(file):
    id = 0

    #we read the lines from the file
    line = file.readline()
    while line != '' and line[0] == '>':
        #we get the title of the sequence here and append it to a global sequence list
        new_sequence = Sequence(id, line[1:-1])
        sequence_list.append(new_sequence)

        #we look for every character between > and the end of the line so that we get every base in there
        #this is then passed to the class that we have declared above for ease of access
        seq_data = ""
        line = file.readline()
        while line != '' and line[0] != '>':
            seq_data += line[:-1]
            line = file.readline()

        new_sequence.seq_data = seq_data
        #we ID each sequence to associate it with a given number
        id += 1
def compute_global_distances():
	i = 0
	j = 0

	while i < len(sequence_list):
		seq_1 = sequence_list[i]

		while j < len(sequence_list):
			if j != i:
				seq_2 = sequence_list[j]
                #we return the matrix calculated below, having associated the sequence from the sequence list with the ID in the class list to get the sequence
				matrix = compute_global_distance_matrix(seq_1, seq_2)

				show_alignment = False
				#as an alignment pair is calculated it gets printed right away, the permission bit is flipped here
				if i % 2 == 0 and j - i == 1:
					show_alignment = True
                #we then calculate the alignment from the matrix above
				calculate_global_distance(seq_1, seq_2, matrix, show_alignment)

			j += 1

		i += 1
		j = i

def compute_global_distance_matrix(seq_1, seq_2):
	#here we begin building our scoring matrix to calculate the alignment, initialize the data from the sequences
	#that we have in the sequence class
	seq_data_1 = ' ' + seq_1.seq_data
	seq_data_2 = ' ' + seq_2.seq_data
	#we initialize a matrix to the size of the sequences
	matrix = numpy.zeros((len(seq_data_2), len(seq_data_1)))

	#for the top row and column, we initialize their values according to the standard scoring system for glbobal alignment
	i = 0
	while i < len(seq_data_1):
		matrix[0][i] = i
		i += 1

	i = 0
	while i < len(seq_data_2):
		matrix[i][0] = i
		i += 1
   
	#we begin to iterate through the matrix to calculate our alignment. as is standard, we try to look for diagonal matches if possible
	i = 1
	j = 1
	while i < len(seq_data_2):
		while j < len(seq_data_1):

			# calc diagonal
			#if characters are aligned, we say its a match and give it a score of 0 so we can minimize it
			#if characters arent aligned we say its a mismatch so that we increase the distance
			diagonal = matrix[i - 1][j - 1]
			if seq_data_1[j] != seq_data_2[i]:
				x = scoringDict[seq_data_1[j]]
				y = scoringDict[seq_data_2[i]]
				diagonal += scoringMatrix[x,y]

			#if its a match we will have a score of the diagonal plus 0 - again so we minimize it
			#if its a mismatch, it will be score of the diagonal plus 1 - increasing the distance
			# calc top
			#same concept for the top, only we specifically use indel here
			x = scoringDict[seq_data_1[j]]
			y = scoringDict[seq_data_2[i]]
			mismatch = scoringMatrix[x,y]
			top = matrix[i - 1][j] + mismatch
      
			# calc left
			#same concept for the left, again only using indel here
			x = scoringDict[seq_data_1[j]]
			y = scoringDict[seq_data_2[i]]
			mismatch = scoringMatrix[x,y]
			left = matrix[i][j - 1] + mismatch
			#we calculate the minmimum to minmize the distance and then put that value in for the matrix
			matrix[i][j] = max(top, left, diagonal)

			j += 1

		i += 1
		j = 1

	return matrix

def calculate_global_distance(seq_1, seq_2, matrix, show_alignment=False):
    #again, we associate the sequences with the sequences in the Sequence class
    seq_data_1 = seq_1.seq_data
    seq_data_2 = seq_2.seq_data
    #we initialize our alignment strings
    align_1 = ""
    align_2 = ""
    #we get the length of the sequences here to properly iterate through the matrix. we will use a backtrace through the matrix
    i = len(seq_data_2) - 1
    j = len(seq_data_1) - 1
    while i > 0 and j > 0:
        #we take the values from three directions, above, left, diagonal. We are looking for the minimum here, and will choose from these values below
        u = matrix[i - 1][j]
        l = matrix[i][j - 1]
        d = matrix[i - 1][j - 1]
        #if the diagonal is at least less than or equal to both the left and top values, we choose the diagonal for our value as it represents a match between the two sequences
        if d >= l and d >= u:
            align_1 = seq_data_1[j] + align_1
            align_2 = seq_data_2[i] + align_2
            i -= 1
            j -= 1
        #if the left side is greater than the diagonal but less then the top, we choose the left side for our value, and it represents a gap in the second sequence
        elif l >= u:
            align_1 = seq_data_1[j] + align_1
            align_2 = '-' + align_2
            j -= 1
        #if the top is smaller then the left and the right, we choose the top and it also represents a gap in the first sequence
        else:
            align_1 = '-' + align_1
            align_2 = seq_data_2[i] + align_2
            i -= 1
    #if we reach the end of the first sequence we append the end of the sequence to the alignment here
    while j >= 0:
        align_1 = seq_data_1[j] + align_1
        j -= 1
    #if we reach the end of the second sequence we append the end of the sequence to the alignment here
    while i >= 0:
        align_2 = seq_data_2[i] + align_2
        i -= 1

    #we compare the lengths of the alignments here - if one is longer than the other, we correct it so that it will align properly
    if len(align_1) > len(align_2):
        align_2 = ((len(align_1) - len(align_2)) * '-') + align_2
    elif len(align_2) > len(align_1):
        align_1 = ((len(align_2) - len(align_1)) * '-') + align_1
    
    distance = 0
    length = len(align_1)

    #if the two sequences do not align, there is either a mismatch or an indel, so we increase the 'distance' by the gap provided on the command prompt
    i = 0
    while i < length:
        if align_1[i] != align_2[i]:
            distance += g

        i += 1
    #we build the distance matrix here. the formula is the mismatches/length, portrayed here as distance/length
    #e compare each sequence to every other sequence, using only the top half of the matrix
    #we do not need to use the full matrix, as the top half/bottom half of the matrix are simply mirrored across the diagonal
    distance_matrix[seq_1.id][seq_2.id] = round(float(distance)/length, 5)
    distance_matrix[seq_2.id][seq_1.id] = round(float(distance) / length, 5)



def print_distance_matrix():
    print "Distance Matrix for all Sequences:"
    print

    line = ""
    #we build up each line in the distance matrix here, and format it according to the specifications that it be within a 0-1 decimal range
    for row in distance_matrix:
        for column in row:
            line += "{0:.5f} ".format(round(column, 5))
        #we print the line, and then clear it for the nextl ine
        print line
        line = ""

if __name__ == "__main__":
	main()