import numpy

DISTANCE_MATRIX = numpy.zeros((5,5))
DISTANCE_MATRIX[0] = [0, 17, 21, 31, 23]
DISTANCE_MATRIX[1] = [17, 0, 30, 34, 21]
DISTANCE_MATRIX[2] = [21, 30, 0, 28, 39]
DISTANCE_MATRIX[3] = [31, 34, 28, 0, 43]
DISTANCE_MATRIX[4] = [23, 21, 39, 43, 0]

SEQUENCE_LIST = ['A', 'B', 'C', 'D', 'E']

SEQUENCE_INDEX = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}


def main():
    cluster_list = []
    for seq in SEQUENCE_LIST:
        cluster_list.append([seq]) # DO NOT WANT TO COPY ENTIRE SEQUENCE IN EACH LIST

    while len(cluster_list) > 1:
        next_c_i, next_c_j, score = find_next_cluster(cluster_list)

        new_cluster_list_item = cluster_list[next_c_i] + cluster_list[next_c_j]

        cluster_list[next_c_i] = new_cluster_list_item
        del cluster_list[next_c_j]


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

    score /= (len(c_i)*len(c_j))

    return score


if __name__ == "__main__":
    main()
