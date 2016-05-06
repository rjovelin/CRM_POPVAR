

def celegans_three_prime_UTR_length(celegans_gff):
    '''
    (file) -> list
    Returns a list with the length of the annotated 3' UTRs in C. elegans gff annotation file
    '''
    #create list to store UTR length
    UTR = []
    # opengff file for reading
    cel = open(celegans_gff, 'r')
    # go through the file, extract UTR length and store in list
    for line in cel:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            if len(line) >= 5:
                if line[1] == 'WormBase':
                    if line[2] == 'three_prime_UTR':
                        UTR_length = (int(line[4]) - int(line[3])) + 1
                        UTR.append(UTR_length)
    cel.close()
    return UTR

def get_percentile(L, percentile):
    '''
    (list) -> num
    Return the value corresponding to the percentile from the list of values L
    Precondition: percentile is not in %
    '''

    # order the list
    L.sort()

    # use the nearest rank method to find the percentile rank
    Q = percentile / 100 * len(L)

    if int(Q+1) >= len(L):
        Qposition = int(Q-1)
    else:
        Qposition = int(Q+1)

    return L[Qposition]
