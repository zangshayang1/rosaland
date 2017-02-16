import sys
import time
import copy

"""
The solution for each question can be found under this class in order
"""
class Bioinfo_StrongHold(object):

    """
    Q1 counting DNA nucleotides
    Input: a file containing a DNA seq terminated by '\n', whose length < 1000nt, which means it can fit in memory.
    Output: four integers indicating the counts of 'A', 'C', 'G', 'T' in this order.
    """
    def count_nucleotides(self, filepath):
        f = open(filepath, 'r')

        # read in all the content.
        # it can be done this way because it fits in the memory.
        seq = f.read().rstrip('\n')

        # init dict data structure this way
        # because input only contains "ATGC" for possibilities.
        dna_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

        for c in seq:
            dna_dict[c] += 1

        return dna_dict['A'], dna_dict['C'], dna_dict['G'], dna_dict['T']


    """
    Q2 transcribing DNA to RNA
    Input: same as Q1
    Output: transcribed RNA corresponding to DNA input
    """
    def transcribeDNA2RNA(self, filepath):
        f = open(filepath, 'r')

        # refer to Q1
        seq = f.read().rstrip('\n')

        # this is the most straightforward way but we are not going to use
        # ---------------------------------------------------------------
        # rna = ""
        # for c in seq:
        #     if c == 'T':
        #         rna += 'U'
        #     else:
        #         rna += c
        # return rna
        # ---------------------------------------------------------------

        # instead we use list comprehension + join() method - refer to string_building_test()
        rna = ''.join(['U' if seq[i] == 'T' else seq[i] for i in range(len(seq))])
        return rna


    """
    Q3 complementing a strand of DNA
    Input: same as Q1
    Output: a complementary strand of the given DNA
    """
    def complementDNA(self, filepath):
        f = open(filepath, 'r')

        # refer to Q1
        seq = f.read().rstrip('\n')

        complementary = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

        # refer to Q2, in case you don't know [1, 2, 3][::-1] -> [3, 2, 1]
        rDNA = ''.join([complementary[seq[i]] for i in range(len(seq))][::-1])
        return rDNA


    """
    Q4 computing GC content
    Input: something like below ...
           NOTE: each line is separated by newline character '\n'
           NOTE: the input contains at most 10 samples with each at most 1kbp, which still fits in the memory
    ------------------------------------------------------------
    >Rosalind_6404
    CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
    TCCCACTAATAATTCTGAGG
    >Rosalind_5959
    CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
    ATATCCATTTGTCAGCAGACACGC
    >Rosalind_0808
    CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
    TGGGAACCTGCGGGCAGTAGGTGGAAT
    ------------------------------------------------------------
    Output:
    ------------------------------------------------------------
    Rosalind_0808 -> highest GC content sample label
    60.919540 -> highest GC content value
    ------------------------------------------------------------
    """
    def computeGC(self, filepath):
        f = open(filepath, 'r')

        # refer to Q1
        data = f.read()

        # NOTE: we left rstrip() this time
        #       becase there are a bunch of '\n' involved in sequences
        #       even when they belongs to the same sample.
        data = data.split('>')[1:] # special case1: ">Haha>Hehe>".split('>') -> ["", "Haha", "Hehe", ""]
        # now data -> ["Rosalind_0808\nXXXXXXX\nXXXXX\n", "Rosalind_5959\nXXXX\n"]
        for i in range(len(data)):
            sample_seq = data[i]
            sample_seq = sample_seq.split('\n')[:-1]
            # sample_seq -> ["Rosalind_0808", "XXXXX", "XXXXX"], why [:-1], refer to special case1
            sample = sample_seq[0] # sample -> "Rosalind_0808"
            seq = ''.join(sample_seq[1:]) # seq -> "XXXXXXXXXXXX"
            data[i] = [sample, seq]
            # store back reconstructed data

            # an unsolved problem for myself
            # -----------------------------------------------------------------
            # sample_seq[:1].append(''.join(sample_seq[1:])) produces None. WHY?
            # -----------------------------------------------------------------


        # the following function computes the GC content of a singe line of sequence.
        # it is nested under the computeGC() which means it can only be called within the scope of computeGC()
        # also because it is not an independent method of the Bioinfo_StrongHold() class
        # no "self" is allowed as the first argument of the function
        # NOTE: the underline '_' in the function name indicates this is a "helper" function called only by computeGC() not by human
        #       while computeGC() is called by human. Put it in terminology, this is a "private" function.
        def _computeGC_percentage(seq):
            gc = 0
            tot = len(seq)
            for c in seq:
                if c == 'G' or c == 'C':
                    gc += 1
            return float(gc) * 100 / tot

        maxSample = None
        maxGC = -1
        for i in range(len(data)):
            sample_seq = data[i]
            gc = _computeGC_percentage(sample_seq[1]) # why using so many functions? so the logic of the code is clear and self-explanatory.
            if gc > maxGC:
                maxGC = gc
                maxSample = sample_seq[0]
        # I didn't use an additional dictionary to store every sample and their GC content?
        # Because it saves "space and time", you don't have to keep the dictionary in memory and no need to loop over the dictionary to find the max
        return maxSample, maxGC


    """
    Q5 count point mutations
    Input: two sequences of the same length
    Output: an integer indicating how many point mutations there are
    """
    def countPointMutation(self, filepath):
        f = open(filepath, 'r')
        content = f.readlines()
        s, t = content[0], content[1]
        count = 0
        for i in range(len(s)):
            if s[i] == t[i]:
                continue
            count += 1
        return count

    """
    Q6 translate RNA to protein
    Input: a string of RNA with length < 10kbp
    Output: a string of protein
    """
    def _load_RNA2protein_table(self):
        # this function will load "dataset/rna2protein.txt" and return a dictionary
        # it is quite ugly due to the given data in a quite ugly format.
        table = {}
        ls_all = []
        f = open('dataset/rna2protein.txt', 'r')
        lines = f.readlines()
        for line in lines:
            ls_line = line.rstrip().split(' ')
            ls_all.extend([item for item in ls_line if item != ''])
        print ls_all
        for i in range(0, len(ls_all), 2):
            tri_key, aa = ls_all[i], ls_all[i + 1]
            table[tri_key] = aa
        return table


    def translateRNA2protein(self, filepath):
        f = open(filepath, 'r')

        # refer to Q1
        rna = f.read().rstrip()

        table = self._load_RNA2protein_table()

        protein = ""
        for i in range(0, len(rna), 3):
            triplet = rna[i: i + 3]
            aa = table[triplet]
            if aa == "Stop":
                if i + 3 < len(rna):
                    raise Exception("ERROR in translateRNA2protein(): Stop condon appears in the middle of the input RNA.")
                return protein
            protein += aa

        raise Exception("ERROR in translateRNA2protein(): Missing stop condon at the end of the input RNA.")


    """
    Q7 find motif in DNA
    Input: a file containing two lines; 1st line is DNA sequence; 2nd line is the motif sequence (the length of motif is less than the length of sequence)
    Output: a sequence of integers indicating where the motif is found in the given DNA. indexing starts from 1 instead of 0.
    """
    def findMotif(self, filepath):
        f = open(filepath, 'r')
        content = f.readlines()
        seq, motif = content[0].rstrip(), content[1].rstrip()

        # init a list to store positions
        positions = []

        # this is a No.1 classic algorithm problem using 2 pointers
        # let's denote the following 2 pointer algorithm as classicAlgo1, it will be referred soon.
        i = 0
        while i < len(seq):
            j = 0 # every time you pick a new starting position, init j to be 0

            # if the first character doesn't match, advance starting position i.
            if seq[i] != motif[j]:
                i += 1
                continue

            # else: go in a while loop to compare characters one by one.
            while i < len(seq) and j < len(motif) and seq[i + j] == motif[j]:
                j += 1

            # if motif found:
            if j >= len(motif):
                positions.append(i + 1)
                i += 1
            # if motif not found:
            elif seq[i + j] != motif[j]:
                # case1: seq[i + j] != motif[j] -> the subsequence starting from seq[i] differs from motif on motif[j]
                i += 1
            else:
                # case2: i >= len(seq) -> end while loop
                break

        return positions



    """
    Q8 overlap graph
    Input: Seq data in FASTA format.
    ------------------------------------------------------------------------
    >Rosalind_0498
    AAATAAA
    >Rosalind_2391
    AAATTTT
    >Rosalind_2323
    TTTTCCC
    >Rosalind_0442
    AAATCCC
    >Rosalind_5013
    GGGTGGG
    ------------------------------------------------------------------------
    Output: list two samples in a line, with the 1st sample's suffix overlap with the 2nd sample's prefix by 3.
    ------------------------------------------------------------------------
    Rosalind_0498 Rosalind_2391
    Rosalind_0498 Rosalind_0442
    Rosalind_2391 Rosalind_2323
    ------------------------------------------------------------------------
    """
    def overlapGraph(self, filepath):
        f = open(filepath, 'r')

        # init a storage for sample and corresponding seq
        fasta = {}

        # because the input is in "ugly" format, we need to reconstruct them in a list before storing them in the above dictionary
        content = []

        # again in real life, it is highly likely that the entire FASTA file doesn't fit in the machine's mmemory
        # in which case, we can't not call f.readlines(), but we can call the following:
        line = f.readline()
        # when the file pointer reaches the end of the file, f.readline() will return an empty string "", which will end the while loop
        while line:
            # clean
            content.append(line.rstrip())
            # read the next line
            line = f.readline()

        # the following algorithm is a variation of the above classicAlgo1
        # ofc, you've encountered similar problem in computeGC(), let's try different ways of doing things to get your brain think.
        i = 0
        while i < len(content):
            sample = content[i][1:] # the first one must be sample name, so grab it.
            seq = "" # start to build sequence
            while i + 1 < len(content) and content[i + 1][0] != '>': # as long as the next element is not sample name
                seq += content[i + 1]
                i += 1

            fasta[sample] = seq # don't forget our goal to build the dictionary

            if i >= len(content):
                break
            else: # another sample name at content[i + 1]
                i += 1

        # now think about the content[] and fasta{}
        # even if you read the input file line by line, both content and fasta basically keep the same amount of data in memory.
        # it doens't make sense because we are assuming this amount of data can't be cached in the memory.
        # the only way to solve this problem is to deal with data in batch and store content and fasta on disk.
        # it gets very complicated here, just keep the question in mind. Maybe you can solve it someday.


        # the idea behind the following algorithm is basically pair-wise comparison.
        # if we have N samples in fasta dictionary, we can compare the 1st sample's suffix with the rest (N - 1) of samples' prefix
        # the above forms (N - 1) comparisons
        # if we do the same for the 2nd sample? again (N - 1) comparisons to be made.
        # THINK about why this algorithm can solve our problem? 
        # so in total, we need to make N(N - 1) comparisons.
        # thus we take the dominant term N^2 and claim the time complexity of this algorithm is O(N^2). This is the so-called bigO notation.
        # big O notation is a measure of how efficient your algorithm is. You need to know this and analyze algorithm like this in the future.
        result = []
        for sample_s in fasta:
            for sample_p in fasta:
                if sample_s == sample_p:
                    continue
                if fasta[sample_s][-3:] != fasta[sample_p][:3]:
                    continue
                result.append((sample_s, sample_p)) # this is where we can use tuple because we don't need to make a change any more

        # output format
        for tup in result:
            print tup[0], tup[1]

        return ;





"""
data path
"""
RELATIVE_PATH = "./dataset/"

q1_data = "countingDNAnucleotides.txt"
q2_data = "transcribingDNAtoRNA.txt"
q3_data = "complementingDNA.txt"
q4_data = "computeGCcontent.txt"
q5_data = "countPointMutations.txt"
q6_data = "translatingRNA2protein.txt"
q7_data = "findMotif_inDNA.txt"
q8_data = "overlapGraph.txt"


"""
run solution code, comment off to see.
"""
so = Bioinfo_StrongHold()
# print so.count_nucleotides(RELATIVE_PATH + q1_data)
# print so.transcribeDNA2RNA(RELATIVE_PATH + q2_data)
# print so.complementDNA(RELATIVE_PATH + q3_data)
# print so.computeGC(RELATIVE_PATH + q4_data)
# print so.countPointMutation(RELATIVE_PATH + q5_data)
# print so.translateRNA2protein(RELATIVE_PATH + q6_data)
# print so.findMotif(RELATIVE_PATH + q7_data)
print so.overlapGraph(RELATIVE_PATH + q8_data)


"""
Tests
"""
class Test(object):

    def string_building_test(self, target_length):
        # the following two are nested functions() under this function - string_building_test()
        # which means the following two can only be called within the scope of this function - string_building_test()

        # since they are nested functions that are belong to string_building_test()
        # they are considered an independent method in class Test(), thus no "self" as the first input argument.
        def build_from_empty_string(target_length):
            start_time = time.time()
            s = ""
            for _ in xrange(target_length):
                s += 'a'
            end_time = time.time()
            print end_time - start_time
            return ;

        def build_from_list_comprehension(target_length):
            start_time = time.time()
            s = ''.join(['a' for _ in xrange(target_length)])
            end_time = time.time()
            print end_time - start_time
            return ;

        build_from_empty_string(target_length)
        build_from_list_comprehension(target_length)
        return ;

"""
Comment off to start test cases
"""
# test = Test()
# test.string_building_test(10000000)
