#Usage: python pretty_bad_simulator.py <fasta.fa> <prob_to_include_a_seq> <subs_prob> <del_prob> <ins_prob>
import sys
import random

file = sys.argv[1]
includeProb=float(sys.argv[2])
subProb=float(sys.argv[3])
delProb=float(sys.argv[4])
insProb=float(sys.argv[5])

with open(sys.argv[1]) as fin:
    for line in fin:
        line = line.rstrip()
        if line[0]==">":
            header = line
        else:
            seq = line

            if random.random()<includeProb: #check if we include this sequence in the output
                #yes
                print header

                #mutate the sequence
                seqMutated = list(seq)
                for i in range(len(seqMutated)):
                    randomNb = random.random()
                    if randomNb<subProb:
                        #do a subs
                        bases=["A", "C", "G", "T"]
                        bases.remove(seqMutated[i])
                        seqMutated[i]=random.sample(bases, 1)[0]
                    elif randomNb<subProb+delProb:
                        #do a del
                        seqMutated[i]=""
                    elif randomNb<subProb+delProb+insProb:
                        #do a insertion
                        bases = ["A", "C", "G", "T"]
                        seqMutated[i]+=random.sample(bases,1)[0]

                #print the mutated sequence
                print "".join(seqMutated)

