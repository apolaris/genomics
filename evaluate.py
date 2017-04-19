def do(seq):
    seqc = []
    col = []
    for i in range(len(seqs[0])):
        for j in range(len(seqs)):
            col.append(seqs[j][i])
        seqc.append(col)
        col = []
    score = 0.0
    full = 0.0
    count = 0
    order = [24, 12, 18, 6, 27, 9, 17, 26, 13, 4, 29, 3, 1, 5, 28, 19, 22, 11, 23, 15, 10, 16, 8, 25, 20, 2, 7, 0, 21, 14]
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            for k in range(len(seq[0])):
                if seq[i][k] != '-':
                    while seqc[count][order[i]] == '_':
                        count += 1
                    if seq[i][k] != seqc[count][order[i]]:
                        print "ERROR"
                        quit()
                    if seq[j][k] == seqc[count][order[j]]:
                        score += 1
                    count += 1
            count = 0
    for line in seqc:
        num = 0
        for c in line:
            if c != '_':
                num += 1
        full += num * (num - 1) / 2.0

    score /= full
    print score


fp = open("419a.txt","r")
fp2 = open("data20.txt","r")
seqa = []
seqb = []
while (1):
    one = []
    a  = fp.readline()
    if not a:
        break
    for i in a:
        one.append(i)
    one.pop()
    seqa.append(one)
while (1):
    one = []
    a  = fp2.readline()
    if not a:
        break
    for i in a:
        one.append(i)
    one.pop()
    seqb.append(one)
fp.close()
fp2.close()
seqs = seqb
do(seqa)
