import multiple

def sc(seq1,seq2):
    s1 = seq1 + []
    s2 = seq2 + []
    l = len(s1)
    i = 0
    score = 0
    
    while(i < len(s1)):
        if(s1[i] == '-' and seq2[i] == '-'):
            s1.pop(i)
            s2.pop(i)
            continue
        i+=1

    if(s1[0] == s2[0]):
        score-=1
    elif (s1[0] == '-' or s2[0] == '-'):
        score += 3
    else:
        score += 2
        
    #print s1,s2
    l = len(s1)
    for i in range(1,l):
        if(s1[i] == s2[i]):
            score-=1
        
        elif(s1[i] == '-'):
            if(s1[i-1] == '-'):
                score += 2
            else:                score += 3
        elif(s2[i] == '-'):
            if(s2[i-1] == '-'):
                score += 2
            else:
                score += 3
        else:
            score += 2
    return score

def evalue(seqs):
    num = len(seqs)
    summ = 0
    score = []
    for i in range(num):
        score.append([0]*num)

    for i in range(num):
        for j in range(i):
            score[i][j] = sc(seqs[i], seqs[j])
            score[j][i] = sc(seqs[i], seqs[j])
            summ += score[i][j]
    print score,"\n\n"
    print summ
    return score, summ



if __name__ == "__main__":
    fp = open("a.txt","r")
    seq = []
    while (1):
        one = []
        a  = fp.readline()
        if not a:
            break
        for i in a:
            one.append(i)
        one.pop()
        seq.append(one)
    print seq[0]
    fp.close()
    evalue(seq[0:6])

    evalue(seq[6:12])

    evalue(seq[12:18])

    evalue(seq[18:24])
