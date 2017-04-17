p = 3.0
def match(a, b):
    global p
    if a == b:
        return -1.0
    else:
        return 2.0

def compare(a, b):
    if a == '-' or b == '-':
        return 0
    elif a == b:
        return 1.0
    else:
        return -2.0


class sequences:
    'sequence list for alignment'

    def __init__(self, seqs, a, b, mis):
        self.seqs = seqs
        self.num = len(seqs)
        self.a = a
        self.b = b
        self.mis = mis

    # do all
    def all(self):
        self.init_matrix()
        self.build()
        self.complete()

    # calculate the difference
    def init_matrix(self):
        self.score = []
        self.ext = []
        for i in range(0, self.num):
            self.score.append([0.0] * self.num)
            self.ext.append([]*self.num)
        for i in range(self.num):
            # num, left, right, seqs, size, distance
            node = [[i], [], [], [self.seqs[i]], 0, -1]
            self.score[i].append(node)
            for j in range(i):
                align = self.pairwise(i, j)
                align2 =[align[1], align[0]]
                leng = len(align[0])
                matc = 0.0
                for i in range(leng):
                    if(align[0][i] == align[1][i]):
                        matc += 1.0
                self.score[i][j] = 1 - matc/leng
                self.score[j][i] = self.score[i][j]
                self.ext[i][j].append(100 * matc/leng)
                self.ext[i][j].append(align)
                kkk1 = []
                kkk2 = []
                for l in range(len(seqs[i])):
                    kkk1.append({})
                for l in range(len(seqs[j])):
                    kkk2.append({})
                numi = -1
                numj = -1
                for l in range(len(align[0])):
                    if(align[0][l] != '-'):
                        numi += 1
                    if(align[1][l] != '-'):
                        numj += 1
                    if(align[0][l] != '-' and align[1][l] != '-'):
                        kkk1[numi][numj] = self.ext[i][j][0]
                        kkk2[numj][numi] = self.ext[i][j][0]
                self.ext[i][j].append(kkk1)
                self.ext[i][j].append(kkk2)
                self.ext[j][i].append(100 * matc/leng)
                self.ext[j][i].append(align2)
        for i in range(self.num):
            for j in range(i):
                for s in range(self.num):
                    if(s == i or s == j):
                        continue
                    ois = self.ext[i][s][1]
                    ojs = self.ext[j][s][1]
                    numiss = -1
                    numisi = -1
                    numjsj = -1
                    numjss = -1
                    numis = -1
                    numjs = -1
                    while(1):
                        while(1):
                            numis+=1
                            if(ois[0][numis] !='-'):
                                numisi += 1
                            if(ois[1][numis] !='-'):
                                numiss += 1
                                break
                        while(1):
                            numjs+=1
                            if(ojs[0][numjs] !='-'):
                                numjsj += 1
                            if(ojs[1][numjs] !='-'):
                                numjss += 1
                                break
                        if(ojs[0][numjs] !='-' and ois[0][numis] !='-'):
                            if (self.ext[i][j][2][numisi].has_key(numjsj)):
                                self.ext[i][j][2][numisi][numjsj] += self.ext[i][s][0]
                            else:
                                self.ext[i][j][2][numisi][numjsj] = self.ext[i][s][0]

                            if (self.ext[i][j][3][numjsj].has_key(numisi)):
                                self.ext[i][j][3][numjsj][numisi] += self.ext[j][s][0]
                            else:
                                self.ext[i][j][3][numjsj][numisi] = self.ext[j][s][0]
                        if(numjss >= len(seqs[s]))
                            break
        for i in range(self.num):
            for j in range(i):
                self.ext[i][j].pop(0)
                self.ext[i][j].pop(0)
                self.ext[j][i].pop(0)
                self.ext[j][i].pop(0)
                self.ext[j][i].append(self.ext[i][j][1])
                self.ext[j][i].append(self.ext[i][j][0])
        for i in self.score:
            print i[0:self.num]

    # compare the two sequences
    def pairwise(self, i, j):
    a = -self.a
    b = -self.b
    s1 = self.seqs[i]
    s2 = self.seqs[j]
    num1 = len(s1)
    num2 = len(s2)
    D1 = [0] * (num1 + 1)
    D2 = [0] * (num1 + 1)
    I1 = [0] * (num1 + 1)
    I2 = [0] * (num1 + 1)
    M1 = [0] * (num1 + 1)
    M2 = [0] * (num1 + 1)
    Path1_d = [[]]  # record the path
    Path1_i = [[]]
    Path1_m = [[]]
    Path2_d = [['I']]
    Path2_i = [['I']]
    Path2_m = [['I']]
    Path_k = [[]]
    K = [0] * (num2 + 1)
    D1[1] = I1[1] = M1[1] = K[1] = b
    key = 2
    for i in range(1, num1 + 1):
        D1[i] = D1[i - 1] + a
        I1[i] = D1[i]
        M1[i] = D1[i]
        Path1_d.append(['D'] * i)
        Path1_i.append(['D'] * i)
        Path1_m.append(['D'] * i)
        Path2_d.append([])
        Path2_i.append([])
        Path2_m.append([])

    for i in range(1, num2 + 1):
        K[i] = K[i - 1] + a
        Path_k.append(['I'] * i)

    for i in range(0, num2):
        if key == 2:
            D2[0] = K[i + 1]
            I2[0] = K[i + 1]
            M2[0] = K[i + 1]
            Path2_m[0] = Path2_i[0] = Path2_d[0] = Path_k[i + 1]
            for j in range(1, num1 + 1):
                score_d = compare(s1[j - 1], s2[i])
            max_d = max(D2[j - 1] + a, I2[j - 1] + b, M2[j - 1] + b)
            max_i = max(I1[j] + a, M1[j] + b, D1[j] + b)
            max_m = max(D1[j - 1], I1[j - 1], M1[j - 1])
            if max_d == D2[j - 1] + a:
                Path2_d[j] = Path2_d[j - 1] + ['D']
            elif max_d == I2[j - 1] + b:
                Path2_d[j] = Path2_i[j - 1] + ['D']
            else:
                Path2_d[j] = Path2_m[j - 1] + ['D']
            if max_i == I1[j] + a:
                Path2_i[j] = Path1_i[j] + ['I']
            elif max_i == D1[j] + b:
                Path2_i[j] = Path1_d[j] + ['I']
            else:
                Path2_i[j] = Path1_m[j] + ['I']
            if max_m == D1[j - 1]:
                Path2_m[j] = Path1_d[j - 1] + ['M']
            elif max_m == I1[j - 1]:
                Path2_m[j] = Path1_i[j - 1] + ['M']
            else:
                Path2_m[j] = Path1_m[j - 1] + ['M']

            D2[j] = max_d
            I2[j] = max_i
            M2[j] = max_m + score_d
            key = 1
        else:
            D1[0] = K[i + 1]
            I1[0] = K[i + 1]
            M1[0] = K[i + 1]
            Path1_m[0] = Path1_i[0] = Path1_d[0] = Path_k[i + 1]
            for j in range(1, num1 + 1):
                score_d = compare(s1[j - 1], s2[i])
            max_d = max(D1[j - 1] + a, I1[j - 1] + b, M1[j - 1] + b)
            max_i = max(I2[j] + a, M2[j] + b, D2[j] + b)
            max_m = max(D2[j - 1], I2[j - 1], M2[j - 1])
            if max_d == D1[j - 1] + a:
                Path2_d[j] = Path2_d[j - 1] + ['D']
            elif max_d == I1[j - 1] + b:
                Path2_d[j] = Path2_i[j - 1] + ['D']
            else:
                Path2_d[j] = Path2_m[j - 1] + ['D']
            if max_i == I2[j] + a:
                Path2_i[j] = Path1_i[j] + ['I']
            elif max_i == D2[j] + b:
                Path2_i[j] = Path1_d[j] + ['I']
            else:
                Path2_i[j] = Path1_m[j] + ['I']
            if max_m == D2[j - 1]:
                Path2_m[j] = Path1_d[j - 1] + ['M']
            elif max_m == I2[j - 1]:
                Path2_m[j] = Path1_i[j - 1] + ['M']
            else:
                Path2_m[j] = Path1_m[j - 1] + ['M']

            D1[j] = max_d
            I1[j] = max_i
            M1[j] = max_m + score_d
            key = 2

    aligned = [[],[]]
    points = [0,0]
    if key == 2:
        max_score = max(M1[num1], D1[num1], I1[num1])
        if M1[num1] == max_score:
            Path = Path1_m[num1]
        elif D1[num1] == max_score:
            Path = Path1_d[num1]
        else:
            Path = Path1_i[num1]
    else:
        max_score = max(M2[num1], D2[num1], I2[num1])
        if M2[num1] == max_score:
            Path = Path2_m[num1]
        elif D2[num1] == max_score:
            Path = Path2_d[num1]
        else:
            Path = Path2_i[num1]
    for c in Path:
        if c == 'M':
            aligned[0].append(s1[points[0]])
            aligned[1].append(s2[points[1]])
            points[0] += 1
            points[1] += 1
        elif c == 'I':
            aligned[0].append('-')
            aligned[1].append(s2[points[1]])
            points[1] += 1
        else:
            aligned[1].append('-')
            aligned[0].append(s2[points[0]])
            points[0] += 1

    return aligned

    def buildtree(self):

        num = self.num
        diff = self.score
        d = 0
        di = -1
        dj = -1
        nj = []
        su = [0]*num
        for i in range(0, num):
            nj.append([0] * num)
        while(num >= 4):
            d = 0
            di = -1
            dj = -1
            for i in range(0, num):
                su[i] = sum(diff[i][:num])
            for i in range(0, num):
                for j in range(0, i):
                    nj[i][j] = (su[i] + su[j]) - (num - 2) * diff[i][j]
                    if(nj[i][j] > d):
                        d = nj[i][j]
                        di = i
                        dj = j

            diff[di][num][5] = 0.5 * diff[di][dj] + (su[di] - su[dj]) / (2 * (num -2))
            diff[dj][num][5] = diff[di][dj] - diff[di][num][5]
            node = [[], diff[di][num], diff[dj][num], [], 0, -1]
            diffk = []

            for i in range(0, num):
                if(i != di and i != dj):
                    #print diffk,i,di,dj,diff[di][i],diff[dj][i],diff[dj][di]

                    diffk.append((diff[di][i] + diff[dj][i] - diff[dj][di]) / 2)
            diffk.append(0)
            diffk.append(node)

            for i in range(0, num):
                if(i != di and i != dj):
                    del diff[i][di]
                    del diff[i][dj]

            del diff[di]
            del diff[dj]
            # now is num-1*num-1 matrix
            for i in range(0, num - 2):
                diff[i].insert(num - 2, diffk[i])
            diff.insert(num - 2, diffk)
            num = num - 1
        if(diff[0][1] < diff[0][2] and diff[0][1] < diff[1][2]):
            di = 0
            dj = 1
            dk = 2
        elif(diff[0][2] < diff[1][2]):
            di = 0
            dj = 2
            dk = 1
        else:
            di = 1
            dj = 2
            dk = 0
        diff[di][3][5]=(0.5*diff[di][dj])
        diff[dj][3][5]=(0.5*diff[di][dj])
        node = [[], diff[di][3], diff[dj][3], [], 0, -1]
        dis = (diff[di][dk] + diff[dj][dk])/2
        diff[dk][3][5]=(dis)
        node[5]=(dis - 0.5*diff[di][dj])
        t = [[],node,diff[dk][3],[],[],-1,-1,0]

        self.tree = t
        print "\n\n"
        #print t
        return t


    def complete(self):
        root = self.tree
        self.dfs(root)
        for seq in root[3]:
            print seq

    def dfs(self, node):
        if not node[1][3]:
            self.dfs(node[1])
        if not node[2][3]:
            self.dfs(node[2])
        seqs = self.multipairwise(node[1], node[2])
        node[3] = seqs
        node[4] = node[1][4] + node[2][4]

    # [num, left, right, seqs, weight, distance, height]
    def multipairwise(self, node1, node2):
        a = -self.a
        b = -self.b
        num1 = len(node1[3][0])
        num2 = len(node2[3][0])
        D1 = [0] * (num1 + 1)
        D2 = [0] * (num1 + 1)
        I1 = [0] * (num1 + 1)
        I2 = [0] * (num1 + 1)
        M1 = [0] * (num1 + 1)
        M2 = [0] * (num1 + 1)
        Path1_d = [[]]  # record the path
        Path1_i = [[]]
        Path1_m = [[]]
        Path2_d = [['I']]
        Path2_i = [['I']]
        Path2_m = [['I']]
        Path_k = [[]]
        K = [0] * (num2 + 1)
        D1[1] = I1[1] = M1[1] = K[1] = b
        key = 2
        num_pair = len(node1[3]) * len(node2[3])
        weight_sum = 0.0
        print node1[4]
        print node2[4]
        for k1 in range(len(node1[3])):
            for k2 in range(len(node2[3])):
                weight_sum += node1[4][k1] * node2[4][k2]
        weight_sum /= num_pair
        for i in range(1, num1 + 1):
            D1[i] = D1[i - 1] + a
            I1[i] = D1[i]
            M1[i] = D1[i]
            Path1_d.append(['D'] * i)
            Path1_i.append(['D'] * i)
            Path1_m.append(['D'] * i)
            Path2_d.append([])
            Path2_i.append([])
            Path2_m.append([])

        for i in range(1, num2 + 1):
            K[i] = K[i - 1] + a
            Path_k.append(['I'] * i)

        for i in range(0, num2):
            if key == 2:
                D2[0] = K[i + 1]
                I2[0] = K[i + 1]
                M2[0] = K[i + 1]
                Path2_m[0] = Path2_i[0] = Path2_d[0] = Path_k[i + 1]
                for j in range(1, num1 + 1):
                    score_d = 0
                    for k1 in range(len(node1[3])):
                        for k2 in range(len(node2[3])):
                            score_d += compare(node1[3][k1][j - 1], node2[3][k2][i]) * node1[4][k1] * node2[4][k2]
                    max_d = max(D2[j - 1] +  a, I2[j - 1] +  b, M2[j - 1] +  b)
                    max_i = max(I1[j] +  a, M1[j] +  b, D1[j] +  b)
                    max_m = max(D1[j - 1], I1[j - 1], M1[j - 1])
                    if max_d == D2[j - 1] +  a:
                        Path2_d[j] = Path2_d[j - 1] + ['D']
                    elif max_d == I2[j - 1] +  b:
                        Path2_d[j] = Path2_i[j - 1] + ['D']
                    else:
                        Path2_d[j] = Path2_m[j - 1] + ['D']
                    if max_i == I1[j] +  a:
                        Path2_i[j] = Path1_i[j] + ['I']
                    elif max_i == D1[j] + b:
                        Path2_i[j] = Path1_d[j] + ['I']
                    else:
                        Path2_i[j] = Path1_m[j] + ['I']
                    if max_m == D1[j - 1]:
                        Path2_m[j] = Path1_d[j - 1] + ['M']
                    elif max_m == I1[j - 1]:
                        Path2_m[j] = Path1_i[j - 1] + ['M']
                    else:
                        Path2_m[j] = Path1_m[j - 1] + ['M']

                    D2[j] = max_d
                    I2[j] = max_i
                    M2[j] = max_m + score_d / num_pair
                    key = 1
            else:
                D1[0] = K[i + 1]
                I1[0] = K[i + 1]
                M1[0] = K[i + 1]
                Path1_m[0] = Path1_i[0] = Path1_d[0] = Path_k[i + 1]
                for j in range(1, num1 + 1):
                    score_d = 0
                    for k1 in range(len(node1[3])):
                        for k2 in range(len(node2[3])):
                            score_d += compare(node1[3][k1][j - 1], node2[3]
                                             [k2][i]) * node1[4][k1] * node2[4][k2]

                    max_d = max(D1[j - 1] +  a, I1[j - 1] +  b, M1[j - 1] +  b)
                    max_i = max(I2[j] +  a, M2[j] +  b, D2[j] +  b)
                    max_m = max(D2[j - 1], I2[j - 1], M2[j - 1])
                    if max_d == D1[j - 1] +  a:
                        Path1_d[j] = Path1_d[j - 1] + ['D']
                    elif max_d == I1[j - 1] +  b:
                        Path1_d[j] = Path1_i[j - 1] + ['D']
                    else:
                        Path1_d[j] = Path1_m[j - 1] + ['D']
                    if max_i == I2[j] +  a:
                        Path1_i[j] = Path2_i[j] + ['I']
                    elif max_i == D2[j] +  b:
                        Path1_i[j] = Path2_d[j] + ['I']
                    else:
                        Path1_i[j] = Path2_m[j] + ['I']
                    if max_m == D2[j - 1]:
                        Path1_m[j] = Path2_d[j - 1] + ['M']
                    elif max_m == I2[j - 1]:
                        Path1_m[j] = Path2_i[j - 1] + ['M']
                    else:
                        Path1_m[j] = Path2_m[j - 1] + ['M']

                    D1[j] = max_d
                    I1[j] = max_i
                    M1[j] = max_m + score_d / num_pair
                    key = 2
        '''need to refine the sequence'''
        aligned = []
        seqnum1 = len(node1[3])
        seqnum2 = len(node2[3])
        points = [0] * (seqnum1 + seqnum2)
        for i in range(seqnum1 + seqnum2):
            aligned.append([])
        if key == 2:
            max_score = max(M1[num1], D1[num1], I1[num1])
            if M1[num1] == max_score:
                Path = Path1_m[num1]
            elif D1[num1] == max_score:
                Path = Path1_d[num1]
            else:
                Path = Path1_i[num1]
        else:
            max_score = max(M2[num1], D2[num1], I2[num1])
            if M2[num1] == max_score:
                Path = Path2_m[num1]
            elif D2[num1] == max_score:
                Path = Path2_d[num1]
            else:
                Path = Path2_i[num1]
        #print "path: "
        #print max_score
        #print Path
        for c in Path:
            if c == 'M':
                for i in range(seqnum1):
                    aligned[i].append(node1[3][i][points[i]])
                    points[i]+=1
                for i in range(seqnum2):
                    aligned[seqnum1 + i].append(node2[3][i][points[seqnum1 + i]])
                    points[seqnum1 + i]+=1
            elif c == 'I':
                for i in range(seqnum1):
                    aligned[i].append('-')
                for i in range(seqnum2):
                    aligned[seqnum1 + i].append(node2[3][i][points[seqnum1 + i]])
                    points[seqnum1 + i]+=1
            else:
                #print points
                for i in range(seqnum1):
                    aligned[i].append(node1[3][i][points[i]])
                    points[i]+=1
                for i in range(seqnum2):
                    aligned[seqnum1 + i].append('-')
        for i in range(seqnum1):
            if points[i] != num1:
                print points[i], num1
                quit()
        for i in range(seqnum2):
            if points[seqnum1 + i] != num2:
                print points[seqnum1 + i], num2
                quit()
        return aligned

# class tree:
if __name__ == "__main__":
    fp = open("data.txt","r")
    seq = []
    while (1):
        one = []
        a  = fp.readline()
        if not a:
            break
        for i in a:
            if i != '_':
                one.append(i)
        one.pop()
        seq.append(one)
    #print seq[0]
    fp.close()
    mul = sequences(seq[0:6], 2.0, 3.0, 1.0)
    mul.init_matrix()
    #for i in mul.score:
        #print i[0:6]
    mul.buildtree()
    mul.complete()
