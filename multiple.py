p = 1.0
def match(a, b):
    global p
    if a == b:
        return 0
    else:
        return p

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
        #print "the num of seqs is "
        #print self.num
        for i in range(0, self.num):
            self.score.append([0.0] * self.num)
        maxx = 0.0
        for i in range(0, self.num):
            # num, left, right, seqs, weight, distance from parent, height,
            # size
            node = [1, [], [], [self.seqs[i]], [], -1, -1, 0]
            self.score[i].append(node)
            for j in range(0, i):
                self.score[i][j] = self.pairwise(i, j)
                self.score[j][i] = self.score[i][j]
                if maxx < self.score[i][j]:
                    maxx = self.score[i][j]
        maxx = 2*maxx
        for i in range ( 0, self.num):
            for j in range(0,i):
                self.score[i][j] /= maxx
                self.score[j][i] /= maxx
        for i in self.score:
            print i[0:self.num]

    # compare the two sequences
    def pairwise(self, i, j):  # penalty for gap is a(k-1)+b
        a = self.a
        b = self.b
        mis = self.mis
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
        K = [0] * (num2 + 1)
        D1[1] = I1[1] = M1[1] = K[1] = b
        key = 2
        for i in range(2, num1 + 1):
            D1[i] = D1[i - 1] + a
            I1[i] = D1[i]
            M1[i] = D1[i]
        for i in range(2, num2 + 1):
            K[i] = K[i-1] + a

        for i in range(0, num2):
            if key == 2:
                D2[0] = K[i + 1]
                I2[0] = K[i + 1]
                M2[0] = K[i + 1]
                for j in range(1, num1 + 1):
                    D2[j] = min(D1[j - 1], I1[j - 1], M1[j - 1]) + \
                        match(s1[j - 1], s2[i])
                    I2[j] = min(I1[j] + a, M1[j] + b, D1[j] + b)
                    M2[j] = min(M2[j - 1] + a, I2[j - 1] + b, D2[j - 1] + b)
                key = 1
            else:
                D1[0] = K[i + 1]
                I1[0] = K[i + 1]
                M1[0] = K[i + 1]
                for j in range(1, num1 + 1):
                    D1[j] = min(D2[j - 1], I2[j - 1], M2[j - 1]) + \
                        match(s1[j - 1], s2[i])
                    I1[j] = min(I2[j] + a, M2[j] + b, D2[j] + b)
                    M1[j] = min(M1[j - 1] + a, I1[j - 1] + b, D1[j - 1] + b)
                key = 2
        if (key == 2):
            return min(M1[num1], D1[num1], I1[num1])
        else:
            return min(M2[num1], D2[num1], I2[num1])

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
            node = [0, diff[di][num], diff[dj][num], [], [], -1, -1, 0]
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
        if(diff[0][1] < diff[0][2]):
            if(diff[0][1] < diff[1][2]):
                di = 0
                dj = 1
                dk = 2
        else:
            if(diff[0][2] < diff[1][2]):
                di = 0
                dj = 2
                dk = 1
        if(diff[1][2]<diff[0][2] and diff[1][2]<diff[0][1]):
            di = 1
            dj = 2
            dk = 0
        diff[di][3][5]=(0.5*diff[di][dj])
        diff[dj][3][5]=(0.5*diff[di][dj])
        node = [0, diff[di][3], diff[dj][3], [],[], -1, -1, 0]
        dis = (diff[di][dk] + diff[dj][dk])/2
        diff[dk][3][5]=(dis)
        node[5]=(dis - 0.5*diff[di][dj])
        t = [0,node,diff[dk][3],[],[],-1,-1,0]






        self.sizeoftree(t)
        self.weighttree(t, 1, -1.0)
        self.tree = t
        print "\n\n"
        #print t
        return t

    def sizeoftree(self, t):
        if(not t[1] and not t[2]):
            t[7] = 1
            return 1
        if(not t[1] and t[2]):
            size = self.sizeoftree(t[2]) + 1
            t[7] = size
            return size
        if(not t[2] and t[1]):
            size = self.sizeoftree(t[1]) + 1
            t[7] = size
            return size
        size = self.sizeoftree(t[1]) + self.sizeoftree(t[2]) + 1
        t[7] = size
        return size

    def weighttree(self, t, h, w):
        t[6] = h
        if(w == -1.0):
            w = 0
        else:
            w += t[5]/t[7]
        if(t[1]):
            self.weighttree(t[1], h + 1, w)
        if(t[2]):
            self.weighttree(t[2], h + 1, w)
        if (not t[1]) and (not t[2]):
            t[4].append(w)
            print t

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
                    max_d = max(D2[j - 1] + weight_sum * a, I2[j - 1] + weight_sum * b, M2[j - 1] + weight_sum * b)
                    max_i = max(I1[j] + weight_sum * a, M1[j] + weight_sum * b, D1[j] + weight_sum * b)
                    max_m = max(D1[j - 1], I1[j - 1], M1[j - 1])
                    if max_d == D2[j - 1] + weight_sum * a:
                        Path2_d[j] = Path2_d[j - 1] + ['D']
                    elif max_d == I2[j - 1] + weight_sum * b:
                        Path2_d[j] = Path2_i[j - 1] + ['D']
                    else:
                        Path2_d[j] = Path2_m[j - 1] + ['D']
                    if max_i == I1[j] + weight_sum * a:
                        Path2_i[j] = Path1_i[j] + ['I']
                    elif max_i == D1[j] + weight_sum *b:
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

                    max_d = max(D1[j - 1] + weight_sum * a, I1[j - 1] + weight_sum * b, M1[j - 1] + weight_sum * b)
                    max_i = max(I2[j] + weight_sum * a, M2[j] + weight_sum * b, D2[j] + weight_sum * b)
                    max_m = max(D2[j - 1], I2[j - 1], M2[j - 1])
                    if max_d == D1[j - 1] + weight_sum * a:
                        Path1_d[j] = Path1_d[j - 1] + ['D']
                    elif max_d == I1[j - 1] + weight_sum * b:
                        Path1_d[j] = Path1_i[j - 1] + ['D']
                    else:
                        Path1_d[j] = Path1_m[j - 1] + ['D']
                    if max_i == I2[j] + weight_sum * a:
                        Path1_i[j] = Path2_i[j] + ['I']
                    elif max_i == D2[j] + weight_sum * b:
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
