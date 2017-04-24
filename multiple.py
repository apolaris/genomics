import time

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
        for i in range(0, self.num):
            self.score.append([0.0] * self.num)
        for i in range(0, self.num):
            # num, left, right, seqs, weight, distance from parent, height,
            # size
            node = [[i], [], [], [self.seqs[i]], [], -1, -1, 0]
            self.score[i].append(node)
            for j in range(0, i):
                self.score[i][j] = self.pairwise(i, j)
                self.score[j][i] = self.score[i][j]
        self.distance = []

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
        Dnum1 = [0.0] * (num1 + 1)
        Dnum2 = [0.0] * (num1 + 1)
        Inum1 = [0.0] * (num1 + 1)
        Inum2 = [0.0] * (num1 + 1)
        Mnum1 = [0.0] * (num1 + 1)
        Mnum2 = [0.0] * (num1 + 1)
        Dn1 = [0.0] * (num1 + 1)
        Dn2 = [0.0] * (num1 + 1)
        In1 = [0.0] * (num1 + 1)
        In2 = [0.0] * (num1 + 1)
        Mn1 = [0.0] * (num1 + 1)
        Mn2 = [0.0] * (num1 + 1)
        K = [0] * (num2 + 1)
        D1[1] = I1[1] = M1[1] = K[1] = b
        key = 2
        kkk = 0
        lll = 0
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
                    if(s1[j - 1]==s2[i]):
                        kkk = -1.0
                        lll = 1.0
                    else:
                        kkk = 2.0
                        lll = 0
                    if(D1[j-1] < I1[j-1] and D1[j-1] < M1[j - 1]):
                        D2[j] = D1[j-1] + kkk
                        Dnum2[j] = Dnum1[j-1] + lll
                        Dn2[j] = Dn1[j-1] + 1.0
                    elif (I1[j-1] < M1[j - 1]):
                        D2[j] = I1[j-1] + kkk
                        Dnum2[j] = Inum1[j-1] + lll
                        Dn2[j] = In1[j-1] + 1.0
                    else:
                        D2[j] = M1[j-1] + kkk
                        Dnum2[j] = Mnum1[j-1] + lll
                        Dn2[j] = Mn1[j-1] + 1.0
                    #########
                    if(I1[j] + a < M1[j] + b and I1[j] + a < D1[j] + b ):
                        I2[j] = I1[j] + a
                        Inum2[j] = Inum1[j]
                        In2[j] = In1[j] + 1.0
                    elif(M1[j]  < D1[j]  ):
                        I2[j] = M1[j] + b
                        Inum2[j] = Mnum1[j]
                        In2[j] = Mn1[j] + 1.0
                    else:
                        I2[j] = D1[j] + b
                        Inum2[j] = Dnum1[j]
                        In2[j] = Dn1[j] + 1.0
                    ##########
                    if(M2[j - 1] + a < I2[j - 1] + b and M2[j - 1] + a < D2[j - 1] + b):
                        M2[j] = M2[j - 1] + a
                        Mnum2[j] = Mnum2[j-1]
                        Mn2[j] = Mn2[j-1] + 1.0
                    if(I2[j - 1] < D2[j - 1] ):
                        M2[j] = I2[j - 1] + a
                        Mnum2[j] = Inum2[j-1]
                        Mn2[j] = In2[j-1] + 1.0
                    else:
                        M2[j] = D2[j - 1] + a
                        Mnum2[j] = Dnum2[j-1]
                        Mn2[j] = Dn2[j-1] + 1.0
                key = 1
            else:
                D1[0] = K[i + 1]
                I1[0] = K[i + 1]
                M1[0] = K[i + 1]
                for j in range(1, num1 + 1):
                    if(s1[j - 1]==s2[i]):
                        kkk = -1.0
                        lll = 1.0
                    else:
                        kkk = 2.0
                        lll = 0
                    if(D2[j-1] < I2[j-1] and D2[j-1] < M2[j - 1]):
                        D1[j] = D2[j-1] + kkk
                        Dnum1[j] = Dnum2[j-1] + lll
                        Dn1[j] = Dn2[j-1] + 1.0
                    elif (I2[j-1] < M2[j - 1]):
                        D1[j] = I2[j-1] + kkk
                        Dnum1[j] = Inum2[j-1] + lll
                        Dn1[j] = In2[j-1] + 1.0
                    else:
                        D1[j] = M2[j-1] + kkk
                        Dnum1[j] = Mnum2[j-1] + lll
                        Dn1[j] = Mn2[j-1] + 1.0
                    #########
                    if(I2[j] + a < M2[j] + b and I2[j] + a < D2[j] + b ):
                        I1[j] = I2[j] + a
                        Inum1[j] = Inum2[j]
                        In1[j] = In2[j] + 1.0
                    elif(M2[j]  < D2[j]  ):
                        I1[j] = M2[j] + b
                        Inum1[j] = Mnum2[j]
                        In1[j] = Mn2[j] + 1.0
                    else:
                        I1[j] = D2[j] + b
                        Inum1[j] = Dnum2[j]
                        In1[j] = Dn2[j] + 1.0
                    ##########
                    if(M1[j - 1] + a < I1[j - 1] + b and M1[j - 1] + a < D1[j - 1] + b):
                        M1[j] = M1[j - 1] + a
                        Mnum1[j] = Mnum1[j-1]
                        Mn1[j] = Mn1[j-1] + 1.0
                    if(I1[j - 1] < D1[j - 1] ):
                        M1[j] = I1[j - 1] + a
                        Mnum1[j] = Inum1[j-1]
                        Mn1[j] = In1[j-1] + 1.0
                    else:
                        M1[j] = D1[j - 1] + a
                        Mnum1[j] = Dnum1[j-1]
                        Mn1[j] = Dn1[j-1] + 1.0
                key = 2
        if (key == 2):
            M2[num1] = M1[num1]
            D2[num1] = D1[num1]
            I2[num1] = I1[num1]
            Mnum2[num1] = Mnum1[num1]
            Dnum2[num1] = Dnum1[num1]
            Inum2[num1] = Inum1[num1]
            Mn2[num1] = Mn1[num1]
            Dn2[num1] = Dn1[num1]
            In2[num1] = In1[num1]
        if(M2[num1] < D2[num1] and M2[num1] < I2[num1]):
            return 1 - Mnum2[num1]/Mn2[num1]
        elif(D2[num1] < I2[num1]):
            return 1 - Dnum2[num1]/Dn2[num1]
        else:
            return 1 - Inum2[num1]/In2[num1]

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
            node = [[], diff[di][num], diff[dj][num], [], [], -1, -1, 0]
            diffk = []

            for i in range(0, num):
                if(i != di and i != dj):
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
        node = [[], diff[di][3], diff[dj][3], [],[], -1, -1, 0]
        dis = (diff[di][dk] + diff[dj][dk])/2
        diff[dk][3][5]=(dis)
        node[5]=(dis - 0.5*diff[di][dj])
        t = [[],node,diff[dk][3],[],[],-1,-1,0]






        self.sizeoftree(t)
        self.weighttree(t, 1, -1.0)
        self.tree = t
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

    def complete(self):
        root = self.tree
        self.dfs(root)
        #print root[0]
        #self.refine()
        '''f = open("424clustal.txt", "w")
        for seq in root[3]:
            print >>f, seq'''

    def dfs(self, node):
        if not node[1][3]:
            self.dfs(node[1])
        if not node[2][3]:
            self.dfs(node[2])
        seqs = self.multipairwise(node[1], node[2])
        node[3] = seqs
        node[4] = node[1][4] + node[2][4]
        node[0] = node[1][0] + node[2][0]

    # [num, left, right, seqs, weight, distance, height]
    '''def multipairwise(self, node1, node2):
        a = -self.a
        b = -self.b
        len1 = len(node1[3][0])
        len2 = len(node2[3][0])
        num_pair = len(node1[3]) * len(node2[3])
        M = []
        D = []
        I = []
        Mp = []
        Dp = []
        Ip = []
        for i in range(len2+1):
            M.append([0]*(len1+1))
            D.append([0]*(len1+1))
            I.append([0]*(len1+1))
            Mp.append([0]*(len1+1))
            Dp.append([0]*(len1+1))
            Ip.append([0]*(len1+1))
        M[0][0] = D[0][0] = I[0][0] = 0
        M[1][0] = D[1][0] = I[1][0] = b
        M[0][1] = D[0][1] = I[0][1] = b
        Mp[1][0] = Dp[1][0] = Ip[1][0] = 1 # 1 is vertical
        Mp[0][1] = Dp[0][1] = Ip[0][1] = -1 # -1 is horizontal
        for i in range(2,len2+1):
            M[i][0] = D[i][0] = I[i][0] = I[i-1][0] + a
            Mp[i][0] = Dp[i][0] = Ip[i][0] = 1 # 1 is vertical
        for i in range(2,len1+1):
            M[0][i] = D[0][i] = I[0][i] = I[0][i-1] + a
            Mp[0][i] = Dp[0][i] = Ip[0][i] = -1 # -1 is horizontal

        weight_sum = 0.0
        for k1 in range(len(node1[3])):
            for k2 in range(len(node2[3])):
                weight_sum += node1[4][k1] * node2[4][k2]
        gap1 = []
        gap2 = []
        for k1 in range(len(node1[3][0])):
            all_gap = True
            for i in range(len(node1[3])):
                if node1[3][i][k1] != '-':
                    all_gap = False
                    break
            gap1.append(all_gap)
        for k1 in range(len(node2[3][0])):
            all_gap = True
            for i in range(len(node2[3])):
                if node2[3][i][k1] != '-':
                    all_gap = False
                    break
            gap2.append(all_gap)
        weight_sum /= num_pair
        pena = weight_sum * a
        penb = weight_sum * b
        for i in range(1, len2 + 1):
            for j in range(1, len1 + 1):
                score_d = 0.0
                for k1 in range(len(node1[3])):
                    for k2 in range(len(node2[3])):
                        score_d += compare(node1[3][k1][j - 1], node2[3][k2][i - 1]) * node1[4][k1] * node2[4][k2]
                score_d /= num_pair
                if M[i-1][j-1] > D[i-1][j-1]:
                    if I[i-1][j-1] > M[i-1][j-1]:
                        M[i][j] = I[i-1][j-1] + score_d
                        Mp[i][j] = 1
                    else:
                        M[i][j] = M[i-1][j-1] + score_d
                else:
                    if I[i-1][j-1] > D[i-1][j-1]:
                        M[i][j] = I[i-1][j-1] + score_d
                        Mp[i][j] = 1
                    else:
                        M[i][j] = D[i-1][j-1] + score_d
                        Mp[i][j] = -1
                if M[i][j-1] > D[i][j-1]:
                    if I[i][j-1] + pena > M[i][j-1] + penb:
                        I[i][j] = I[i][j-1] + pena
                        Ip[i][j] = 1
                    else:
                        I[i][j] =  M[i][j-1] + penb
                else:
                    if I[i][j-1] + pena > D[i][j-1] + penb:
                        I[i][j] = I[i][j-1] + pena
                        Ip[i][j] = 1
                    else:
                        I[i][j] = D[i][j-1] + penb
                        Ip[i][j] = -1
                if M[i-1][j] > I[i-1][j] :
                    if D[i-1][j] + pena > M[i-1][j] + penb:
                        D[i][j] = D[i-1][j] + pena
                        Dp[i][j] = -1
                    else:
                        D[i][j] = M[i-1][j] + penb
                else:
                    if D[i-1][j] + pena > I[i-1][j] + penb:
                        D[i][j] = D[i-1][j] + pena
                        Dp[i][j] = -1
                    else:
                        D[i][j] = I[i-1][j] + penb
                        Dp[i][j] = 1
        if M[len2][len1] > I[len2][len1]:
            if M[len2][len1] > D[len2][len1]:
                key = 0
            else:
                key = -1
        else:
            if I[len2][len1] > D[len2][len1]:
                key = 1
            else:
                key = -1
        #need to refine the sequence
        aligned = []
        seqnum1 = len(node1[3])
        seqnum2 = len(node2[3])
        points = [0] * (seqnum1 + seqnum2)
        for i in range(seqnum1 + seqnum2):
            aligned.append([])
        i = len2
        j = len1
        while i != 0 and j != 0:
            if key == 0:
                key = Mp[i][j]
                i -= 1
                j -= 1
                for k in range(seqnum1):
                    aligned[k].append(node1[3][k][j])
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append(node2[3][k][i])
            elif key == 1:
                key = Ip[i][j]
                j -= 1
                for k in range(seqnum1):
                    aligned[k].append(node1[3][k][j])
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append('-')
            else:
                key = Dp[i][j]
                i -= 1
                for k in range(seqnum1):
                    aligned[k].append('-')
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append(node2[3][k][i])
        while i != 0:
            i -= 1
            for k in range(seqnum1):
                aligned[k].append('-')
            for k in range(seqnum2):
                aligned[seqnum1 + k].append(node2[3][k][i])
        while j != 0:
            j -= 1
            for k in range(seqnum1):
                aligned[k].append(node1[3][k][j])
            for k in range(seqnum2):
                aligned[seqnum1 + k].append('-')
        for i in range(seqnum1):
            aligned[i] = aligned[i][::-1]
        for i in range(seqnum2):
            aligned[seqnum1 + i] = aligned[seqnum1 + i][::-1]
        return aligned'''

    def multipairwise(self, node1, node2):
        a = -self.a
        b = -self.b
        len1 = len(node1[3][0])
        len2 = len(node2[3][0])
        num_pair = len(node1[3]) * len(node2[3])
        M = []
        D = []
        I = []
        Mp = []
        Dp = []
        Ip = []
        for i in range(len2+1):
            M.append([0]*(len1+1))
            D.append([0]*(len1+1))
            I.append([0]*(len1+1))
            Mp.append([0]*(len1+1))
            Dp.append([0]*(len1+1))
            Ip.append([0]*(len1+1))
        M[0][0] = D[0][0] = I[0][0] = 0
        M[1][0] = D[1][0] = I[1][0] = b
        M[0][1] = D[0][1] = I[0][1] = b
        Mp[1][0] = Dp[1][0] = Ip[1][0] = 1 # 1 is vertical
        Mp[0][1] = Dp[0][1] = Ip[0][1] = -1 # -1 is horizontal
        for i in range(2,len2+1):
            M[i][0] = D[i][0] = I[i][0] = I[i-1][0] + a
            Mp[i][0] = Dp[i][0] = Ip[i][0] = 1 # 1 is vertical
        for i in range(2,len1+1):
            M[0][i] = D[0][i] = I[0][i] = I[0][i-1] + a
            Mp[0][i] = Dp[0][i] = Ip[0][i] = -1 # -1 is horizontal

        weight_sum = 0.0
        for k1 in range(len(node1[3])):
            for k2 in range(len(node2[3])):
                weight_sum += node1[4][k1] * node2[4][k2]
        gap1 = []
        gap2 = []
        for k1 in range(len(node1[3][0])):
            all_gap = True
            for i in range(len(node1[3])):
                if node1[3][i][k1] != '-':
                    all_gap = False
                    break
            gap1.append(all_gap)
        for k1 in range(len(node2[3][0])):
            all_gap = True
            for i in range(len(node2[3])):
                if node2[3][i][k1] != '-':
                    all_gap = False
                    break
            gap2.append(all_gap)
        weight_sum /= num_pair
        pena = weight_sum * a
        penb = weight_sum * b
        for i in range(1, len2 + 1):
            for j in range(1, len1 + 1):
                score_d = 0.0
                for k1 in range(len(node1[3])):
                    for k2 in range(len(node2[3])):
                        score_d += compare(node1[3][k1][j - 1], node2[3][k2][i - 1]) * node1[4][k1] * node2[4][k2]
                score_d /= num_pair
                if M[i-1][j-1] > D[i-1][j-1]:
                    if I[i-1][j-1] > M[i-1][j-1]:
                        M[i][j] = I[i-1][j-1] + score_d
                        Mp[i][j] = 1
                    else:
                        M[i][j] = M[i-1][j-1] + score_d
                else:
                    if I[i-1][j-1] > D[i-1][j-1]:
                        M[i][j] = I[i-1][j-1] + score_d
                        Mp[i][j] = 1
                    else:
                        M[i][j] = D[i-1][j-1] + score_d
                        Mp[i][j] = -1
                if gap1[j-1]:
                    if M[i][j-1] > D[i][j-1]:
                        if I[i][j-1] > M[i][j-1]:
                            I[i][j] = I[i][j-1]
                            Ip[i][j] = 1
                        else:
                            I[i][j] =  M[i][j-1]
                    else:
                        if I[i][j-1] > D[i][j-1]:
                            I[i][j] = I[i][j-1]
                            Ip[i][j] = 1
                        else:
                            I[i][j] = D[i][j-1]
                            Ip[i][j] = -1
                else:
                    if M[i][j-1] > D[i][j-1]:
                        if I[i][j-1] + pena > M[i][j-1] + penb:
                            I[i][j] = I[i][j-1] + pena
                            Ip[i][j] = 1
                        else:
                            I[i][j] =  M[i][j-1] + penb
                    else:
                        if I[i][j-1] + pena > D[i][j-1] + penb:
                            I[i][j] = I[i][j-1] + pena
                            Ip[i][j] = 1
                        else:
                            I[i][j] = D[i][j-1] + penb
                            Ip[i][j] = -1
                if gap2[i-1]:
                    if M[i-1][j] > I[i-1][j] :
                        if D[i-1][j] > M[i-1][j]:
                            D[i][j] = D[i-1][j]
                            Dp[i][j] = -1
                        else:
                            D[i][j] = M[i-1][j]
                    else:
                        if D[i-1][j] + pena > I[i-1][j]:
                            D[i][j] = D[i-1][j]
                            Dp[i][j] = -1
                        else:
                            D[i][j] = I[i-1][j]
                            Dp[i][j] = 1
                else:
                    if M[i-1][j] > I[i-1][j] :
                        if D[i-1][j] + pena > M[i-1][j] + penb:
                            D[i][j] = D[i-1][j] + pena
                            Dp[i][j] = -1
                        else:
                            D[i][j] = M[i-1][j] + penb
                    else:
                        if D[i-1][j] + pena > I[i-1][j] + penb:
                            D[i][j] = D[i-1][j] + pena
                            Dp[i][j] = -1
                        else:
                            D[i][j] = I[i-1][j] + penb
                            Dp[i][j] = 1
        if M[len2][len1] > I[len2][len1]:
            if M[len2][len1] > D[len2][len1]:
                key = 0
            else:
                key = -1
        else:
            if I[len2][len1] > D[len2][len1]:
                key = 1
            else:
                key = -1
        #need to refine the sequence
        aligned = []
        seqnum1 = len(node1[3])
        seqnum2 = len(node2[3])
        points = [0] * (seqnum1 + seqnum2)
        for i in range(seqnum1 + seqnum2):
            aligned.append([])
        i = len2
        j = len1
        while i != 0 and j != 0:
            if key == 0:
                key = Mp[i][j]
                i -= 1
                j -= 1
                for k in range(seqnum1):
                    aligned[k].append(node1[3][k][j])
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append(node2[3][k][i])
            elif key == 1:
                key = Ip[i][j]
                j -= 1
                for k in range(seqnum1):
                    aligned[k].append(node1[3][k][j])
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append('-')
            else:
                key = Dp[i][j]
                i -= 1
                for k in range(seqnum1):
                    aligned[k].append('-')
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append(node2[3][k][i])
        while i != 0:
            i -= 1
            for k in range(seqnum1):
                aligned[k].append('-')
            for k in range(seqnum2):
                aligned[seqnum1 + k].append(node2[3][k][i])
        while j != 0:
            j -= 1
            for k in range(seqnum1):
                aligned[k].append(node1[3][k][j])
            for k in range(seqnum2):
                aligned[seqnum1 + k].append('-')
        for i in range(seqnum1):
            aligned[i] = aligned[i][::-1]
        for i in range(seqnum2):
            aligned[seqnum1 + i] = aligned[seqnum1 + i][::-1]
        return aligned

    def sortdistance(self, node):
        if node[1]:
            self.sortdistance(node[1])
        if node[2]:
            self.sortdistance(node[2])
        index = 0
        for i in range(len(self.distance)):
            if (self.distance[i][5] < node[5]):
                break
            else :
                index += 1
        self.distance.insert(index, node)

    def sc(self, seq1,seq2):
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

    def evalue(self, seqs):
        num = len(seqs)
        summ = 0
        score = []
        for i in range(num):
            score.append([0]*num)

        for i in range(num):
            for j in range(i):
                score[i][j] = self.sc(seqs[i], seqs[j])
                score[j][i] = self.sc(seqs[i], seqs[j])
                summ += score[i][j]
        #print score,"\n\n"
        return summ

    def refine(self, times):
        root = self.tree
        self.sortdistance(root)
        for i in range(len(self.distance)):
            print self.distance[i][0], self.distance[i][5]
        max_score = self.evalue(root[3])
        print max_score
        for i in range(times):
            node1 = self.distance[i]
            node2 = [[], [], [], [], [], -1, -1, 0]
            for j in range(len(root[0])):
                if j not in node1[0]:
                    node2[0].append(j)
                    node2[3].append(root[3][root[0].index(j)])
                    node2[4].append(root[4][root[0].index(j)])
            aligned = self.multipairwise(node1, node2)
            score = self.evalue(aligned)
            if (score < max_score):
                max_score = score
                root[3] = aligned
                root[0] = node1[0] + node2[0]
                root[4] = node1[4] + node2[4]
            print max_score
        f = open("424clustal2.txt", "w")
        for seq in root[3]:
            print >>f, seq


# class tree:
if __name__ == "__main__":
    fp = open("good3.txt","r")
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
    fp.close()
    mul = sequences(seq[0:20], 2.0, 3.0, 1.0)
    t1 = time.clock()
    mul.init_matrix()

    t2 = time.clock()
    mul.buildtree()

    t3 = time.clock()
    mul.complete()

    t4 = time.clock()
    mul.refine(10)
    t5 = time.clock()

    print "%f, %f, %f, %f" % (t2-t1, t3-t2, t4-t3, t5-t4)
