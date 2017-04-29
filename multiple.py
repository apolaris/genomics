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
        self.score2 = []
        for i in range(0, self.num):
            self.score.append([0.0] * self.num)
            self.score2.append([0.0] * self.num)
        for i in range(0, self.num):
            # num, left, right, seqs, weight, distance from parent, height,
            # size
            node = [[i], [], [], [self.seqs[i]], [], -1, -1, 0]
            self.score[i].append(node)
            for j in range(0, i):
                print i,j
                self.score[i][j] = self.pairwise2(i, j)
                self.score2[i][j] = self.score2[j][i] = self.score[i][j]
                self.score[j][i] = self.score[i][j]
        self.distance = []

    # compare the two sequences
    def pairwise(self, i, j):
        s1 = self.seqs[i]
        s2 = self.seqs[j]
        len1 = len(s1)
        len2 = len(s2)
        M = []
        D = []
        I = []
        Mp = []
        Dp = []
        Ip = []
        a = self.a
        b = self.b
        for i in range(len2+1):
            M.append([0]*(len1+1))
            D.append([0]*(len1+1))
            I.append([0]*(len1+1))
            Mp.append([0]*(len1+1))
            Dp.append([0]*(len1+1))
            Ip.append([0]*(len1+1))
        M[0][0] = D[0][0] = I[0][0] = 0
        M[1][0] = D[1][0] = I[1][0] = 3
        M[0][1] = D[0][1] = I[0][1] = 3
        Mp[1][0] = Dp[1][0] = Ip[1][0] = 1 # 1 is vertical
        Mp[0][1] = Dp[0][1] = Ip[0][1] = -1 # -1 is horizontal
        for i in range(2,len2+1):
            M[i][0] = D[i][0] = I[i][0] = I[i-1][0] + 2
            Mp[i][0] = Dp[i][0] = Ip[i][0] = 1
        for i in range(2,len1+1):
            M[0][i] = D[0][i] = I[0][i] = I[0][i-1] + 2
            Mp[0][i] = Dp[0][i] = Ip[0][i] = -1
        for i in range(1,len2+1):
            for j in range(1,len1+1):
                if(M[i-1][j-1] <= D[i-1][j-1]):
                    if(I[i-1][j-1] < M[i-1][j-1]):
                        M[i][j] = I[i-1][j-1] + match(s1[j-1], s2[i-1])
                        Mp[i][j] = 1
                    else:
                        M[i][j] = M[i-1][j-1] + match(s1[j-1], s2[i-1])
                else:
                    if(I[i-1][j-1] < D[i-1][j-1]):
                        M[i][j] = I[i-1][j-1] + match(s1[j-1], s2[i-1])
                        Mp[i][j] = 1
                    else:
                        M[i][j] = D[i-1][j-1] + match(s1[j-1], s2[i-1])
                        Mp[i][j] = -1
                if(M[i][j-1] <= D[i][j-1]):
                    if(I[i][j-1] + a < M[i][j-1] + b):
                        I[i][j] = I[i][j-1] + a
                        Ip[i][j] = 1
                    else:
                        I[i][j] =  M[i][j-1] + b
                else:
                    if(I[i][j-1] + a < D[i][j-1] + b):
                        I[i][j] = I[i][j-1] + a
                        Ip[i][j] = 1
                    else:
                        I[i][j] = D[i][j-1] + b
                        Ip[i][j] = -1
                if(M[i-1][j] <= I[i-1][j]):
                    if(D[i-1][j] + a < M[i-1][j] + b):
                        D[i][j] = D[i-1][j] + a
                        Dp[i][j] = -1
                    else:
                        D[i][j] = M[i-1][j] + b
                else:
                    if(D[i-1][j] + a < I[i-1][j] + b):
                        D[i][j] = D[i-1][j] + a
                        Dp[i][j] = -1
                    else:
                        D[i][j] = I[i-1][j] + b
                        Dp[i][j] = 1
        i = len2
        j = len1
        if(M[len2][len1] < I[len2][len1]):
            if(M[len2][len1] < D[len2][len1]):
                key = 0
            else:
                key = -1
        else:
            if(I[len2][len1] < D[len2][len1]):
                key = 1
            else:
                key = -1
        align = [[],[]]
        while(i != 0 and j != 0):
            if(key == 0):
                key = Mp[i][j]
                i-=1
                j-=1
                align[0].append(s1[j])
                align[1].append(s2[i])
            elif(key == 1 ):
                key = Ip[i][j]
                j-=1
                align[0].append(s1[j])
                align[1].append('-')
            else:
                key = Dp[i][j]
                i-=1
                align[0].append('-')
                align[1].append(s2[i])
        while(i!=0):
            i-=1
            align[0].append('-')
            align[1].append(s2[i])
        while(j!=0):
            j-=1
            align[0].append(s1[j])
            align[1].append('-')
        leng = len(align[0])
        same = 0.0
        for i in range(leng):
            if(align[0][i] == align[1][i]):
                same += 1
        return same/leng

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
            d = -1
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
            #print diff[di][dj],num,di,dj
            diff[di][num][5] = 0.5 * diff[di][dj] + (su[di] - su[dj]) / (2 * (num -2))
            diff[dj][num][5] = diff[di][dj] - diff[di][num][5]
            #print su[di] - su[dj], diff[di][num][5], diff[dj][num][5]
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
        print root[4]
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
                weight_sum += self.score2[node1[0][k1]][node2[0][k2]] * 0.5 + 0.5
                #weight_sum += node1[4][k1] * node2[4][k2]
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
        #pena = weight_sum * a
        #penb = weight_sum * b
        pena = a
        penb = b
        for i in range(1, len2 + 1):
            for j in range(1, len1 + 1):
                score_d = 0.0
                for k1 in range(len(node1[3])):
                    for k2 in range(len(node2[3])):
                        #score_d += compare(node1[3][k1][j - 1], node2[3][k2][i - 1]) * node1[4][k1] * node2[4][k2]
                        #score_d += compare(node1[3][k1][j - 1], node2[3][k2][i - 1]) * (self.score2[node1[0][k1]][node2[0][k2]] * 0.5 + 0.5)
                        score_d += compare(node1[3][k1][j - 1], node2[3][k2][i - 1])
                score_d /= num_pair
                if M[i-1][j-1] >= D[i-1][j-1]:
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
                if M[i][j-1] >= D[i][j-1]:
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
                if M[i-1][j] >= I[i-1][j] :
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
        if M[len2][len1] >= I[len2][len1]:
            if M[len2][len1] >= D[len2][len1]:
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
        return aligned'''

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
        stop = 0
        for i in range(len(self.distance)):
            print self.distance[i][0], self.distance[i][5]
        max_score = self.evalue(root[3])
        print max_score
        for i in range(len(self.seqs)):
            node1 = self.distance[i]
            node2 = [[], [], [], [], [], -1, -1, 0]
            if stop >= 5:
                break
            if node1 == root[1] or node1 == root[2]:
                continue
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
                stop = 0
            else:
                stop += 1
            print max_score
        index = 0
        for k1 in range(len(root[3][0])):
            all_gap = True
            for i in range(len(root[3])):
                if root[3][i][index] != '-':
                    all_gap = False
                    break
            if all_gap:
                for i in range(len(root[3])):
                    del root[3][i][index]
            else:
                index += 1
        f = open("428clustal4.txt", "w")
        for seq in root[3]:
            print >>f, seq

    def pairwise2(self, i, j):
        s1 = self.seqs[i]
        s2 = self.seqs[j]
        len1 = len(s1)
        len2 = len(s2)
        M = []
        D = []
        I = []
        Mp = []
        Dp = []
        Ip = []
        a = self.a
        b = self.b
        for i in range(len2+1):
            M.append([100000000]*(len1+1))
            D.append([100000000]*(len1+1))
            I.append([100000000]*(len1+1))
            Mp.append([0]*(len1+1))
            Dp.append([0]*(len1+1))
            Ip.append([0]*(len1+1))
        M[0][0] = D[0][0] = I[0][0] = 0
        M[1][0] = D[1][0] = I[1][0] = 3
        M[0][1] = D[0][1] = I[0][1] = 3
        Mp[1][0] = Dp[1][0] = Ip[1][0] = 1 # 1 is vertical
        Mp[0][1] = Dp[0][1] = Ip[0][1] = -1 # -1 is horizontal
        for i in range(2,len2+1):
            M[i][0] = D[i][0] = I[i][0] = I[i-1][0] + 2
            Mp[i][0] = Dp[i][0] = Ip[i][0] = 1
        for i in range(2,len1+1):
            M[0][i] = D[0][i] = I[0][i] = I[0][i-1] + 2
            Mp[0][i] = Dp[0][i] = Ip[0][i] = -1
        quart = len1/4
        begin = 0
        end = 0
        for i in range(1,len2+1):
            begin = i - quart
            end = i + quart
            if(begin <= 0):
                begin = 1
            if(end > len1):
                end = len1+1
            for j in range(begin,end):
                if(M[i-1][j-1] <= D[i-1][j-1]):
                    if(I[i-1][j-1] < M[i-1][j-1]):
                        M[i][j] = I[i-1][j-1] + match(s1[j-1], s2[i-1])
                        Mp[i][j] = 1
                    else:
                        M[i][j] = M[i-1][j-1] + match(s1[j-1], s2[i-1])
                else:
                    if(I[i-1][j-1] < D[i-1][j-1]):
                        M[i][j] = I[i-1][j-1] + match(s1[j-1], s2[i-1])
                        Mp[i][j] = 1
                    else:
                        M[i][j] = D[i-1][j-1] + match(s1[j-1], s2[i-1])
                        Mp[i][j] = -1
                if(M[i][j-1] <= D[i][j-1]):
                    if(I[i][j-1] + a < M[i][j-1] + b):
                        I[i][j] = I[i][j-1] + a
                        Ip[i][j] = 1
                    else:
                        I[i][j] =  M[i][j-1] + b
                else:
                    if(I[i][j-1] + a < D[i][j-1] + b):
                        I[i][j] = I[i][j-1] + a
                        Ip[i][j] = 1
                    else:
                        I[i][j] = D[i][j-1] + b
                        Ip[i][j] = -1
                if(M[i-1][j] <= I[i-1][j]):
                    if(D[i-1][j] + a < M[i-1][j] + b):
                        D[i][j] = D[i-1][j] + a
                        Dp[i][j] = -1
                    else:
                        D[i][j] = M[i-1][j] + b
                else:
                    if(D[i-1][j] + a < I[i-1][j] + b):
                        D[i][j] = D[i-1][j] + a
                        Dp[i][j] = -1
                    else:
                        D[i][j] = I[i-1][j] + b
                        Dp[i][j] = 1
        i = len2
        j = len1
        if(M[len2][len1] < I[len2][len1]):
            if(M[len2][len1] < D[len2][len1]):
                key = 0
            else:
                key = -1
        else:
            if(I[len2][len1] < D[len2][len1]):
                key = 1
            else:
                key = -1
        align = [[],[]]
        while(i != 0 and j != 0):
            if(key == 0):
                key = Mp[i][j]
                i-=1
                j-=1
                align[0].append(s1[j])
                align[1].append(s2[i])
            elif(key == 1 ):
                key = Ip[i][j]
                j-=1
                align[0].append(s1[j])
                align[1].append('-')
            else:
                key = Dp[i][j]
                i-=1
                align[0].append('-')
                align[1].append(s2[i])
        while(i!=0):
            i-=1
            align[0].append('-')
            align[1].append(s2[i])
        while(j!=0):
            j-=1
            align[0].append(s1[j])
            align[1].append('-')
        leng = len(align[0])
        same = 0.0
        for i in range(leng):
            if(align[0][i] == align[1][i]):
                same += 1
        return same/leng

    def multipairwise2(self, node1, node2):
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
            M.append([-10000000]*(len1+1))
            D.append([-10000000]*(len1+1))
            I.append([-10000000]*(len1+1))
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
        begin = 0
        end = 1
        quart = len1/4

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
        #pena = weight_sum * a
        #penb = weight_sum * b
        pena = a
        penb = b
        for i in range(1, len2 + 1):
            begin = i - quart
            end = i + quart
            if(begin <= 0):
                begin = 1
            if(end > len1):
                end = len1+1
            for j in range(begin, end):
                score_d = 0.0
                for k1 in range(len(node1[3])):
                    for k2 in range(len(node2[3])):
                        #score_d += compare(node1[3][k1][j - 1], node2[3][k2][i - 1]) * node1[4][k1] * node2[4][k2]
                        score_d += compare(node1[3][k1][j - 1], node2[3][k2][i - 1])
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
        return aligned


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
    mul = sequences(seq[0:50], 2.0, 3.0, 1.0)
    t1 = time.clock()
    mul.init_matrix()

    t2 = time.clock()
    mul.buildtree()

    t3 = time.clock()
    mul.complete()

    t4 = time.clock()
    mul.refine(50)
    t5 = time.clock()

    print "%f, %f, %f, %f" % (t2-t1, t3-t2, t4-t3, t5-t4)
