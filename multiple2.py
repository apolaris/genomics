import time

class sequences:
    'sequence list for alignment'

    def __init__(self, seqs, a, b, match, mismatch, write_name):
        self.seqs = seqs
        self.num = len(seqs)
        self.a = a
        self.b = b
        self.mat = match
        self.mismat = mismatch
        self.write_name = write_name

    def match(self, a, b):
        if a == b:
            return self.mat
        else:
            return self.mismat

    def compare(self, a, b):
        if a == '-' or b == '-':
            return 0
        elif a == b:
            return -self.mat
        else:
            return -self.mismat

    # do all
    def all(self):
        print "begin init"
        self.init_matrix()
        print "begin build"
        self.build()
        print "begin align"
        self.complete()

    # calculate the difference
    def init_matrix(self):
        self.score = []
        self.ext = []
        self.distance = []
        for i in range(0, self.num):
            self.score.append([0.0] * self.num)
            tmplist = []
            for j in range(self.num):
                tmplist.append([])
            self.ext.append(tmplist)

        for i in range(self.num):
            # num, left, right, seqs, size, distance, realindex
            node = [[i], [], [], [self.seqs[i]], 0, -1,[range(len(self.seqs[i]))]]
            self.score[i].append(node)
            for j in range(i):
                align = self.pairwise2(i, j)
                align2 =[align[1], align[0]]
                leng = len(align[0])
                matc = 0.0
                for k in range(leng):
                    if(align[0][k] == align[1][k]):
                        matc += 1.0
                self.score[i][j] = 1 - matc/leng
                self.score[j][i] = self.score[i][j]
                self.ext[i][j].append(100 * matc/leng)
                self.ext[i][j].append(align)
                kkk1 = []
                kkk2 = []
                for l in range(len(self.seqs[i])):
                    kkk1.append({})
                for l in range(len(self.seqs[j])):
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
        t2 = time.clock()
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
                    min_ext = min(self.ext[i][s][0], self.ext[j][s][0])
                    lens = len(self.seqs[s])-1
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
                                self.ext[i][j][2][numisi][numjsj] += min_ext
                            else:
                                self.ext[i][j][2][numisi][numjsj] = min_ext

                            if (self.ext[i][j][3][numjsj].has_key(numisi)):
                                self.ext[i][j][3][numjsj][numisi] += min_ext
                            else:
                                self.ext[i][j][3][numjsj][numisi] = min_ext
                        if(numiss >= lens):
                            break
        for i in range(self.num):
            for j in range(i):
                self.ext[i][j].pop(0)
                self.ext[i][j].pop(0)
                self.ext[j][i].pop(0)
                self.ext[j][i].pop(0)
                self.ext[j][i].append(self.ext[i][j][1])
                self.ext[j][i].append(self.ext[i][j][0])

        t3 = time.clock()
        print "extension time: ", t3-t2

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

            diff[di][num][5] = 0.5 * diff[di][dj] + (su[di] - su[dj])/ (2 * (num -2))
            diff[dj][num][5] = diff[di][dj] - diff[di][num][5]
            node = [[], diff[di][num], diff[dj][num], [], 0, -1, []]
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
        node = [[], diff[di][3], diff[dj][3], [], 0, -1, []]
        dis = (diff[di][dk] + diff[dj][dk])/2
        diff[dk][3][5]=(dis)
        node[5]=(dis - 0.5*diff[di][dj])
        t = [[],node,diff[dk][3],[],[],-1,-1,0]

        self.tree = t
        return t


    def complete(self):
        root = self.tree
        self.dfs(root)


    def dfs(self, node):
        if not node[1][3]:
            self.dfs(node[1])
        if not node[2][3]:
            self.dfs(node[2])
        node[3], node[6] = self.multipairwise2(node[1], node[2])
        node[0] = node[1][0] + node[2][0]
        node[4] = node[1][4] + node[2][4]

    def showExt (self, extlist, i, j):
        if i == -1:
            return 0
        else:
            dic = extlist[i]
            if dic.has_key(j):
                return dic[j]
            else:
                return 0

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
            score += self.mat
        elif (s1[0] == '-' or s2[0] == '-'):
            score += self.b
        else:
            score += self.a

        l = len(s1)
        for i in range(1,l):
            if(s1[i] == s2[i]):
                score += self.mismat

            elif(s1[i] == '-'):
                if(s1[i-1] == '-'):
                    score += self.a
                else:
                    score += self.b
            elif(s2[i] == '-'):
                if(s2[i-1] == '-'):
                    score += self.a
                else:
                    score += self.b
            else:
                score += self.mismat
        return score

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
        return summ

    def refine(self, times):
        # num, left, right, seqs, size, distance, realindex
        root = self.tree
        self.sortdistance(root)
        stop = 0
        t1 = time.clock()
        max_score = self.evalue(root[3])
        print "the SP score:", max_score
        for i in range(times):
            node1 = self.distance[i]
            node2 = [[], [], [], [], -1, -1, []]
            if node1 == root[1] or node1 == root[2]:
                continue
            for j in range(len(root[0])):
                if j not in node1[0]:
                    node2[0].append(j)
                    node2[3].append(root[3][root[0].index(j)])
                    node2[6].append(root[6][root[0].index(j)])
            aligned, realIndex = self.multipairwise2(node1, node2)
            score = self.evalue(aligned)
            if (score < max_score):
                max_score = score
                root[3] = aligned
                root[0] = node1[0] + node2[0]
                root[6] = realIndex
                stop = 0
            else:
                stop += 1
            print "the SP score:", max_score
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
        f = open(self.write_name, "w")
        for seq in root[3]:
            for c in seq:
                f.write(c)
            f.write('\n')


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
        slope = (float(len1))/(len2)
        for i in range(1,len2+1):
            i_slope = int(i*slope)
            begin = i_slope - quart
            end = i_slope + quart
            if(begin <= 0):
                begin = 1
            if(end > len1):
                end = len1+1
            for j in range(begin,end):
                if(M[i-1][j-1] <= D[i-1][j-1]):
                    if(I[i-1][j-1] < M[i-1][j-1]):
                        M[i][j] = I[i-1][j-1] + self.match(s1[j-1], s2[i-1])
                        Mp[i][j] = 1
                    else:
                        M[i][j] = M[i-1][j-1] + self.match(s1[j-1], s2[i-1])
                else:
                    if(I[i-1][j-1] < D[i-1][j-1]):
                        M[i][j] = I[i-1][j-1] + self.match(s1[j-1], s2[i-1])
                        Mp[i][j] = 1
                    else:
                        M[i][j] = D[i-1][j-1] + self.match(s1[j-1], s2[i-1])
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
        align[0] = align[0][::-1]
        align[1] = align[1][::-1]
        return align

    def multipairwise2(self, node1, node2):
        a = 0.0
        b = 0.0
        Ext = self.ext
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
            M.append([-1]*(len1+1))
            D.append([-1]*(len1+1))
            I.append([-1]*(len1+1))
            Mp.append([0]*(len1+1))
            Dp.append([0]*(len1+1))
            Ip.append([0]*(len1+1))

        for i in range(1,len2+1):
            M[i][0] = D[i][0] = I[i][0] = 0
            Mp[i][0] = Dp[i][0] = Ip[i][0] = 1 # 1 is vertical
        for i in range(1,len1+1):
            M[0][i] = D[0][i] = I[0][i] = 0
            Mp[0][i] = Dp[0][i] = Ip[0][i] = -1 # -1 is horizontal
        begin = 0
        end = 1
        quart = len1/4
        slope = (float(len1))/(len2)
        for i in range(1, len2 + 1):
            i_slope = int(i*slope)
            begin = i_slope - quart
            end = i_slope + quart
            if(begin <= 0):
                begin = 1
            if(end > len1):
                end = len1+1
            for j in range(begin, end):
                score_d = 0.0
                for k1 in range(len(node1[3])):
                    for k2 in range(len(node2[3])):
                        extlist = Ext[node1[0][k1]][node2[0][k2]][0]
                        score_d += self.showExt(extlist, node1[6][k1][j - 1], node2[6][k2][i - 1])
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
                if M[i-1][j] >= I[i-1][j] :
                    if D[i-1][j] > M[i-1][j] :
                        D[i][j] = D[i-1][j]
                        Dp[i][j] = -1
                    else:
                        D[i][j] = M[i-1][j]
                else:
                    if D[i-1][j] > I[i-1][j]:
                        D[i][j] = D[i-1][j]
                        Dp[i][j] = -1
                    else:
                        D[i][j] = I[i-1][j]
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
        '''need to refine the sequence'''
        aligned = []
        realIndex = []
        seqnum1 = len(node1[3])
        seqnum2 = len(node2[3])
        for i in range(seqnum1 + seqnum2):
            aligned.append([])
            realIndex.append([])
        i = len2
        j = len1
        while i != 0 and j != 0:
            if key == 0:
                key = Mp[i][j]
                i -= 1
                j -= 1
                for k in range(seqnum1):
                    aligned[k].append(node1[3][k][j])
                    realIndex[k].append(node1[6][k][j])
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append(node2[3][k][i])
                    realIndex[seqnum1 + k].append(node2[6][k][i])
            elif key == 1:
                key = Ip[i][j]
                j -= 1
                for k in range(seqnum1):
                    aligned[k].append(node1[3][k][j])
                    realIndex[k].append(node1[6][k][j])
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append('-')
                    realIndex[seqnum1 + k].append(-1)
            else:
                key = Dp[i][j]
                i -= 1
                for k in range(seqnum1):
                    aligned[k].append('-')
                    realIndex[k].append(-1)
                for k in range(seqnum2):
                    aligned[seqnum1 + k].append(node2[3][k][i])
                    realIndex[seqnum1 + k].append(node2[6][k][i])
        while i != 0:
            i -= 1
            for k in range(seqnum1):
                aligned[k].append('-')
                realIndex[k].append(-1)
            for k in range(seqnum2):
                aligned[seqnum1 + k].append(node2[3][k][i])
                realIndex[seqnum1 + k].append(node2[6][k][i])
        while j != 0:
            j -= 1
            for k in range(seqnum1):
                aligned[k].append(node1[3][k][j])
                realIndex[k].append(node1[6][k][j])
            for k in range(seqnum2):
                aligned[seqnum1 + k].append('-')
                realIndex[seqnum1 + k].append(-1)
        for i in range(seqnum1):
            aligned[i] = aligned[i][::-1]
            realIndex[i] = realIndex[i][::-1]
        for i in range(seqnum2):
            aligned[seqnum1 + i] = aligned[seqnum1 + i][::-1]
            realIndex[seqnum1 + i] = realIndex[seqnum1 + i][::-1]
        return aligned, realIndex
