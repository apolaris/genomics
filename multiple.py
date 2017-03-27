def match (a, b)
    global p
    return (a==b)? 0:p



def pariwise (a, b)
    global gap
    global mis
    alen = len(a)
    blen = len(b)
    col1 = [0]*(alen+1);
    col2 = [0]*(alen+1);
    col1[0] = 0
    for i in range(1,alen+1)
        col1[i] = col1[i-1] - gap;
    now = 2
    for i in range(0,blen)
        if(now == 1)
            col1 = col2[0] - gap
            for j in range(1,alen+1)
                col1[j] = max(col1[j-1] - gap, col2[j] - gap, col2[j-1] + match(a[j-1], b[i]))
            now = 2
        else if(now == 2)
            col2 = col1[0] - gap
            for j in range(1,alen+1)
                col2[j] = max(col2[j-1] - gap, col1[j] - gap, col1[j-1] + match(a[j-1], b[i]))
            now = 1
    if(now == 1)
        return col2[alen]
    else
        return col1[alen]



def build_tree (deff, num)
#deff is the difference between the sequences
    score = [[0]*num]*num
    for i in range(0,num)
        for j in range(0,num)
            score[i][j] += pairwise(seqs[i], seqs[j])





























class sequences:
    'sequence list for alignment'
    def __init__(self, seqs, a, b, mis):
        self.seqs = seqs
        self.num = len(seqs)
        self.a = a
        self.b = b
        self.mis = mis

    #do all
    def all (self):
        self.init_matrix()
        self.build()
        self.complete()

    #calculate the difference
    def init_matrix (self):
        self.score = [[0]*self.num]*self.num
        for i in range(0,self.num):
            node = [1,[],[],seqs[i],-1,-1,-1] #num, left, right, seqs, weight, difstance from parent, height
            self.score[i].append(node)
            for j in range(0,i):
                self.score[i][j] = self.pairwise(i,j)
                self.score[j][i] = self.score[i][j]


    #compare the two sequences
    def pairwise (self, i, j):#penalty for gap is a(k-1)+b
        a = self.a
        b = self.b
        mis = mis
        s1 = self.seqs[i]
        s2 = self.seqs[j]
        num1 = len(s1)
        num2 = len(s2)
        D1 = [0]*(num1+1)
        D2 = [0]*(num1+1)
        I1 = [0]*(num1+1)
        I2 = [0]*(num1+1)
        M1 = [0]*(num1+1)
        M2 = [0]*(num1+1)
        K = [0]*(num1+1)
        D1[1] = I1[1] = M1[1] = K[1] = b
        key = 2
        D2[0] = b-a
        for i in range(2,num1+1):
            D1[i] = D1[i-1] + a
            I1[i] = D1[i]
            M1[i] = D1[i]
            K[i] = D1[i]

        for i in range(0,num2):
            if k==2:
                D2[0] = K[i+1]
                I2[0] = K[i+1]
                M2[0] = K[i+1]
                for j in range(1,num1+1):
                    D2[j] = min(D1[j-1], I1[j-1], M1[j-1]) + match(s1[j-1], s2[i-1])
                    I2[j] = min(I1[j]+a, M1[j]+b, D1[j]+b)
                    M2[j] = min(M2[j-1]+a, I2[j-1]+b, D2[j-1]+b)
            else:
                D1[0] = K[i+1]
                I1[0] = K[i+1]
                M1[0] = K[i+1]
                for j in range(1,num1+1):
                    D1[j] = min(D2[j-1], I2[j-1], M2[j-1]) + match(s1[j-1], s2[i-1])
                    I1[j] = min(I2[j]+a, M2[j]+b, D2[j]+b)
                    M1[j] = min(M1[j-1]+a, I1[j-1]+b, D1[j-1]+b)
        if (k==2):
            return min(M1[num1], D1[num1], I1[num1])
        else
            return min(M2[num1], D2[num1], I2[num1])


    def buildtree (self)

        num = self.num
        diff = self.score1
        d = 0
        di = -1
        dj = -1
        nj = [[0]*num]*num
        while(num>=2):
            for i in range(0,num):
                su[i] = sum(diff[i][:num])
            for i in range(0,num):
                for j in range(0,i):
                    nj[i][j] = (su[i] + su[j]) - (num-2)*diff[i][j]
                    if(nj[i][j] > d):
                        d = nj[i][j]
                        di = i
                        dj = j

            diff[i][num][3] = 0.5*diff[i][j]+(su[i] - su[j])/(2*(num-2))
            diff[j][num][3] = diff[i][j] - diff[i][num][3]
            node = [0,diff[i][num],diff[j][num],-1,-1,-1]
            diffk = []
            for i in range(0,num):
                if(i!=di && i != dj):
                    diffk.append((diff[di][i] + diff[dj][i]- diff[dj][di])/2)
            diffk.append(0)
            diffk.append(node)

            for i in range(0,num):
                if(i!=di && i != dj):
                    del diff[i][di]
                    del diff[i][dj]

            del diff[di]
            del diff[dj]
            #now is num-1*num-1 matrix
            for i in range(0,num-2):
                diff[i].insert(num-2,diffk[i])
            diff.insert(num-2, diffk)
            num=num-1
        self.weighttree(diff[1],1)
        self.tree = diff[1]
        return diff[1]
    def weighttree (self,t,h)
        t[6] = h
        if(!t[1]):
            weighttree(t[1],h+1)
        if(!t[2]):
            weighttree(t[2],h+2)


#class tree:

