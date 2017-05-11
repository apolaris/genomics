from sys import argv

if __name__ == "__main__":
    if len(argv) != 11:
        print "python main.py <read_file_name> <write_file_name> <opengap_penalty> <extendgap_penalty> <match_score> <mismatch _score> <seq_num> <option> <refine_time> <input_type>"
        quit()
    read_name = argv[1]
    write_name = argv[2]
    open_gap = float(argv[3])
    extension_gap = float(argv[4])
    match_score = float(argv[5])
    mismatch_score = float(argv[6])
    seq_num = int(argv[7])
    option = int(argv[8])
    refine_time = int(argv[9])
    input_type = int(argv[10])
    seq = []
    if input_type == 0:
        fp = open(read_name,"r")
        one = []
        count = 0
        flag = False
        while (1):
            a  = fp.read(1)
            if not a:
                break
            elif a == '>':
                if one :
                    seq.append(one)
                one = []
                count = 0
                flag = False
            elif a != '\n' and flag:
                one.append(a)
                count += 1
            elif a == '\n' and not flag:
                flag = True
        '''for line in seq:
            for c in line:
                fp2.write(c)
            fp2.write("\n")'''
        fp.close()
    elif input_type == 1:
        fp = open(read_name,"r")
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
    else:
        print "the input type does not exist!"
        quit()
    if option == 0:
        from multiple import *
    elif option == 1:
        from multiple2 import *
    elif option == 2:
        from multiple3 import *
    else:
        print "the option does not exist!"
        quit()
    #[seqs, extend gap, open gap, match, mismatch
    if seq_num == 0:
        seq_num = len(seq)
    mul = sequences(seq[0:seq_num], extension_gap, open_gap, match_score, mismatch_score, write_name)
    t1 = time.clock()
    mul.init_matrix()

    t2 = time.clock()
    mul.buildtree()

    t3 = time.clock()
    mul.complete()

    t4 = time.clock()
    mul.refine(refine_time)
    t5 = time.clock()

    print "calculate distance matrix: %f s" % (t2 - t1)
    print "build tree: %f s" % (t3 - t2)
    print "progressive alignment: %f s" % (t4 - t3)
    print "refinement: %f s" % (t5 - t4)
