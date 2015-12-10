#!env python

import operator
import numpy as np

from ViennaWrappers import runRNAFold, runRNAEval


def displaySecStr(l,n):
    exp = ["."]*n
    for left,right in l :
        exp[left] = "("
        exp[right]= ")"
    return "".join(exp)


def parseSecStr(struct):
    n = len(struct)
    stack = []
    l     = []
    for i,c in enumerate(struct) :
        if c == "(" :
            stack.append(i)
        elif c == ")" :
            right = i
            left  = stack.pop()
            l.append([left,right])
    l = sorted(l, key=operator.itemgetter(0))
    return (l,n)


def countSecStr_recursive(r, debug=False,theta=1):
    bp = [("G","C"),("C","G"),("G","U"),("U","G"),("A","U"),("U","A")]
    n = len(r)
    m = np.ones((n,n))
    lm= np.zeros((n,n))
    def Nij(r,m,lm,i,j):
        def call(i,j) :
            if i >= j-1 :
                return 1
            if not lm[i,j] :
                m[i,j] = Nij(r,m,lm,i,j)
                lm[i,j] += 1
            return m[i,j]
        m[i,j] = call(i+1,j)
        for k in range(i+theta+1,j+1) :
            if debug or (r[i],r[k]) in bp or (r[k],r[j]) in bp :
                m[i,j] += call(i+1,k-1) * call(k+1,j)
        return m[i,j]
    return int(Nij(r,m,lm,0,n-1))


def countSecStr_prog_dyna(r, debug=False,theta=1) :
    bp = [("G","C"),("C","G"),("G","U"),("U","G"),("A","U"),("U","A")]
    n  = len(r)
    if n < 3 :
        return 1
    m  = np.eye(n) + np.eye(n,k=1)
    for i in range(n-2,-1,-1) :
        for j in range(i+2,n) :
            m[i,j] = m[i+1,j]
            for k in range(i+theta+1,j+1) :
                if debug or (r[i],r[k]) in bp or (r[k],r[j]) in bp :
                    if k >= j :
                        m[i,j] += m[i+1,k-1]
                    else :
                        m[i,j] += m[i+1,k-1] * m[k+1,j]
    return int(m[0,-1])


def E(n) :
    gi = {"A":0, "U":1, "G":2, "C":3}
    g  = np.zeros((4,4))
    g[gi["A"],gi["U"]] = -1
    g[gi["U"],gi["A"]] = -1
    g[gi["G"],gi["C"]] = -1
    g[gi["C"],gi["G"]] = -1
    g[gi["G"],gi["U"]] = -1
    g[gi["U"],gi["G"]] = -1
    return (g,gi)


def fillMatrix_recursive(r,theta=1) :
    n = len(r)
    m = np.zeros((n,n))
    lm= np.zeros((n,n))
    (g,gi) = E(n)
    def Nij(r,m,lm,i,j):
        def call(i,j) :
            if i >= j - 1 :
                return 0
            if lm[i,j] <= 0 :
                lm[i,j] += 1
                m[i,j] = Nij(r,m,lm,i,j)
            return m[i,j]
        comp = [ call(i+1,j) ]
        for k in range(i+theta+1,j+1) :
            e = g[gi[r[i]],gi[r[k]]]
            comp.append(e + call(i+1,k-1) + call(k+1,j))
        m[i,j] = min(comp)
        lm[i,j] += 1
        return m[i,j]
    [ Nij(r,m,lm,0,l) for l in range(2,n)]
    return m.tolist()


def fillMatrix_prog_dyna(r,theta=1) :
    n  = len(r)
    (g,gi) = E(n)
    if n < 3 :
        return np.zeros((n,n)).tolist()
    m  = np.zeros((n,n))
    for i in range(n-2,-1,-1) :
        for j in range(i+2,n) :
            comp = [ m[i+1,j] ] 
            for k in range(i+theta+1,j+1) :
                e = g[gi[r[i]],gi[r[k]]]
                if k >= j :
                    comp.append(m[i+1,k-1] + e)
                else :
                    comp.append(m[i+1,k-1] + m[k+1,j] + e)
            m[i,j] = min(comp)
    return m.tolist()

    
def traceback(tab):
    ss = []
    def tb(i,j):
        # print tab[i,j],tab[i+1,j],tab[i,j-1],tab[i+1,j-1]
        if j <= i :
            return
        elif tab[i,j] + 0. == tab[i+1,j] + 0. :
            tb(i+1,j)
        elif tab[i,j] + 0. == tab[i,j-1] + 0. :
            tb(i,j-1)
        elif tab[i,j] + 0. == tab[i+1,j-1] + -1. :
            ss.append((i,j))
            tb(i+1,j-1)
        else :
            for k in range(i+1,j-1) :
                if tab[i,j] + 0. == tab[i,k-1] + tab[k+1,j-1] + 0. :
                    tb(i,k-1)
                    tb(k+1,j-1)
    tb(0,len(tab)-1)
    return ss

def nussinov(seq,theta=1,fillMatrix=fillMatrix_prog_dyna) :
    tab = fillMatrix(seq,theta=theta)
    ss = traceback(np.array(tab))
    struct = displaySecStr(ss,len(seq))
    energy = runRNAEval(seq,struct)
    return (struct,energy)

def test_suite():

    print " #### warm-up test : #### "
    print displaySecStr([(2,8),(3,7)],10)

    print "\n #### parsing secondary structures test : #### "
    print parseSecStr("((.(...)))")

    def countSecStr_test(countSecStr) :
        print countSecStr_recursive("A")
        print countSecStr_recursive("CAG")
        print countSecStr_recursive("CAGU")
        print countSecStr_recursive("AAAA",debug=True)
        print countSecStr_recursive("AAAAAAAA",debug=True)
        print [countSecStr_recursive("A"*i,debug=True) for i in range(1,13)]

    print "\n #### recursive counting compatible structures test : #### "
    countSecStr_test(countSecStr_recursive)
    
    print "\n #### dynamic programming counting compatible structures test : #### " 
    countSecStr_test(countSecStr_prog_dyna)
    
    print "\n #### recursive nussinov test : #### "
    print np.array(fillMatrix_recursive("CCCCUUUUGGGGG",3))
    
    print "\n #### dynamic programming nussinov test : #### "
    print np.array(fillMatrix_prog_dyna("CCCCUUUUGGGGG",3))

    print "\n #### Nussinov Traceback test : #### "
    print traceback(np.array(fillMatrix_prog_dyna("CCCCUUUUGGGGG",3)))
    print displaySecStr(traceback(np.array(fillMatrix_prog_dyna("CCCCUUUUGGGGG",3))),13)
    print nussinov("CCCCUUUUGGGGG",theta=3)

#### model discrepancies ####

def delta(seq) :
    (struct_1,energy_1) = runRNAFold(seq)
    (struct_2,energy_2) = nussinov(seq)
    # print (struct_1,energy_1), (struct_2,energy_2)
    return energy_2 - energy_1


def compareSS(s1,s2):
    count = 0
    for e in s1 :
        if e in s2 :
            count += 1
    return count


def benchmark(known,structure_predict):
    (seq_k,stru_k) = known
    correct = []
    for e in seq_k :
        if e in structure_predict :
            correct.append(e)
    return correct
    
            
if __name__ == "__main__" :
    pass
    # test_suite()
    
