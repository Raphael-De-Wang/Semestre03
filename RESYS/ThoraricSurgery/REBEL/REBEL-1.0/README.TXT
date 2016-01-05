Some Instructions for Using REBEL Version 1.0

I. POLICY OF USAGE

You may use the software freely for academic, non-profit purposes only. 
Do not modify the source code nor reuse any pieces of it. 
Do not redistribute. 

The purpose of these conditions is to ensure correct and coherent referencing to
the original software. The source code may not be understandable to others than
the author.

It would be appreciated if you would send a note at
mikko.koivisto(at)cs.helsinki.fi stating that you have downloaded the software
and your intended usage.

Please send bug reports to mikko.koivisto(at)cs.helsinki.fi.


II. DESCRIPTION

Implements most of the algorithms presented in 

[1] Mikko Koivisto and Kismat Sood: Exact Bayesian structure discovery in
Bayesian networks. Journal of Machine Learning Research, 5: 549-573, 2004.

[2] Mikko Koivisto: Advances in exact Bayesian structure discovery in Bayesian
networks. In: R. Dechter and T. Richardson (eds.), Proc. of the 22nd Conference
on Uncertainty in Artificial Intelligence (UAI 2006), pp. 241-248, AUAI Press,
2006.

Please cite these references when publishing results obtained using this
software.

In the current version: 
(1) Layering structure (described in [1]) NOT implemented. 
(2) Finding an optimal network structure (described in [1]) NOT implemented.
(3) Other than a default prior NOT supported.


III. COMPILE (USING, E.G., g++)

 g++ -Wall -O2 -o rebel main.cc


IV. RUN

Example:

 rebel -m 4 -d data/n20-m4-r2-i0-s10000.dat -u 500  > test.res

Options: 

 -d Data file:
                                     data/n20-m4-r2-i0-s10000.dat
 -m Maximum indegree:
                                                                4
 -u Maximum number of data records read:
                                                              500

V. SEE RESULTS

File "test.res" will contain the posterior probability of each edge:

n2 -> n1 0.000504846
n3 -> n1 0.000512037
n4 -> n1 0.000627938
n5 -> n1 0.00142807
n6 -> n1 0.000219346
.
.
.
n17 -> n20 0.000705036
n18 -> n20 0.00132344
n19 -> n20 0.000626858


