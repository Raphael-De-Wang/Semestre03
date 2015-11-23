#!env Rscript


# 3.1.2 Application with igraph
ER.graph <- erdos.renyi.game(10000,(10/10000))

ER.degrees <- degree.distribution(ER.graph)
plot(ER.degrees,type='h')

ER.transitivity <- transitivity(ER.graph,type="local")
ER.argmax10ind <- order(ER.transitivity, decreasing=TRUE)[1:10]
ER.argmax10val <- sort(ER.transitivity, decreasing=TRUE)[1:10]

ER.averagePathLen <- average.path.length(ER.graph)


# 3.2.2 Application 

AB.graph <- barabasi.game(10000, power=0.3, m=6)

AB.degrees <- degree.distribution(AB.graph)
plot(AB.degrees,type='h')

AB.powlaw <- power.law.fit(AB.degrees)

AB.transitivity <- transitivity(AB.graph,type="local")
AB.argmax10ind <- order(AB.transitivity, decreasing=TRUE)[1:10]
AB.argmax10val <- sort(AB.transitivity, decreasing=TRUE)[1:10]

AB.averagePathLen <- average.path.length(AB.graph)



