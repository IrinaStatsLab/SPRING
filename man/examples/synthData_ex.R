# goal is to generate synthetic data with a prescribed graph structure.
# load real data "QMP" in SPRING package.
data(QMP)
set.seed(12345) # set the seed number for make_graph part.
p1 = ncol(QMP) # the number of nodes.
e1 = 2*p1 # the number of edges is set as twice the number of nodes.
gtype = "cluster"
# available types in SpiecEasi: "band", "cluster", "scale_free", "erdos_renyi", "hub", "block".
graph_p1 <- SpiecEasi::make_graph(gtype, p1, e1) # adjacency matrix. 1: edge, 0: no edge.
Prec1  <- SpiecEasi::graph2prec(graph_p1) # precision matrix. inverse of covariance.
Cor1   <- cov2cor(SpiecEasi::prec2cov(Prec1)) # correlation matrix.

X1_count <- synthData_from_ecdf(QMP, Sigma = Cor1, n = 100)
# generate data of size n by p.
# p = ncol(Cor1) = ncol(QMP) should hold.
# need to specify sample size n.




