%% Bipartite Graph
% Make a random MxN adjacency matrix
m = 1000;
n = 1000;
p = 0.3;
a = rand(n,n) < p;
% Expand out to symmetric (M+N)x(M+N) matrix
big_a = [zeros(m,m), a;
         a', zeros(n,n)];     
g = graph(big_a);

% Grab adjacency list
adjList = g.Edges;
adjList(:,[2]) = []; 

% Write to text file
writetable(adjList, 'bipartite.txt', 'Delimiter', ' ', 'WriteVariableNames', 0)

%% Undirected Graph
% Choose number of nodes and probability of edge
n = 1000; 
p = 0.3;
a = rand(n,n) < p;
a = triu(a,1);
a = a + a';
g = graph(a);

% Grab adjaacency list
adjList = g.Edges;
adjList(:,[2]) = [];

% Write to text file
writetable(adjList, 'graph.txt', 'Delimiter', ' ', 'WriteVariableNames', 0)

%% Plot
h = plot(g);
% Make it pretty
h.XData(1:m) = 1;
h.XData((m+1):end) = 2;
h.YData(1:m) = linspace(0,1,m);
h.YData((m+1):end) = linspace(0,1,n);

%% Plotting Data

% Changing |U|
t = [0.004884, 0.003156, 0.004021, 0.002878, 0.002757, 0.002291, 0.002199
0.00299
0.003345
0.002951];
u = [100	1000	2500	5000	10000];
boxplot(t,u)

%% Finding fit of plot

x = [100, 1000, 5000, 10000, 20000, 30000, 40000, 50000];
Waxman = [0.1719285	0.172152	0.2218442	0.2733775	0.3825667	0.4633624	0.49554	0.4716769];
Barabasi = [0.1557554	0.161688	0.2169493	0.2541176	0.362568	0.398463	0.4370652	0.4429345];


plot(x, Waxman, 'o-') %% Cubic

plot(x, Barabasi, 'o-') %% Cubic