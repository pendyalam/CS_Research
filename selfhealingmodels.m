%% Beta vs. Vertices 
    % How the fairness index changes as you add more nodes?
    
    % Xn is the average number of packets 
    
   % B = ((n^2)*(n+1^2))/(4*(n - (1/6)*n*(n+1)*(2n*1))); %Simplified everything
    
    
    % Using connectivity matrix (adjacency matrix) and another matrix that
    % keeps the state of number of packets for each node. 
    
    % Packet state matrix will be used to calculate Jain index - x(i) =
    % number of packets
    
    
%% Jain Fairness Index - holding x constant and changing n (n is the same as v) 
    % Plug this into dv/dt 

%% What is A? How fast you want the network to grow?

%% Simulation of network

% Creating graph

A = zeros(4,4); %Start off with a 4x4 matrix
A = diag([1 1 1],-1)+diag([1 1 1 1],0)+diag([1 1 1],1)+diag(1,-3)+diag(1,3); %Set connections

G = graph(A);
G.Nodes.Names = {'1' '2' '3' '4'}';
G.Nodes.ProcTime = [3 2 6 4]';
G.Nodes.numPackets = [0 0 0 0]';
G.Nodes.probability = [1 0.5 0.5 1]';
plot(G)
hold on
% First packet

s = randi([1 length(nodeVector)],1,1);      %rand source
d = randi([1 length(nodeVector)],1,1);      %rand dest
path = shortestpath(G, s, d);
p = packets(s, d, path);
packetVector = [p];
pi = 2;                                     %packet array index
disp(packetVector(1))

events = 10;   %Number of events


%% For loop
for i = 1:events
    
    % Create a new packet every 5 events and add to the packet vector
    if rem(i, 5) == 0
        disp(num2str(1))
        s = randi([1 numnodes(G)],1,1);      %rand source
        d = randi([1 numnodes(G)],1,1);      %rand dest
        path = shortestpath(G, s, d);               %create path from s to d
        p = packets(s, d, path);                    %create packet object
        packetVector(pi) = p;                       %add packet object to array
        
        disp(p.currentLoc)
        G.Nodes.numPackets(p.source) = G.Nodes.numPackets(p.source) + 1;    %Increment number of packets
        
        disp(packetVector(pi))
        pi = pi + 1;
        
        G.Nodes
    end
    
    % Move packets every 3 events
    if rem(i, 3) == 0
        for j = 1:length(packetVector)
            disp(num2str(2))
            p = packetVector(j);
            cur = p.currentLoc;
            if (p.path(cur) == p.dest && p.time == G.Nodes.ProcTime(p.dest))        %If packet is at dest and is done processing, remove the packet and decrement numPackets in node
                disp(num2str(3))
                G.Nodes.numPackets(p.dest) = G.Nodes.numPackets(p.dest) - 1;
                packetVector(j) = [];
            elseif (p.path(cur) == p.dest && p.time ~= G.Nodes.ProcTime(p.dest))    %If packet is a dest but has not finished processing, add to time   
                disp(num2str(4))
                p.time = p.time + 1;
            else                                                                    %If packet is not at dest, decrement numPacket at current node, increment curLoc, and increment numPacket at newNode
                disp(num2str(4))
                G.Nodes.numPackets(p.path(cur)) = G.Nodes.numPackets(p.path(cur)) - 1;
                p.currentLoc = cur + 1;
                cur = p.currentLoc;
                G.Nodes.numPackets(p.path(cur)) = G.Nodes.numPackets(p.path(cur)) + 1;
            end
            disp([num2str(p.path(cur)), num2str(p.time)])
            G.Nodes
        end
    end
end

   %% 
    % Increment packet number in corresponding node
    for j = 1:length(packetVector)
        p = packetVector(j);
        disp(p)
        for k = 1:numnodes(G)
            if p.currentLoc == nodeVector(k).nodeNum
                nodeVector(k).numPackets = nodeVector(k).numPackets + 1;
            end
            disp(nodeVector(k))
        end 
    end
    plot(G)
%% Random rate of sending number of packets

% Fairness index - rank each node by congestion + connectivity of the nodes
    % Notion of centrality - page rank (net rank)
    % Eigenvector of the matrix ***
    % See how adding nodes and egdes affects centrality
    % To test: shortest path (ECMP)
    % ADD SEED FOR REPRODUCIBILITY
    
    