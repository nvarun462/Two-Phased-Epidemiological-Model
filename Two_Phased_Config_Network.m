%This code creates a configuration network in which nodes(humans) with
%varying levels of contacts(degree) are put under SIR dynamics during a
%quarantine (Q1). The quarantine is then loosened(Q2) after certain criteria are met
%and the overall spread is analyzed.

%clear all

%CLEARING SPECIFIC VARIABLES FROM PREVIOUS RUNS
 clear e;
 clear e2;
 clear c;
 clear c_new;
 clear edge;
 clear edge1;
 clear allEdges;
 clear q;
 clear k;

N = 10000; % Number of Nodes
A = 9000; x = 2; %A Nodes each with x edges
B = 1000; y = 50; %B Nodes each with y edges

if ((A+B)~= N) %Confirms Correct Node Count
    disp 'Error: All types of Nodes should add to N';
    return;
end
S = A*x + B*y; %Total number of Stubs(two stubs connect to make a contact pair)

if (mod(S,2)==1) 
    disp 'Error: Total Stubs must be even';
    return;
end

%Creates an Matrix to keep track of which edge belongs to which node.
%The designation of edge E can be found with e(E,:)
k = 1;
for i=1:A
    for j=1:x
        q = strlength(['N',num2str(i),'E',num2str(j)]);
        e(k,1:q)= ['N',num2str(i),'E',num2str(j)];
        k= k+1;
    end
end
for i=(A+1):(B+A)
    for j=1:y
        q = strlength(['N',num2str(i),'E',num2str(j)]);
        e(k,1:q)= ['N',num2str(i),'E',num2str(j)];
        k= k+1;
    end
end



%Randomly Pairing Edges Together
Rand = randperm(S);  %Generates a random permutation of all Stubs, with every 2 representing a pair
c(1:2,1:(S/2)) = 0;  %Matrix of all Connections between 2 Stubs made by segmenting Rand into pieces of length 2
for i = 1:(S/2)  %Populates c matrix
    c(1,i) = Rand((2*i)-1);
    c(2,i) = Rand(2*i);
end 

edge(1:2, 1:(S/2)) = 0;
for i = 1:S
    edge(i) = getNode(e(c(i),:));
end

OutA(1:x,1:A)= 0; %Output Matrix for A showing the connections of each Node
for i = 1:(A*x) %Populates OutA
I = find(c == i);
if (mod(I,2) == 0) %Depending on if I is bottom/top row, the corresponding pair is different
    OutA(i) = getNode(e(c((I-1)),:)); %Find the designation of the corresponding edge and use getNode
else
    OutA(i) = getNode(e(c((I+1)),:));
end
end  

OutB(1:y,1:B)= 0; %Output Matrix for B showing the connections of each Node
for i = ((A*x)+1):((B*y)+(A*x)) %Populates OutB
I = find(c == i);
if (mod(I,2) == 0)
    OutB(i-((A*x))) = getNode(e(c((I-1)),:));
else
    OutB(i-((A*x))) = getNode(e(c((I+1)),:));
end
end  


% Modeling SIR on the network
transmission = 0.01; % probability to transmit the disease per day (one meeting)
duration = 15; % disease duration in days
tmax1 = 100; %Max days for Q1
tmax2 = 160; %Max days for Q2
for ttt = 1:(tmax1+tmax2)
    susceptible(ttt) = 0;
    infected(ttt) = 0;
    recovered(ttt) = 0;
end
for ttt = 1:(tmax1+tmax2)
    susceptiblehigh(ttt) = 0;
    infectedhigh(ttt) = 0;
    recoveredhigh(ttt) = 0;
end
%1 - susceptible, 2 - infected, 3 - recovered
%initial conditions
for i=1:N,
    nodes(i) = 1;
end
index = floor(N*rand)+1; % randomly infect a single node
nodes(index) = 2;
% start the dynamics
internalclock(1:N) = 0; % represents the internal clock of infected

flagg = 0; %Flagg determines if Q1 should shift to Q2

t=0;
while (flagg < 1)
    % single day dynamics
    for i=1:S/2
        indexedge = floor((S/2)*rand)+1; % randomly choose an edge
        n1 = edge(1,indexedge);
        n2 = edge(2,indexedge);
        if ((nodes(n1)==1)&&(nodes(n2)==2)&&(rand < transmission))
            nodes(n1) = 2;
            internalclock(n1) = 1;
        end
        if ((nodes(n2)==1)&&(nodes(n1)==2)&&(rand < transmission))
            nodes(n2) = 2;
            internalclock(n2) = 1;
        end
    end
    % updating the recovery
    for i=1:N
        if (nodes(i)==2)
            internalclock(i) = internalclock(i) + 1;
            if (internalclock(i)>duration)
                nodes(i) = 3;
            end
        end
    end
                
    t=t+1;
    
    sus = 0; infec = 0; rec = 0;
    sushigh = 0; infechigh = 0; rechigh = 0;
    % some counting
    for i=1:A
        if (nodes(i)==1) sus = sus + 1;
        elseif (nodes(i)==2) infec = infec + 1;
        else             rec = rec +1;
        end
    end
    for i=(A+1):N
        if (nodes(i)==1) sushigh = sushigh + 1;
        elseif (nodes(i)==2) infechigh = infechigh + 1;
        else             rechigh = rechigh +1;
        end
    end
    susceptible(t) = sus+sushigh;
    infected(t) = infec+infechigh;
    recovered(t) = rec+rechigh;
    susceptiblehigh(t) = sushigh;
    infectedhigh(t) = infechigh;
    recoveredhigh(t) = rechigh;
    
    if ((t==duration*2+1) && (infected(duration*2+1) == 0)) %(If everyone that the first infected person infects dies without infecting anyone else, run is treated as a "fail"
        disp 'Failure: Epidemic Did Not Occur';
        pass = 0;
        return;
    end
    
    if ((t==tmax1)||((infected(t)/N<0.002)&&((t>1)&& (infected(t)/N < infected(t-1)/N)))) %Criteria to determine if flagg is triggered
        flagg = 2;
        tShift = t;
    end
 
end % finishes the while loop

%New Code to add edges
A2 = 2000; %number of A nodes that increase connections
addedEdges = 23; %number of new edges each A2 node gets
defaultEdges = 4; %A nodes not selected will still get defaultEdges more edges

if ((A2)> A) %Confirms A2 is Less than A
    disp 'Error: A2 is too high';
    return;
end

nodeSelector = randperm(A);
updatedNodes = nodeSelector(1,1:A2); %selects A2 nodes at random from A

%Creates an Matrix to keep track of which new edge belongs to which node.
%The designation of edge E can be found with e(E,:)
k = 1;
for i=1:A
if (ismember(i,updatedNodes))
    for j=1:(addedEdges)
        q = strlength(['N',num2str(i),'E',num2str(j)]);
        e2(k,1:q)= ['N',num2str(i),'E',num2str(j)];
        k= k+1;
    end
else
    for j=1:defaultEdges 
        q = strlength(['N',num2str(i),'E',num2str(j)]);
        e2(k,1:q)= ['N',num2str(i),'E',num2str(j)];
        k= k+1;
    end
end
end

S_new = (A2*addedEdges) + ((A-A2)*defaultEdges); %all new Stubs
Rand_new = randperm(S_new); 
c_new(1:2,1:(S_new/2)) = 0;  %Matrix of all Connections between 2 new Stubs made by segmenting Rand into pieces of length 2
for i = 1:(S_new/2)  %Populates c matrix
    c_new(1,i) = Rand_new((2*i)-1);
    c_new(2,i) = Rand_new(2*i);
end 

edge1(1:2, 1:(S_new/2)) = 0; 
for i = 1:S_new
    edge1(i) = getNode(e2(c_new(i),:));
end

allEdges = [edge edge1]; %one big array with all old and new connections


% Phase two

while (t<(tmax1+tmax2))
    % single day dynamics
    for i=1:(S+S_new)/2
        indexedge = floor(((S+S_new)/2)*rand)+1; % randomly choose an edge
        n1 = allEdges(1,indexedge);
        n2 = allEdges(2,indexedge);
        if ((nodes(n1)==1)&&(nodes(n2)==2)&&(rand < transmission))
            nodes(n1) = 2;
            internalclock(n1) = 1;
        end
        if ((nodes(n2)==1)&&(nodes(n1)==2)&&(rand < transmission))
            nodes(n2) = 2;
            internalclock(n2) = 1;
        end
    end
    % updating the recovery
    for i=1:N
        if (nodes(i)==2)
            internalclock(i) = internalclock(i) + 1;
            if (internalclock(i)>duration)
                nodes(i) = 3;
            end
        end
    end
                
    t=t+1;
    
    sus = 0; infec = 0; rec = 0;
    sushigh = 0; infechigh = 0; rechigh = 0;
    % some counting
    for i=1:A
        if (nodes(i)==1) sus = sus + 1;
        elseif (nodes(i)==2) infec = infec + 1;
        else             rec = rec +1;
        end
    end
    for i=(A+1):N
        if (nodes(i)==1) sushigh = sushigh + 1;
        elseif (nodes(i)==2) infechigh = infechigh + 1;
        else             rechigh = rechigh +1;
        end
    end
    susceptible(t) = sus+sushigh;
    infected(t) = infec+infechigh;
    recovered(t) = rec+rechigh;
    susceptiblehigh(t) = sushigh;
    infectedhigh(t) = infechigh;
    recoveredhigh(t) = rechigh;
    
    if ((t==duration*2+1) && (infected(duration*2+1) == 0))
        disp 'Failure: Epidemic Did Not Occur';
        pass = 0;
        return;
    end
    
end % finishes the while loop
pass = 1; %for use in MC Simulations Code

function Nreturn = getNode(Esel) %Parses the Designation to get just the Node number
P = find(Esel == 'E') - 1;
Nreturn = str2double(Esel(:,2:P));
end