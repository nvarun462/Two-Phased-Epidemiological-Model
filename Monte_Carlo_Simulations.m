%This code will run Two_Phased_Config_Network multiple times and automatically
%average only the successful runs and plot/save the data

tic %Starts Timer

%Same variables as Two_Phased_Config_Network - ENSURE CONTINUITY
N = 10000;
tmax1 = 100;
tmax2 = 160;
tmax = tmax1+tmax2;

%Pre-Allocates Matrices
 for ttt = 1:tmax
     susoverall(ttt) = 0;
     infoverall(ttt) = 0;
     recoverall(ttt) = 0;
     susoverallhigh(ttt) = 0;
     infoverallhigh(ttt) = 0;
     recoverallhigh(ttt) = 0;
 end
 
simulations = 0;
pass = 0;
while (simulations < 5) %Keeps running until desired "good" runs are reached
    Two_Phased_Config_Network;
    if (pass==0) %Skips the failed runs
        continue; 
    end
    pass = 0;
    simulations = simulations + 1
    susoverall(1:tmax) = susoverall(1:tmax) + susceptible(1:tmax);
    infoverall(1:tmax) = infoverall(1:tmax) + infected(1:tmax);
    recoverall(1:tmax) = recoverall(1:tmax) + recovered(1:tmax);
    susoverallhigh(1:tmax) = susoverallhigh(1:tmax) + susceptiblehigh(1:tmax);
    infoverallhigh(1:tmax) = infoverallhigh(1:tmax) + infectedhigh(1:tmax);
    recoverallhigh(1:tmax) = recoverallhigh(1:tmax) + recoveredhigh(1:tmax);
end
plot(susoverall/simulations/N, 'b');
save('Data.mat','susoverall', 'infoverall', 'recoverall', 'susoverallhigh', 'infoverallhigh', 'recoverallhigh');
toc %Ends Timer