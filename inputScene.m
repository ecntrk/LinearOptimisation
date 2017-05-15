%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[att]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inputScene
%Takes the input variables in the scenario
%   manually change the values here.


%%Time tick limit
tick = 30; % this means every discretised tick is for 30 minutes.
%Dr. Doan's original code had 20 mins but for this example, 30 min is
%better choice
%
%The reason to put this here is, once use put minutes as time inputs, 
%they don't ahve to change everything if they try to lower the 
%discreet time limit. just make tick any other value!
%
%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% locations
%%%%%%%%%%%%%%%%%%%%%
global N; %number of locations
global epsilon_i; %nearby locations for each i. This is a vector defining j.
global p_ij; %preference of nearby loc j given loc i. ii is pref 1. ij >1
global r_ij; %response time of j to i. ii = 0.

N = 8;

epsilon_i = zeros(N, N-1);
%making the adjacent cities list epsilon_i
for count = 1:N
    temp  = 1:N;
    tempp = temp(find(temp~=count));
    epsilon_i(count,:) = tempp;
end

p_ij = zeros(N,N);
%            A B C D E F G H   from table 2.
p_ij(1,:) = [1,2,3,4,6,8,5,7];
p_ij(2,:) = [4,1,3,2,5,8,6,7];
p_ij(3,:) = [6,3,1,2,7,8,4,5];
p_ij(4,:) = [7,3,2,1,4,8,5,6];
p_ij(5,:) = [8,6,4,3,1,2,5,7];
p_ij(6,:) = [8,7,4,6,3,1,2,5];
p_ij(7,:) = [8,7,4,5,6,2,1,3];
p_ij(8,:) = [8,6,3,5,7,4,2,1];


ar_ij = zeros(N,N);
%            A   B   C   D   E   F   G   H   from table 1, picture.
ar_ij(1,:) = [20,150,230,250,410,560,390,420]; %A
ar_ij(2,:) = [150,25,120,100,260,410,280,280]; %B
ar_ij(3,:) = [230,120,25, 80,240,270,160,160]; %C
ar_ij(4,:) = [250,100,80, 20,160,310,240,240]; %D
ar_ij(5,:) = [410,260,240,160,25,150,260,380]; %E
ar_ij(6,:) = [560,410,270,310,150,25,110,230]; %F
ar_ij(7,:) = [390,280,160,240,260,110,20,120]; %G
ar_ij(8,:) = [420,280,160,240,380,230,120,25]; %H


r_ij = ceil (ar_ij/tick); %so now, r_ij has the discreetised times

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% type of equipments
%%%%%%%%%%%%%%%%%%%%%
global R; %set of equipment types
global alphaR; %number of specialists
global betaR; %number of supporting generalist
global C_r; % current capacity of type r equipment

R = 3; %IRU, USAR, HVP

alphaR = zeros(1,R); 
alphaR(:) = [10, 10, 2]; % actual data from table 3

betaR = zeros(1,R); 
betaR(:) = [20, 25, 5]; %actual data from table 3

C_r = zeros(1,R);
C_r = [72, 20, 46]; %from table 3

%There are penalties for equipment efficiency. 
global phi_r
global theta;
global mu;

phi_r = zeros(1,R);
phi_r(:) = 1; %please change this if you have data about efficiency

theta = 1;
mu = 1;


%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% setup Time
%%%%%%%%%%%%%%%%%%%%%
global sr ;  %setup time for all r.
global tij; %Travel time from j to i.
global Td; %deployment time.

sr = zeros(1,R); 
asr = [60, 30, 60]; %setup times from table 3
sr = ceil (asr/tick);

tij = r_ij; %because t_ij is travel times and r_ij is transfer time of equipments from j to i.

Td = 300/tick; %becaus eit is mentioned, firefighters need 5 hours break and we know shift time = break time 


%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Disaster types
%%%%%%%%%%%%%%%%%%%%%
global K; %disaster types
global R_k; %equipments needed in case of k.
global T_k; %time points for k where next resource is necessary

K = 2; % terrorist attack and major flood. page 12.
R_k(1,:) = [1,2]; %k = 1, bomb, need IRU, USAR, 1,2
R_k(2,:) = [2,3]; %k = 2, flood, need USAR, HVP. 2,3

aT_k{1} = [120, 180, 240, 600]; %(2, 3, 4, 10 hrs). 30 mints per tick.
aT_k{2} = [120, 240, 600]; %(2 hr, 4 hr, 10 hr.) 30 mints per tick.
%will move to cell later. now -100 is end point make even indices!
%correction: made it cell anywway!
for count = 1:K
    T_k{count} = ceil(aT_k{count}/tick);
end


%number of possible disasters (l) is 
global L;
L = N*K; %each city can have k disaster types

%each l from L has il and kl, location and disaster types.
%trivial to find them

%this here maps the il to a corresponding l value
global whatIL; %il is the location of disaster l.
global whatKL; %kl is disaster type of disaster l.

whatIL = zeros(1,L);
for count = 1:L
    whatIL(count) = mod(count,N);
end
whatIL(whatIL==0) = N;

whatKL = zeros(1,L);
for count = 1:L
    whatKL(count) = mod(count,K);
end
whatKL(whatKL==0) = K;


%R_k_l is the set of type r equipment. this is basically R_k where k = k_l
%T_k_l is the set of time points. this is basically T_k where k = k_l

%next is d_lrt which is requirement of r at time t in cae of l
%this is a cell array. 1st dimension is l, for each l, we have a 2d cell
%for t (1st row) and r (2:last rows). 
%TO access d_lrt, access the cell l, then traverse row 1 for a t match.
%then go to r+1 row for the r. 


global d_lrt;

d_lrt = cell (2,L);
d_lrt{1,1} = 1; %this means at city A, disaster 1 (bomb) happening
d_lrt{2,1} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              10, 0, 0, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)

d_lrt{1,2} = 2; %city B, disaster 1 (bomb) happening
d_lrt{2,2} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              0, 0, 10, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)

d_lrt{1,3} = 3; %city C, disaster 1 (bomb) happening
d_lrt{2,3} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              0, 0, 10, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)


d_lrt{1,4} = 4; %city D, disaster 1 (bomb) happening
d_lrt{2,4} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              0, 0, 10, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)


d_lrt{1,5} = 5; %city E, disaster 1 (bomb) happening
d_lrt{2,5} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              0, 0, 10, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)


d_lrt{1,6} = 6; %city F, disaster 1 (bomb) happening
d_lrt{2,6} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              0, 0, 10, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)


d_lrt{1,7} = 7; %city G, disaster 1 (bomb) happening
d_lrt{2,7} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              0, 0, 10, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)


d_lrt{1,8} = 8; %city H, disaster 1 (bomb) happening
d_lrt{2,8} = [T_k{1}; %time ticks (2, 3, 4, 10 hours)
              0, 10, 0, 4; % IRU units (300 casualty per unit)
              0, 0, 10, 5; %USAR units
              0, 0, 0 ,0]; %HVP (not needed at all)

   
d_lrt{1,9} = 9; %city A, disaster 2 (flood) happening
d_lrt{2,9} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              10, 0, 5; %USAR units
              10, 0 ,0]; %HVP (not needed at all)
          
d_lrt{1,10} = 10; %city B, disaster 2 (flood) happening
d_lrt{2,10} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              0,10, 5; %USAR units
              3, 0 ,0]; %HVP (not needed at all)          

          
d_lrt{1,11} = 11; %city C, disaster 2 (flood) happening
d_lrt{2,11} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              0,10, 5; %USAR units
              3, 0 ,0]; %HVP (not needed at all)   
          
d_lrt{1,12} = 12; %city D, disaster 2 (flood) happening
d_lrt{2,12} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              0,10, 5; %USAR units
              5, 0 ,0]; %HVP (not needed at all)          

          
d_lrt{1,13} = 13; %city E, disaster 2 (flood) happening
d_lrt{2,13} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              0,10, 5; %USAR units
              3, 0 ,0]; %HVP (not needed at all) 
          
d_lrt{1,14} = 14; %city F, disaster 2 (flood) happening
d_lrt{2,14} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              0,10, 5; %USAR units
              3, 0 ,0]; %HVP (not needed at all)          

          
d_lrt{1,15} = 15; %city G, disaster 2 (flood) happening
d_lrt{2,15} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              0,10, 5; %USAR units
              3, 0 ,0]; %HVP (not needed at all) 
          
d_lrt{1,16} = 16; %city H, disaster 2 (flood) happening
d_lrt{2,16} = [T_k{2}; %time ticks (2, 4, 10 hours)
              0, 0, 0; % IRU units (300 casualty per unit)
              0,10, 5; %USAR units
              5, 0 ,0]; %HVP (not needed at all) 



%Simultaneously there's  d0_lt          
%We do not have any time requirement on the generalists.
%in case, we do, this can be easily changed to be like the d_lrt

global d0_lt; %l is only index
%       A      B    C     D   E    F      G    H
d0_lt = [200, 130, 130, 150, 130 , 130 , 150, 130, % for k = 1
        150,  80,  80,   100, 80,  80,   100, 80];  % for k = 2
          
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Scenarios
%%%%%%%%%%%%%%%%%%%%%
global S; %set of scenarios

%for now, lets consider only 2 scenario with 
%s1 = bomb in A, B, C. 
%s2 = flood in B,C,D.
S = 3;


%this here maps the l to a corresponding s value
global L_s; %simultaneous disasters l's are happening in each s from S.
%this has to be cell because variable numebr of locaitons possible
% L_s{1} = [1,2,3]; %these are l values
% L_s{2} = [11,12,13];

%creating all possible combination of scenarios of same type
tempV = 1:8;
temp1 = nchoosek(tempV, 3);%generating combinaitons of 8C3
tempV = 9:16;
temp2 = nchoosek(tempV, 3);%generating combinaitons of 8C3
temp = [temp1;temp2];
[rowcomb,~] = size(temp);
for iter = 1:rowcomb
L_s{iter} = temp(iter,:);
end



global Ts;
% considering 30 mins as 1 time tick, two s from S has 10 
%hours of operations.
Ts = zeros(1,S);
Ts(:) = 600 / tick; %Ts is anyway 10 hours for evey situation.
%here's provision if the situation changes anytime in future!

%There is the probability of having a particular diaster s!
global q_s;
q_s = zeros(1,S);
q_s(:) = 1/S; %right now, making it same probability. 
%chnage it if you may for each s

%keyboard
end

