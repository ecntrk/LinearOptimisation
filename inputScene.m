function [ ] = inputScene ( )
%Takes the input variables in the scenario
%   manually change the values here.


%%Time tick limit
tick = 30; % this means every discretised tick is for 30 minutes.
%Dr. Doan's original code had 20 mins but for this example, 30 min is
%better choice
%
%The reason to put this here is, once use put minutes as time inputs, 
%they don't ahve to change everything if they try to lower the 
%discreettime limit. just make tick any other value!
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

epsilon_i = [8,8,8,8,8,8,8,8];

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


r_ij = ceil (ar_ij/tick);

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

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% setup Time
%%%%%%%%%%%%%%%%%%%%%
global sr ;  %setup time for all r.
global tij; %Travel time from j to i.

sr = zeros(1,R); 
asr = [60, 30, 60]; %setup times from table 3
sr = ceil (asr/tick);

tij = r_ij; %because t_ij is travel times and r_ij is transfer time of equipments from j to i.



%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Disaster types
%%%%%%%%%%%%%%%%%%%%%
global K; %disaster types
global R_k; %equipments needed in case of k.
global T_k; %time points for k where next resource is necessary

K = 2; % terrorist attack and major flood. page 12.
R_k(1,:) = [1,2]; %k = 1, bomb, need IRU, USAR, 1,2
R_k(2,:) = [2,3]; %k = 2, flood, need USAR, HVP. 2,3

aT_k(1,:) = [120, 180, 240, 600] %(2, 3, 4, 10 hrs). 30 mints per tick.
aT_k(2,:) = [120, 240, 600, -100] %(2 hr, 4 hr, 10 hr.) 30 mints per tick.
%will move to cell later. now -100 is end point make even indices!
T_k = ceil (aT_k/tick);


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
%this is a 3D cell array. 1st dimension is l, for each l, we have a 2d cell
%for t (row 1) and r (rest rows).

global d_lrt;




%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Scenarios
%%%%%%%%%%%%%%%%%%%%%
global S; %set of scenarios

%for now, lets consider only 2 scenario with 
%s1 = bomb in A, B, C. 
%s2 = flood in B,C,D.
S = 2;


%this here maps the l to a corresponding s value
global whatL; %simultaneous disasters l's are happening in each s from S.
whatL = zeros(S,3);
whatL(1,:) = [1,2,3]; %these are l values
whatL(1,:) = [11,12,13];

global Ts;
% considering 30 mins as 1 time tick, two s from S has 10 
%hours of operations.
Ts = 600 / tick; %Ts is anyway 10 hours for evey situation.

%here is a potential bug. if you don't update the j_ and only whatJ, the
%indexes the resolver finds will be wrong.




end

