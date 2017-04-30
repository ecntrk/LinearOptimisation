function [ ] = inputScene ( )
%Takes the input variables in the scenario
%   manually change the values here.

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% locations
%%%%%%%%%%%%%%%%%%%%%
global N; %number of locations
global epsilon_i; %nearby locations for each i. This is a vector defining j.
global p_ij; %preference of nearby loc j given loc i. ii is pref 1. ij >1
global r_ij; %response time of j to i. ii = 0.

N = 8;
epsilon_i = [8,8,8,8,8,8,8,8];
p_ij = zeros(8,8);
%            A B C D E F G H
p_ij(1,:) = [1,2,3,4,6,8,5,7];
p_ij(2,:) = [4,1,3,2,5,8,6,7];
p_ij(3,:) = [4,3,6,2,7,8,4,5];
p_ij(4,:) = [4,1,3,2,5,8,6,7];
p_ij(5,:) = [4,3,6,2,7,8,4,5];
p_ij(6,:) = [4,1,3,2,5,8,6,7];
p_ij(7,:) = [4,3,6,2,7,8,4,5];
p_ij(8,:) = [4,3,6,2,7,8,4,5];


%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% type of equipments
%%%%%%%%%%%%%%%%%%%%%
global R; %set of equipment types
global betaR; %number of supporting generalist
global alphaR; %number of specialists
global C_r; % current capacity of type r equipment

R = 5;

betaR = zeros(1,R); 
betaR(:) = 5; %change this with actual data

alphaR = zeros(1,R); 
alphaR(:) = 5; %change this with actual data

C_r = zeros(1,R);


%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% setup Time
%%%%%%%%%%%%%%%%%%%%%
global sr ;  %setup time for all r.
global tij;
sr = zeros(1,R); 
tij = 0;


%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Disaster types
%%%%%%%%%%%%%%%%%%%%%
global K; %disaster types
global R_k; %equipments needed in case of k.
global T_k; %time points for k where next resource is necessary


%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Scenarios
%%%%%%%%%%%%%%%%%%%%%
global S; %set of scenarios

%this here maps the l to a corresponding s value
global whatL; %simultaneous disasters l's are happening in each s from S.
whatL = zeros(1,maxV_(5));
whatL(:) = maxV_(5);
%here is a potential bug. if you don't update the j_ and only whatJ, the
%indexes the resolver finds will be wrong.

%this here maps the il to a corresponding l value
global whatIL; %il is the location of disaster l.
global whatKL; %kl is disaster type of disaster l.
whatIL = zeros(1,maxV_(3));
whatIL(:) = 1;

whatKL = zeros(1,maxV_(3));
whatKL(:) = 1;
%here is a potential bug. if you don't update the j_ and only whatJ, the
%indexes the resolver finds will be wrong.

%R_k_l is the set of type r equipment. this is basically R_k where k = k_l
%T_k_l is the set of time points. this is basically T_k where k = k_l



end

