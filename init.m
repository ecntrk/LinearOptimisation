%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[att]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init ()
%Initialisation of the various variables.
%User needs to change the numbers manually in case of any changes
%Var Order: s, l, r, i, j, t
%Decision Var Order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z, 

%inputScene();

%%%
% globals from input
%%%
global N; 
%global epsilon_i; 
global R; 
global S;
global L;
global Ts;
%global L_s;
%global T_k;

%Define all the maximum values of the variables. For now, only 3!
i_ = N;
j_ = N-1; %changed to dependant in generateindices
l_ = L;
r_ = R;
s_ = S; 
t_ = max(Ts);

global varLen; %The number of variables i,j,l,r,s,t
varLen = 6;

global maxV_;
maxV_ = zeros(1,varLen);

%the order is according to dependencies. refer to documentaiton sec 2..
maxV_(1) = s_; %because s is top of th efood chain! 
maxV_(2) = l_; 
maxV_(3) = r_;
maxV_(4) = i_;
maxV_(5) = j_;
maxV_(6) = t_;


%%%%%%%%% OBSOLETE
%this here maps the j to a corresponding i value
% global whatJ;
% whatJ = zeros(1,maxV_(1));
% whatJ(:) = maxV_(2);
%here is a potential bug. if you don't update the j_ and only whatJ, the
%indexes the resolver finds will be wrong.


% this is needed to minimize footpront of DVar u.
% Currently commented out for simplicity
% calculating max set length of Tk
% max_T_k = 0;
% for count = 1:length(T_k)
%     a = length(T_k{count});
%     if(max_T_k < a)
%         max_T_k = a;
%     end
% end
% 
% % calculating max set length of Ls
% max_L_s = 0;
% for count = 1:length(L_s)
%     a = length(L_s{count});
%     if(max_L_s < a)
%         max_L_s = a;
%     end
% end

global vecLen;
vecLen = l_*r_*s_*t_ + i_*r_*s_*t_ + 2*i_*r_ + i_*j_*s_*t_ ...
    + 2*i_*j_*r_*s_*t_ + i_*s_*t_ + 3*i_*r_*s_*t_% + l_*t_ + l_*r_*t_;

%Ranges for decisionvariables in the vector. 11 element for 11 variables in
%the predetermined order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z

global dVarRanges;
dVarRanges = zeros(1,11);

dVarRanges(1) = l_*r_*s_*t_;            % for u (because only u is different)
dVarRanges(2) = dVarRanges(1)+i_*r_*s_*t_ ;     % for v
dVarRanges(3) = dVarRanges(2)+i_*r_;            % for w
dVarRanges(4) = dVarRanges(3)+i_*r_;            % for wbar
dVarRanges(5) = dVarRanges(4)+i_*j_*s_*t_;      % for x0
dVarRanges(6) = dVarRanges(5)+i_*j_*r_*s_*t_;   % for x
dVarRanges(7) = dVarRanges(6)+i_*j_*r_*s_*t_;   % for xbar
dVarRanges(8) = dVarRanges(7)+i_*s_*t_;         % for y0
dVarRanges(9) = dVarRanges(8)+i_*r_*s_*t_;      % for y
dVarRanges(10) = dVarRanges(9)+i_*r_*s_*t_;     % for ybar
dVarRanges(11) = dVarRanges(10)+i_*r_*s_*t_;    % for z

disp(dVarRanges);

%debug purposes
%s, l, r, i, j, t
% resolvePos(1, [0,0,1,2,1,2])
% resolvePos(2, [1,0,0,1,1,1])
% resolvePos(3, [1,0,0,1,0,0])
% resolvePos(4, [1,0,0,1,0,0])
% resolvePos(5, [1,1,0,0,1,1])
% resolvePos(6, [1,1,0,1,1,1])
% resolvePos(7, [1,1,0,1,1,1])
% resolvePos(8, [2,0,0,0,1,1])
% resolvePos(9, [1,0,0,1,1,1])
% resolvePos(10, [1,0,0,1,1,1])
% resolvePos(11, [1,0,0,1,1,1])
%Vec_ = construct();
%resolvepos is working perfectly!
end