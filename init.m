%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init ()
%Initialisation of the various variables.
%User needs to change the numbers manually in case of any changes
%Var Order: i, j, l, r, s, t
%Decision Var Order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z


%Define all the maximum values of the variables. For now, only 3!
i_ = 3;
j_ = 3; %change to dependant later
l_ = 3;
r_ = 3;
s_ = 3; %change to dependant later
t_ = 3;

global varLen; %The number of variables i,j,l,r,s,t
varLen = 6;

global maxV_;
maxV_ = zeros(6);
maxV_(1) = i_;maxV_(2) = j_;maxV_(3) = l_;maxV_(4) = r_;maxV_(5) = s_;maxV_(6) = t_;

global vecLen;
vecLen = l_*r_*s_*t_ + i_*r_*s_*t_ + 2*i_*r_ + i_*j_*s_*t_ ...
    + 2*i_*j_*r_*s_*t_ + i_*s_*t_ + 3*i_*r_*s_*t_;

%Ranges for decisionvariables in the vector. 11 element for 11 variables in
%the predetermined order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z
global dVarRanges;
dVarRanges = zeros(11);
dVarRanges(1) = l_*r_*s_*t_;                    % for u
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

%i, j, l, r, s, t
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
end