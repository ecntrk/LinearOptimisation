%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vector ] = findVariables( index )
%FINDVARIABLES Summary of this function goes here
%   Given and index, it returns a vector of 7 elements. 1st element is the
%   decision variable number. Next 6 elements are indexes
% Var Order: s, l, r, i, j, t
% Decision Var Order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z,
% [6, 2,0,1,2,3,4] means variable x with indices i= 2, j = 3, r = 1, s = 2,
% t = 4.

global varLen;
global dVarRanges; global maxV_;
global S; global L; global R; global Ts;

varP = 1;
%             s, l, r, i, j, t
indexGuide = [1, 1, 1, 0, 0, 1;
              1, 0, 1, 1, 0, 1;
              0, 0, 1, 1, 0, 0;
              0, 0, 1, 1, 0, 0;
              1, 0, 0, 1, 1, 1;
              1, 0, 1, 1, 1, 1;
              1, 0, 1, 1, 1, 1;
              1, 0, 0, 1, 0, 1;
              1, 0, 1, 1, 0, 1;
              1, 0, 1, 1, 0, 1;
              1, 0, 1, 1, 0, 1];



for iter = 1:11
   index = index - dVarRanges(iter);
   if (index <= 0)
        index = index+dVarRanges(iter);
        varP = iter;
        break;
   end   
end


vector = [varP];
end

