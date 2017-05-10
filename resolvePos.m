%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos] = resolvePos(dVar, indices) 
%Finds the position in the big Vector of DVar Coefficients
%input: dVar, indices
%Pass the var dVar number and the slrijt positions in the "indices". Don't
%omit 0 for indices that were not used.
%Return value integer pointing to the position in vec
%
%Returns -1 if there is something wrong.

%First let's find the position form the indices
global varLen; 
global dVarRanges; global maxV_;
global S; global L; global R; global Ts;
pos = 1;
interval = 1;

if(dVar == 1)
    interval = dVarRanges(1);
    a = [S,L,R,0,0,Ts];
    b = 0;  
elseif(dVar>1 && dVar<= length(dVarRanges))
    a = maxV_;
    b = dVarRanges(dVar-1);
    interval = dVarRanges(dVar)-dVarRanges(dVar-1);
elseif (dVar<1 || dVar > length(dVarRanges))
    pos = -1; %that serves as error code.
    return;
end

for count= varLen:-1:1
        if(indices(count) ~= 0)
            interval = interval./a(count);
            pos = pos + interval*(indices(count)-1);
        end
end

pos =  b + pos;


end