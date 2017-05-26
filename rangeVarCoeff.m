%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm] = rangeVarCoeff(dVar, indicesRange, tStart, coeff)
%function [vec_] = rangeVarCoeff(dVar, indicesRange, tStart, coeff)
%Finds the coefficient vectors for a given DVar and a given indices range.

%global vecLen;

%calculating how many different indices for dVar to produce
% nzi = nonzeros(indicesRange);
% nzi(end) = nzi(end)-tStart+1;
% totalRange = prod(nzi)

indicesArr = generateIndices(indicesRange, tStart);
[tot ,~] = size(indicesArr);

%vec_ = zeros(totalRange,vecLen);
spm = zeros(tot,3); %sparse matrix
%ind = indices; ind(ind~=0)=1;

for iter = 1:tot
    la = indicesArr(iter,:);
    %aa = singleVarCoeff(dVar, la, coeff);
    pos = resolvePos(dVar, la);
    if (pos == -1)
        display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
        return;
    end
     %display(aa(1:10));
    %vec_(iter,:) = aa;
    spm(iter,:) = [iter,pos,coeff];
end
end

%%%%%%%%%%%%%%%%%%%%%%%
function [vec_] = singleVarCoeff(dVar, indices, coeff)
%Finds the coefficient vectors for a given DVar and a given indices 

global vecLen;
vec_ = zeros(1,vecLen);
pos = resolvePos(dVar, indices);
if (pos == -1)
    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
    return;
end
vec_(pos) = coeff;

end