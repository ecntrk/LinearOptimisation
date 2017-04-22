%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
% Finds the coefficient vectors for a given DVar and a given indices  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos1, pos2] = doubleVarCoeff(dVar, indices)
%Finds the coefficient vectors for a given DVar and a given indices 

%for t
pos1 = resolvePos(dVar, indices);
if (pos1 == -1)
    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
    return;
end
indices(6) = indices(6)+1;

%for t+1
pos2 = resolvePos(dVar, indices);
if (pos2 == -1)
    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
    return;
end

end


