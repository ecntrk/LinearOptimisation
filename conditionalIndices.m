%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[att]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [arr] = conditionalIndices(vec)
%generates variable length indices combinations
%sequentially. Can't be made parallel unfortunately right now.
%
%Input: vector[s,l,r,i,j,t]. put 1 in the indices wanted. 0 otherwise.
%t=0, no t. t=1 refers to Tkl. t = 2 refers to full iteration of Ts.

%S is always there, i,j will never be there.
global S;
global L_s;
global whatKL; %determines k from l.
global R_k;
global T_k;
global Ts;
%arr = [0,0,0,0,0,0];
count = 1;
%display (S);
for Scn = 1:S 
    %display(Scn);
    Lrange = L_s{Scn};
    for Lcn = Lrange;
        kk = whatKL(Lcn);
        Rrange = R_k(kk,:);
        if(vec(3) == 0)
            Rrange = [0];
        end
        for Rcn = Rrange
            if (vec(6) == 0)
                Trange = [0];
            elseif(vec(6) == 1)
                Trange = T_k{kk}; %T_k is cell
            else
                Trange = 1:Ts(Scn);
            end 
            
            for Tcn = Trange
                temp = [Scn, Lcn, Rcn, 0, 0, Tcn];
                arr(count,:) = temp;
                count = count + 1;
            end
        end
        
    end
end

%display(indArr);
 
end


