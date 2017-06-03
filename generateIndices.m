%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[att]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indArr] = generateIndices(iR, tSt)
%Clever trick to generate variable length indices combinations
%sequentially. Can't be made parallel unfortunately right now.
%
%Input: indices range values at correct places, tStart  is the starting of
%variable t (because of conditions like t+1>sr+tij

%extracting non zero elements form range
nzi = nonzeros(iR);

%making the starting positions
iRs = iR;
iRs(iRs>0)=1;
iRs(end) = tSt;
%guide index for position of the nonzero index elements
gInd = zeros(1, nnz(iR));

count = 1;
for c=1:length(iR)
    if(iR(c) >0)
        gInd(count) = c;
        count = count+1;
    end
end

%This many combinations to produce
nzi(end) = nzi(end)-tSt+1;
tot = prod(nzi);
%initialising the output array of indices combinations
indArr = zeros(tot,length(iR)); 

%display (gInd); %debug stuff

%Here goes the actual logic
%interval = 1;
interval = tot;
%we iterate coloumnwise
for c = 1:length(gInd)
    %display(c);
    maxVal = iR(gInd(c));
    minVal  = iRs(gInd(c));
    %gInd(c)
%    if(c<1)
        %interval = interval * (iR(gInd(c-1)));%-iRs(gInd(c-1)))
        interval = interval/maxVal;
 %   end
    cc =1;
    %iterate val with 
    asd = tot/((maxVal-minVal+1)*interval);
    for rowSeg = 1:asd
        for val = minVal:maxVal
           for in = 1:interval
               indArr(cc,gInd(c)) = val;
               cc=cc+1; %faster than calculating index with rowSeg.
           end
        end
    end
    
end
%display(indArr);
%eliminating indices where i = j. simpler than iterating in range.
if (iR(4) ~= 0 && iR(5) ~= 0)
[rr,~] = size(indArr);
    for count = 1:rr
        if (indArr(count,4) == indArr(count,5)) 
            indArr(count,:) = [];
        end
    end
end

end


