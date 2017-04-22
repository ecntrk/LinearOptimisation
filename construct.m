%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_] = construct()
%Constructs coefficient vectors for a group of linear functions.
%Var Order: i, j, l, r, s, t
%Decision Var Order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z


%initialising all global vars
init();

%vec_ = rangeVarCoeff(1, [0,0,3,3,3,3], 1);
%a= generateIndices([2,0,3,5,2,3,4]);

%constructing equations one by one:
%eq1();
%eq2();
%eq3();
vec_ = eq4();
end

function eq1()
%constructign equation 10.a
%y0 (i,s, t=1) = 0

global maxV_;
rangeVarCoeff(8,[maxV_(1),0,0,0,maxV_(5),1], 1, 1);
end


function eq2()
%constructign equation 10.b
%y (i,r,s, t=1) = 0

global maxV_;
rangeVarCoeff(9,[maxV_(1),0,0,maxV_(4),maxV_(5),1], 1, 1);
end



function eq3()
%constructign equation 10.c
%ybar (i,r,s, t=1) = 0

global maxV_;
rangeVarCoeff(10,[maxV_(1),0,0,maxV_(4),maxV_(5), 1],1,1);
end

function [vec_] = eq4()
%constructign equation 11
%y (i,r,s, t+1) - y (i,r,s, t+1) - sum(x ijrst) = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t+1>sr+tij).
global maxV_; global vecLen; global whatJ;

sr = 0; tij =0;
tCondition = (1+sr+tij)+1;
indArr = generateIndices([maxV_(1),0,0,maxV_(4),maxV_(5),maxV_(6)], 1); 

%making the coeff vector
vec_ = zeros(length(indArr),vecLen);

%disp(indArr);
for iter = 1:length(indArr)
temp = zeros(1,vecLen);
%calculating positions for y irs(t), y irs(t+1)
[p1, p2] = doubleVarCoeff(9, indArr(iter,:));
temp(p1) = 1;
temp(p2) = -1;

%calculaitng the positions in the sum
%only then t condition is met!
if(indArr(iter,6) >= tCondition)
    %do the actual sum
    for j = 1:whatJ(indArr(iter,1)) %calculates j from i conditionally
        arr = indArr(iter,:);
        arr(2) = j; %iterating over j
        pos = resolvePos(6, arr);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
            return;
        end
        temp(pos) = 1;
    end
end

vec_(iter,:) = temp;

end
end


