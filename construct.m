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
eq1();
eq2();
eq3();

end

function eq1()
%constructign equation 10.a
%y0 (i,s, t=1) = 0

global maxV_;
rangeVarCoeff(8,[maxV_(1),0,0,0,maxV_(5), 1],1);
end


function eq2()
%constructign equation 10.b
%y (i,r,s, t=1) = 0

global maxV_;
rangeVarCoeff(9,[maxV_(1),0,0,maxV_(4),maxV_(5), 1],1);
end



function eq3()
%constructign equation 10.c
%ybar (i,r,s, t=1) = 0

global maxV_;
rangeVarCoeff(10,[maxV_(1),0,0,maxV_(4),maxV_(5), 1],1);
end

function eq4()
%constructign equation 11
%y (i,r,s, t+1) - y (i,r,s, t+1) - sum(x ijrst) = 0

global maxV_;
rangeVarCoeff(10,[maxV_(1),0,0,maxV_(4),maxV_(5), 1],1);
end
