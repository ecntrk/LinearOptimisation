%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allvar, u, v, w, wR, x0, x, xR, y0, y, yR, z] = optimise(loadfile)
% This function generates the coefficient vectors and optimises it by CPLEX.
% Usage [x] = optimise(loadFile); loadfile is a mat file for a variable
% environment. If you don't have it, put noting. It will take values from 
% the inputScene() function.
%
% Tip: Make sure you have your input right at inputScene. You can also 
% load an environment my a mat file. 

global maxV_;

if (nargin == 0)
    inputScene();
else
    load(loadfile);
end

init();
[a,b,f, ae,be,ane,bne] = construct();

X = cplexmilp(f,ane,bne,[a{1};a{2};a{3};a{4};a{6};a{7};a{8}], [b{1};b{2};b{3};b{4};b{6};b{7};b{8}]);

 x1(1,:) = find(X);
 x1(2,:) = X(X~=0);
 
 allvar = x1';

%doing this for each decision variable
iter = 1;
%for u
u = zeros(maxV_(1),maxV_(2),maxV_(3),maxV_(6));
for t = 1:maxV_(6)
   for r = 1: maxV_(3)
    for l = 1: maxV_(2)
        for s = 1: maxV_(1)
            u(s,l,r,t) = X(iter);
            iter = iter+1;
        end
    end
   end
end
display(iter);
%for v
v = zeros(maxV_(1),maxV_(3),maxV_(4),maxV_(6));
for t = 1:maxV_(6)
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
        for s = 1: maxV_(1)
            v(s,r,i,t) = X(iter);
            iter = iter+1;
        end
    end
   end
end
display(iter);

%for w
w = zeros(maxV_(3),maxV_(4));
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
            w(r,i) = X(iter);
            iter = iter+1;
    end
   end
display(iter);

%for w_bar
wR = zeros(maxV_(3),maxV_(4));
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
            wR(r,i) = X(iter);
            iter = iter+1;
    end
   end
   
display(iter);
   
%for x_0 (s,i,j,t)
x0 = zeros(maxV_(1),maxV_(4),maxV_(5),maxV_(6));
for t = 1:maxV_(6)
   for j = 1: maxV_(5)
    for i = 1: maxV_(4)
        for s = 1: maxV_(1)
            x0(s,i,j,t) = X(iter);
            iter = iter+1;
        end
    end
   end
end
display(iter);
   
%for x (srijt)
x = zeros(maxV_(1),maxV_(3),maxV_(4),maxV_(5),maxV_(6));
for t = 1:maxV_(6)
   for j = 1: maxV_(5)
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
        for s = 1: maxV_(1)
            x(s,r,i,j,t) = X(iter);
            iter = iter+1;
        end
    end
   end
   end
end
  
display(iter);

%for x_bar (srijt)
xR = zeros(maxV_(1),maxV_(3),maxV_(4),maxV_(5),maxV_(6));
for t = 1:maxV_(6)
   for j = 1: maxV_(5)
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
        for s = 1: maxV_(1)
            xR(s,r,i,j,t) = X(iter);
            iter = iter+1;
        end
    end
   end
   end
end
display(iter);

%for y0
y0 = zeros(maxV_(1),maxV_(4),maxV_(6));
for t = 1:maxV_(6)
   for i = 1: maxV_(4)
        for s = 1: maxV_(1)
            y0(s,i,t) = X(iter);
            iter = iter+1;
        end
   end
end
display(iter);

%for y
y = zeros(maxV_(1),maxV_(3),maxV_(4),maxV_(6));
for t = 1:maxV_(6)
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
        for s = 1: maxV_(1)
            y(s,r,i,t) = X(iter);
            iter = iter+1;
        end
    end
   end
end
display(iter);


%for y_bar
yR = zeros(maxV_(1),maxV_(3),maxV_(4),maxV_(6));
for t = 1:maxV_(6)
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
        for s = 1: maxV_(1)
            yR(s,r,i,t) = X(iter);
            iter = iter+1;
        end
    end
   end
end
display(iter);

%for z
z = zeros(maxV_(1),maxV_(3),maxV_(4),maxV_(6));
for t = 1:maxV_(6)
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
        for s = 1: maxV_(1)
            z(s,r,i,t) = X(iter);
            iter = iter+1;
        end
    end
   end
end
display(iter);



end
