%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ u, v, w, wR, x0, x, xR, y0, y, yR, z] = optimise(loadfile)
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
tic;
disp('generating coefficient for S,L,R,I,J,T: ');disp(maxV_)
%[f, ae, be, ane, bne] = construct();
[f, A, lhs, rhs, lb, ub] = constructSPM();
disp('coeff matrix generation complete. Optimising now...')
toc
%starting optimisaiton
tic;
%defining the Cplex object
cObj = Cplex;

cObj.Model.name = 'disaster1';
cObj.Model.sense = 'minimize';
cObj.Model.obj = f;
cObj.Model.A = A;
cObj. Model.lhs = lhs;
cObj.Model.rhs = rhs;
cObj.Model.lb = lb;
cObj.Model.ub = ub;
cObj.writeModel('disaster1.lp');

%This is a mixed integer proble and CPLEX automatically detects it (becuase
%we have intergers and no quadratic term)
options = cplexoptimset('Display', 'on');
ansObj = cObj.solve();

X = ansObj.x;

%X = cplexmilp(f,ane,bne,ae,be);





disp('optimisation complete')
toc
% 
%  x1(1,:) = find(X);
%  x1(2,:) = X(X~=0);
%  
%  allvar = x1';

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
%disp('number of u: ',iter);
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

%disp('number of v: ',iter);

%for w
w = zeros(maxV_(3),maxV_(4));
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
            w(r,i) = X(iter);
            iter = iter+1;
    end
   end
%disp(['number of w: ',iter]);

%for w_bar
wR = zeros(maxV_(3),maxV_(4));
   for i = 1: maxV_(4)
    for r = 1: maxV_(3)
            wR(r,i) = X(iter);
            iter = iter+1;
    end
   end
   
%disp('number of w_bar: ',iter);
   
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

%disp('number of x0: ',iter);
   
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
  
%disp('number of x: ',iter);

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

%disp('number of x_bar: ',iter);

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

%disp('number of y0: ',iter);

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
%disp('number of y: ',iter);


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
%disp('number of y_bar: ',iter);

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
%disp('number of z: ',iter);



end
