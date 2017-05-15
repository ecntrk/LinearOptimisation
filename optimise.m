%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[at]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = optimise(loadfile)
% This function generates the coefficient vectors and optimises it by CPLEX.
% Usage [x] = optimise(loadFile); loadfile is a mat file for a variable
% environment. If you don't have it, put noting. It will take values from 
% the inputScene() function.
%
% Tip: Make sure you have your input right at inputScene. You can also 
% load an environment my a mat file. 

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

x = x1';
end
