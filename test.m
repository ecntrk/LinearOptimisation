cObj = Cplex;


a = full(A);
r = full(rhs);

cObj.Model.name = 'disaster1';
%cObj.Model.sense = 'minimize';
cObj.Model.obj = f';
cObj.Model.lb = lb;
cObj.Model.ub = ub;
cObj.Model.A = A;
cObj.Model.lhs = lhs;
cObj.Model.rhs = rhs;

cObj.writeModel('disaster1.lp');

options = cplexoptimset('Display', 'on');
%X = cObj.solve();



%    % Initialize the CPLEX object
%    cplex = Cplex('dispconflict');
%    cplex.Model.sense = 'maximize';
%    
%    % Now populate the problem with the data
%    % Use arrays to populate the model
%       cplex.Model.obj   = [ 1;   2;   3; 1];
%       cplex.Model.lb    = [ 0;   0;   0; 2];
%       cplex.Model.ub    = [40; inf; inf; 3];
%       cplex.Model.ctype = 'CCCI';
%       cplex.Model.A =  [-1  1  1  10;
%          1 -3  1   0;
%          0  1  0  -3.5;
%          0  1 1  0];
%       cplex.Model.lhs = [-inf; -inf; 0; -inf];
%       cplex.Model.rhs = [  20;   30; 0; 0];
%       cplex.writeModel('test.lp');