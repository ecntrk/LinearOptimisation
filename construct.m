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
%eq11();
%eq12();
%eq13();
%eq14();
vec_ = eq15();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq1()
%constructign equation 10.a
%y0 (i,s, t=1) = 0

    global maxV_;
    rangeVarCoeff(8,[maxV_(1),0,0,0,maxV_(5),1], 1, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq2()
%constructign equation 10.b
%y (i,r,s, t=1) = 0

    global maxV_;
    rangeVarCoeff(9,[maxV_(1),0,0,maxV_(4),maxV_(5),1], 1, 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq3()
%constructign equation 10.c
%ybar (i,r,s, t=1) = 0

    global maxV_;
    rangeVarCoeff(10,[maxV_(1),0,0,maxV_(4),maxV_(5), 1],1,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_] = eq11()
%constructign equation 11
%y (i,r,s, t+1) - y (i,r,s, t+1) - sum(x ijrst) = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t+1>sr+tij).
    global maxV_; global vecLen; global whatJ;

    global sr ; global tij;
    tCondition = (1+sr+tij)+1;
    indArr = generateIndices([maxV_(1),0,0,maxV_(4),maxV_(5),maxV_(6)], 1); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        %calculating positions for y irs(t), y irs(t+1)
        [p1, p2] = doubleVarCoeff(9, indArr(iter,:));
        temp(p1) = -1;
        temp(p2) = 1;

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
                temp(pos) = -1;
            end
        end

        vec_(iter,:) = temp;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_] = eq12()
%constructign equation 12
%y (i,r,s, t+1) - y (i,r,s, t+1) - sum(x_0 ijst) - sum (sum (betaR * x (ijrst) ))= 0

%generate indices combinations for i,s,t and then iterate j,r for sum(x).
%put combination (t+1>tij).
    global maxV_; global vecLen; global whatJ;

    global tij;
    tCondition = (1+tij)+1;
    indArr = generateIndices([maxV_(1),0,0,0,maxV_(5),maxV_(6)], 1);

    %having the coefficients for x ijrs(t)
    global betaR;

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    %building the coeff vector for each of i,s,t
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        %calculating positions for y_0 is(t), y_0 is(t+1)
        [p1, p2] = doubleVarCoeff(8, indArr(iter,:));
        temp(p1) = -1;
        temp(p2) = 1;

        %calculaitng the positions in the sum
        %only then t condition is met!
        if(indArr(iter,6) >= tCondition)
            %do the actual sum over j
            for j = 1:whatJ(indArr(iter,1)) %calculates j from i conditionally
                arr = indArr(iter,:);

                %ifirst, single sub over j for x_0.
                arr(2) = j; %iterating over j 

                pos = resolvePos(5, arr);
                if (pos == -1)
                    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
                    return;
                end
                temp(pos) = -1;

                %now the double sum of j and r for x
                for r = 1:maxV_(4)
                    arr(4) = r; %iterating over r
                    pos = resolvePos(6, arr);
                    if (pos == -1)
                        display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
                        return;
                    end
                    temp(pos) = -betaR(r);
                end

            end
        end

        vec_(iter,:) = temp;

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_] = eq13()
%constructign equation 13
%yBar (i,r,s, t+1) - yBar (i,r,s, t+1) - sum(xBar (ijrst) + sum( alphaR * x(ijrst)) ) = 0

%generate indices combinations for i,r,s,t and then iterate r for sum(xbar + alpa * x).
%put combination (t+1>tij).

    global maxV_; global vecLen; global whatJ;

    global tij;
    tCondition = (1+tij)+1;
    indArr = generateIndices([maxV_(1),0,0,maxV_(4),maxV_(5),maxV_(6)], 1);

    %having the coefficients for x ijrs(t)
    global alphaR;

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    %building the coeff vector for each of i,s,t
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        %calculating positions for yBar irs(t), yBar irs(t+1)
        [p1, p2] = doubleVarCoeff(10, indArr(iter,:));
        temp(p1) = -1;
        temp(p2) = 1;

        %calculaitng the positions in the sum
        %only then t condition is met!
        if(indArr(iter,6) >= tCondition)
            %do the actual sum over j
            for j = 1:whatJ(indArr(iter,1)) %calculates j from i conditionally
                arr = indArr(iter,:);

                %ifirst, single sub over j for xBar.
                arr(2) = j; %iterating over j 

                pos = resolvePos(7, arr);
                if (pos == -1)
                    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
                    return;
                end
                temp(pos) = -1;

                %for x
                pos = resolvePos(6, arr);
                if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
                return;
                end
                %indArr(iter,4) is the current value of r
                temp(pos) = -alphaR(indArr(iter,4)); 

            end
        end

        vec_(iter,:) = temp;

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_] = eq14()
%constructign equation 14
%z (i,r,s, t) - y (i,r,s, t+1) - sum(x ijrst) = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t+1>sr+tij).
    global maxV_; global vecLen; global whatJ;

    global sr ; global tij;
    tCondition = (1+sr+tij)+1;
    indArr = generateIndices([maxV_(1),0,0,maxV_(4),maxV_(5),maxV_(6)], 1); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        %calculating position for z irs(t)
        pos = resolvePos(11, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
            return;
        end
        temp(pos) = 1;
        
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
                temp(pos) = -1;
            end
        end

        vec_(iter,:) = temp;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_] = eq15()
%constructign equation 15
%v (i,r,s, t) - sum z (i,r,s, t cond) - sum z (i,r,s, t cond)  = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t+1>sr+tij).
    global maxV_; global vecLen;

    Td = 2;

    indArr = generateIndices([maxV_(1),0,0,maxV_(4),maxV_(5),maxV_(6)], 1); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        %calculating position for v irs(t)
        pos = resolvePos(2, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
            return;
        end
        temp(pos) = 1;
        
        %calculaitng the positions in the conditional sum
        %tau = 0 .. min(t, Td)
        t = indArr(iter,6);
        for tau = 0:t-1 %calculates j from i conditionally
            arr = indArr(iter,:);
            arr(6) = t-tau; %iterating over (t-tau)
            pos = resolvePos(6, arr);
            if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indices(:)) ));
                return;
            end
            %here comes the condition
            if (tau <= min(t,Td)-1)
                temp(pos) = -1;
            else
                temp(pos) = -2;
            end
        end

        vec_(iter,:) = temp;

    end
end
