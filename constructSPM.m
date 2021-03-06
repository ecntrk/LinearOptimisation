%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[att]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [f, aeq, beq, aneq, bneq] = constructSPM()

function [f, A, lhs, rhs, lb, ub] = constructSPM()
%Constructs coefficient vectors for a group of linear functions.
%Var Order: i, j, l, r, s, t
%Decision Var Order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z


%initialising all global vars
%inputScene(); %for debug
%init(); %for debug

global vecLen;
global dVarRanges;

global Naeq;
global Naneq;
Naeq = 0;
Naneq = 0;

%constructs the main equation to be minimised
%f = sparse=(eqnF());
f = eqnF();


% %constructing equations one by one:

%
%block for the = equations
%
spm = eq1();
spmat = spm;
tt = Naeq

spm = eq2();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
tt = Naeq+tt

spm = eq3();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
tt = Naeq+tt

spm = eq11();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
tt = Naeq+tt

spm = eq12();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
tt = Naeq+tt

spm = eq13();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
tt = Naeq+tt

spm = eq14();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
tt = Naeq+tt

spm = eq15();
spm(:,1) = spm(:,1) + tt;
spmat = [spmat;spm];
tt = Naeq+tt


%Starting the >= equations


tt1 = tt %to identify number of >= equaitons

[spm, bspm] = eq16();
spm(:,1) = spm(:,1) + tt; 
spm(:,3) = spm(:,3)*-1; 
bspm(:,3) = bspm(:,3)*-1; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = bspm;
tt = Naneq+tt

[spm, bspm] = eq17();
spm(:,1) = spm(:,1) + tt; 
spm(:,3) = spm(:,3)*-1; 
bspm(:,3) = bspm(:,3)*-1; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = [bspmat;bspm];
tt = Naneq+tt

[spm, bspm] = eq18();
spm(:,1) = spm(:,1) + tt; 
spm(:,3) = spm(:,3)*-1; 
bspm(:,3) = bspm(:,3)*-1; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = [bspmat;bspm];
tt = Naneq+tt


[spm, bspm] = eq19();
spm(:,1) = spm(:,1) + tt; 
spm(:,3) = spm(:,3)*-1; 
bspm(:,3) = bspm(:,3)*-1; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = [bspmat;bspm];
tt = Naneq+tt

[spm, bspm] = eq20();
spm(:,1) = spm(:,1) + tt; 
spm(:,3) = spm(:,3)*-1; 
bspm(:,3) = bspm(:,3)*-1; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = [bspmat;bspm];
tt = Naneq+tt

%moving onto <=


tt2 = tt;


[spm, bspm] = eq21();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = [bspmat;bspm];
tt = Naneq+tt


[spm, bspm] = eq22();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = [bspmat;bspm];
tt = Naneq+tt

[spm, bspm] = eq23();
spm(:,1) = spm(:,1) + tt; 
spmat = [spmat;spm];
bspm(:,1) = bspm(:,1) + tt; 
bspmat = [bspmat;bspm];
tt = Naneq+tt



%
%Equation construction complete. moving to last adjustments

lb = -Inf(vecLen,1);
ub = Inf(vecLen,1);
%here, w, wbar, x0, x, xbar all are positive integers. so putting lb = 1
%for them
lb(dVarRanges(2)+1 : dVarRanges(7), 1) = 0;


temp = zeros(tt1,1); 
temp1 = -Inf(tt-tt1,1);
lhs = [temp;temp1];
rhs = [temp;bspmat(:,3)];




%lhs = sparse(lhs(:,1), lhs(:,2), lhs(:,3), tt,1);
%rhs = sparse(rhs(:,1), rhs(:,2), rhs(:,3), tt,1);



A = sparse(spmat(:,1)',spmat(:,2)',spmat(:,3)',tt,vecLen);




end



%
% Starting function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [svec] = eq1()
%constructing equation 10.a
%y0 (i,s, t=1) = 0

    global maxV_;
    svec = rangeVarCoeff(8,[maxV_(1),0,0,maxV_(4),0,1], 1, 1);
    %[r,~] = size(vec);
    %svec2 = zeros(r,1); %RHS is zero
    global Naeq;
    [r,~] = size(svec);
    Naeq = r;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [svec]=eq2()
%constructign equation 10.b
%y (i,r,s, t=1) = 0

    global maxV_;
    svec = rangeVarCoeff(9,[maxV_(1),0,maxV_(3),maxV_(4),0,1], 1, 1);
    
    
    global Naeq;
    [r,~] = size(svec);
    Naeq = r;
    %vec2 = zeros(r,1); %RHS is zero
    

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [svec] = eq3()
%constructign equation 10.c
%ybar (i,r,s, t=1) = 0

    global maxV_;
    svec = rangeVarCoeff(10,[maxV_(1),0,maxV_(3),maxV_(4),0, 1],1,1);
    
    %[r,~] = size(vec);
    global Naeq;
    [r,~] = size(svec);
    Naeq = r;

    %vec2 = zeros(r,1); %RHS is zero
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm] = eq11()
%constructign equation 11
%y (i,r,s, t+1) - y (i,r,s, t) - sum(x ijrst) = 0

%generate indices combinations for i,r,s,t = (Ts-1) and then iterate j for sum(x).

    global maxV_; 
    %global vecLen; 

    global sr ; global tij; global epsilon_i;
    %tCondition = (1+sr(r)+tij(i,j)+1;
    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)-1], 1); 


    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %total coeff on one row
    cN = 1 + 1 + maxV_(5);
    spm = zeros(rN*cN,3);
    num = 1;
    
    global Naeq;
    Naeq = rN;

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);
        %calculating positions for y irs(t), y irs(t+1)
        [p1, p2] = doubleVarCoeff(9, indArr(iter,:));
        %temp(p1) = -1;
        %temp(p2) = 1;
        spm(num,:)=[iter,p1,-1]; num= num+1;
        spm(num,:)=[iter,p2,1]; num= num+1;
        
        %calculaitng the positions in the sum
            for j = epsilon_i(indArr(iter,4),:) %calculates j from i conditionally
                    
                %put combination (t+1>sr+tij). 
                condition = sr(indArr(iter,3)) + tij(indArr(iter,4),j);
                if ((indArr(iter,6)+1) > condition)         %only then t condition is met!

                    %do the actual sum
                    arr = indArr(iter,:);
                    arr(5) = j; %iterating over j
                    arr(6) = arr(6) + 1 - condition;
                    pos = resolvePos(6, arr);
                    if (pos == -1)
                        display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                        return;
                    end
                    %temp(pos) = -1;
                    spm(num,:)=[iter,pos,-1]; num= num+1;

                end

            end
        %vec_(iter,:) = temp;

    end
    
    %vec2 = zeros(length(indArr),1); %RHS is zero
    %
    spm( ~any(spm,2), : ) = [];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm] = eq12()
%constructign equation 12
%y_0 (i,s, t+1) - y_0 (i,s, t) - sum(x_0 ijst) - sum (sum (betaR * x (ijrst) ))= 0

%generate indices combinations for i,s,t and then iterate j,r for sum(x).
%put combination (t+1>tij).
    global maxV_; 
    %global vecLen; 
    global epsilon_i;
    global tij;

    %having the coefficients for x ijrs(t)
    global betaR;

    indArr = generateIndices([maxV_(1),0,0,maxV_(4),0,maxV_(6)-1], 1);

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %total coeff on one row
    cN = 1 + 1 + maxV_(5) + maxV_(5)*maxV_(3);
    spm = zeros(rN*cN,3);
    num = 1;

    global Naeq;
    Naeq = rN;

    %disp(indArr);
    %building the coeff vector for each of i,s,t
    for iter = 1:rN
        %temp = zeros(1,vecLen);
        %calculating positions for y_0 is(t), y_0 is(t+1)
        [p1, p2] = doubleVarCoeff(8, indArr(iter,:));
        %temp(p1) = -1;
        %temp(p2) = 1;
                    
        spm(num,:)=[iter,p1,-1]; num= num+1;
        spm(num,:)=[iter,p2,1]; num= num+1;

        %calculaitng the positions in the sum
        %only then t condition is met!

            %do the actual sum over j
            for j = epsilon_i(indArr(iter,4), :) %calculates j from i conditionally
                %display(j);
                %put combination (t+1>tij). 
                condition = tij(indArr(iter,4),j);
                if ((indArr(iter,6)+1) > condition)         %only then t condition is met!

                    arr = indArr(iter,:);

                    %ifirst, single sum over j for x_0.
                    arr(5) = j; %iterating over j 
                    %updating t
                    arr(6) = arr(6)+1-condition;
                    pos = resolvePos(5, arr);
                    if (pos == -1)
                        display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                        return;
                    end
                    %temp(pos) = -1;
                    spm(num,:)=[iter,pos,-1]; num= num+1;

                    %now the double sum of j and r for x
                    for r = 1:maxV_(3)
                        arr(3) = r; %iterating over r
                        pos = resolvePos(6, arr);
                        if (pos == -1)
                            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                            return;
                        end
                        %temp(pos) = -betaR(r);
                        spm(num,:)=[iter,pos, -betaR(r)]; num= num+1;

                    end
                end
            end

        %vec_(iter,:) = temp;

    end
    
    %vec2 = zeros(length(indArr),1); %RHS is zero
    % stripping off extra rows with zeros 
    spm( ~any(spm,2), : ) = [];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm] = eq13()
%constructign equation 13
%yBar (i,r,s, t+1) - yBar (i,r,s, t) - sum(xBar (ijrst) + sum( alphaR * x(ijrst)) ) = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(xbar + alpa * x).
%put combination (t+1>tij).

    global maxV_; 
    %global vecLen; 
    global epsilon_i;
    global tij;

    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)-1], 1); %t = 1..Ts-1

    %having the coefficients for x ijrs(t)
    global alphaR;

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %total coeff on one row
    cN = 1 + 1 + maxV_(5) + maxV_(5);
    spm = zeros(rN*cN,3);
    num = 1;

    global Naeq;
    Naeq = rN;

    %disp(indArr);
    %building the coeff vector for each of i,s,t
    for iter = 1:rN
        %temp = zeros(1,vecLen);
        %calculating positions for yBar irs(t), yBar irs(t+1)
        [p1, p2] = doubleVarCoeff(10, indArr(iter,:));
        %temp(p1) = -1;
        %temp(p2) = 1;
        
        spm(num,:)=[iter,p1,-1]; num= num+1;
        spm(num,:)=[iter,p2,1]; num= num+1;

        %calculaitng the positions in the sum
        %only then t condition is met!
        %do the actual sum over j
        for j = epsilon_i(indArr(iter,4),:) %calculates j from i conditionally
             %put combination (t+1>tij). 
            condition = tij(indArr(iter,4),j);
            if ((indArr(iter,6)+1) > condition)         %only then t condition is met!

                arr = indArr(iter,:);
                arr(5) = j; %iterating over j 
                %updating t
                arr(6) = arr(6)+1-condition;

                %ifirst, single sub over j for xBar.
                pos = resolvePos(7, arr);
                if (pos == -1)
                    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                    return;
                end
                %temp(pos) = -1;
               spm(num,:)=[iter,pos,-1]; num= num+1;

                %for x
                pos = resolvePos(6, arr);
                if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                return;
                end
                %indArr(iter,4) is the current value of r
                %temp(pos) = -alphaR(indArr(iter,3)); 
                spm(num,:)=[iter,pos,-alphaR(indArr(iter,3))]; num= num+1;

            end

        end

        %vec_(iter,:) = temp;

    end
    
    %vec2 = zeros(length(indArr),1); %RHS is zero

    % stripping off extra rows with zeros 
    spm( ~any(spm,2), : ) = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm] = eq14()
%constructign equation 14
%z (i,r,s, t)  - sum(x ijrst) = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t>sr+tij), not t+1.
    global maxV_; 
    %global vecLen; 
    global epsilon_i;

    global sr ; 
    global tij;
    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)], 1); %t = 1..Ts

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %total coeff on one row
    cN = 1 +  maxV_(5);
    spm = zeros(rN*cN,3);
    num = 1;

    global Naeq;
    Naeq = rN;

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);
        %calculating position for z irs(t)
        pos = resolvePos(11, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;

        %calculaitng the positions in the sum
        %only then t condition is met!
        
        %do the actual sum
        for j = epsilon_i(indArr(iter,4),:) %calculates j from i conditionally
            %put combination (t+1>sr+tij). 
            condition = sr(indArr(iter,3)) + tij(indArr(iter,4),j);
            if ((indArr(iter,6)) > condition)         %only then t condition is met!
                arr = indArr(iter,:);
                arr(5) = j; %iterating over j
                arr(6) = arr(6)-condition; %updating t
                %display(arr);
                pos = resolvePos(6, arr);
                if (pos == -1)
                    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                    return;
                end
                %display (pos);
                %temp(pos) = -1;
                spm(num,:)=[iter,pos,-1]; num= num+1;

            end
        end

        %vec_(iter,:) = temp;

    end
    %vec2 = zeros(length(indArr),1); %RHS is zero
    
    % stripping off extra rows with zeros
    spm( ~any(spm,2), : ) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm] = eq15()
%constructign equation 15
%v (i,r,s, t) - sum z (i,r,s, t cond) - sum z (i,r,s, t cond)  = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t+1>sr+tij).
    global maxV_; 
    %global vecLen;
    global Td;


    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)], 1); 

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %total coeff on one row
    cN = 1 + maxV_(6)+ maxV_(6);
    spm = zeros(rN*cN,3);
    num = 1;

    global Naeq;
    Naeq = rN;

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);
        %calculating position for v irs(t)
        pos = resolvePos(2, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;

        %calculaitng the positions in the conditional sum of z
        %tau = 0 .. min(t, Td)
        t = indArr(iter,6);
        for tau = 0:t-1 %calculates j from i conditionally
            arr = indArr(iter,:);
            arr(6) = t-tau; %iterating over (t-tau)
            pos = resolvePos(11, arr);
            if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                return;
            end
            %here comes the condition
            if (tau <= min(t,Td)-1)
                %temp(pos) = -1;
                spm(num,:)=[iter,pos,-1]; num= num+1;

            else
                %temp(pos) = -2;
                spm(num,:)=[iter,pos,-2]; num= num+1;

            end
        end

        %vec_(iter,:) = temp;

    end
    
    %vec2 = zeros(length(indArr),1); %RHS is zero
    
    % stripping off extra rows with zeros
    spm( ~any(spm,2), : ) = [];
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm, bspm] = eq16()
%constructign equation 16
%y_0 (il, t) >= d_0 (l,t) 

%generate indices combinations for l,s,t
%although l depends on s, it'll be dealt in conditionalIndices()
%it  will generate the correct number of l for each s! Don't worry!
%
%il is a constant depnding on l so don't need to iterate i.


    %global vecLen; 
    global whatIL; 
    global d0_lt;
    %global maxV_;

    indArr = conditionalIndices([1,1,0,0,0,1]); %unlike generateindices, you put 1 on values you want.

 

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN = 1 ;
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    bnum = 1;
        
    global Naneq;
    Naneq = rN;
    
    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);
        
        %fill it with il wrt l
        ind = indArr(iter,:);
        ind(4) = whatIL(ind(2));
        ind(2) = 0; %we don't need l now
        
        %calculating position for y_0 (il,s,t)
        pos = resolvePos(8, ind);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',ind(:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;

        
        %last resort for eqn 12

        %This line is causing conflict with Eqn 12. 
        %adding generalists (d0_lt) at every time tick is not feasibly
        %Assuming generalists as cumulative solves the problem
        %so, for l = 1, d0_lt = 0,0,0,200.
        %so, for l = 9, d0_lt = 0,0,150.
        
        if (ind(6) == 20) %only at the last time tick (t=20, add generalists
            bspm(bnum,:)=[iter,1,d0_lt(indArr(iter,2))]; bnum= bnum+1;
        else %all the other time ticks add 0.
            bspm(bnum,:)=[iter,1,0]; bnum= bnum+1;
        end
                
        
        
        %lets deal with the d (l,t)
        %vec2(iter) = d0_lt(indArr(iter,2)); %possible change when there more t.
        
        %vec_(iter,:) = temp;

         
    end
    
    % stripping off extra rows with zeros
    spm( ~any(spm,2), : ) = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm, bspm] = eq17()
%constructign equation 17
%y_bar (il,r,s,t) - v (il,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    %global vecLen; 
    global whatIL;
    
    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,2]); % t=2 means it'll iterate over Ts not T_kl

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN = 1 + 1;
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    %bnum = 1;

    global Naneq;
    Naneq = rN;
    

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);
        
        %fill it with il wrt l
        ind = indArr(iter,:);
        ind(4) = whatIL(ind(2));
        ind(2) = 0;
        
        %calculating position for y_bar (il,r,s,t)
        pos = resolvePos(10, ind);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',ind(:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;

        %V(il, r,s,t)        
        pos = resolvePos(2, ind);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',ind(:)) ));
            return;
        end
        %temp(pos) = -1;
        spm(num,:)=[iter,pos,-1]; num= num+1;

        %vec_(iter,:) = temp;
        bspm(iter,:)=[iter,1, 0];

    end
    
    %vec2 = zeros(length(indArr),1); %RHS is zero
    
    spm( ~any(spm,2), : ) = [];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm,bspm] = eq18()
%constructign equation 18
%u (l,r,s,t) - d(l,r,t) + y (il,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global d_lrt; 
    %global vecLen; 
    global whatIL; 
    global whatKL;

    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,1]); 

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN = 1 + 1;
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    bnum = 1;

    global Naneq;
    Naneq = rN;
    

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);

        %pos for u
        pos = resolvePos(1, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indArr(:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;

        %fill it with il wrt l
        ind = indArr(iter,:);
        ind(4) = whatIL(ind(2));
        l = ind(2);
        ind(2) = 0;
        
        %calculating position for y (il,r,s,t)
        pos = resolvePos(9, ind);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',ind(:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;

        %d(l, r,t)        
        k = whatKL(l);
        temp_v =  d_lrt{2,k};
        tRange = temp_v(1,:);
        %display(ind(6));
        col = 1;
        for ct = tRange
            if(ct == ind(6))
               break;
            end
            col = col+1;
        end
        
        %vec2(iter) = temp_v(indArr(iter,3)+1, col); %d_lrt, r+1 and mathcing t. see inputscene for clue.
        
        %vec_(iter,:) = temp;

        bspm(bnum,:)=[iter,1,temp_v(indArr(iter,3)+1, col)]; bnum= bnum+1;

    end
    
    spm( ~any(spm,2), : ) = [];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm, bspm] = eq19()
%constructign equation 19
%u (l,r,s,t) - 1/t *d(l,r,t) + sum y (il,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global d_lrt; 
    %global vecLen; 
    global whatIL; 
    global whatKL;
    global maxV_;
    
    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,1]); 

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN = 1 + maxV_(6);
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    bnum = 1;

    global Naneq;
    Naneq = rN;
    

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);

        %pos for u
        pos = resolvePos(1, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;
        
        %y sum (il,r,s,t)
        %fill it with il wrt l
        ind = indArr(iter,:);
        ind(4) = whatIL(ind(2));
        l = ind(2);
        ind(2) = 0;
        
        %calculating position for y (il,r,s,t)
        %fix sum over tau.
        %tau = 0 .. min(t, Td)
        t = indArr(iter,6);
        for tau = 1:t % tau iterates
            arr = ind;
            arr(6) = tau; %iterating over (tau)
            pos = resolvePos(9, arr);
            if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                return;
            end
            %temp(pos) = 1/ind(6); %1/t
            spm(num,:)=[iter,pos,1/ind(6)]; num= num+1;

        end
        
        %d(l, r,t)
        k = whatKL(l);
        temp_v =  d_lrt{2,k}; 
        tRange = temp_v(1,:);
        col = 1;
        for ct = tRange
            if(ct == ind(6))
               break;
            end
            col = col+1;
        end
        
        %vec2(iter) = temp_v(indArr(iter,3)+1, col) / ind(6); %d_lrt* 1/t, r+1 and mathcing t. see inputscene for clue.
        bspm(bnum,:)=[iter,1,temp_v(indArr(iter,3)+1, col) / ind(6)]; bnum= bnum+1;

        %update entire vector
        %vec_(iter,:) = temp;

    end
    
    spm( ~any(spm,2), : ) = [];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm, bspm] = eq20()
%constructign equation 20
%u (l,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


     %global vecLen; global whatIL;


    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,1]); 

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN = 1 ;
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    %bnum = 1;

    global Naneq;
    Naneq = rN;
    

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);

        %pos for u (l,r,s,t)
        pos = resolvePos(1, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indArr(iter,:)) ));
            return;
        end
        %temp(pos) = 1;
        spm(num,:)=[iter,pos,1]; num= num+1;
        
        bspm(iter,:)=[iter,1, 0];


        %vec_(iter,:) = temp;

    end
    
    %vec2 = zeros(length(vec_),1); %RHS is zero
    spm( ~any(spm,2), : ) = [];

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm,bspm] = eq21()
%constructign equation 21
% sum w (i,r)  = Cr

%generate indices combinations for r
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global maxV_; 
    %global vecLen; 
    global C_r;

    indArr = generateIndices([0,0,maxV_(3),0,0,0], 1);

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN =  maxV_(4);
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    bnum = 1;


    global Naneq;
    Naneq = rN;
    
%disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);

        %pos for sum w (i,r,) over i
        arr = indArr(iter,:);
        for i = 1:maxV_(4);
            arr(4) = i;
            pos = resolvePos(3, arr);
            if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                return;
            end
            %temp(pos) = 1;
            spm(num,:)=[iter,pos,1]; num= num+1;

        end
        
        %vec2(iter) = C_r(arr(3));
        bspm(iter,:)=[iter,1,C_r(arr(3))]; bnum= bnum+1;

        %vec_(iter,:) = temp;

    end
        spm( ~any(spm,2), : ) = [];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm,bspm] = eq22()
%constructign equation 22
% double sum x (J, I, r, s,t)  <= w(i,r)

%generate indices combinations for s,r,i
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global maxV_; 
    %global vecLen; 
    global Ts; 
    global epsilon_i;

    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,0], 1);

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN = 1 + maxV_(5)*maxV_(6);
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    %bnum = 1;

    global Naneq;
    Naneq = rN;
    
    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);

        %pos for doublesum x (s,r,j, i,t) over Ts, j
        arr = indArr(iter,:);
        i = arr(4);
        for t = 1:Ts(maxV_(1));            
            for j = epsilon_i(i,:)
                arr(4) = j; %because it has j in place of i.
                arr(5) = i; %because it has i in place of j.
                arr(6) = t;
                pos = resolvePos(6, arr);
                if (pos == -1)
                    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                    return;
                end
                %temp(pos) = 1;
                spm(num,:)=[iter,pos,1]; num= num+1;

            end
        end
        
        %for w(i,r)
        arr(:) = 0;
        arr(3) = indArr(iter,3);
        arr(4) = i;
        
        pos = resolvePos(3, arr);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        %temp(pos) = -1;
        spm(num,:)=[iter,pos,-1]; num= num+1;
       
        bspm(iter,:)=[iter,1, 0];

        %vec_(iter,:) = temp;

    end
    
    spm( ~any(spm,2), : ) = [];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spm,bspm] = eq23()
%constructign equation 23
% double sum x_bar (J, I, r, s,t) + alpha_r * x (J, I, r, s,t) <= w_bar(i,r)

%generate indices combinations for s,r,i
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global maxV_; 
    %global vecLen; 
    global Ts; 
    global epsilon_i;
    global alphaR;

    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,0], 1);

    %making the coeff vector
    [rN , ~] = size(indArr);
    %vec_ = zeros(length(indArr),vecLen);
    %vec2 = zeros(length(indArr),1); %for RHS
    %total coeff on one row
    cN = 1 + maxV_(5)*maxV_(6) + maxV_(5)*maxV_(6);
    spm = zeros(rN*cN,3);
    bspm = zeros(rN,3);
    num = 1;
    %bnum = 1;

    global Naneq;
    Naneq = rN;
    

    %disp(indArr);
    for iter = 1:rN
        %temp = zeros(1,vecLen);

        %pos for doublesum x (s,r,j, i,t) over Ts, j
        arr = indArr(iter,:);
        i = arr(4);
        for t = 1:Ts(maxV_(1));            
            for j = epsilon_i(i,:)
                arr(4) = j; %because it has j in place of i.
                arr(5) = i; %because it has i in place of j.
                arr(6) = t;
                pos = resolvePos(7, arr);
                if (pos == -1)
                    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                    return;
                end
                %temp(pos) = 1;
                spm(num,:)=[iter,pos,1]; num= num+1;

                %for x(jirst)
                pos = resolvePos(6, arr);
                if (pos == -1)
                    display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                    return;
                end
                %temp(pos) = alphaR(arr(3));
                spm(num,:)=[iter,pos,alphaR(arr(3))]; num= num+1;

            end
        end
        
        %for w(i,r)
        arr(:) = 0;
        arr(3) = indArr(iter,3);
        arr(4) = i;
        
        pos = resolvePos(4, arr);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        %temp(pos) = -1;
        spm(num,:)=[iter,pos,-1]; num= num+1;
        bspm(iter,:)=[iter,1,0];

        %vec_(iter,:) = temp;

    end
        spm( ~any(spm,2), : ) = [];

end


%%%%%%
%End of all equaitons
%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     To Be minimised
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vec] = eqnF()
global vecLen;
vec = zeros(1,vecLen);
global S; 
global L_s;
global phi_r;
global theta;
global mu;
global N;
global Ts;
global epsilon_i;
global R;
global R_k;
global T_k;
global whatKL;
global p_ij;
global q_s;


%hardcoding the function to minimise
for s = 1:S
    for l = L_s{s}
        kk = whatKL(l);
        Rrange = R_k(kk,:);
        Trange = T_k{kk}; %T_k is cell

        for r = Rrange
            for t = Trange
                %pos for u (l,r,s,t)
                pos = resolvePos(1, [s,l,r,0,0,t]);
                asd = q_s(s)*phi_r(r)*(Ts(s)-t+1);
                vec(1,pos) = asd;
            end
        end
    end
end

%2nd term
for i = 1:N
    for r = 1:R
        pos = resolvePos(4, [0,0,r,i,0,0]);
        vec(1,pos) = theta;
    end
end

for s = 1:S
    for t = 1:Ts(s)
        for i = 1:N
            for j = epsilon_i(i,:)
                %pos for x0 (s i j t)
                pos = resolvePos(5, [s,0,0,i,j,t]);
                vec(1,pos) = mu*p_ij(i,j);
                
                %pos for x (s r i j t)
                for r = 1:R
                    pos = resolvePos(6, [s,0,r,i,j,t]);
                    vec(1,pos) = mu*p_ij(i,j);
                end  
                
            end
        end
    end
end

%vec(:) = vec(:)*-1;
end

