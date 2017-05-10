%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Debmalya Sinha. debmalya.01[att]gmail.com
%Copyleft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, aeq, beq, aineq, bineq] = construct()
%Constructs coefficient vectors for a group of linear functions.
%Var Order: i, j, l, r, s, t
%Decision Var Order: u, v, w, wbar, x0, x, xbar, y0, y, ybar, z

%inputScene();

%initialising all global vars
init();

%vec_ = rangeVarCoeff(1, [0,0,3,3,3,3], 1);
%a= generateIndices([2,0,3,5,2,3,4]);
f = eqnF();
%constructing equations one by one:
[a{1}, b{1}] = eq1();
[a{2}, b{2}] = eq2();
[a{3}, b{3}] = eq3();
[a{4}, b{4}] = eq11();
[a{5}, b{5}] = eq12();
[a{6}, b{6}] = eq13();
[a{7}, b{7}] = eq14();
[a{8}, b{8}] = eq15();
[a{9}, b{9}] = eq16();
[a{10}, b{10}] = eq17();
[a{11}, b{11}] = eq18();
[a{12}, b{12}] = eq19();
[a{13}, b{13}] = eq20();
[a{14}, b{10}] = eq21();

aeq = a{1}; beq = b{1};
for i = 2:8
    temp = [aeq;a{i}];
    aeq = temp; 
    temp = [beq;b{i}];
    beq = temp; 
end

aineq = a{9}; bineq = b{9};
for i = 10:13
    temp = [aineq;a{i}];
    aineq = temp; 
    temp = [bineq;b{i}];
    bineq = temp; 
end
%because both are >= and cplex wants <=
aineq = aineq*-1; 
bineq = bineq*-1;

%the last ones. No need to change sign. these are <=
aineq = [aineq;a{14}]; bineq = [bineq;b{14}];
%aineq = [aineq;a{15}]; bineq = [bineq;b{15}];
%aineq = [aineq;a{14}]; bineq = [bineq;b{16}];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec,vec2] = eq1()
%constructing equation 10.a
%y0 (i,s, t=1) = 0

    global maxV_;
    vec = rangeVarCoeff(8,[maxV_(1),0,0,maxV_(4),0,1], 1, 1);

    vec2 = zeros(length(vec),1); %RHS is zero

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec, vec2]=eq2()
%constructign equation 10.b
%y (i,r,s, t=1) = 0

    global maxV_;
    vec = rangeVarCoeff(9,[maxV_(1),0,maxV_(3),maxV_(4),0,1], 1, 1);
    
    vec2 = zeros(length(vec),1); %RHS is zero

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec,vec2] = eq3()
%constructign equation 10.c
%ybar (i,r,s, t=1) = 0

    global maxV_;
    vec = rangeVarCoeff(10,[maxV_(1),0,maxV_(3),maxV_(4),0, 1],1,1);
    
    vec2 = zeros(length(vec),1); %RHS is zero
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_, vec2] = eq11()
%constructign equation 11
%y (i,r,s, t+1) - y (i,r,s, t) - sum(x ijrst) = 0

%generate indices combinations for i,r,s,t = (Ts-1) and then iterate j for sum(x).

    global maxV_; global vecLen; 

    global sr ; global tij; global epsilon_i;
    %tCondition = (1+sr(r)+tij(i,j)+1;
    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)-1], 1); 


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
            for j = epsilon_i(indArr(iter,4)) %calculates j from i conditionally
                    
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
                    temp(pos) = -1;
                
                end

            end
        vec_(iter,:) = temp;

    end
    
    vec2 = zeros(length(vec_),1); %RHS is zero

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_, vec2] = eq12()
%constructign equation 12
%y (i,r,s, t+1) - y (i,r,s, t) - sum(x_0 ijst) - sum (sum (betaR * x (ijrst) ))= 0

%generate indices combinations for i,s,t and then iterate j,r for sum(x).
%put combination (t+1>tij).
    global maxV_; global vecLen; global epsilon_i;
    global tij;

    %having the coefficients for x ijrs(t)
    global betaR;

    indArr = generateIndices([maxV_(1),0,0,maxV_(4),0,maxV_(6)-1], 1);

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

            %do the actual sum over j
            for j = epsilon_i(indArr(iter,4)) %calculates j from i conditionally
                
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
                    temp(pos) = -1;

                    %now the double sum of j and r for x
                    for r = 1:maxV_(3)
                        arr(3) = r; %iterating over r
                        pos = resolvePos(6, arr);
                        if (pos == -1)
                            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                            return;
                        end
                        temp(pos) = -betaR(r);
                    end
                end
            end

        vec_(iter,:) = temp;

    end
    
    vec2 = zeros(length(vec_),1); %RHS is zero

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_, vec2] = eq13()
%constructign equation 13
%yBar (i,r,s, t+1) - yBar (i,r,s, t) - sum(xBar (ijrst) + sum( alphaR * x(ijrst)) ) = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(xbar + alpa * x).
%put combination (t+1>tij).

    global maxV_; global vecLen; global epsilon_i;
    global tij;

    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)-1], 1); %t = 1..Ts-1

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
        %do the actual sum over j
        for j = epsilon_i(indArr(iter,4)) %calculates j from i conditionally
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
                temp(pos) = -1;

                %for x
                pos = resolvePos(6, arr);
                if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                return;
                end
                %indArr(iter,4) is the current value of r
                temp(pos) = -alphaR(indArr(iter,3)); 
            end

        end

        vec_(iter,:) = temp;

    end
    
    vec2 = zeros(length(vec_),1); %RHS is zero

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_, vec2] = eq14()
%constructign equation 14
%z (i,r,s, t)  - sum(x ijrst) = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t>sr+tij), not t+1.
    global maxV_; global vecLen; global epsilon_i;

    global sr ; global tij;
    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)], 1); %t = 1..Ts

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        %calculating position for z irs(t)
        pos = resolvePos(11, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        temp(pos) = 1;
        
        %calculaitng the positions in the sum
        %only then t condition is met!
        
        %do the actual sum
        for j = epsilon_i(indArr(iter,4)) %calculates j from i conditionally
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
                temp(pos) = -1;
            end
        end

        vec_(iter,:) = temp;

    end
    vec2 = zeros(length(vec_),1); %RHS is zero

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_,vec2] = eq15()
%constructign equation 15
%v (i,r,s, t) - sum z (i,r,s, t cond) - sum z (i,r,s, t cond)  = 0

%generate indices combinations for i,r,s,t and then iterate j for sum(x).
%put combination (t+1>sr+tij).
    global maxV_; global vecLen, global Td;


    indArr = generateIndices([maxV_(1),0,maxV_(3),maxV_(4),0,maxV_(6)], 1); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        %calculating position for v irs(t)
        pos = resolvePos(2, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        temp(pos) = 1;
        
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
                temp(pos) = -1;
            else
                temp(pos) = -2;
            end
        end

        vec_(iter,:) = temp;

    end
    
    vec2 = zeros(length(vec_),1); %RHS is zero
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_,vec2] = eq16()
%constructign equation 16
%y_0 (il, t) - d_0 (l,t)  = 0

%generate indices combinations for l,s,t
%although l depends on s, it'll be dealt in conditionalIndices()
%it  will generate the correct number of l for each s! Don't worry!
%
%il is a constant depnding on l so don't need to iterate i.


    global vecLen; global whatIL; global d0_lt;

    indArr = conditionalIndices([1,1,0,0,0,1]); %unlike generateindices, you put 1 on values you want.

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);
    vec2 = zeros(length(indArr),1); %for RHS
    
    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        
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
        temp(pos) = 1;
        
        %lets deal with the d (l,t)
        vec2(iter) = d0_lt(indArr(iter,2)); %possible change when there more t.
        
        vec_(iter,:) = temp;

    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_, vec2] = eq17()
%constructign equation 17
%y_bar (il,r,s,t) - v (il,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global vecLen; global whatIL;
    
    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,2]); % t=2 means it'll iterate over Ts not T_kl

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);
        
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
        temp(pos) = 1;
        
        %V(il, r,s,t)        
        pos = resolvePos(2, ind);
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',ind(:)) ));
            return;
        end
        temp(pos) = -1;

        vec_(iter,:) = temp;

    end
    
    vec2 = zeros(length(vec_),1); %RHS is zero

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_,vec2] = eq18()
%constructign equation 18
%u (l,r,s,t) - d(l,r,t) + y (il,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global d_lrt; global vecLen; global whatIL; global whatKL;

    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,1]); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);
    vec2 = zeros(length(indArr),1);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);

        %pos for u
        pos = resolvePos(1, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indArr(:)) ));
            return;
        end
        temp(pos) = 1;
        
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
        temp(pos) = 1;
        
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
        
        vec2(iter) = temp_v(indArr(iter,3)+1, col); %d_lrt, r+1 and mathcing t. see inputscene for clue.
        
        vec_(iter,:) = temp;

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_,vec2] = eq19()
%constructign equation 19
%u (l,r,s,t) - 1/t *d(l,r,t) + sum y (il,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global d_lrt; global vecLen; global whatIL; global whatKL;

    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,1]); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);
    vec2 = zeros(length(indArr),1);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);

        %pos for u
        pos = resolvePos(1, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
            return;
        end
        temp(pos) = 1;
        
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
            arr = indArr(iter,:);
            arr(6) = tau; %iterating over (tau)
            pos = resolvePos(9, arr);
            if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                return;
            end
            temp(pos) = 1/ind(6); %1/t
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
        
        vec2(iter) = temp_v(indArr(iter,3)+1, col) / ind(6); %d_lrt* 1/t, r+1 and mathcing t. see inputscene for clue.

        %update entire vector
        vec_(iter,:) = temp;

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_,vec2] = eq20()
%constructign equation 20
%u (l,r,s,t)  = 0

%generate indices combinations for l,r,s,t
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


     global vecLen; global whatIL;


    %unlike generateindices, you put 1 on values you want.
    indArr = conditionalIndices([1,1,1,0,0,1]); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);

        %pos for u (l,r,s,t)
        pos = resolvePos(1, indArr(iter,:));
        if (pos == -1)
            display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',indArr(iter,:)) ));
            return;
        end
        temp(pos) = 1;
        

        vec_(iter,:) = temp;

    end
    
    vec2 = zeros(length(vec_),1); %RHS is zero

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_,vec2] = eq21()
%constructign equation 21
% sum w (i,r)  = Cr

%generate indices combinations for r
%The upper bounds like Rkl is dependent and can be dealt easily.
% We also need to iterate l because il depends on l.
%il is a constant depnding on l so don't need to iterate i.


    global maxV_; global vecLen; global C_r;

    indArr = generateIndices([0,0,maxV_(3),0,0,0], 1); 

    %making the coeff vector
    vec_ = zeros(length(indArr),vecLen);
    vec2 = zeros(length(indArr),1);

    %disp(indArr);
    for iter = 1:length(indArr)
        temp = zeros(1,vecLen);

        %pos for sum w (i,r,) over i
        for i = 1:maxV_(4);
            arr = indArr(iter,:);
            arr(4) = i;
            pos = resolvePos(3, arr);
            if (pos == -1)
                display(strcat(sprintf('Error: position not found for dVar:%d and indices:',dVar), sprintf(' %d',arr(:)) ));
                return;
            end
            temp(pos) = 1;
        end
        
        vec2 = C_r(arr(3));

        vec_(iter,:) = temp;

    end
end
















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


end

