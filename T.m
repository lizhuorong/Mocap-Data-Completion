close all, clear all
% combined 
A = dlmread('Maju Undur.txt');
A01 = dlmread('Melingkar.txt');
A02 = dlmread('Maju Undur.txt');
A05 = dlmread('Setempat.txt');%260
A06 = dlmread('Langkah_Maju_Undur.txt');%460
AN = [A01(1:500,:);A02(1:200,:);A05(1:200,:);A05(151:250,:);A06(1:400,:);A06(261:460,:)];%COMBINED %1600

% A = dlmread('Melingkar.txt');
% AN = A(1:300,:);% Maju Undur  %3D % 300 in 439
% A = dlmread('fastsong7_Take_001.txt');% man_cut_7

% AN = A(1:400,:);%melingkar
% AN = [A(251:450,:);A(551:850,:)];% man_cut_7  %2D % 700 in 856

meanAN = repmat(sum(AN,1),length(AN),1)/length(AN);
B = AN - meanAN;  % only use ground truth mean for test, may further try A1's mean...

[mae_PLOS1_T,mae_xyz_T,recovery] = compute_T(A,B); % missing joint
% [mae_PLOS1_T1,mae_xyz_T1] = compute_T1(A,B);


function [mae_PLOS1,mae_xyz] = compute_T1(A,B)
B = B';
M = B*B';
[u,d,v]=svd(M);

mean1 = repmat(sum(A(401:end,:),1),101,1)/101;

for joint=1:size(M)/3
    joint1 = (joint-1)*3+1;
    joint2 = joint*3;

A1 = (A(401:end,:) - mean1)';
GT = A(401:end,:)';
GT = GT(joint1:joint2,:);
if joint1==1
    A10 = A1(joint2+1:end,:);
    alpha = u(joint2+1:end,:)'*inv(u(joint2+1:end,:)*u(joint2+1:end,:)')*A10;
elseif joint2==size(M)/3
    A10 = A1(1:end-3,:);
    alpha = u(1:end-3,:)'*inv(u(1:end-3,:)*u(1:end-3,:)')*A10;
else
    A10 = A1([1:joint1-1,joint2+1:end],:);
    alpha = u([1:joint1-1,joint2+1:end],:)'*inv(u([1:joint1-1,joint2+1:end],:)*u([1:joint1-1,joint2+1:end],:)')*A10;
end

%alpha = u([1:69,73:end],:)'*inv(u([1:69,73:end],:)*u([1:69,73:end],:)')*A10;

A11 = u(joint1:joint2,:)*alpha;
Astar = A11 + mean1(:,joint1:joint2)';
[n,m] = size(Astar);

mae_PLOS1(joint) = sum(sqrt(sum((Astar-GT).*(Astar-GT),1)))/m;
mae_xyz(:,joint) = sum(abs(Astar-GT),2)/m;
end
end


function [mae_PLOS1,mae_xyz,recovery_ele] = compute_T(A,B)
B = B';
M = B*B';
[u,d,v]=svd(M);
tt=19;% num of joint
recovery_ele =repmat(A(301:400,:),tt,1);% maju % only replace the zero-elements

for joint=1:size(M)/3
    joint1 = (joint-1)*3+1;
    joint2 = joint*3;
Bs = [];
B0 = [];
u0 = [];

% recovery_mat = repmat(A(301:400,:),tt,1);% Astar
for i=1:16
    tmp = B(:,(i-1)*100+1:i*100);
    Bs = [Bs,tmp];
    if joint1==1
        tmp = tmp(joint2+1:end,:);
    elseif joint2==size(M)/3
        tmp = tmp(1:end-3,:);
    else
        tmp = tmp([1:joint1-1,joint2+1:end],:);
    end
    B0 = [B0,tmp];
    [ut,dt,vt] = svd(tmp*tmp');
    u0 = [u0,ut];
    
end
        
alpha = [];
alpha0 = [];
for i=1:16
    alpha =  [alpha, u'                    *Bs(:,(i-1)*100+1:i*100)];
    alpha0 = [alpha0,u0(:,(i-1)*54+1:i*54)'*B0(:,(i-1)*100+1:i*100)];% 3D: 54;2D:72
end

T = (alpha*alpha0')*inv(alpha0*alpha0');
% GT mean
mean1 = repmat(sum(A(301:400,:),1),100,1)/100;

A1 = (A(301:400,:) - mean1)';
GT = A(301:400,:)';
GT = GT(joint1:joint2,:);%3x100
if joint1==1
    A10 = A1(joint2+1:end,:);
elseif joint2==size(M)/3
    A10 = A1(1:end-3,:);
else
    A10 = A1([1:joint1-1,joint2+1:end],:);
end

[u10,d10,v10] = svd(A10*A10');

alpha1 = u10'*A10;

A11 = u(joint1:joint2,:)*T*alpha1;
Astar = A11 + mean1(:,joint1:joint2)';%3x100
tmp1 = Astar';%100x3
recovery_ele((joint-1)*100+1:(joint-1)*100+100,joint1:joint2) = tmp1;
[n,m] = size(Astar);
mae_PLOS1(joint) = sum(sqrt(sum((Astar-GT).*(Astar-GT),1)))/m;
mae_xyz(:,joint) = sum(abs(Astar-GT),2)/m;
end
end




