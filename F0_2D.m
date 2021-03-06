close all, clear all
traindata  = dlmread('outputChaimue.data');% 900
% traindata = dlmread('outputman_cut_7.data');% 500
A=[];
for i = 1:25
    A = [A,traindata(:,(i-1)*3+1:(i-1)*3+2)];
end  

AN = [A(1951:2650,:);A(2751:2950,:)];% chaimue_exp2 %900
% AN = [A(251:450,:);A(551:850,:)];% man_cut_7   % 500
% meanA = repmat(sum(A,1),length(A),1)/length(A);
meanAN = repmat(sum(AN,1),length(AN),1)/length(AN);
B = AN - meanAN; 

[mae_F0,mae_xyz,rate_F0,recovery]= compute_F0(A,B);% final

function [mae_PLOS1,mae_xyz,rate,recovery] = compute_F0(A,B) 
% B = [B(1:100,:),B(101:200,:),B(201:300,:),B(301:400,:),B(401:500,:)]';
B = [B(1:100,:),B(101:200,:),B(201:300,:),B(301:400,:),B(401:500,:),...
    B(501:600,:),B(601:700,:),B(701:800,:),B(801:900,:)]';% Chaimue

M = B'*B;
[u,d,v]=svd(M);
% u = u(:,1:57);
% [m,n]=size(A(1:100,:));
m = 100;
n =size(A,2);
u = u(:,1:n);
% s = 0;
% while s<15
%     mask = dlmread('1_25.data'); 
%     mask = rand(n,m);
%     mask(find(mask>0.75))=1;
%     mask(find(mask<1))=0;
%     s= min([min(sum(mask,1)),min(sum(mask,2))]);
% end
% rate = sum(sum(mask))/(n*m);
mask_clip = dlmread('0924_interval25.txt');
[mask,gap,num] = gen_mask(mask_clip);
rate = 1-sum(sum(mask))/(n*m);
B0 = [];
u0 = [];
for i=1:9
    tmp = B((i-1)*n+1:i*n,:);
    tmp = tmp.*mask;
    B0 = [B0;tmp];
    [ut,dt,vt] = svd(tmp'*tmp);
    u0 = [u0;ut(:,1:n)];
end
        
alpha = [];
alpha0 = [];
for i=1:9
    alpha = [alpha;B((i-1)*n+1:i*n,:)*u];
    alpha0 = [alpha0;B0((i-1)*n+1:i*n,:)*u0((i-1)*m+1:i*m,:)];
end

F = inv(alpha0'*alpha0)*(alpha0'*alpha);

% A1 MEAN
a=A(2651:2750,:).*mask';% Chaimue
% a=A(451:550,:).*mask';% man_cut_7
meanvec= sum(a)./sum(a~=0);
meanvec(isnan(meanvec)) = 0;
mean1 = repmat(meanvec,100,1);


GT = A(2651:2750,:)';% Chaimue
% GT = A(451:550,:)';% man_cut_7

A1 = (GT'-mean1)';
A10 = A1.*mask;

[u10,d10,v10] = svd(A10'*A10);
u10 = u10(:,1:n);
alpha1 = A10*u10;
A11 = alpha1*F*u';
Astar = A11 + mean1';

temppp=Astar';
% a1=GT.* mask;
% a2=Astar.*(~mask);
% temppp = GT.* mask + Astar.*(~mask);
% temppp = temppp';


recovery =  [A(1:2650,:);temppp;A(2751:end,:)];% Chaimue
% recovery =  [A(1:450,:);temppp;A(551:end,:)];%man_cut_7
% recovery = [];

% 0928
% [n,m] = size(Astar);
error = ((Astar-GT).*(Astar-GT)).*(~mask);
mae_PLOS1 = 0;
[n,m1] = size(gap);
for i=1:n
    tmp = gap(i,:);
    mae_PLOS1 = mae_PLOS1+sum(sqrt(sum(error(tmp(1):tmp(1)+2,tmp(2):tmp(3)),1)))/(tmp(3)-tmp(2)+1);
end
mae_PLOS1 = mae_PLOS1/n;
mae_xyz = sum(abs(Astar-GT).*(~mask),2)./sum(~mask,2);
end

function [mask,gap,num] = gen_mask(mask_clip)
% s = 0;
% n=n/3;
% while s<15
%     mask = rand(n,m);
%     mask(find(mask>threshold))=1;
%     mask(find(mask<1))=0;
%     s= min(sum(mask,2));
% end
temp = mask_clip;
mask = mask_clip';
% temp = mask;
% mask = mask';
mask = [mask(:),mask(:)]';
tmp = [];
n=75;
m=100;
n1=n/3;
for i=1:n1
    tmp = [tmp;mask(:,(i-1)*m+1:i*m)]; 
end
mask = tmp;
num = sum(sum(~temp));

gap = []; j1=0; j2=0;
for i=1:n1
    for j=1:m
        if temp(i,j)==0
            if j1==0
                j1=j;
            else
                j2=j;
            end
        else
            if j1>0
                if j2==0
                    j2=j1;
                end
                gap = [gap;[i,j1,j2]];
                j1=0; j2=0;
            end
        end
    end
    if j1>0
        j2 = m;
        gap = [gap;[i,j1,j2]];
        j1=0; j2=0;
    end
end
end
