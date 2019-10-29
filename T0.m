close all, clear all

% A = dlmread('Melingkar.txt');

% combined Malaysia male dance and test on Maju
A = dlmread('Maju Undur.txt');
A01 = dlmread('Melingkar.txt');
A02 = dlmread('Maju Undur.txt');
A05 = dlmread('Setempat.txt');%260
A06 = dlmread('Langkah_Maju_Undur.txt');%460
AN = [A01(1:500,:);A02(1:200,:);A05(1:200,:);A05(151:250,:);...
    A06(1:400,:);A06(261:460,:)];%COMBINED DATA:1+2+5+6 %1600

% % combined fast_song and test on man_cut_7
% A = dlmread('.\data\mat\fastsong7_Take_001.txt');% man_cut_7
% A01 = dlmread('.\data\mat\fastsong1_Take_001.txt');%630
% A02 = dlmread('.\data\mat\fastsong2_Take_001.txt');%301
% A03 = dlmread('.\data\mat\fastsong3_Take_001.txt');%644
% A04 = dlmread('.\data\mat\fastsong4_Take_001.txt');%396
% A05 = dlmread('.\data\mat\fastsong5_Take_001.txt');%500
% A06 = dlmread('.\data\mat\fastsong6_Take_001.txt');%619
% A07 = dlmread('.\data\mat\fastsong7_Take_001.txt');%856
% A08 = dlmread('.\data\mat\fastsong8_Take_001.txt');%807
% A09 = dlmread('.\data\mat\fastsong9_Take_001.txt');%194
% AN = [A02(1:300,:);A03(1:600,:);A04(81:380,:);A05(1:500,:);...
%     A06(1:600,:);A08(1:800,:);A07(251:450,:);A07(551:850,:)];%COMBINED DATA %3600

meanAN = repmat(sum(AN,1),length(AN),1)/length(AN);
B = AN - meanAN; 
[mae_T0,mae_xyz,rate,recovery] = compute_T0(A,B);
% final

function [mae_PLOS1,mae_xyz,rate,recovery] = compute_T0(A,B)
B = B';
M = B*B';
[u,d,v]=svd(M);
[n,m] = size(B(:,1:100));
% s = 0;
% while s<15
%     mask = dlmread('1_25.data');
%     mask = rand(n,m);
%     mask(find(mask>0.75))=1;
%     mask(find(mask<1))=0;
%     s= min([min(sum(mask,1)),min(sum(mask,2))]);
% end
% rate = sum(sum(mask))/(n*m);

mask_clip = dlmread('0924_interval10.txt');
[mask,gap,num] = gen_mask(mask_clip);
rate = 1-sum(sum(mask))/(n*m);

Bs = [];
B0 = [];
u0 = [];

for i=1:16
    tmp = B(:,(i-1)*100+1:i*100);
    Bs = [Bs,tmp];
    tmp = tmp.*mask;
    B0 = [B0,tmp];
    [ut,dt,vt] = svd(tmp*tmp');
    u0 = [u0,ut];
end
        
alpha = [];
alpha0 = [];
for i=1:16
    alpha =  [alpha, u'                  *Bs(:,(i-1)*100+1:i*100)];
    alpha0 = [alpha0,u0(:,(i-1)*n+1:i*n)'*B0(:,(i-1)*100+1:i*100)];
end

T = (alpha*alpha0')*inv(alpha0*alpha0');

% man_cut_7
% GT MEAN
% GT = A(451:550,:)';
% mean1=repmat(sum(A(451:550,:),1),100,1)/100;
% A1 MEAN
% a=A(451:550,:).*mask';% man_cut_7
% GT = A(451:550,:)';% man_cut_7
a=A(301:400,:).*mask';% 
GT = A(301:400,:)';% 

meanvec= sum(a)./sum(a~=0);
meanvec(isnan(meanvec)) = 0;
mean1 = repmat(meanvec,100,1);
A1 = (GT'-mean1)';
A10 = A1.*mask;
[u10,d10,v10] = svd(A10*A10');
alpha1 = u10'*A10;

A11 = u*T*alpha1;
Astar = A11 + mean1';

temppp=Astar';
% a1=GT.* mask;
% a2=Astar.*(~mask);
% temppp = GT.* mask + Astar.*(~mask);
% temppp = temppp';

% recovery =  [A(1:200,:);temppp;A(301:end,:)];
recovery =  [A(1:300,:);temppp;A(401:end,:)];
% recovery =  [A(1:450,:);temppp;A(551:end,:)];%man_cut_7
% recovery = [];

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
mask = mask_clip(1:19,:)';
% temp = mask;
% mask = mask';
mask = [mask(:),mask(:),mask(:)]';
tmp = [];
n=57;
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

