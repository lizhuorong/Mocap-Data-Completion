close all, clear all

% traindata  = dlmread('outputChaimue.data');% 900
traindata = dlmread('outputman_cut_7.data');% 500
% traindata = dlmread('test.data');% 400
A=[];
for i = 1:25
    A = [A,traindata(:,(i-1)*3+1:(i-1)*3+2)];
end  

% AN = [A(1951:2650,:);A(2751:2950,:)];% chaimue_exp2 %900
AN = [A(251:450,:);A(551:850,:)];% man_cut_7  % 500
% AN = A(1:400,:);% melingkar   % 400
meanAN = repmat(sum(AN,1),length(AN),1)/length(AN);
B = AN - meanAN; 

[mae_PLOS1_joint_F1,mae_PLOS1_frame_F1,mae_xyz_F1,A_GAP_ms_frame] = compute_F1(A,B);% final

% final
function [mae_PLOS1_joint,mae_PLOS1_frame,mae_xyz,A_GAP] = compute_F1(A,B)
% B = [B(1:100,:),B(101:200,:),B(201:300,:),B(301:400,:),B(401:500,:),B(501:600,:),B(601:700,:),B(701:800,:),B(801:900,:)]';
B = [B(1:100,:),B(101:200,:),B(201:300,:),B(301:400,:),B(401:500,:)]';% man_cut_7
% B = [B(1:100,:),B(101:200,:),B(201:300,:),B(301:400,:)]';% melingkar

tt = 50;% maximum missing frames
m = 100;
n =size(A,2);
M = B'*B;
[u,d,v]=svd(M);
u = u(:,1:n);

A_GAP =repmat(A(451:550,:),tt,1);% man_cut_7 
% A_GAP =repmat(A(2651:2750,:),tt,1);% Chaimue
% A_GAP =repmat(A(401:500,:),tt,1);% melingkar

% GT mean
% mean1 = repmat(sum(A(451:550,:),1),m,1)/m;% man_cut_7
% mean1 = repmat(sum(A(2651:2730,:),1),80,1)/80;% Chaimue_test2

% frame mean
% mean2 = repmat(A(2700,:),m,1);%Chaimue, incease missing frames from 2701 to 2750
mean2 = repmat(A(470,:),m,1);% man_cut_7 
% mean2 = repmat(A(420,:),m,1);%

mean1 = mean2;

mae_PLOS1_frame = [];
mae_xyz = [];

% for frame=1:tt/2 % every two frames
for frame=1:tt % every single frame

% A1 = (A(2651:2750,:) - mean1)';%Chaimue
% GT = A(2651:2750,:)';%Chaimue
A1 = (A(451:550,:) - mean1)';% man_cut_7 
GT = A(451:550,:)';% man_cut_7 
% A1 = (A(401:500,:) - mean1)';% melingkar
% GT = A(401:500,:)';% melingkar 

% % every 2 frames from 51
% GT = GT(:,51:50+frame*2);
% A10 = A1(:,[1:50,50+frame*2+1:end]);
% alpha = A10*u([1:50,50+frame*2+1:end],:)*inv(u([1:50,50+frame*2+1:end],:)'*u([1:50,50+frame*2+1:end],:));
% A11 = alpha*u(51:50+frame*2,:)';
% Astar = A11 + mean1(51:50+frame*2,:)';% 50*ms_frame

% every single frame:from21 % man_cut_7 and, Melingkar
GT = GT(:,21:20+frame);
A10 = A1(:,[1:20,20+frame+1:end]);
alpha = A10*u([1:20,20+frame+1:end],:)*inv(u([1:20,20+frame+1:end],:)'*u([1:20,20+frame+1:end],:));
A11 = alpha*u(21:20+frame,:)';
Astar = A11 + mean1(21:20+frame,:)';% 57*ms_frame

% % every single frame:from 51 % chaimue
% GT = GT(:,51:50+frame);
% A10 = A1(:,[1:50,50+frame+1:end]);
% alpha = A10*u([1:50,50+frame+1:end],:)*inv(u([1:50,50+frame+1:end],:)'*u([1:50,50+frame+1:end],:));
% A11 = alpha*u(51:50+frame,:)';
% Astar = A11 + mean1(51:50+frame,:)';% 57*ms_frame


% for rendering
% A_GAP((frame-1)*100+51:(frame-1)*100+50+frame*2,:)=Astar';%filling gaps gradually % every 2 frames
A_GAP((frame-1)*100+21:(frame-1)*100+20+frame,:)=Astar';%filling gaps gradually % every single frames
% A_GAP((frame-1)*100+51:(frame-1)*100+50+frame,:)=Astar';%filling gaps gradually % every single frames


[n,ms_frame] = size(Astar);% n: 50, m: ms frame
tmp = (Astar-GT).*(Astar-GT);

mae_PLOS1 = 0;
for i=1:2:n
%     mae_PLOS1_joint(i) = sum(sqrt(sum(tmp(i:i+2,:),1)))/m;%YU
%     mae_PLOS1 = mae_PLOS1 + sum(sqrt(sum(tmp(i:i+2,:),1)))/m;%YU
    mae_PLOS1_joint(frame,i) = sum(sqrt(sum(tmp(i:i+1,:),1)))/ms_frame;% calculate mae as per every single joint
    mae_PLOS1 = mae_PLOS1 + sum(sqrt(sum(tmp(i:i+1,:),1)))/ms_frame;
    
end
% mae_PLOS1 = mae_PLOS1/(n/3);%YU
% mae_xyz = sum(abs(Astar-GT),2)/m;%YU
tmp1 = mae_PLOS1/(n/2);
mae_PLOS1_frame = [mae_PLOS1_frame,tmp1];% calculate the mae as per joint

tmp2 = sum(abs(Astar-GT),2)/ms_frame;
mae_xyz = [mae_xyz,tmp2];% num of joints x num of times
end
mae_PLOS1_joint = mae_PLOS1_joint(:,[1:2:n]);
mae_xyz = mae_xyz';

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

