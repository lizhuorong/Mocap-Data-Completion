close all, clear all,
addpath('./data/aniage');
addpath('./data/Malaysia_male');
addpath('./mask/V2');
addpath('./mat');

figure_tra=0;
 
missing_mask_file = '1_100f_100g.txt';

datafile_list={
    'test_Les.txt',...
    'test_Ragam Tiga Ting Ting.txt',...
    'test_Ragam_Dua.txt',...
    'test_Ragam_Satu_Asas.txt',...
    'test_Angkat kaki ke hadapan dan undur ke belakang_R.txt',...
    'test_Angkat Kaki Ke Sisi Kanan Dan Kiri.txt',...
    'test_Angkat kaki ke sisi kiri.txt',...
    'test_Pergerakan ragam  kunang-kunang mabuk.txt',...
    'test_Tarancang_Melompat.txt',...
    'test_Terancang Setempat.txt',...
   'test_Langkah Acah.txt',...
    'test_Langkah_Maju_Undur.txt',...
    'test_Maju Undur.txt',...
    'test_Setempat.txt',...
    'test_Banga_Alid_Depan_Take_001.txt',...
    'test_bunga_lapan_Take_001.txt',...
    'test_bunga_les_Take_001.txt',...
    'test_Taksim.txt',...
    'test_Melingkar_03.txt',...
    'test_Langkah_Maju_Undur_02.txt',...
     };

MAE_T1=[];
avg_mae_T1=0;

MAE_T2=[];
avg_mae_T2=0;

MAE_T3=[];
avg_mae_T3=0;

MAE_plos1=[];
avg_mae_plos1=0;

Rec_T1=[];
GT_list=[];
% read_s_e=-1;
 read_s_e=[1,130];
% data_length=-1;
data_length=[1,130];
test_times=size(datafile_list,2);
for k=1:test_times
%% load testing data
formatSpec1='Test %d \n';
fprintf(formatSpec1,k);
dataFile =datafile_list{k};
original = read_matrix_data(dataFile, read_s_e);
[original.data,repeatID]=checkrepeat(original.data);
% GT_list=[GT_list; original.data(101:200,:)];
[incomplete, mask] = prepare_data(original, missing_mask_file, data_length);
ColumnsWithGaps = find(any(isnan(incomplete.data),1));
% mask=1-dlmread('2d_3col.txt');
show_figure = 0;
if show_figure == 1
    figure(2);imagesc(mask);
    title('Simulated gaps');
    xlabel('Markers in 3D data (interleaved x,y,z)');
%     xlabel('Markers in 2D data (interleaved x,y)');
    ylabel('Frames');
end
mask=1-mask;%framex15
  
%% preprocess
Data_gaps_plos=incomplete.data;

trainingdata=load_trainingdata(101);
Data_gaps_ours=[trainingdata;incomplete.data];
[M_zeros,N_no_gaps,N_zeros,mean_N_no_gaps,stdev_N_no_gaps,mean_x,...
    mean_y,mean_z,weightvector,after_subt_mean] = PredictMissingMarkers_preprocess(Data_gaps_ours);
AN = N_no_gaps;
AN0 = N_zeros;
A1 = M_zeros;

%% reconstruction
Alength=size(N_no_gaps,1);
GapfilledDataSet_T1 = incomplete.data;
GapfilledDataSet_T2 = incomplete.data;
GapfilledDataSet_T3 = incomplete.data;
IndexesWithGaps = find(isnan(incomplete.data));

% PLOS1
% ReconstructedData_PLOS = PLOS1(AN_plos,AN0_plos,A1_plos);
 [ReconstructedData_PLOS] = PredictMissingMarkers(incomplete.data,'Algorithm',2);
 
%T1
[UN,UN0,T1,r] = computingT1(AN,AN0,Alength);  % Ai*UN*T=Ai0*UN0
Astar_T1=A1*UN0*T1*UN';
Astar_T1 = repmat(mean_N_no_gaps,size(Astar_T1,1),1)...
    + Astar_T1.*repmat(stdev_N_no_gaps,size(Astar_T1 ,1),1)./...
    repmat(reshape([1 1 1]'*weightvector,1,[]),size(M_zeros ,1),1);

Reconstruction_T1 = after_subt_mean;
for j = ColumnsWithGaps 
    Reconstruction_T1(:,j) = Astar_T1(:,j);
end
columns=size(Reconstruction_T1,2);
Reconstruction_T1(:,1:3:columns) = ...
    Reconstruction_T1(:,1:3:columns) + ...
    repmat(mean_x,1,columns/3);
Reconstruction_T1(:,2:3:columns) = ...
    Reconstruction_T1(:,2:3:columns) + ...
    repmat(mean_y,1,columns/3);
Reconstruction_T1(:,3:3:columns) = ...
    Reconstruction_T1(:,3:3:columns) + ...
    repmat(mean_z,1,columns/3);
% Reconstruction_T1=Reconstruction_T1(end-199:end,:);
Reconstruction_T1=Reconstruction_T1(end-data_length(2)+1:end,:);
GapfilledDataSet_T1(IndexesWithGaps) = Reconstruction_T1(IndexesWithGaps);     
Astar_T1=GapfilledDataSet_T1;
% error1 = abs(Astar_T1 -ReconstructedData_PLOS);
% Rec_T1=[Rec_T1;Astar_T1(101:200,:)];

%T2
[UN,UN0,T2,r] = computingT2(AN,AN0,Alength);  % Ai*UN*T=Ai0*UN0
Astar_T2=A1*UN0*T2*UN';
Astar_T2 = repmat(mean_N_no_gaps,size(Astar_T2,1),1)...
    + Astar_T2.*repmat(stdev_N_no_gaps,size(Astar_T2 ,1),1)./...
    repmat(reshape([1 1 1]'*weightvector,1,[]),size(M_zeros ,1),1);

Reconstruction_T2 = after_subt_mean;
for j = ColumnsWithGaps 
    Reconstruction_T2(:,j) = Astar_T2(:,j);
end
columns=size(Reconstruction_T2,2);
Reconstruction_T2(:,1:3:columns) = ...
    Reconstruction_T2(:,1:3:columns) + ...
    repmat(mean_x,1,columns/3);
Reconstruction_T2(:,2:3:columns) = ...
    Reconstruction_T2(:,2:3:columns) + ...
    repmat(mean_y,1,columns/3);
Reconstruction_T2(:,3:3:columns) = ...
    Reconstruction_T2(:,3:3:columns) + ...
    repmat(mean_z,1,columns/3);
% Reconstruction_T1=Reconstruction_T1(end-199:end,:);
Reconstruction_T2=Reconstruction_T2(end-data_length(2)+1:end,:);
GapfilledDataSet_T2(IndexesWithGaps) = Reconstruction_T2(IndexesWithGaps);     
Astar_T2=GapfilledDataSet_T2;
% error1 = abs(Astar_T1 -ReconstructedData_PLOS);
% Rec_T1=[Rec_T1;Astar_T1(101:200,:)];

%T3
[UN,UN0,T3,r] = computingT3(AN,AN0,Alength);  % UN=UN0*T
Astar_T3=A1*UN0*T3*UN';
Astar_T3 = repmat(mean_N_no_gaps,size(Astar_T3,1),1)...
    + Astar_T3.*repmat(stdev_N_no_gaps,size(Astar_T3 ,1),1)./...
    repmat(reshape([1 1 1]'*weightvector,1,[]),size(M_zeros ,1),1);

Reconstruction_T3 = after_subt_mean;
for j = ColumnsWithGaps 
    Reconstruction_T3(:,j) = Astar_T3(:,j);
end
Reconstruction_T3(:,1:3:columns) = ...
    Reconstruction_T3(:,1:3:columns) + ...
    repmat(mean_x,1,columns/3);
Reconstruction_T3(:,2:3:columns) = ...
    Reconstruction_T3(:,2:3:columns) + ...
    repmat(mean_y,1,columns/3);
Reconstruction_T3(:,3:3:columns) = ...
    Reconstruction_T3(:,3:3:columns) + ...
    repmat(mean_z,1,columns/3);
% Reconstruction_T3=Reconstruction_T3(end-199:end,:);
Reconstruction_T3=Reconstruction_T3(end-data_length(2)+1:end,:);
GapfilledDataSet_T3(IndexesWithGaps) = Reconstruction_T3(IndexesWithGaps);     
Astar_T3=GapfilledDataSet_T3;

% error2 = abs(Astar_T3-Astar_T1);
% error3 = abs(Astar_T3-ReconstructedData_PLOS);

% 
% [UN,UN0,T2,r] = computingT2(AN,AN0,Alength);  % UN'*Ai=T*U_i'*Ai0
% 
% [u,d,v]=svd(A1,'econ');
% u=u(:,1:r);
% Astar=UN*T1*u'*A1;
%% evaluation
[mae_PLOS1]=calMAE(ReconstructedData_PLOS,original.data,mask);
MAE_plos1=[MAE_plos1,mae_PLOS1];
avg_mae_plos1=avg_mae_plos1+mae_PLOS1/test_times;

[mae_T1]=calMAE(Astar_T1,original.data,mask);
MAE_T1=[MAE_T1,mae_T1];
avg_mae_T1=avg_mae_T1+mae_T1/test_times;

[mae_T2]=calMAE(Astar_T2,original.data,mask);
MAE_T2=[MAE_T2,mae_T2];
avg_mae_T2=avg_mae_T2+mae_T2/test_times;

[mae_T3]=calMAE(Astar_T3,original.data,mask);
MAE_T3=[MAE_T3,mae_T3];
avg_mae_T3=avg_mae_T3+mae_T3/test_times;
% save('comparison_results.mat','MAE_plos1','MAE_T1','MAE_T3','avg_mae_plos1','avg_mae_T1','avg_mae_T3');



end

if figure_tra==1
    figure
    for j=1:20
    subplot(4,5,j)
    x=linspace(1,100);
    y1 = GT_list((j-1)*100+1:j*100,ColumnsWithGaps);
    y2 = Rec_T1((j-1)*100+1:j*100,ColumnsWithGaps);
    p=plot(x,y1,'k--','LineWidth',1);
    hold on
    plot(x,y2,'LineWidth',1)
    % plot(x,y2,x,y1,'k--')
    str=['Test sample ', num2str(j)];
    % str=['Predicted trajectories of test sample ', num2str(k)];
    title(str);
    xlabel('Frames');
    ylabel('Markers coordinate');
    % hold off
    end
    legend(p,'Ground true')
end

function ReconstructedData = PLOS1(N_no_gaps,N_zeros,M_zeros)

    [PC_vectors_no_gaps,sqrtEigvals_no_gaps] = PCA(N_no_gaps);

    [PC_vectors_zeros,sqrtEigvals_zeros] = PCA(N_zeros);

    % Select the number of PV-vectors to include in the analysis
    n_eig = settingrank(sqrtEigvals_no_gaps,sqrtEigvals_zeros);

    PC_vectors_no_gaps = PC_vectors_no_gaps(:,1:n_eig);
    PC_vectors_zeros = PC_vectors_zeros(:,1:n_eig);
    % Calculate Transformation Matrix
    T = PC_vectors_no_gaps'*PC_vectors_zeros;
    
    % Transform Data first into incomplete-, then into full-PC basis system.
    ReconstructedData = M_zeros*PC_vectors_zeros*T*PC_vectors_no_gaps';
end

function r=settingrank(varargin)
MinCumSV=0.99;
num = length(varargin);
r=[];
for i=1:num
    eigv=varargin{i};
    r = [r;find(cumsum(eigv)>=MinCumSV*sum(eigv),1,'first')];
%        find(cumsum(eigv2)>=MinCumSV*sum(eigv2),1,'first')];
end
    r = max(r);
end

function [PC,sqrtEV] = PCA(Data)
    [N,M] = size(Data);
    Y = Data / sqrt(N-1);
    [U,sqrtEV,PC] = svd(Y,'econ');
    sqrtEV = diag(sqrtEV);
end

function [UN,UN0,T,r]=computingT3(AN,AN0,Alength)% UN=UN0*T

AN=AN/sqrt(Alength-1);
AN0=AN0/sqrt(Alength-1);

K = size(AN,1)/Alength;
[u0,d0,v0] = svd(AN0,'econ');
[u,d,v] = svd(AN,'econ');

r = settingrank(diag(d),diag(d0));

UN = v(:,1:r);
UN0 = v0(:,1:r);

T = UN0'*UN;
% T = UN'*UN0;%li
end

function [UN,UN0,T,r] = computingT1(AN,AN0,Alength)% Ai*UN*T=Ai0*UN0
K = size(AN,1)/Alength;
AN = AN/sqrt(Alength-1);
AN0= AN0/sqrt(Alength-1);

[u0,d0,v0] = svd(AN0,'econ');
[u,d,v] = svd(AN,'econ');

tmp = 0;
for i=1:K
    Ai = AN((i-1)*Alength+1:i*Alength,:);
    A0i = AN0((i-1)*Alength+1:i*Alength,:);
    tmp = tmp + v0'*A0i'*A0i*v0;
end
dtmp=svd(tmp);
r = settingrank(diag(d),diag(d0),dtmp);
%r = rank(tmp,0.001);
UN = v(:,1:r);
UN0 = v0(:,1:r);

tmp = 0;
tmpN = 0;
for i=1:K
    Ai = AN((i-1)*Alength+1:i*Alength,:);
    A0i = AN0((i-1)*Alength+1:i*Alength,:);
    tmp = tmp+UN0'*A0i'*A0i*UN0;
    tmpN = tmpN+UN0'*A0i'*Ai*UN;
end
T = tmpN*inv(tmp);
end

function [UN,UN0,T,r] = computingT2(AN,AN0,Alength)
K = size(AN,1)/Alength;
AN = AN/sqrt(Alength-1);
AN0= AN0/sqrt(Alength-1);

[u0,d0,v0] = svd(AN0,'econ');
[u,d,v] = svd(AN,'econ');

tmp = 0;
for i=1:K
    Ai0 = AN0((i-1)*Alength+1:i*Alength,:);
    [ui,di,vi]=svd(Ai0,'econ');
    tmp = tmp + vi'*Ai0'*Ai0*vi;
end
dtmp=svd(tmp);
r = settingrank(diag(d),diag(d0),dtmp);

UN = v(:,1:r);
UN0 = v0(:,1:r);

tmp = 0;
tmpN = 0;
for i=1:K
    Ai = AN((i-1)*Alength+1:i*Alength,:);
    Ai0 = AN0((i-1)*Alength+1:i*Alength,:);
    [ui,di,vi]=svd(Ai0,'econ');
    vi = vi(:,1:r);
    tmp = tmp+vi'*Ai0'*Ai0*vi;
    tmpN = tmpN+vi'*Ai0'*Ai*UN;
end
T = tmpN*inv(tmp);
end

function [mae_val]=calMAE(Astar,GT,mask) 
temp=mask';% mask 100x15
mm=size(mask,1);
n1=size(mask,2);
mask = [mask(:),mask(:),mask(:)]';%3x1500
tmppp = [];
for i=1:n1
    tmppp = [tmppp;mask(:,(i-1)*mm+1:i*mm)]; 
end
mask = tmppp;%45x100
% num = sum(sum(~temp));

gap = []; j1=0; j2=0;
for i=1:n1
    for j=1:mm
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
        j2 = mm;
        gap = [gap;[i,j1,j2]];
        j1=0; j2=0;
    end
end

if size((Astar),1)~=45 
    Astar=Astar';  
end    
if size((GT),1)~=45 
    GT=GT';
end   
error = ((Astar-GT).*(Astar-GT)).*(~mask);%45x100
mae_val = 0;
[n,m1] = size(gap);

for i=1:n
    tmp = gap(i,:);
    mae_val = mae_val+sum(sqrt(sum(error(tmp(1)*3-2:tmp(1)*3,tmp(2):tmp(3)),1)))/(tmp(3)-tmp(2)+1);
end
AE=1;%Absolute error
if AE==1
    mae_val = mae_val/1;
else
    mae_val = mae_val/n;
end
% mae_xyz = sum(abs((Astar-GT)).*(~mask),2)./sum(~mask,2);
end

function [distArray] = distance2marker(MarkerData,colWithGaps)
    [n,m] = size(MarkerData);
    MarkerWithGaps = colWithGaps(3:3:end)/3;
    nMarkerWidthGaps = length(MarkerWithGaps);
    MarkerData = reshape(MarkerData',3,m/3,n);

    distArray = nan(nMarkerWidthGaps,m/3,n);
    for i=1:n
        distArray(:,:,i) =pdist2(MarkerData(:,MarkerWithGaps,i)',...
            MarkerData(:,:,i)','euclidean');
    end
    distArray = nanmean(distArray,3);
end

function [AN]=load_trainingdata(train_idx)
A01 = dlmread('Les.txt');
A02 = dlmread('Ragam Tiga Ting Ting.txt');
A03 = dlmread('Ragam_Dua.txt');
A04 = dlmread('Ragam_Satu_Asas.txt');
A05 = dlmread('Angkat kaki ke hadapan dan undur ke belakang.txt');
A06 = dlmread('Angkat Kaki Ke Sisi Kanan Dan Kiri.txt');
A07 = dlmread('Angkat kaki ke sisi kiri.txt');
A08 = dlmread('Pergerakan ragam  kunang-kunang mabuk.txt');
A09 = dlmread('Pergerakan Setempat.txt');
A10= dlmread('Tarancang_Melompat.txt');
A11= dlmread('Terancang Setempat.txt');
A12= dlmread('Langkah Acah.txt');
A13= dlmread('Langkah_Maju_Undur.txt');
A14= dlmread('Maju Undur.txt');
A15= dlmread('Melingkar.txt');
A16= dlmread('Setempat.txt');
A17= dlmread('Banga_Alid_Depan_Take_001.txt');
A18= dlmread('bunga_lapan_Take_001.txt');
A19= dlmread('bunga_les_Take_001.txt');
A20= dlmread('Taksim.txt');   
training_data=[...
    A01(train_idx:end,:);A02(train_idx:end,:);A03(train_idx:end,:);...
    A04(train_idx:end,:);A05(train_idx:end,:);A06(train_idx:end,:);A07(train_idx:end,:);...
    A08(train_idx:end,:);A09(train_idx:end,:);A10(train_idx:end,:);A11(train_idx:end,:);...
    A12(train_idx:end,:);A13(train_idx:end,:);A14(train_idx:end,:);A15(train_idx:end,:);...
    A16(train_idx:end,:);A17(train_idx:end,:);A18(train_idx:end,:);A19(train_idx:end,:);...
    A20(train_idx:end,:);...
    ];
AN=training_data;
[AN,repeatID] = checkrepeat(AN);

end
function [incomplete, missing_mask] = prepare_data(original, missing_mask_file, data_length)
% missing_mask = readmatrix(missing_mask_file);
missing_mask = importdata(missing_mask_file);
if data_length ~= -1
    missing_mask = missing_mask(data_length(1):data_length(2), :);
end
missing_mask = 1 - missing_mask;
incomplete = original;
colNum=size(missing_mask,2);
full_missing_mask=zeros(size(missing_mask,1),colNum*3); 
for i=1:colNum
    full_missing_mask(:,i*3-2:i*3)=repmat(missing_mask(:,i),1,3);
end
full_missing_mask = logical(full_missing_mask);
incomplete.data(full_missing_mask)=nan;
end
function [B,repeatID]=checkrepeat(B)
[m,Jn]=size(B);
n = Jn/3;
repeatID = [];
for i=1:n-1
    for j=i+1:n
        tmp = sum(sum(abs(B(:,(i-1)*3+1:i*3) - B(:,(j-1)*3+1:j*3))));
        if tmp==0
            repeatID = [repeatID,j];
        end
    end
end

% if ~isempty(repeatID)
%     for i=length(repeatID):-1:1
%         B(:,(repeatID(i)-1)*3+1:repeatID(i)*3) = [];
%     end
% end
if ~isempty(repeatID)
    tmp = [];
    for i=1:length(repeatID)
        tmp = [tmp,(repeatID(i)-1)*3+1:repeatID(i)*3];
    end
    B(:,tmp) = [];
end
end

function [original] = read_matrix_data(data_file, data_length)
original = mcread('format.c3d');
data = importdata(data_file);
original.data = data;
original.filename = data_file;
original.nFrames = size(data, 1);
original.nMarkers = size(data, 2) / 3;
original.markerName = original.markerName(1: original.nMarkers);
original.freq = 30;
if data_length ~= -1
original.data = original.data(data_length(1):data_length(2), :);
original.nFrames = data_length(2)-data_length(1)+1;
original.other.residualerror = original.other.residualerror(data_length(1):data_length(2), :);
end
end

