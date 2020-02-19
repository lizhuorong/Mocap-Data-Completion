function [mae_val]=mae(Astar,GT,mask) 
temp=mask';% mask framesx15
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
mae_val = mae_val/n;
% mae_xyz = sum(abs((Astar-GT)).*(~mask),2)./sum(~mask,2);
end