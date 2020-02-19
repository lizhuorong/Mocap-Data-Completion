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




