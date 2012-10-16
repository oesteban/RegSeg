function [nii1, nii2, stats] = discriminant_image(coeff_file,pointsfile,nii,maskfile, norm,save)

[coeffs, headertext] = xlsread(coeff_file);
selected_vols = headertext(4:(end-2),1);
volumes = zeros(size(selected_vols));

[~,~,ext] = fileparts(pointsfile);

if (ext=='.xls')
    % XLS Read (only windows)
    [points, vol_names] = xlsread(pointsfile);
else
    % CSV Read
    points = csvread(pointsfile, 1, 0);
    fid = fopen( pointsfile, 'r' );
    vol_names = regexp( fgetl(fid), '\,', 'split' );
    fclose(fid);
end

vol_names = vol_names(1,2:end);

for i = 1:numel(selected_vols)
    for j = 1:numel(vol_names)
        if strcmp(selected_vols{i},vol_names{j})
            volumes(i) = j;
        end
    end
end

nii.img = double(nii.img);
F1 = zeros(size(nii.img(:,:,:,1))); 
F2 = F1;
for i = 1:numel(volumes)
    F1 = coeffs(i,1).*nii.img(:,:,:,volumes(i)) + F1;
    F2 = coeffs(i,2).*nii.img(:,:,:,volumes(i)) + F2;
end
F1 = F1 + coeffs(numel(volumes)+1,1);
F2 = F2 + coeffs(numel(volumes)+1,2);

% Images normalization-----------------------------------------------------
maxIntensity = norm;
if norm ~= 0
    min_F1 = min(F1(:));
    F1 = F1 - min_F1;
    max_F1 = max(F1(:));
    F1 = F1/(max_F1);
    F1 = F1.*maxIntensity;
    
    min_F2 = min(F2(:));
    F2 = F2 - min_F2;
    max_F2 = max(F2(:));
    F2 = F2/max_F2;
    F2 = F2.*maxIntensity;
end

% masking
mask = load_nii( maskfile );
brain_mask = double(mask.img);
brain_mask(brain_mask~=0)=1;
F1 = F1.*brain_mask;
F2 = F2.*brain_mask;

% reorientation
F1 = flipdim(flipdim(F1,1),1);
F2 = flipdim(flipdim(F2,1),1);


% saving images
nii1 = nii;
nii2 = nii;
nii1.img = uint16(F1);
nii2.img = uint16(F2);
nii1.img = (F1);
nii2.img = (F2);
nii1.hdr.dime.dim(5) = 1;
nii2.hdr.dime.dim(5) = 1;
nii1.hdr.dime.dim(1) = 3;
nii2.hdr.dime.dim(1) = 3;
%nii1.hdr.dime.datatype = 512;
%nii2.hdr.dime.datatype = 512;

if save == 1
    save_nii(nii1,'F1composed.nii')
    save_nii(nii2,'F2composed.nii')
end
status = view_nii(nii1);
status = view_nii(nii2);

% -------------------------------------------------------------------------
% Statistics computation---------------------------------------------------

volumes = volumes + 1;
indx = points(:,1);
classes = unique(indx);

%stats = cell(numel(classes)+1,5);
%stats(1,:) = {'Class', 'mean(F1)', 'mean(F2)', 'std(F1)', 'std(F2)'};
%stats(2:end,1) = num2cell(classes);

% Points in the new discriminant space

P1 = zeros(size(points(:,1))); 
P2 = P1;
for i = 1:numel(volumes)
    P1 = coeffs(i,1).*points(:,volumes(i)) + P1;
    P2 = coeffs(i,2).*points(:,volumes(i)) + P2;
end
P1 = P1 + coeffs(i+1,1);
P2 = P2 + coeffs(i+1,2);

if norm ~= 0
    % Apply normalization to the new points
    P1 = P1 - min_F1;
    P1 = P1/max_F1;
    P1 = P1.*maxIntensity;

    P2 = P2 - min_F2;
    P2 = P2/max_F2;
    P2 = P2.*maxIntensity;
end

statistics = zeros( numel(classes), 6 );

for k = 1:numel(classes)
    statistics(classes(k),1) = mean(P1(indx == classes(k)));
    statistics(classes(k),2) = mean(P2(indx == classes(k)));    
    covariance = cov(P1(indx == classes(k)), P2(indx == classes(k)) );
    statistics(classes(k),3:end) = reshape(covariance,1,4);
end

%stats(2:end,2:3) = num2cell(means_F12);
%stats(2:end,4:5) = num2cell(std_F12);
stats = num2cell(statistics);



% Saving the statistics data
if save == 1
    if norm ~= 0,
        name = 'Stats_normalized.csv';
    else
        name = 'Stats_in_LDA_space.csv';
    end
    csvwrite(name,stats)
end