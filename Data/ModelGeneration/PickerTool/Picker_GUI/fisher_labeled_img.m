function [nii1, nii2, nii3] = fisher_labeled_img(nii, excelfile, pointsfile, mask)

[coeffs, headertext] = xlsread(excelfile);

names = headertext(2:(end-1),1);

[points, vol_names] = xlsread(pointsfile);
vol_names = vol_names(2:end);

volumes = zeros(size(names));
for i = 1:numel(names)
    for j = 1:numel(vol_names)
        if strcmp(names{i},vol_names{j})
            volumes(i) = j;
        end
    end
end

imgs = double(nii.img);

F1 = zeros(size(imgs(:,:,:,1)));
F2 = F1;
F3 = F1;

for m = 1:numel(volumes)

    F1 = imgs(:,:,:,volumes(m)).*coeffs(m,1) + F1;
    F2 = imgs(:,:,:,volumes(m)).*coeffs(m,2) + F2;
    F3 = imgs(:,:,:,volumes(m)).*coeffs(m,3) + F3;

end
    
F1 = F1 + coeffs(m+1,1);
F2 = F2 + coeffs(m+1,2);
F3 = F3 + coeffs(m+1,3);

roi1 = (F1>F2).*(F1>F3);
roi2 = (F2>F1).*(F2>F3);
roi3 = (F3>F1).*(F3>F2);

nii1 = nii;
nii2 = nii;
nii3 = nii;

nii1.img = roi1;
nii2.img = roi2;
nii3.img = roi3;

nii1.hdr.dime.dim(1) = 3;
nii2.hdr.dime.dim(1) = 3;
nii3.hdr.dime.dim(1) = 3;
nii1.hdr.dime.dim(5) = 1;
nii2.hdr.dime.dim(5) = 1;
nii3.hdr.dime.dim(5) = 1;

% Mask and labeled image
%mask = load_nii('DTI1_mask.nii');
% mask = load_nii('DSI_mask.nii');

mask = double(mask.img);

out = mask.*(roi1)*2 + mask.*(roi2)*6 + mask.*(roi3)*9;

out = flipdim(flipdim(out,1),2);

nii_out = nii;
nii_out.img = out;
nii_out.hdr.dime.dim(1) = 3;
nii_out.hdr.dime.dim(5) = 1;
status = view_nii(nii_out);

% out = (roi1)*3 + (roi2)*6.75 + (roi3)*9;
% nii_out = nii;
% nii_out.img = out;
% nii_out.hdr.dime.dim(1) = 3;
% nii_out.hdr.dime.dim(5) = 1;
% status = view_nii(nii_out);


% Delete Background from roi2 ---------------------------------------------

% im1 = zeros(size(roi2));
% im2 = zeros(size(roi2));
% im3 = zeros(size(roi2));
% 
% for i = 1:size(roi2,1)
%     im1(i,:,:) = imclearborder(squeeze(roi2(i,:,:)),8);
% end
% 
% for i = 1:size(roi2,2)
%     im2(:,i,:) = imclearborder(squeeze(roi2(:,i,:)),8);
% end
% 
% for i = 1:size(roi2,3)
%     im3(:,:,i) = imclearborder(roi2(:,:,i),8);
% end
% 
% out = im1 | im2 | im3;
% nii_out = nii;
% nii_out.img = im1;
% nii_out.hdr.dime.dim(1) = 3;
% nii_out.hdr.dime.dim(5) = 1;
% status = view_nii(nii_out);
% 
% 
% nii_out = nii;
% nii_out.img = imclearborder(roi2, 6);
% nii_out.hdr.dime.dim(1) = 3;
% nii_out.hdr.dime.dim(5) = 1;
% status = view_nii(nii_out);
