function nii_out = add_img_series(file_serie, imgs2add)

% The first parameter "file_serie" is the DTI series images and the second
% parameter "imgs2add" is the image to add at the end of the DTI series.
% Both parameters are strings and the files must be in the path. The
% parameter "imgs2add" may also be the path where are the images to add, in
% this case is possible to add more than only one image. The function will
% add all the image volumes contained in that folder. (Be careful, the
% images must be 3D, not 4D!)


[path, name ,ext] = fileparts(imgs2add); 

if strcmp(ext,'.nii') || strcmp(ext,'.gz')
    
    if ~isstruct(file_serie)
        nii_in = load_nii(file_serie);
    else
        nii_in = file_serie;
    end
    nii2add = load_nii(imgs2add);

    img1 = nii_in.img;
    img2 = nii2add.img;
    
    if numel(img2) ~= numel(img1)
        err = MException('SizeChk:Missmatch', ...
        'Adding images with different size is not supported');
        throw (err)
    end

    [a b c d] = size(img1);
    new_img = zeros(a,b,c,d + 1);

    new_img(:,:,:,1:end-1) = img1;
    new_img(:,:,:,end) = img2;

    nii_out = nii_in;
    nii_out.img = new_img;
    nii_out.hdr.dime.dim(5) = nii_out.hdr.dime.dim(5) + 1;
    nii_out.original.hdr.dime.dim(5) = nii_out.original.hdr.dime.dim(5) + 1;

elseif isempty(ext)
    nii_out = load_nii(file_serie);
    dirData = dir(imgs2add);      %# Get the data for the current directory
    dirIndex = [dirData.isdir];  %# Find the index for directories
    fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
    if ~isempty(fileList)
        fileList = cellfun(@(x) fullfile(imgs2add,x),fileList,'UniformOutput',false);
    end
    for k = 1 : numel(fileList)
        fileName = fileList(k);
        if ( ~strcmp(fileName, '.' ) && ~strcmp(fileName, '..' ) )
            nii_out = add_img_series(nii_out,  char(fileName) );
            fprintf('File  "%s"  added.\n', fileName )
        end
    end
    
end