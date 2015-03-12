function nii_out = add_imgs( fnames )
nii_out = '';
for cell = fnames
    fname = char(cell);
    [path, name ,ext] = fileparts(fname);
    
    if strcmp(ext,'.nii') || strcmp(ext,'.gz')
        if ~isstruct( nii_out )
            im1 = load_nii( fname );
            nii_out = im1;
        else
            nii2add = load_nii(fname);
            img = nii2add.img;
            prev_img = nii_out.img;

            if numel(img) ~= numel(im1.img)
                err = MException('SizeChk:Missmatch', ...
                'Adding images with different size is not supported');
                throw (err)
            end
            
            [a b c d] = size(prev_img);
            new_img = zeros(a,b,c,d + 1);
            new_img(:,:,:,1:end-1) = prev_img;
            new_img(:,:,:,end) = img;

            nii_out.img = new_img;
            nii_out.hdr.dime.dim(5) = nii_out.hdr.dime.dim(5) + 1;
            nii_out.original.hdr.dime.dim(5) = nii_out.original.hdr.dime.dim(5) + 1;
        end

    end
end