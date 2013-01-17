function VoxIntensity = add_new_voxel(nii_view, VoxIntensity, img, mask, imgvalue)


handles = get(nii_view.handles.hButton,'Children');
roi = findobj(handles,'Value',1);
name = get(roi,'Tag');
idxClass = str2double(name(6));

points = get(nii_view.handles.vox_select,'Data');
nameRoi = get(nii_view.handles.vox_select,'ColumnName');
npoints = size(points,1);
vox_stats = get(nii_view.handles.vox_stats,'Data');
nVol = size(VoxIntensity,2)-1;

if (VoxIntensity{idxClass, 1}<npoints) 
    
    set(nii_view.handles.hDel,'Enable','on')
    
    k = 1; % Loop to find the next gap to fill
    while ~isempty(points{k,idxClass})&&(k<numel(points(:,idxClass))), 
        k = k+1; 
    end
    points{k,idxClass} = imgvalue;
    
    statsClass = cell(nVol, 2);
    
    for iScan = 1:nVol
        img_neighbour = mask.*double(squeeze(img(:,:,:, iScan)));
        imgpoints = img_neighbour(mask~=0);
        VoxIntensity{idxClass,iScan+1} = [VoxIntensity{idxClass,iScan+1}; imgpoints];
        
        statsClass{iScan, 1} = mean(VoxIntensity{idxClass,iScan+1});
        statsClass{iScan, 2} = std(VoxIntensity{idxClass,iScan+1});
    end
    VoxIntensity{idxClass,1} = VoxIntensity{idxClass,1}+1;
    
    x0 = (idxClass - 1)*nVol + 1;
    x1 = nVol*idxClass;
    vox_stats(x0:x1,3:4) = statsClass;
    set(nii_view.handles.vox_select,'Data',points);
    set(nii_view.handles.vox_stats,'Data',vox_stats);
   
    fprintf('Region = %s ...... [%d %d %d] ...... Intensity = %d  \n',...
              nameRoi{idxClass}, nii_view.slices.sag,nii_view.slices.cor,nii_view.slices.axi,...
              imgvalue), 
          
end
              