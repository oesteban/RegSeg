function delete_last_point(fig)

VoxIntensity = getappdata(fig,'VoxIntensity');
nii_view = getappdata(fig,'nii_view');
radioButtons = get(nii_view.handles.hButton,'Children');
idxHandle = findobj(radioButtons,'Value',1);
name = get(idxHandle,'Tag');
idxClass = str2double(name(6));
nVol = size(VoxIntensity,2)-1;


if VoxIntensity{idxClass,1} ~= 0 
    
    % VoxIntensity data ---------------------------------------------------
    neighborhood = 19;
       
    for i = 2:(nVol+1)

        actual_points = VoxIntensity{idxClass,i};
        new_points = actual_points(1:(end - neighborhood));
        if isempty(new_points), new_points = []; end
        VoxIntensity{idxClass,i} = new_points;

    end
    VoxIntensity{idxClass,1} = VoxIntensity{idxClass,1} - 1;
    setappdata(fig,'VoxIntensity',VoxIntensity);
    
    % Statistics table ----------------------------------------------------
    statsClass = cell(nVol, 2);
    
    if VoxIntensity{idxClass,1} ~= 0
        for i = 1:nVol
            statsClass{i, 1} = mean(VoxIntensity{idxClass,i + 1});
            statsClass{i, 2} = std(VoxIntensity{idxClass,i + 1});
        end
    end
    x0 = (idxClass - 1)*nVol + 1;
    x1 = nVol*idxClass;
    vox_stats = get(nii_view.handles.vox_stats,'Data');
    vox_stats(x0:x1,3:4) = statsClass;
    set(nii_view.handles.vox_stats,'Data',vox_stats);
    
    % Points table --------------------------------------------------------
    points = get(nii_view.handles.vox_select,'Data');

    k = 1; % Loop to find the next gap to fill
    while ~isempty(points{k,idxClass})&&(k<numel(points(:,idxClass))), 
        k = k+1; 
    end
    if isempty(points{k,idxClass})
        points{k-1,idxClass} = [];    
    else
        points{k,idxClass} = [];
    end
    
    set(nii_view.handles.vox_select,'Data',points);
    
    % Disable the delete button if applies
    if sum(cell2mat(VoxIntensity(:,1))) == 0
        set(nii_view.handles.hDel,'Enable','off');
    end
    
end