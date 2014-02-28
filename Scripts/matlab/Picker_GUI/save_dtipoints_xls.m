function save_dtipoints_xls(fig)

warning off

if iscell(fig)
    VoxIntensity = fig;
elseif ishandle(fig)
    VoxIntensity = getappdata(fig, 'VoxIntensity');
end

[filename pathname] = uiputfile({'*.xls';'*.xlsx'},'Save as');

if ~isempty(VoxIntensity)&&(~isnumeric(filename))
    disp('Exporting data to excel...')
  
    nClasses = size(VoxIntensity,1);
    roiPoints = zeros(nClasses, 1);
    for i = 1:nClasses
        roiPoints(i) = numel(VoxIntensity{i,2});
    end
    
    npoints = sum(roiPoints);
    nVolumes = size(VoxIntensity,2) - 1;
    xlsmatrix = cell(npoints + 1, nVolumes + 1);
    xlsmatrix{1,1} = 'Matter';

    roiTag = [];
    for i = 1:nClasses
        roiTag_i = i*ones(roiPoints(i), 1);
        roiTag = [roiTag; roiTag_i];
    end
    xlsmatrix(2:end,1) = num2cell(roiTag);

    for i = 2:nVolumes+1
        xlsmatrix{1,i} = sprintf('V%d',i-2);
        
        roiPoints = [];
        for j = 1:nClasses
            roiPoints_j = VoxIntensity{j,i};
            roiPoints = [roiPoints; roiPoints_j];
        end
        xlsmatrix(2:end,i) = num2cell(roiPoints);
    end
    xlsmatrix{1,2} = 'B0';
    filename = fullfile(pathname,filename);
    save_data(filename, xlsmatrix);
    disp('done.')
end
