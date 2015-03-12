function save_dtipoints_mat(fig)

VoxIntensity = getappdata(fig, 'VoxIntensity');
% dat = datestr(clock);
% dat(12) = '_';
% dat(15) = '_';
% dat(18) = '_';
% name = sprintf('DTIpoints_%s',dat);

if ~isempty(VoxIntensity)
    disp('Exporting data to matfile...')
    save DTIpoints VoxIntensity 
    disp('Done.')
end

