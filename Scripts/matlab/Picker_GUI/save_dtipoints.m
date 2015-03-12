function save_dtipoints(fig)

VoxIntensity = getappdata(fig, 'VoxIntensity');
% dat = datestr(clock);
% dat(12) = '_';
% dat(15) = '_';
% dat(18) = '_';
% name = sprintf('DTIpoints_%s',dat);
save DTIpoints VoxIntensity 