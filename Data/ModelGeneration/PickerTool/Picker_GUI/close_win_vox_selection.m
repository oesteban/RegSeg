function close_win_vox_selection(fig)

if ishandle(fig)
    nii_view = getappdata(fig,'nii_view');
    nii_view.handles = rmfield(nii_view.handles,'vox_select');
    setappdata(fig,'nii_view',nii_view);
    
    set(fig,'units','char')
    Pos = get(fig,'Position');
    PosMain = Pos + [-(52/2)   0   0   0];
    set(fig,'Position',PosMain);

end