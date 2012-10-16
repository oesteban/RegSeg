function roi_seg_menu(fig,regions,npix)

nii_view = getappdata(fig,'nii_view');

icon_path = strcat (fileparts(which( 'roi_seg_menu')),'/icon/');



if ~isfield(nii_view.handles,'vox_select')

    set(fig,'units','char')
    Pos = get(fig,'Position');
    PosMain = Pos + [(52/2)   0   0   0];
    set(fig,'Position',PosMain);

    % Main figure ---------------------------------------------------------
    
    PosWin = [11 Pos(2) 52 Pos(4)];
    hSelect = figure('unit','char','position',PosWin,...
         'menubar', 'none','name','Voxel Intensities','numbertitle','off','resize','on',...
         'DeleteFcn',sprintf('close_win_vox_selection(%d)',fig));
     
    % Points to be acquired -----------------------------------------------
    numPoints = npix; 
    nclass = numel(regions);
    points = cell(numPoints,nclass);

    if numPoints > 8 && numPoints <= 12 
        d_off = 1.3; % to increase table size from 8 rows to the maximum
    elseif (numPoints <=8) || (numPoints > 12)
        d_off = 0; % for number of pixels less than 8 -> do not change position
    end
    tablePosition = [0.25,   (38.6 - d_off.*(numPoints - 8)),  51.4,   (12 + d_off.*(numPoints - 8))];

    nii_view.handles.vox_select = uitable('Parent',hSelect,'Units','char',...
                'Position', tablePosition, 'Data', points,... 
                'ColumnName', regions,...
                'Columnwidth',num2cell(repmat(45,1,nclass)));
            
    % Buttongroup ---------------------------------------------------------
    hButton = uibuttongroup('Parent',hSelect,'visible','on','Units','char',...
        'Position',[6.3   50.6   45    1.5],'BorderType','beveledout');
    
    enButton = {'off' 'off' 'off' 'off' 'off'};
    for k = 1:numel(regions)
        enButton{k} = 'on';
    end
    % Create three radio buttons in the button group.
    class1 = uicontrol('Style','Radio','Tag','class1',...
        'Units','normalized','pos',[0.07 0.1 0.1 0.85],'parent',hButton,'Enable',enButton{1});
    class2 = uicontrol('Style','Radio','Tag','class2',...
        'Units','normalized','pos',[0.27 0.1 0.1 0.85],'parent',hButton,'Enable',enButton{2});
    class3 = uicontrol('Style','Radio','Tag','class3',...
        'Units','normalized','pos',[0.47 0.1 0.1 0.85],'parent',hButton,'Enable',enButton{3});
    class4 = uicontrol('Style','Radio','Tag','class4',...
        'Units','normalized','pos',[0.67 0.1 0.1 0.85],'parent',hButton,'Enable',enButton{4});
    class5 = uicontrol('Style','Radio','Tag','class5',...
        'Units','normalized','pos',[0.87 0.1 0.1 0.85],'parent',hButton,'Enable',enButton{5});
    
    nii_view.handles.hButton = hButton;
         
    % "Save as" panel -----------------------------------------------------
    
    savePanelPos = [27,  (34.5 - d_off.*(numPoints - 8)),  23.4,  3.7];

    hSavePanel = uipanel('Parent',hSelect,'title','Save as...','Backgroundcolor',[0.8 0.8 0.8],...
        'Units','Char','Position',savePanelPos);

    icRGB = read_icon_menu( strcat( icon_path, 'matlab.png' ));
    hSave1 = uicontrol('Parent',hSavePanel,'Style','pushbutton','CData',icRGB,...
        'Units','normalized','Position',[0.04 0.08 0.35 0.95],'Callback',sprintf('save_dtipoints_mat(%d)',fig));
    
    icRGB = read_icon_menu(strcat( icon_path, 'excel.png'));
    hSave2 = uicontrol('Parent',hSavePanel,'Style','pushbutton','CData',icRGB,...
        'Units','normalized','Position',[0.6 0.08 0.35 0.95],'Callback',sprintf('save_dtipoints_xls(%d)',fig));
    
    % Delete button -------------------------------------------------------
    
    delButtonPos = [6.76,   (36 - d_off.*(numPoints - 8)),    5.72,    2.11];
    
    icRGB = imread(strcat( icon_path, 'delete.png'));
    nii_view.handles.hDel = uicontrol('Parent',hSelect,'Style','pushbutton',...
        'CData',icRGB,'Units','Char','Enable','off','Callback',sprintf('delete_last_point(%d)',fig),...
        'Position',delButtonPos,'BackgroundColor',[0.8 0.8 0.8]);
    
    % Mean and standard deviation -----------------------------------------
    nVolumes = nii_view.numscan;
    statsTable = cell(nVolumes*numel(regions),4);
    
    for j = 1:nclass
        nameVolumes = cell(nVolumes,2);
        nameVolumes{1,1} = regions{j};
        nameVolumes{1,2} = 'B0';
        
        for i = 2:nVolumes
            nameVolumes{i,1} = regions{j};
            nameVolumes{i,2} = sprintf('V%d',i-1);
        end
        statsTable( (nVolumes*(j-1) + 1) : nVolumes*j, 1:2) = nameVolumes;
    end
    
    % Create statistics table ---------------------------------------------
    statsPos = [1.0400    3.2238   50.0240   (30.1692 - d_off.*(numPoints - 8))];
%     statsPos = statsPos - [0 0 0 10];
    nii_view.handles.vox_stats = uitable('Parent',hSelect,'Units','Char','Position',...
                statsPos, 'Data', statsTable,'Tag','statistics',... 
                'ColumnName', {'Matter','Volume','Mean', 'Std.Dev.'},...
                'Columnformat',{'char' 'char' 'bank' 'bank'},...
                'Columnwidth',{51 44 47 49});     
    
    icRGB = read_icon_menu(strcat(icon_path , 'export.jpg'));
    hExportTable = uicontrol('Parent',hSelect,'Style','pushbutton','CData',icRGB,...
        'Units','char','Position',[41.6    0.5854    5.72    2.1074],...
        'Callback',sprintf('export_excel_table(%d)',hSelect));
    
    uicontrol('Parent',hSelect,'Style','text','String','Export to excel file...',...
        'Backgroundcolor',[0.8 0.8 0.8],'Units','normalized','Position',[0.38 0.021 0.42 0.02]);        
            
    setappdata(fig,'nii_view',nii_view);    
        
    % Variable VoxIntensity -----------------------------------------------
    nImg = nii_view.numscan;
    VoxIntensity = cell(nclass,nImg+1);
    for n = 1:nclass,   VoxIntensity{n,1} = 0;    end
    
    setappdata(fig, 'VoxIntensity', VoxIntensity);    
end


function icRGB = read_icon_menu(name_icon)

    clear icRGB
    [ic,~] = imread(name_icon);
    tamX = 22; tamY = 22;
    icR = imresize(ic(:,:,1), [tamX tamY]);
    icG = imresize(ic(:,:,2), [tamX tamY]);
    icB = imresize(ic(:,:,3), [tamX tamY]);
    icRGB(:,:,1) = icR;
    icRGB(:,:,2) = icG;
    icRGB(:,:,3) = icB;
    
return;