function varargout = dti_segmentation(varargin)

% if nargin == 0,  help atlas3dgui; return;  end

% execute callback function then return;
if ischar(varargin{1}) && ~isempty(findstr(varargin{1},'Callback')),
  if nargout
    [varargout{1:nargout}] = feval(varargin{:});
  else
    feval(varargin{:});
  end
  return;
end


ANA = {};  PERMUTE_VEC = [];
%     ANA = anaload(Ses,varargin{2}); 
    
namefile = varargin{1};
spm_defaults;
V = spm_vol(namefile);
vol = spm_read_vols(V);
vol = abs(vol);
flip = [2 3];
for N = 1:length(flip),
    vol = flipdim(vol,flip(N));
end
ANA.dat = vol;
ANA.dim = V.dim;
ANA.ds = abs(V.private.mat0([1 6 11]));
if nargin > 2,  PERMUTE_VEC = varargin{3};  end

roiFile = load(varargin{2});
% Roi = ATLAS_ROI;
if isfield(roiFile,'RoiMatrix')
    Roi = roiFile.RoiMatrix;
    RoiNames = Roi(1,:);
elseif isfield(roiFile,'Atlas_mdeftinj')
    Roi = roiFile.Atlas_mdeftinj;
    RoiNames = Roi.roinames;
    Roi.fileName = varargin{2};
end

if isempty(ANA),
  fprintf('\n%s ERROR: no way to get anatomy data.\n',mfilename);
  return;
end

if ndims(ANA.dat) > 3,  ANA.dat = squeeze(ANA.dat);  end

% do permutation, if given
if ~isempty(PERMUTE_VEC),
  ANA.dat = permute(ANA.dat,PERMUTE_VEC);
end

ANA.dat = double(ANA.dat); 

if ~isfield(ANA,'session') || isempty(ANA.session),
  ANA.session = 'unknown';
end
if ~isfield(ANA,'grpname') || isempty(ANA.grpname),
  ANA.grpname = 'unknown';
end

if isempty(PERMUTE_VEC), permuteV = [1 2 3];
else 
    permuteV = PERMUTE_VEC;
end

% GET SCREEN SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oldunits = get(0,'units');
set(0,'units','char');
SZscreen = get(0,'ScreenSize');
set(0,'units',oldunits);
% scrW = SZscreen(3);  scrH = SZscreen(4);
% 
% figW = 232; figH = 51;
% figX = 12.8;  figY = scrH-figH-5.75;

% [figX figY figW figH]
firstScreen = [0.6000    4.6923  254.8000   51.1538];
secondScreen = [16.4000    4.6923  223.8000   51.1538];

% color = [0.99 0.90 0.70];

color = [0.9300    0.8500    0.7100];
% CREATE A MAIN FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hMain = figure(...
    'Name',sprintf('%s: Visualize Atlas over experimental Image',mfilename),...
    'NumberTitle','off', 'toolbar','figure',...
    'Tag','main', 'units','char', 'pos',secondScreen,...
    'HandleVisibility','on', 'Resize','on',...
    'DoubleBuffer','on', 'BackingStore','on', 'Visible','on',...
    'DefaultAxesFontSize',12, 'Color', color,...
    'DefaultAxesFontName', 'Comic Sans MS',...
    'Pointer', 'watch',...
    'DefaultAxesfontweight','bold');


set(hMain,'toolBar','none');
set(hMain,'menubar','none');

H = 2; XSZ = 55; YSZ = 20;
XINI = 4;

%POSITION = [xinicial yinicial ancho alto];


PanelAna = uipanel(...
    'Parent',hMain,...'Title','Anatomical Image','FontWeight','bold','FontSize',11,...
    'BorderType','beveledin ',...
    'BackgroundColor',get(hMain,'Color')-[0.045 0.06 0.06],...
    'Units','char','Position',[XINI-0.5 H-0.5 31 49]);
TitleAna = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI 48.7 30 1.5],...
    'String','Anatomical Image','FontWeight','bold','FontSize',12,...
    'HorizontalAlignment','center',...
    ...'FontName', 'Corbel',...'FontSize',9,... 
    'BackgroundColor',get(PanelAna,'BackgroundColor'));

if exist(namefile,'file')
    FileTxt = uicontrol(...
        'Parent',hMain,'Style','Text',...
        'Units','char','Position',[XINI+1.5 45 27 1.3],...
        'String',sprintf('%s',namefile),'FontSize',9,...
        'HorizontalAlignment','left',...
        ...'FontName', 'Corbel',...'FontSize',9,... 
        'BackgroundColor',color);
else
    SesTxt = uicontrol(...
        'Parent',hMain,'Style','Text',...
        'Units','char','Position',[XINI+1.5 46.5 27 1.3],...
        'String',sprintf('SES: %s',varargin{1}),...
        'HorizontalAlignment','left',...
        'FontName', 'Comic Sans MS','FontSize',9,... 
        'Tag','SesTxt','BackgroundColor',color);
    GrpTxt = uicontrol(...
        'Parent',hMain,'Style','Text',...
        'Units','char','Position',[XINI+1.5 45 27 1.3],...
        'String',sprintf('GRP:  %s',varargin{2}),...
        'HorizontalAlignment','left',...
        'FontName', 'Comic Sans MS','FontSize',9,... 
        'Tag','GrpTxt','BackgroundColor',color);
end
DimTxt = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI+1.5 43.5 27 1.3],...
    'String',sprintf('SIZE: %dx%dx%d',ANA.dim),...
    'HorizontalAlignment','left',...
    'FontName', 'Comic Sans MS','FontSize',9,... 
    'Tag','DimTxt','BackgroundColor',color);
ResTxt = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI+1.5 42 27 1.3],...
    'String',sprintf('RES: %gx%gx%g',ANA.ds),...
    'HorizontalAlignment','left',...
    'FontName', 'Comic Sans MS','FontSize',9,... 
    'Tag','ResTxt','BackgroundColor',color);
PermuteTxt = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI+1.5 40.5 27 1.3],...
    'String',sprintf('PERMUTE: [%g %g %g]',permuteV),...
    'HorizontalAlignment','left',...
    'FontName', 'Comic Sans MS','FontSize',9,... 
    'Tag','ResTxt','BackgroundColor',color);

infopanel = uipanel(...
    'Parent',hMain,...
    'BorderType','etchedin ',...
    'BackgroundColor',color,...
    'Units','char','Position',[XINI+0.5 40 29 8]);

% BREGMA/VOXEL COORDINATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hSel = uibuttongroup('Parent',hMain,'visible','on','Units','char',...
    'Position',[XINI 35.3 30 4],'BackgroundColor',get(PanelAna,'BackgroundColor'),...
    'BorderType','none ');

SelVox = uicontrol('Style','Radio','String','Voxels','tag','SelVox',...
    'Background',get(PanelAna,'BackgroundColor'),'FontWeight','bold','FontSize',8,...
    'Units','char','pos',[0.5 2.7 12 1.2],'parent',hSel,'HandleVisibility','off');
SelBreg = uicontrol('Style','Radio','String','Milimeters','tag','SelBreg',...
    'Background',get(PanelAna,'BackgroundColor'),'FontWeight','bold','FontSize',8,...
    'Units','char','pos',[12.8 2.7 16 1.2],'parent',hSel,'HandleVisibility','off');

set(hSel,'SelectionChangeFcn','atlas3dgui(''Main_Callback'',gcbo,''change-coord'',guidata(gcbo))');


GridCheck = uicontrol(...
    'Parent',hMain,'Style','Checkbox',...
    'Units','char','Position',[XINI+1 35.7 12 1.2],...
    'Tag','GridCheck','Value',0,...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''grid-check'',guidata(gcbo))',...
    'String','Grid','FontWeight','bold',...
    'TooltipString','Grid on/off','BackgroundColor',get(PanelAna,'BackgroundColor'));

% CHECK BOX FOR "show-mri" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ShowMRICheck = uicontrol(...
    'Parent',hMain,'Style','Checkbox',...
    'Units','char','Position',[XINI+13 35.7 14 1.2],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''showmri'',guidata(gcbo))',...
    'Tag','ShowMRICheck','Value',1,...
    'String','Anatomy','FontWeight','bold',...
    'TooltipString','Show MRI planes','BackgroundColor',get(PanelAna,'BackgroundColor'));


% VIEW MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ViewModeCmb = uicontrol(...
    'Parent',hMain,'Style','Popupmenu',...
    'Units','char','Position',[XINI H+9.5 30 1.5],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''view-mode'',guidata(gcbo))',...
    'String',{'Orthogonal','Lightbox-Coronal','Lightbox-Sagittal','Lightbox-Transverse'},...
    'Tag','ViewModeCmb','Value',1,...
    'TooltipString','Select the view mode',...
    'FontWeight','bold');
ViewPageList = uicontrol(...
    'Parent',hMain,'Style','Listbox',...
    'Units','char','Position',[XINI H 30 9],...
    'String',{'page1','page2','page3','page4'},...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''view-page'',guidata(gcbo))',...
    'HorizontalAlignment','left',...
    'FontName','Comic Sans MS','FontSize',10,...
    'Tag','ViewPageList','Background','white','visible','off');

% AXES FOR COLORBAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=14;
ColorbarAxs = axes(...
    'Parent',hMain,'Tag','ColorbarAxs',...
    'units','char','Position',[XINI+1 H XSZ*0.1 YSZ],...
    'FontSize',8,...
    'Box','off','YAxisLocation','right','XTickLabel',{},'XTick',[]);
ColorbarMinEdt = uicontrol(...
    'Parent',hMain,'Style','Edit',...
    'Units','char','Position',[XINI+14.5 H 12 1.5],...
    'Callback','atlas3dgui(''Plot_Callback'',gcbo,[],[])',...
    'String','','Tag','ColorbarMinEdt',...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''update-clim'',guidata(gcbo))',...
    'HorizontalAlignment','center',...
    'TooltipString','set colorbar minimum',...
    'FontWeight','Bold');
ColorbarMaxEdt = uicontrol(...
    'Parent',hMain,'Style','Edit',...
    'Units','char','Position',[XINI+14.5 H+YSZ-1.5 12 1.5],...
    'Callback','atlas3dgui(''Plot_Callback'',gcbo,[],[])',...
    'String','','Tag','ColorbarMaxEdt',...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''update-clim'',guidata(gcbo))',...
    'HorizontalAlignment','center',...
    'TooltipString','set colorbar maximum',...
    'FontWeight','Bold');


% GAMMA SETTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI+14.5, H+YSZ/2+5 15.5 1.25],...
    'String','Gamma: ','FontWeight','bold',...
    'HorizontalAlignment','left',...
    'BackgroundColor',get(PanelAna,'BackgroundColor'));
GammaEdt = uicontrol(...
    'Parent',hMain,'Style','Edit',...
    'Units','char','Position',[XINI+14.5, H+YSZ/2+3.5 10 1.5],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''set-gamma'',guidata(gcbo))',...
    'String','1.8','Tag','GammaEdt',...
    'HorizontalAlignment','center',...
    'TooltipString','set a gamma value for image',...
    'FontWeight','bold');


% CHECK BOX FOR X,Y,Z direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XReverseCheck = uicontrol(...
    'Parent',hMain,'Style','Checkbox',...
    'Units','char','Position',[XINI+14.5 H+YSZ/2 15.5 1.5],...
    'Tag','XReverseCheck','Value',0,...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''dir-reverse'',guidata(gcbo))',...
    'String','X-Reverse','FontWeight','bold',...
    'TooltipString','Xdir reverse','BackgroundColor',get(PanelAna,'BackgroundColor'));
YReverseCheck = uicontrol(...
    'Parent',hMain,'Style','Checkbox',...
    'Units','char','Position',[XINI+14.5 H+YSZ/2-2.5 15.5 1.5],...
    'Tag','YReverseCheck','Value',0,...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''dir-reverse'',guidata(gcbo))',...
    'String','Y-Reverse','FontWeight','bold',...
    'TooltipString','Ydir reverse','BackgroundColor',get(PanelAna,'BackgroundColor'));
ZReverseCheck = uicontrol(...
    'Parent',hMain,'Style','Checkbox',...
    'Units','char','Position',[XINI+14.5 H+YSZ/2-2.5*2 15.5 1.5],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''dir-reverse'',guidata(gcbo))',...
    'Tag','ZReverseCheck','Value',0,...
    'String','Z-Reverse','FontWeight','bold',...
    'TooltipString','Zdir reverse','BackgroundColor',get(PanelAna,'BackgroundColor'));


% CHECK BOX FOR "cross-hair" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CrosshairCheck = uicontrol(...
    'Parent',hMain,'Style','Checkbox',...
    'Units','char','Position',[XINI+14.5 H+YSZ/2-7.5 15.5 1.5],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''crosshair'',guidata(gcbo))',...
    'Tag','CrosshairCheck','Value',1,...
    'String','Crosshair','FontWeight','bold',...
    'TooltipString','Show a crosshair','BackgroundColor',get(PanelAna,'BackgroundColor'));


% AXES for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AXES FOR LIGHT BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=2;

VisualPanel = uipanel(...
    'Parent',hMain,...
    'BorderType','etchedin ',...
    'BackgroundColor','none',...
    'Units','char','Position',[XINI+30+6-4 H-1.5 XSZ*2+12+6 YSZ*2+6.5+3.5]);

LightiboxAxs = axes(...
    'Parent',hMain,'Tag','LightboxAxs',...
    'Units','char','Position',[XINI+30+6 H XSZ*2+12 YSZ*2+6.5],...
    'Box','off','color','black','Visible','off');


% AXES FOR ORTHOGONAL VIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 XSZ = 55; 
YSZ = 20; 
XINI = XINI+30+6;
offs = 4;
TriplotAxs = axes(...
    'Parent',hMain,'Tag','TriplotAxs',...
    'Units','char','Position',[XINI+10+XSZ-offs H XSZ+offs YSZ+offs],...
    'Box','off','color','white');

TransverseTxt = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI H+YSZ+0.3 23 1.3],...
    'String','Transverse (X-Y)',...'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Tag','TransverseTxt',...
    'FontName','Comic Sans MS','FontSize',10,...
    'BackgroundColor',get(hMain,'Color'));
TransverseEdt = uicontrol(...
    'Parent',hMain,'Style','Edit',...
    'Units','char','Position',[XINI+22 H+YSZ+0.2 8 1.5],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''edit-transverse'',guidata(gcbo))',...
    'String','','Tag','TransverseEdt',...
    'HorizontalAlignment','center',...
    'TooltipString','set transverse slice',...
    'FontWeight','Bold');
TransverseSldr = uicontrol(...
    'Parent',hMain,'Style','slider',...
    'Units','char','Position',[XINI+XSZ*0.6 H+YSZ+0.2 XSZ*0.4 1.2],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''slider-transverse'',guidata(gcbo))',...
    'Tag','TransverseSldr','SliderStep',[1 4],...
    'TooltipString','transverse slice');
TransverseAxs = axes(...
    'Parent',hMain,'Tag','TransverseAxs',...
    'Units','char','Position',[XINI H XSZ YSZ],...
    'Box','off','Color','black','XTickLabel',{},'XTick',[],'YTickLabel',{},'YTick',[]);

H = 28;

SagitalTxt = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI+10+XSZ H+YSZ+0.3 20 1.35],...
    'String','Sagittal (Y-Z)',...'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Tag','SagitalTxt',...
    'FontName','Comic Sans MS','FontSize',10,...
    'BackgroundColor',get(hMain,'Color'));
SagitalEdt = uicontrol(...
    'Parent',hMain,'Style','Edit',...
    'Units','char','Position',[XINI+10+XSZ+22 H+YSZ+0.2 8 1.5],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''edit-sagital'',guidata(gcbo))',...
    'String','','Tag','SagitalEdt',...
    'HorizontalAlignment','center',...
    'TooltipString','set sagital slice',...
    'FontWeight','Bold');
SagitalSldr = uicontrol(...
    'Parent',hMain,'Style','slider',...
    'Units','char','Position',[XINI+10+XSZ*1.6 H+YSZ+0.2 XSZ*0.4 1.2],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''slider-sagital'',guidata(gcbo))',...
    'Tag','SagitalSldr','SliderStep',[1 4],...
    'TooltipString','sagital slice');
SagitalAxs = axes(...
    'Parent',hMain,'Tag','SagitalAxs',...
    'Units','char','Position',[XINI+10+XSZ H XSZ YSZ],...
    'Box','off','Color','black','XTickLabel',{},'XTick',[],'YTickLabel',{},'YTick',[]);

CoronalTxt = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI H+YSZ+0.3 30 1.3],...
    'String','Coronal (X-Z)',...'FontWeight','bold',...
    'HorizontalAlignment','left',...
    'FontName','Comic Sans MS','FontSize',10,...
    'Tag','CoronalTxt',...
    'BackgroundColor',get(hMain,'Color'));
CoronalEdt = uicontrol(...
    'Parent',hMain,'Style','Edit',...
    'Units','char','Position',[XINI+22 H+YSZ+0.2 8 1.5],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''edit-coronal'',guidata(gcbo))',...
    'String','','Tag','CoronalEdt',...
    'HorizontalAlignment','center',...
    'TooltipString','set coronal slice',...
    'FontWeight','Bold');
CoronalSldr = uicontrol(...
    'Parent',hMain,'Style','slider',...
    'Units','char','Position',[XINI+XSZ*0.6 H+YSZ+0.2 XSZ*0.4 1.2],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''slider-coronal'',guidata(gcbo))',...
    'Tag','CoronalSldr','SliderStep',[1 4],...
    'TooltipString','coronal slice');
CoronalAxs = axes(...
    'Parent',hMain,'Tag','CoronalAxs',...
    'Units','char','Position',[XINI H XSZ YSZ],...
    'Box','off','Color','black','XTickLabel',{},'XTick',[],'YTickLabel',{},'YTick',[]);


%--------------------------------------------------------------------------
% LIST WITH ALL THE ROIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = 28 + 20;
% XDSP=XDSP+XSZ+7;
XINI = XINI+XSZ*2+12+4;

PanelAtlas = uipanel(...
    'Parent',hMain,...'Title','Atlas','FontName','Arial','FontWeight','bold','FontSize',11,...
    'BorderType','beveledin ',...
    'BackgroundColor',get(PanelAna,'BackgroundColor'),...
    'Units','char','Position',[XINI-0.5 H-15.5 55 18]);
titleAtlas = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI H+0.8 53 1.4],...
    'String','Atlas','FontWeight','bold','FontSize',12,...
    'HorizontalAlignment','center',...
    ...'FontName', 'Corbel',...'FontSize',9,... 
    'BackgroundColor',get(PanelAna,'BackgroundColor'));

RoiList = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI+1 H-1 25 1.5],...
    'String','Brain ROIs',...
    'HorizontalAlignment','left',...
    'FontName','Arial','FontSize',10,...
    'Tag','RoiList',...
    'BackgroundColor',get(PanelAna,'BackgroundColor'));
AllRoiList = uicontrol(...
    'Parent',hMain,'Style','Listbox',...
    'Units','char','Position',[XINI+1 H-13 22 12],...
    'String',RoiNames,...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''AllRoiList'',guidata(gcbo))',...
    'HorizontalAlignment','left',...
    'FontName','Comic Sans MS','FontSize',10,...
    'Tag','AllRoiList','Background','white');


% BUTTON TRANSFER ROI
ButtonRoi = uicontrol(...
    'Parent',hMain,'Style','pushbutton',...
    'Units','char','Position',[XINI+23 H-7 8 1.5],...
    'String','<>','FontWeight','bold',...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''ButtonRoi'',guidata(gcbo))',...
    'HorizontalAlignment','left',...
    'FontSize', 11, 'FontName', 'Arial',...
    'Tag','ButtonRoi',...
    'BackgroundColor',get(PanelAna,'BackgroundColor'));


% LIST WITH THE ROIS TO SHOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RoiShow = uicontrol(...
    'Parent',hMain,'Style','Text',...
    'Units','char','Position',[XINI+23+8 H-1 20 1.5],...
    'String','Selected ROI',...
    'HorizontalAlignment','left',...
    'FontName','Arial','FontSize',10,...
    'Tag','RoiShow',...
    'BackgroundColor',get(PanelAna,'BackgroundColor'));
ShowRoiList = uicontrol(...
    'Parent',hMain,'Style','Listbox',...
    'Units','char','Position',[XINI+23+8 H-13 22 12],...
    'String',{},'Max',2,'Min',0,...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''showroilist'',guidata(gcbo))',...
    'HorizontalAlignment','left',...
    'FontName','Comic Sans MS','FontSize',9,...
    'Tag','ShowRoiList','Background','white');

% CHECK BOX FOR "show-atlas" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ShowAtlasCheck = uicontrol(...
    'Parent',hMain,'Style','Checkbox',...
    'Units','char','Position',[XINI+1 H-15 10 1.5],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''showatlas'',guidata(gcbo))',...
    'Tag','ShowAtlasCheck','Value',1,...
    'String','Atlas','FontWeight','bold',...
    'TooltipString','Show ROIs in 3D','BackgroundColor',get(PanelAna,'BackgroundColor'));


%--------------------------------------------------------------------------
%%% STAT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PanelStat = uipanel(...
    'Parent',hMain,...'Title','Atlas','FontName','Arial','FontWeight','bold','FontSize',11,...
    'BorderType','beveledin','tag','PanelStat',...
    'BackgroundColor',get(PanelAna,'BackgroundColor'),...
    'Units','char','Position',[XINI-0.5 1.5 55 30]);

H=3;XINI=1;
titleStat = uicontrol(...
    'Parent',PanelStat,'Style','Text',...
    'Units','char','Position',[0 28.3 53 1.4],...
    'String','Statistical Analysis','FontWeight','bold','FontSize',12,...
    'HorizontalAlignment','center',...
    ...'FontName', 'Corbel',...'FontSize',9,... 
    'BackgroundColor',get(PanelAna,'BackgroundColor'));

[ic,~] = imread('open.jpg'); 
icR = imresize(ic(:,:,1),[24 24]);
icG = imresize(ic(:,:,2),[24 24]);
icB = imresize(ic(:,:,3),[24 24]);
icRGB(:,:,1) = icR;
icRGB(:,:,2) = icG;
icRGB(:,:,3) = icB;
StatmapReadBtn = uicontrol(...
    'Parent',PanelStat,'Style','PushButton',...
    'Units','char','Position',[XINI H+23.5-1.7 6 2.3],...18 2],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''browse-statmap'',guidata(gcbo))',...
    'Tag','StatmapReadBtn','String','','CData',icRGB,...
    'TooltipString','Browse a mapfile','FontWeight','Bold');
clear icRGB;

StatmapEdt = uicontrol(...
    'Parent',PanelStat,'Style','Edit',...
    'Units','char','Position',[XINI+6.5 H+23.85-1.7 46 1.5],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''edit-statmap'',guidata(gcbo))',...
    'String','','Tag','StatmapEdt','Fontsize',6.5,...
    'HorizontalAlignment','left',...
    'TooltipString','Statistical mapfile',...
    'FontWeight','bold');

% CHECK BOX FOR "statistical analysis" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StatMapCheck = uicontrol(...
    'Parent',PanelStat,'Style','Checkbox',...
    'Units','char','Position',[XINI 22.5-1 20 1.3],...
    'Callback','atlas3dgui(''OrthoView_Callback'',gcbo,''showstatmaps'',guidata(gcbo))',...
    'Tag','StatMapCheck','Value',0,'enable','off',...
    'String','Statistic Maps','FontWeight','bold',...
    'TooltipString','Statistical maps for analysis','BackgroundColor',get(PanelStat,'BackgroundColor'));

% P-value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = H - 22;
AlphaTxt = uicontrol(...
    'Parent',PanelStat,'Style','Text',...
    'Units','char','Position',[XINI H+15.5-1 10 1.5],...
    'String','Alpha:','FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Tag','PvalueTxt',...
    'BackgroundColor',get(PanelStat,'BackgroundColor'));
AlphaEdt = uicontrol(...
    'Parent',PanelStat,'Style','Edit',...
    'Units','char','Position',[XINI+8 H+15.6-1 12 1.5],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''update-cmap'',guidata(gcbo))',...
    'String','0.000001','Tag','AlphaEdt',...
    'HorizontalAlignment','left',...
    'TooltipString','alpha for significance level',...
    'FontWeight','Bold');


panelrun = uibuttongroup('Parent',PanelStat,'visible','on','Units','char',...
    'Position',[XINI+21.1 H+15.2-1 31.8 6+1.5],'title','Run analysis through ROIs',...
    'BackgroundColor',get(PanelStat,'BackgroundColor'),...
    'BorderType','line');

runSelect = uicontrol('Style','Radio','String','selected in ROI list','tag','runSelect',...
    'Background',get(PanelStat,'BackgroundColor'),...
    'Units','char','pos',[0.5 2.9+1.5 25 1.2],'parent',panelrun,'value',1,'HandleVisibility','off');
runAll = uicontrol('Style','Radio','String','percentage threshold','tag','runAll',...
    'Background',get(PanelStat,'BackgroundColor'),...
    'Units','char','pos',[0.5 1+1.5 25 1.2],'parent',panelrun,'HandleVisibility','off');
connectEdt = uicontrol(...
    'Parent',panelrun,'Style','Edit',...
    'Units','char','Position',[25 0.9+1.5 6 1.5],...
    'String','...%','Tag','connectEdt',...
    'HorizontalAlignment','center',...
    'TooltipString','select percetange of connectivity',...
    'FontWeight','Bold');
runMinVox = uicontrol('Style','Radio','String','voxels threshold','tag','runMinVox',...
    'Background',get(PanelStat,'BackgroundColor'),...
    'Units','char','pos',[0.5 0.6 25 1.2],'parent',panelrun,'HandleVisibility','off');
MinVoxEdt = uicontrol(...
    'Parent',panelrun,'Style','Edit',...
    'Units','char','Position',[25 0.5 6 1.5],...
    'String','...','Tag','MinVoxEdt',...
    'HorizontalAlignment','center',...
    'TooltipString','select min number of connectivity voxels',...
    'FontWeight','Bold');

[ic,~] = imread('play.png'); 
icR = imresize(ic(:,:,1),[26 26]);
icG = imresize(ic(:,:,2),[26 26]);
icB = imresize(ic(:,:,3),[26 26]);
icRGB(:,:,1) = icR;
icRGB(:,:,2) = icG;
icRGB(:,:,3) = icB;
RunStat = uicontrol(...
    'Parent',PanelStat,'Style','PushButton',...
    'Units','char','Position',[XINI+45.1 H+9-1 6 2.3],...
    ...'Callback','atlas3dgui(''Main_Callback'',gcbo,''run-analysis'',guidata(gcbo))',...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''run-analysis'',guidata(gcbo))',...
    'Tag','RunStat','CData',icRGB,...
    'TooltipString','Run statistical analysis','FontWeight','Bold');
clear icRGB
runTxt = uicontrol(...
    'Parent',PanelStat,'Style','Text',...
    'Units','char','Position',[XINI+36.1 H+9.5-1 8 1.3],...
    'String','Start...','FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Tag','runTxt',...
    'BackgroundColor',get(PanelStat,'BackgroundColor'));

[ic,~] = imread('excel.jpg'); 
icR = imresize(ic(:,:,1),[26 26]);
icG = imresize(ic(:,:,2),[26 26]);
icB = imresize(ic(:,:,3),[26 26]);
icRGB(:,:,1) = icR;
icRGB(:,:,2) = icG;
icRGB(:,:,3) = icB;
ExportBtn = uicontrol(...
    'Parent',PanelStat,'Style','PushButton',...
    'Units','char','Position',[XINI+45.1 H+6-1 6 2.3],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''export-data'',guidata(gcbo))',...
    'Tag','ExportBtn','CData',icRGB,...
    'TooltipString','Export data','FontWeight','Bold');
clear icRGB

expTxt = uicontrol(...
    'Parent',PanelStat,'Style','Text',...
    'Units','char','Position',[XINI+33.1 H+6.5-1 11 1.3],...
    'String','Export to','FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Tag','expTxt',...
    'BackgroundColor',get(PanelStat,'BackgroundColor'));

[ic,~] = imread('analysis.jpg'); 
icR = imresize(ic(:,:,1),[26 26]);
icG = imresize(ic(:,:,2),[26 26]);
icB = imresize(ic(:,:,3),[26 26]);
icRGB(:,:,1) = icR;
icRGB(:,:,2) = icG;
icRGB(:,:,3) = icB;
ViewBtn = uicontrol(...
    'Parent',PanelStat,'Style','PushButton',...
    'Units','char','Position',[XINI+45.1 H+3-1 6 2.3],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''view-results'',guidata(gcbo))',...
    'Tag','ViewBtn','CData',icRGB,'enable','off',...
    'TooltipString','View analysis results','FontWeight','Bold');
clear icRGB

viewTxt = uicontrol(...
    'Parent',PanelStat,'Style','Text',...
    'Units','char','Position',[XINI+29.2 H+3.5-1 16 1.3],...
    'String','View results','FontWeight','bold',...
    'HorizontalAlignment','left',...
    'Tag','expTxt',...
    'BackgroundColor',get(PanelStat,'BackgroundColor'));


% AXES FOR STATISTICAL COLORBAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = 1;
StatColorbarAxs = axes(...
    'Parent',PanelStat,'Tag','StatColorbarAxs','Visible','on',...
    'units','char','Position',[XINI+0.5 H-0.5 6 15],...
    'FontSize',8,...
    'Box','off','YAxisLocation','right','XTickLabel',{},'XTick',[]);
StatColorbarMinEdt = uicontrol(...
    'Parent',PanelStat,'Style','Edit',...
    'Units','char','Position',[XINI+9 H-0.5 10 1.5],...
    'Callback','atlas3dgui(''Plot_Callback'',gcbo,[],[])',...
    'String','','Tag','StatColorbarMinEdt',...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''update-cmap'',guidata(gcbo))',...
    'HorizontalAlignment','center',...
    'TooltipString','set colorbar minimum',...
    'FontWeight','Bold');
StatColorbarMaxEdt = uicontrol(...
    'Parent',PanelStat,'Style','Edit',...
    'Units','char','Position',[XINI+9 H+15-1.5-0.5 10 1.5],...
    'String','','Tag','StatColorbarMaxEdt',...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''update-cmap'',guidata(gcbo))',...
    'HorizontalAlignment','center',...
    'TooltipString','set colorbar maximum',...
    'FontWeight','Bold');
% H = H+17-1.5;
uicontrol(...
    'Parent',PanelStat,'Style','Text','Enable','on',...
    'Units','char','Position',[XINI+10 H+15-4-1-0.5 12 1.5],...
    'String','ColorMap: ','FontWeight','bold',...
    'HorizontalAlignment','left',...
    'BackgroundColor',get(PanelStat,'BackgroundColor'));
cmaps = {'autumn','spring','winter','cool','jet','hot','hsv','copper'};

StatColormapCmb = uicontrol(...
    'Parent',PanelStat,'Style','Popupmenu','Enable','on',...
    'Units','char','Position',[XINI+10 H+15-4-2.5-0.5 13 1.5],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''update-cmap'',guidata(gcbo))',...
    'String',cmaps,'Value',1,'Tag','StatColormapCmb',...
    'TooltipString','Select colormap',...
    'FontWeight','bold');
clear cmaps idx;

% GAMMA STAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uicontrol(...
    'Parent',PanelStat,'Style','Text',...
    'Units','char','Position',[XINI+10 H+15-6-3-0.5 12 1.3],...
    'String','Gamma: ','FontWeight','bold',...
    'HorizontalAlignment','left',...
    'BackgroundColor',get(PanelStat,'BackgroundColor'));
StatGammaEdt = uicontrol(...
    'Parent',PanelStat,'Style','Edit',...
    'Units','char','Position',[XINI+10 H+15-6-4.5-0.5 13 1.5],...
    'Callback','atlas3dgui(''Main_Callback'',gcbo,''update-cmap'',guidata(gcbo))',...
    'String','1.0','Tag','StatGammaEdt',...
    'HorizontalAlignment','center',...
    'TooltipString','Set a gamma value for the statistical map',...
    'FontWeight','bold');


% get widgets handles at this moment
HANDLES = findobj(hMain);


% INITIALIZE THE APPLICATION
setappdata(hMain,'ANA',ANA);
setappdata(hMain,'colorGUI',color);

set(PanelStat,'visible','on');
setappdata(hMain,'Roi',Roi);
setappdata(hMain,'PERMUTE_VEC',PERMUTE_VEC);
Main_Callback(SagitalAxs,'init');
set(hMain,'visible','on');


% NOW SET "UNITS" OF ALL WIDGETS AS "NORMALIZED".
HANDLES = HANDLES(find(HANDLES ~= hMain));
set(HANDLES,'units','normalized');

set(hMain, 'Pointer', 'arrow')
% RETURNS THE WINDOW HANDLE IF REQUIRED.
if nargout,
  varargout{1} = hMain;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Main_Callback(hObject,eventdata,handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wgts = guihandles(hObject);
ANA  = getappdata(wgts.main,'ANA');

switch lower(eventdata),
 case {'init'}
  Roi  = getappdata(wgts.main,'Roi');
  PERMUTE_VEC = getappdata(wgts.main,'PERMUTE_VEC');

  STATMINV = -10;  STATMAXV = 10;
  set(wgts.StatColorbarMinEdt,'string',sprintf('%.1f',STATMINV));
  set(wgts.StatColorbarMaxEdt,'string',sprintf('%.1f',STATMAXV));
  setappdata(wgts.main,'STATMINV',STATMINV);
  setappdata(wgts.main,'STATMAXV',STATMAXV);
  
  MINV = 0;  MAXV = round(double(max(ANA.dat(:))) * 0.8/1000)*1000;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  if isfield(Roi,'roi')
        %create progress win rois
        [hWinPro, LoadingTxt, hBar, RegionsTxt] = createProgressWin...
            (numel(Roi.roinames),getappdata(wgts.main,'colorGUI'));
        drawnow
        set(LoadingTxt,'String','Creating ROI vector...');
        drawnow
        %create RoiMatrix vector
        RoiNames = Roi.roinames;
        RoiMatrix = cell(11,numel(RoiNames));
        nx = ANA.dim(1); ny = ANA.dim(2); nz = ANA.dim(3);
        COLORS = 'crgbkmy';
        j = 1;      
        for k = 1:numel(Roi.roi)
            RoiEl = Roi.roi{k};
            RoiMatrix{4,j}(:,:,RoiEl.slice) = RoiEl.mask; 

            if k == numel(Roi.roi) %brain
                RoiMatrix{1,j} = RoiEl.name;
                RoiMatrix{2,j} = RoiEl.px;
                if RoiEl.slice~= nz,  RoiMatrix{4,j}(:,:,nz) = false(nx,ny); end
                RoiMatrix{3,j} = RoiEl.py; 
                if ~isempty(PERMUTE_VEC),
                    RoiMatrix{4,j} = permute(RoiMatrix{4,j}, PERMUTE_VEC);
                end
                RoiMatrix{5,j} = [1 0.7 1];
                RoiMatrix{6,j} = false;
                continue,
            end    

            if strcmp(Roi.roi{k+1}.name,RoiEl.name), continue, end

            RoiMatrix{1,j} = RoiEl.name;
            RoiMatrix{2,j} = RoiEl.px;
            RoiMatrix{3,j} = RoiEl.py;    
            if RoiEl.slice~= nz,  RoiMatrix{4,j}(:,:,nz) = false(nx,ny); end
            
            if ~isempty(PERMUTE_VEC),
                RoiMatrix{4,j} = permute(RoiMatrix{4,j}, PERMUTE_VEC);
            end
            cidx = find(strcmpi(RoiNames,RoiEl.name));
            if isempty(cidx),  cidx = 1;  end
            cidx = mod(cidx(1),length(COLORS)) + 1;
            
            if strcmp(COLORS(cidx),'k')
                RoiMatrix{5,j} = [0.3 0.3 0.3];
            else
                RoiMatrix{5,j} = COLORS(cidx);
            end
            RoiMatrix{6,j} = false;
            j = j+1;
            
        end
        
        set(LoadingTxt,'String','Creating ROI vector...');
        drawnow
        win.hWinPro = hWinPro;
        win.hBar = hBar;
        win.RegionsTxt = RegionsTxt;
        win.LoadingTxt = LoadingTxt;
        setappdata(wgts.main,'hWin',win)
        setappdata(wgts.main,'RoiMatrix',RoiMatrix);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if MAXV == 0,  MAXV = 100;  end
  set(wgts.ColorbarMinEdt,'string',sprintf('%.1f',MINV));
  set(wgts.ColorbarMaxEdt,'string',sprintf('%.1f',MAXV));
  setappdata(wgts.main,'MINV',MINV);
  setappdata(wgts.main,'MAXV',MAXV);
  
  % initialize view
  if nargin < 3,
    OrthoView_Callback(hObject(1),'init');
    LightboxView_Callback(hObject(1),'init');
  else
    OrthoView_Callback(hObject(1),'init',handles);
    LightboxView_Callback(hObject(1),'init',handles);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  RoiMatrix = getappdata(wgts.main,'RoiMatrix');
  kref1 = 1;
  kref2 = 1;
  RoiRef1 = 'SHi';
  RoiRef2 = 'M1';
  for k=1:length(RoiMatrix)
    
      if strcmp(RoiMatrix{1,k},RoiRef1)
          kref1 = k; 
      elseif  strcmp(RoiMatrix{1,k},RoiRef2)
          kref2 = k; 
      end
      
  end
  bregma = findBregma(RoiMatrix{4,kref1},RoiMatrix{4,kref2});
%   bregma = [0 0 0];
  setappdata(wgts.main,'bregma',bregma);  
  CoordAxs.TriX = get(wgts.TriplotAxs,'XTick');
  CoordAxs.TriY = get(wgts.TriplotAxs,'YTick');
  CoordAxs.TriZ = get(wgts.TriplotAxs,'ZTick');
  CoordAxs.CorX = get(wgts.CoronalAxs,'XTick');
  CoordAxs.CorY = get(wgts.CoronalAxs,'YTick');
  CoordAxs.SagX = get(wgts.SagitalAxs,'XTick');
  CoordAxs.SagY = get(wgts.SagitalAxs,'YTick');
  CoordAxs.TranX = get(wgts.TransverseAxs,'XTick');
  CoordAxs.TranY = get(wgts.TransverseAxs,'YTick');
  setappdata(wgts.main, 'CoordAxs', CoordAxs);
  
 case {'update-clim'}
  MINV = str2double(get(wgts.ColorbarMinEdt,'String'));
  if isempty(MINV),
    MINV = getappdata(wgts.main,'MINV');
    set(wgts.ColorbarMinEdt,'String',sprintf('%.1f',MINV));
  end
  MAXV = str2double(get(wgts.ColorbarMaxEdt,'String'));
  if isempty(MAXV),
    MAXV = getappdata(wgts.main,'MAXV');
    set(wgts.ColorbarMaxEdt,'String',sprintf('%.1f',MAXV));
  end
  % update tick for colorbar
  GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
  if ~isempty(GRAHANDLE),
    ydat = [0:255]/255 * (MAXV - MINV) + MINV;
    set(GRAHANDLE.colorbar,'ydata',ydat);
    set(wgts.ColorbarAxs,'ylim',[MINV MAXV]);
  end
  % update color for images
  haxs = [wgts.SagitalAxs, wgts.CoronalAxs, wgts.TransverseAxs, wgts.TriplotAxs, wgts.LightboxAxs];
  set(haxs,'clim',[MINV MAXV]);
  setappdata(wgts.main,'MINV',MINV);
  setappdata(wgts.main,'MAXV',MAXV);
  
 case {'view-mode'}
  ViewMode = get(wgts.ViewModeCmb,'String');
  ViewMode = ViewMode{get(wgts.ViewModeCmb,'Value')};
  hL = [wgts.LightboxAxs];
  hO = [wgts.CoronalTxt, wgts.CoronalEdt, wgts.CoronalSldr, wgts.CoronalAxs,...
        wgts.SagitalTxt, wgts.SagitalEdt, wgts.SagitalSldr, wgts.SagitalAxs,...
        wgts.TransverseTxt, wgts.TransverseEdt, wgts.TransverseSldr, wgts.TransverseAxs,...
        wgts.CrosshairCheck, wgts.ButtonRoi, wgts.ShowAtlasCheck, wgts.ShowMRICheck];
  
  if strcmpi(ViewMode,'Orthogonal'),
    set(wgts.ViewPageList,'visible','off');
    set(hL,'visible','off');
    set(findobj(hL),'visible','off');
    set(hO,'visible','on');
    set(findobj([wgts.TriplotAxs],'Type','axes'),'visible','on');
    hImg = findobj([wgts.CoronalAxs, wgts.SagitalAxs, wgts.TransverseAxs],'Type','image');
    set(hImg,'visible','on');
    hslice = findobj([wgts.CoronalAxs, wgts.SagitalAxs, wgts.TransverseAxs],'tag','roi');
    set(hslice,'visible','on');
    hline = findobj([wgts.CoronalAxs, wgts.SagitalAxs, wgts.TransverseAxs],'tag','line');
    set(hline,'visible','on');
    if get(wgts.ShowMRICheck,'value') == 1
        htri = findobj([wgts.TriplotAxs],'Type','surface');
        set(htri,'visible','on');
    end
    if get(wgts.ShowAtlasCheck,'value') == 1
        hRois = findobj([wgts.TriplotAxs],'Type','patch');
        set(hRois,'visible','on');
    end
  else
    set(wgts.ViewPageList,'visible','on');
    set(hL,'visible','on');
    set(findobj(hL),'visible','on');
    set(hO,'visible','off');
    h = findobj([wgts.CoronalAxs, wgts.SagitalAxs, wgts.TransverseAxs, wgts.TriplotAxs]);
    set(h,'visible','off');
    LightboxView_Callback(hObject,'init',[]);
    LightboxView_Callback(hObject,'redraw',[]);
  end

 case {'view-page'}
  ViewMode = get(wgts.ViewModeCmb,'String');
  ViewMode = ViewMode{get(wgts.ViewModeCmb,'Value')};
  if ~isempty(strfind(ViewMode,'Lightbox')),
    LightboxView_Callback(hObject,'redraw',[]);
  end
  
 case {'dir-reverse'}
  ViewMode = get(wgts.ViewModeCmb,'String');
  ViewMode = ViewMode{get(wgts.ViewModeCmb,'Value')};
  if ~isempty(strfind(ViewMode,'lightbox')),
    LightboxView_Callback(hObject,'redraw',[]);
  else
    OrthoView_Callback(hObject,'dir-reverse',[]);
  end
 
 case {'change-coord'}
    if get(wgts.SelBreg,'Value')
        subBregmaCoord(wgts, getappdata(wgts.main,'bregma'));
    elseif get(wgts.SelVox,'Value')   
        CoordAxs = getappdata(wgts.main,'CoordAxs');
        CorX = CoordAxs.CorX;
        CorY = CoordAxs.CorY;
        SagX = CoordAxs.SagX;
        SagY = CoordAxs.SagY;
        TranX = CoordAxs.TranX;
        TranY = CoordAxs.TranY;
        set(wgts.CoronalAxs,'XTick',CorX,'XTicklabel',CorX,...
            'YTick',CorY,'YTicklabel',CorY);
        set(wgts.SagitalAxs,'XTick',SagX,'XTicklabel',SagX,...
            'YTick',SagY,'YTicklabel',SagY);
        set(wgts.TransverseAxs,'XTick',TranX,'XTicklabel',TranX,...
            'YTick',TranY,'YTicklabel',TranY);
        set(wgts.TriplotAxs,'XTick',CoordAxs.TriX,'XTicklabel',CoordAxs.TriX,...
            'YTick',CoordAxs.TriY,'YTicklabel',CoordAxs.TriY,'ZTick',...
            CoordAxs.TriZ,'ZTicklabel',CoordAxs.TriZ);
    end
   
 case {'grid-check'}
     axes(wgts.SagitalAxs), grid
     axes(wgts.CoronalAxs), grid
     axes(wgts.TransverseAxs), grid
     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 case {'buttonroi'}
      
    ViewMode = get(wgts.ViewModeCmb,'String');
    ViewMode = ViewMode{get(wgts.ViewModeCmb,'Value')};
    
    if strcmp(ViewMode,'Orthogonal') 
        set(wgts.ButtonRoi,'Enable','off')
        RoiMatrix  = getappdata(wgts.main,'RoiMatrix');
        Roi1Indx = get(wgts.AllRoiList,'Value');
        Roi2Indx = get(wgts.ShowRoiList,'Value');
        if ~isempty(Roi1Indx)
            LeftList = get(wgts.AllRoiList,'String');
            SelRoi = LeftList{Roi1Indx};
        elseif ~isempty(Roi2Indx)
            RightList = get(wgts.ShowRoiList,'String');
            SelRoi = RightList{Roi2Indx};
        else
            set(wgts.ButtonRoi,'Enable','on')
            return;
        end
        
        ridx = find(strcmp(SelRoi,RoiMatrix(1,:)));
        ridx = ridx(1);                             %multiples nombres!!!!
        RoiMatrix{6,ridx} =~ RoiMatrix{6,ridx}; 
        setappdata(wgts.main,'RoiMatrix',RoiMatrix);
        left = ~cell2mat(RoiMatrix(6,:));
        right = cell2mat(RoiMatrix(6,:));
        axes(wgts.TriplotAxs);

        if ~isempty(Roi1Indx) %Complete RoiList -> Selected Roi's

            if get(wgts.ShowAtlasCheck,'value'), show = 'on'; 
            else show = 'off';end

                if isempty(RoiMatrix{8,ridx})
                    RoiMatrix{8,ridx} = roi3dMatrix(RoiMatrix(:,ridx), 0, show,'roi', mfilename);
                else
                    set(RoiMatrix{8,ridx},'visible',show)
                end
                iX = str2double(get(wgts.SagitalEdt,'String'));
                axes(wgts.SagitalAxs); 
                hold on;
                RoiMatrix(:,ridx) = plotRoiMatrix(RoiMatrix(:,ridx), iX, 'sagital',show, mfilename);
                hold off; 
                iY = str2double(get(wgts.CoronalEdt,'String'));
                axes(wgts.CoronalAxs); 
                hold on; 
                RoiMatrix(:,ridx) = plotRoiMatrix(RoiMatrix(:,ridx), iY, 'coronal',show, mfilename);
                hold off;
                iZ = str2double(get(wgts.TransverseEdt,'String'));
                axes(wgts.TransverseAxs); 
                hold on; 
                RoiMatrix(:,ridx) = plotRoiMatrix(RoiMatrix(:,ridx), iZ, 'transverse',show, mfilename);
                hold off;
                setappdata(wgts.main,'RoiMatrix',RoiMatrix);
%                  subRedrawCrosshair(wgts,ANA)
            if numel(LeftList) > Roi1Indx
                set(wgts.AllRoiList,'string',RoiMatrix(1,left)); 
            elseif numel(LeftList) == Roi1Indx
                if Roi1Indx == 1
                     set(wgts.AllRoiList,'Value',[]);
                else
                     set(wgts.AllRoiList,'Value',(Roi1Indx-1));   
                end
                set(wgts.AllRoiList,'string',RoiMatrix(1,left)); 
            end
            set(wgts.ShowRoiList,'string',RoiMatrix(1,right));

        elseif ~isempty(Roi2Indx) %Selected Roi's -> Complete RoiList 

            set(RoiMatrix{8,ridx},'visible','off')
            delete(RoiMatrix{8,ridx});
            RoiMatrix{8,ridx}=[];
            set(RoiMatrix{9,ridx},'XData',[],'YData',[])
            delete(RoiMatrix{9,ridx});
            RoiMatrix{9,ridx}=[];
            set(RoiMatrix{10,ridx},'XData',[],'YData',[])
            delete(RoiMatrix{10,ridx});
            RoiMatrix{10,ridx}=[];
            set(RoiMatrix{11,ridx},'XData',[],'YData',[])
            delete(RoiMatrix{11,ridx});
            RoiMatrix{11,ridx}=[];
            setappdata(wgts.main,'RoiMatrix',RoiMatrix);
            
            TXTHANDLE = getappdata(wgts.main,'TXTHANDLE');
            showTextRoiMatrix('sagital', 0, TXTHANDLE, 0);
            showTextRoiMatrix('coronal', 0, TXTHANDLE, 0);
            showTextRoiMatrix('transverse', 0, TXTHANDLE, 0);
            if numel(RightList) > Roi2Indx
                set(wgts.ShowRoiList,'string',RoiMatrix(1,right));

            elseif numel(RightList) == Roi2Indx
                if Roi2Indx == 1
                     set(wgts.ShowRoiList,'Value',[]);
                else
                     set(wgts.ShowRoiList,'Value',(Roi2Indx-1));   
                end
                set(wgts.ShowRoiList,'string',RoiMatrix(1,right)); 
            end
            set(wgts.AllRoiList,'string',RoiMatrix(1,left));
        end
        set(wgts.ButtonRoi,'Enable','on')
    end
       
  case {'allroilist'}
   set(wgts.ShowRoiList,'Max',2,'Min',0,'Value',[]);
   
  case {'showroilist'}
   set(wgts.AllRoiList,'Max',2,'Min',0,'Value',[]);
   
  case {'export-data'}
    test = getappdata(wgts.main,'test');
    if ~isempty(test)
        sprintf('exporting data...')
        if isfield(test,'tbase')
            v={'','Voxels','% of Brain','Average Intensity','Std.Deviation', ...
            'Mn Voxels','Connectivity in %','Average Significance','Std. Deviation'};
        else
            v={'','Voxels','% of Brain','Average Intensity','Std.Deviation'};
        end
        
        namefile = sprintf('analysis_%s_%s',date);
%         infomat = cell(8,2);
%         infomat{:,2} = ['',test.tbase,test.tsel,test.df,test.smooth,test.smooth_hsize,...
%                         test.norm,test.date,test.time];
        matcell = [v;test.roi];
        xlswrite(namefile, matcell, 'Statistical Analysis')
    end
    
  case{'run-analysis'}
    STATS = getappdata(wgts.main,'STATS');
    ANA  = getappdata(wgts.main,'ANA');
    RoiMatrix = getappdata(wgts.main,'RoiMatrix');
    activeRoi = RoiMatrix(:,cell2mat(RoiMatrix(6,:)));
    set(wgts.main,'Pointer','watch')
    
    if ~isempty(STATS)
        alpha = str2double(get(wgts.AlphaEdt,'String'));
        statmat = (STATS.p<alpha).*STATS.dat;
        brainRoi = RoiMatrix{4,end};
        brainRoi = brainRoi(brainRoi==1);
        connect_th = str2double(get(wgts.connectEdt,'String'));
        minvox_th = str2double(get(wgts.MinVoxEdt,'String'));
        if get(wgts.runSelect,'value')&&~isempty(activeRoi)
            testRoi = cell(size(activeRoi,2),8);
            for k = 1:size(activeRoi,2)
                vol = activeRoi{4,k}.*ANA.dat;
                vol = vol(vol~=0);
                volstat = activeRoi{4,k}.*statmat;
                volstat = volstat(volstat~=0);
                if ~isempty(vol)
                testRoi{k,1} = activeRoi{1,k};
                testRoi{k,2} = numel(vol(:));
                testRoi{k,3} = numel(vol(:))/numel(brainRoi).*100;
                testRoi{k,4} = mean(vol(:));
                testRoi{k,5} = std(vol(:));
                testRoi{k,6} = numel(volstat);
                testRoi{k,7} = numel(volstat)/numel(vol(:)).*100;
                testRoi{k,8} = mean(volstat(:));
                testRoi{k,9} = std(volstat(:));
                end
            end
            test.roi = testRoi;
            test.tbase = STATS.tbase;
            test.tsel = STATS.tsel;
            test.df = STATS.df;
            test.smooth = STATS.flags.smooth;
            test.smooth_hsize = STATS.flags.smooth_hsize;
            test.norm = STATS.flags.normalize;
            test.type = 'roilist';
            test.date = date;
            test.time = datestr(now,'HH:MM:SS');
            setappdata(wgts.main,'test',test);
            dispResults(test,get(wgts.main,'Color'),'roilist');

        elseif get(wgts.runAll,'value')&&~isnan(connect_th)
            j=1;
            for k = 1:size(RoiMatrix,2)
                volstat = RoiMatrix{4,k}.*statmat;
                volstat = volstat(volstat~=0);
                vol = RoiMatrix{4,k}.*ANA.dat;
                vol = vol(vol~=0);
                ratio = numel(volstat(:))./numel(vol).*100;
                if ratio < connect_th, continue, end
                if ~isempty(vol)
                    testRoi{j,1} = RoiMatrix{1,k};
                    testRoi{j,2} = numel(vol(:));
                    testRoi{j,3} = numel(vol(:))/numel(brainRoi).*100;
                    testRoi{j,4} = mean(vol(:));
                    testRoi{j,5} = std(vol(:));
                    testRoi{j,6} = numel(volstat);
                    testRoi{j,7} = numel(volstat)/numel(vol(:)).*100;
                    testRoi{j,8} = mean(volstat(:));
                    testRoi{j,9} = std(volstat(:));
                    j = j+1;
                end
            end
            test.roi = testRoi;
            test.tbase = STATS.tbase;
            test.tsel = STATS.tsel;
            test.df = STATS.df;
            test.smooth = STATS.flags.smooth;
            test.smooth_hsize = STATS.flags.smooth_hsize;
            test.norm = STATS.flags.normalize;
            test.type = 'percent';
            test.date = date;
            test.time = datestr(now,'HH:MM:SS');
            setappdata(wgts.main,'test',test);
            dispResults(test,get(wgts.main,'Color'),'percent');
            set(wgts.ViewBtn,'enable','on');
        elseif get(wgts.runMinVox,'value')&&~isnan(minvox_th)
            j=1;
            for k = 1:size(RoiMatrix,2)-1
                volstat = RoiMatrix{4,k}.*statmat;
                volstat = volstat(volstat~=0);
                vol = RoiMatrix{4,k}.*ANA.dat;
                vol = vol(vol~=0);
                if numel(volstat) < minvox_th, continue, end
                if ~isempty(vol)
                    testRoi{j,1} = RoiMatrix{1,k};
                    testRoi{j,2} = numel(vol(:));
                    testRoi{j,3} = numel(vol(:))/numel(brainRoi).*100;
                    testRoi{j,4} = mean(vol(:));
                    testRoi{j,5} = std(vol(:));
                    testRoi{j,6} = numel(volstat);
                    testRoi{j,7} = numel(volstat)/numel(vol(:)).*100;
                    testRoi{j,8} = mean(volstat(:));
                    testRoi{j,9} = std(volstat(:));
                    j = j+1;
                end
            end
            test.roi = testRoi;
            test.tbase = STATS.tbase;
            test.tsel = STATS.tsel;
            test.df = STATS.df;
            test.smooth = STATS.flags.smooth;
            test.smooth_hsize = STATS.flags.smooth_hsize;
            test.norm = STATS.flags.normalize;
            test.type = 'minvox';
            test.date = date;
            test.time = datestr(now,'HH:MM:SS');
            setappdata(wgts.main,'test',test);
            dispResults(test,get(wgts.main,'Color'),'minvox');
            set(wgts.ViewBtn,'enable','on');
        end
    elseif ~isempty(activeRoi)
        brainRoi = RoiMatrix{4,end};
        brainRoi = brainRoi(brainRoi==1);
        if get(wgts.runSelect,'value')&&~isempty(activeRoi)
            testRoi = cell(size(activeRoi,2),5);
            for k = 1:size(activeRoi,2)
                vol = activeRoi{4,k}.*ANA.dat;
                vol = vol(vol~=0);
                if ~isempty(vol)
                testRoi{k,1} = activeRoi{1,k};
                testRoi{k,2} = numel(vol(:));
                testRoi{k,3} = numel(vol(:))/numel(brainRoi).*100;
                testRoi{k,4} = mean(vol(:));
                testRoi{k,5} = std(vol(:));
                end
            end
            test.roi = testRoi;
            test.date = date;
            test.time = datestr(now,'HH:MM:SS');
            test.type = 'anatomical';
            setappdata(wgts.main,'test',test);
            dispResults(test,get(wgts.main,'Color'),'anatomical');
        end
        set(wgts.ViewBtn,'enable','on');
    end
    set(wgts.main,'Pointer','arrow');
    
  case {'view-results'}
      test = getappdata(wgts.main,'test');
      if ~isempty(test)
          dispResults(test,get(wgts.main,'Color'),test.type);
      end
        
  case{'browse-statmap'}
      mapfile = get(wgts.StatmapEdt,'String');
      [mapfile, pathname] = uigetfile(...
          {'*.mat', 'Mat-files (*.mat)'}, 'Pick a statistial map file',mapfile);
      if ~isequal(mapfile,0) && ~isequal(pathname,0),
        fullpathname = fullfile(pathname,mapfile);
        set(wgts.StatmapEdt,'String',fullpathname);
        Main_Callback(hObject,'edit-statmap',[]);
      end
      
 case {'edit-statmap','select-roi'}
      
      mapfile = get(wgts.StatmapEdt,'String');
      if ~exist(mapfile,'file'),
        setappdata(wgts.main,'STATS',[]);
        return
      end
      STATS = load(mapfile);
      fname = fieldnames(STATS);
      STATS = STATS.(fname{1});
      if iscell(STATS),  STATS = STATS{1};  end
      
      if isfield(STATS,'dat')&&isfield(STATS,'p')
      STATMINV = 0;
      STATMAXV = round(double(max(STATS.dat(:))) * 0.35/10)*10;
      permute_vec = getappdata(wgts.main,'PERMUTE_VEC');
      if ~isempty(permute_vec),
        STATS.dat = permute(STATS.dat, permute_vec);
        STATS.p   = permute(STATS.p,   permute_vec);
      end
      mask = (ANA.dat~=0);
      STATS.dat = STATS.dat.*mask;
      STATS.p = STATS.p.*mask;
      setappdata(wgts.main,'STATS',STATS);
      set(wgts.StatColorbarMinEdt,'string',sprintf('%.1f',STATMINV));
      set(wgts.StatColorbarMaxEdt,'string',sprintf('%.1f',STATMAXV));
      Main_Callback(hObject,'update-cmap',[]);
      setappdata(wgts.main,'STATMINV',STATMINV);
      setappdata(wgts.main,'STATMAXV',STATMAXV);
      set(wgts.StatMapCheck,'Enable','on');
      set(wgts.StatMapCheck,'value',1);

      OrthoView_Callback(hObject,'showstatmaps',[]);
      end
      
 case {'update-cmap'}
      GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
      STATMINV = str2double(get(wgts.StatColorbarMinEdt,'String'));
      if isempty(STATMINV),
        STATMINV = getappdata(wgts.main,'STATMINV');
        set(wgts.StatColorbarMinEdt,'String',sprintf('%.1f',STATMINV));
      end
      STATMAXV = str2double(get(wgts.StatColorbarMaxEdt,'String'));
      if isempty(STATMAXV),
        STATMAXV = getappdata(wgts.main,'STATMAXV');
        set(wgts.StatColorbarMaxEdt,'String',sprintf('%.1f',STATMAXV));
      end
      setappdata(wgts.main,'STATMINV',STATMINV);
      setappdata(wgts.main,'STATMAXV',STATMAXV);
      statcmap = subGetColormap(wgts,(STATMAXV-STATMINV));
      axes(wgts.StatColorbarAxs), imagesc(statcmap);
      set(wgts.StatColorbarAxs,'tag','StatColorbarAxs');	% set this again, some will reset.
      set(wgts.StatColorbarAxs,'xtick',[],'ytick',[]);
      setappdata(wgts.main,'CMAP',statcmap);

      if isfield(GRAHANDLE,'triStat')
        set(wgts.main,'Pointer','watch')          
        axes(wgts.TriplotAxs)        
        STATS = getappdata(wgts.main,'STATS');
        delete(GRAHANDLE.triStat);
        alpha = str2double(get(wgts.AlphaEdt,'string'));
        statmat = (STATS.p<alpha).*STATS.dat;
        cmapstr = get(wgts.StatColormapCmb,'String');
        cmapstr = cmapstr{get(wgts.StatColormapCmb,'value')};
        GRAHANDLE.triStat = statmap3d(statmat, 0, 'off', cmapstr, mfilename);
        setappdata(wgts.main,'GRAHANDLE',GRAHANDLE)  
        set(wgts.main,'Pointer','arrow')       
      end
      if get(wgts.StatMapCheck,'Value')
        set(GRAHANDLE.triStat,'visible','on');
        OrthoView_Callback(hObject,'slider-sagital',[]);
        OrthoView_Callback(hObject,'slider-coronal',[]);
        OrthoView_Callback(hObject,'slider-transverse',[]);
      end

 otherwise
end
  
return;


       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to handle orthogonal view
function OrthoView_Callback(hObject,eventdata,handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wgts = guihandles(get(hObject,'Parent'));
ANA  = getappdata(wgts.main,'ANA');
STATS  = getappdata(wgts.main,'STATS');
MINV = getappdata(wgts.main,'MINV');
MAXV = getappdata(wgts.main,'MAXV');
STATMINV = getappdata(wgts.main,'STATMINV');
STATMAXV = getappdata(wgts.main,'STATMAXV');

switch lower(eventdata),
 case {'init'}
  
  RoiMatrix = getappdata(wgts.main,'RoiMatrix');
  Roi = getappdata(wgts.main,'Roi');

%   iX = 1;  iY = 1;  iZ = 1;  
  nX = size(ANA.dat,1);  nY = size(ANA.dat,2);  nZ = size(ANA.dat,3);  
  iX = round(nX/2);  iY = round(nY/2);  iZ = round(nZ/2);
  
  statcmap = subGetColormap(wgts,(STATMAXV-STATMINV));
  axes(wgts.StatColorbarAxs), imagesc(statcmap);
  set(wgts.StatColorbarAxs,'tag','StatColorbarAxs');	% set this again, some will reset.
  set(wgts.StatColorbarAxs,'xtick',[],'ytick',[]);
  setappdata(wgts.main,'CMAP',statcmap);
  
  % set slider edit value
  set(wgts.SagitalEdt,   'String', sprintf('%d',iX));
  set(wgts.CoronalEdt,   'String', sprintf('%d',iY));
  set(wgts.TransverseEdt,'String', sprintf('%d',iZ));
  % set slider, add +0.01 to prevent error.
  set(wgts.SagitalSldr,   'Min',1,'Max',nX+0.01,'Value',iX);
  set(wgts.CoronalSldr,   'Min',1,'Max',nY+0.01,'Value',iY);
  set(wgts.TransverseSldr,'Min',1,'Max',nZ+0.01,'Value',iZ);
  % set slider step, it is normalized from 0 to 1, not min/max
  set(wgts.SagitalSldr,   'SliderStep',[1, 2]/max(1,nX));
  set(wgts.CoronalSldr,   'SliderStep',[1, 2]/max(1,nY));
  set(wgts.TransverseSldr,'SliderStep',[1, 2]/max(1,nZ));
  
  cmap = gray(256);
  gammav = str2double(get(wgts.GammaEdt,'String'));
  if ~isempty(gammav),
    cmap = cmap.^(1/gammav);
  end
  
%   AXISCOLOR = [0.8 0.2 0.8];
  AXISCOLOR = [0 0.2 0.8];

% now draw images
  axes(wgts.SagitalAxs);     colormap(cmap);
  hSag = imagesc(1:nY,1:nZ,squeeze(ANA.dat(iX,:,:))');
  axis image;
  set(hSag,'ButtonDownFcn',...
      'atlas3dgui(''OrthoView_Callback'',gcbo,''button-sagital'',guidata(gcbo))');
  set(wgts.SagitalAxs,'tag','SagitalAxs');	% set this again, some will reset.
  SagText = text(0,0,'','fontsize',10,'tag','roi','fontweight','bold');
  
  axes(wgts.CoronalAxs);     colormap(cmap);
  hCor = imagesc(1:nX,1:nZ,squeeze(ANA.dat(:,iY,:))');
  axis image;
  set(hCor,'ButtonDownFcn',...
      'atlas3dgui(''OrthoView_Callback'',gcbo,''button-coronal'',guidata(gcbo))');
  set(wgts.CoronalAxs,'tag','CoronalAxs');  % set this again, some will reset.
  CorText = text(0,0,'','fontsize',10,'tag','roi','fontweight','bold');
  
  axes(wgts.TransverseAxs);  colormap(cmap);
  hTra = imagesc(1:nX,1:nY,squeeze(ANA.dat(:,:,iZ))');
  axis image;
  set(hTra,'ButtonDownFcn',...
      'atlas3dgui(''OrthoView_Callback'',gcbo,''button-transverse'',guidata(gcbo))');
  set(wgts.TransverseAxs,'tag','TransverseAxs');	% set this again, some will reset.
  TranText = text(0,0,'','fontsize',10,'tag','roi','fontweight','bold');
  
  % now draw a color bar
  axes(wgts.ColorbarAxs);
  ydat = [0:255]/255 * (MAXV - MINV) + MINV;
  hColorbar = imagesc(1,ydat,[0:255]'); colormap(cmap);
  set(wgts.ColorbarAxs,'Tag','ColorbarAxs');  % set this again, some will reset.
  set(wgts.ColorbarAxs,'ylim',[MINV MAXV],'xcolor',AXISCOLOR,'ycolor',AXISCOLOR,...
                    'YAxisLocation','right','XTickLabel',{},'XTick',[],'Ydir','normal');
  
  haxs = [wgts.SagitalAxs, wgts.CoronalAxs, wgts.TransverseAxs];
  set(haxs,'fontsize',8,'xcolor',AXISCOLOR,'ycolor',AXISCOLOR);

  GRAHANDLE.sagital    = hSag;
  GRAHANDLE.coronal    = hCor;
  GRAHANDLE.transverse = hTra;
  GRAHANDLE.colorbar   = hColorbar;
%   GRAHANDLE.Statcolorbar   = hStatColorbar;

  TXTHANDLE.sagital    = SagText;
  TXTHANDLE.coronal    = CorText;
  TXTHANDLE.transverse = TranText;
  
  % draw crosshair(s)
  axes(wgts.SagitalAxs);
  hSagV = line([iY iY],[ 1 nZ],'color','b','tag','line');
  hSagH = line([ 1 nY],[iZ iZ],'color','b','tag','line');
  set([hSagV hSagH],...
      'ButtonDownFcn','atlas3dgui(''OrthoView_Callback'',gcbo,''button-sagital'',guidata(gcbo))');
  axes(wgts.CoronalAxs);
  hCorV = line([iX iX],[ 1 nZ],'color','b','tag','line');
  hCorH = line([ 1 nX],[iZ iZ],'color','b','tag','line');
  set([hCorV hCorH],...
      'ButtonDownFcn','atlas3dgui(''OrthoView_Callback'',gcbo,''button-coronal'',guidata(gcbo))');
  axes(wgts.TransverseAxs);
  hTraV = line([iX iX],[ 1 nY],'color','b','tag','line');
  hTraH = line([ 1 nX],[iY iY],'color','b','tag','line');
  set([hTraV hTraH],...
      'ButtonDownFcn','atlas3dgui(''OrthoView_Callback'',gcbo,''button-transverse'',guidata(gcbo))');
  if get(wgts.CrosshairCheck,'Value') == 0,
    set([hSagV hSagH hCorV hCorH hTraV hTraH],'visible','off');
  end
  GRAHANDLE.sagitalV    = hSagV;
  GRAHANDLE.sagitalH    = hSagH;
  GRAHANDLE.coronalV    = hCorV;
  GRAHANDLE.coronalH    = hCorH;
  GRAHANDLE.transverseV = hTraV;
  GRAHANDLE.transverseH = hTraH;
  
  % tri-plot
  axes(wgts.TriplotAxs);
  tmpv = squeeze(ANA.dat(iX,:,:));
  [xi,yi,zi] = meshgrid(iX,1:nY,1:nZ);
  hSag = surface(...
      'xdata',reshape(xi,[nY,nZ]),'ydata',reshape(yi,[nY,nZ]),'zdata',reshape(zi,[nY,nZ]),...
      'cdata',tmpv,...
      'facecolor','texturemap','edgecolor','none',...
      'CDataMapping','scaled','linestyle','none');
  tmpv = squeeze(ANA.dat(:,iY,:));
  [xi,yi,zi] = meshgrid(1:nX,iY,1:nZ);
  hCor = surface(...
      'xdata',reshape(xi,[nX,nZ]),'ydata',reshape(yi,[nX,nZ]),'zdata',reshape(zi,[nX,nZ]),...
      'cdata',tmpv,...
    'facecolor','texturemap','edgecolor','none',...
      'CDataMapping','scaled','linestyle','none');
  tmpv = squeeze(ANA.dat(:,:,iZ));
  [xi,yi,zi] = meshgrid(1:nX,1:nY,iZ);
  hTra = surface(...
      'xdata',1:nX,'ydata',1:nY,'zdata',reshape(zi,[nY,nX]),...
      'cdata',tmpv',...
      'facecolor','texturemap','edgecolor','none',...
      'CDataMapping','scaled','linestyle','none');
 
  set(gca,'Tag','TriplotAxs');
  set(gca,'clim',[MINV MAXV],'fontsize',8,...
          'xlim',[1 nX],'ylim',[1 nY],'zlim',[1 nZ],'zdir','reverse');
  view(320,36);grid on;  axis equal;
  xlabel('X'); ylabel('Y');  zlabel('Z');
  set(wgts.TriplotAxs,'clim',[MINV MAXV],'fontsize',8,'xcolor',AXISCOLOR,'ycolor',AXISCOLOR,'zcolor',AXISCOLOR);

  set([gca hSag hCor hTra],...
      'ButtonDownFcn','atlas3dgui(''OrthoView_Callback'',gcbo,''button-triplot'',guidata(gcbo))');
  TriText = text(0,0,0,'','fontsize',10,'tag','roi','fontweight','bold');
  
  GRAHANDLE.triSagital    = hSag;
  GRAHANDLE.triCoronal    = hCor;
  GRAHANDLE.triTransverse = hTra;
  
  TXTHANDLE.triText       = TriText;
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %load rois in a matrix

  if isfield(Roi,'roi')
      hWin = getappdata(wgts.main, 'hWin');  
      figure(hWin.hWinPro)
      set(hWin.LoadingTxt,'String','Loading Atlas regions...');
      for i=1:size(RoiMatrix,2)
        VOI = RoiMatrix{4,i};
        VOI = permute(VOI,[2 1 3]);
        p2=isosurface(VOI,.4);
        RoiMatrix{7,i} = p2;
        set(hWin.hBar,'XData',[0 i i 0]);
        set(hWin.RegionsTxt,'String',sprintf('%s...\n%g of %g',RoiMatrix{1,i},i,size(RoiMatrix,2)));
        drawnow
      end
      set(hWin.LoadingTxt,'String','Loading Atlas regions...done');
      drawnow
      roidir = which(Roi.fileName);
      [pathstr name ext] = fileparts(roidir);
      matfile = [pathstr '\RoiMatrix' Roi.fileName(4:end-4)];
      set(hWin.LoadingTxt,'String','Saving as a matfile in:');
      set(hWin.RegionsTxt,'String',[matfile '.mat' '...']);
      drawnow
      save(matfile,'RoiMatrix');
      set(hWin.LoadingTxt,'String','Saving as a matfile...done');
      fprintf('Saved as: %s\n',[matfile '.mat'])
      drawnow
      pause(1)
      delete(hWin.hWinPro);
  else
      RoiMatrix = Roi;
  end

  setappdata(wgts.main,'RoiMatrix',RoiMatrix)
  setappdata(wgts.main,'GRAHANDLE',GRAHANDLE);
  setappdata(wgts.main,'TXTHANDLE',TXTHANDLE);
  OrthoView_Callback(hObject,'dir-reverse',[]);

 case {'slider-sagital'}
  GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
  TXTHANDLE = getappdata(wgts.main,'TXTHANDLE');
  RoiMatrix = getappdata(wgts.main,'RoiMatrix');
  showTextRoiMatrix('sagital', 0, TXTHANDLE, 0);
  gammav = str2double(get(wgts.GammaEdt,'String'));
  gamma_stat = str2double(get(wgts.StatGammaEdt,'String'));
  cmapstr = get(wgts.StatColormapCmb,'String');
  cmapstr = cmapstr{get(wgts.StatColormapCmb,'value')};
  
  if ~isempty(GRAHANDLE)
    iX = round(get(wgts.SagitalSldr,'Value'));
    
    if get(wgts.StatMapCheck,'Value')
        alpha = str2double(get(wgts.AlphaEdt,'string'));
        statmat = (STATS.p<alpha).*STATS.dat;
        anaimg = squeeze(ANA.dat(iX,:,:))';
        statimg = squeeze(statmat(iX,:,:))';
        imgrgb = subFuseImage(anaimg,statimg,MINV,MAXV,STATMINV,STATMAXV,...
                            gammav,gamma_stat,cmapstr);
        set(GRAHANDLE.sagital,'Cdata',imgrgb);
    else
        set(GRAHANDLE.sagital,'cdata',squeeze(ANA.dat(iX,:,:))');
    end
    set(GRAHANDLE.coronalV,   'xdata',[iX iX]);
    set(GRAHANDLE.transverseV,'xdata',[iX iX]);
    set(wgts.SagitalEdt,'String',sprintf('%d',iX));
    xdata = get(GRAHANDLE.triSagital,'xdata');
    xdata(:) = iX;
    set(GRAHANDLE.triSagital,'xdata',xdata,'cdata',squeeze(ANA.dat(iX,:,:)));
    axes(wgts.SagitalAxs); 
    hold on; 
    
    if get(wgts.ShowAtlasCheck,'Value')
        RoiMatrix = plotRoiMatrix(RoiMatrix, iX, 'sagital','on', mfilename);
    elseif get(wgts.StatMapCheck,'Value')
        RoiMatrix = plotRoiMatrix(RoiMatrix, iX, 'sagital','off', mfilename);
    end
    hold off;    
    setappdata(wgts.main,'RoiMatrix',RoiMatrix);
    
  end
  
  
 case {'slider-coronal'}
  GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
  TXTHANDLE = getappdata(wgts.main,'TXTHANDLE');
  RoiMatrix = getappdata(wgts.main,'RoiMatrix');
  showTextRoiMatrix('coronal', 0, TXTHANDLE, 0);
  gammav = str2double(get(wgts.GammaEdt,'String'));
  gamma_stat = str2double(get(wgts.StatGammaEdt,'String'));
  cmapstr = get(wgts.StatColormapCmb,'String');
  cmapstr = cmapstr{get(wgts.StatColormapCmb,'value')};
  
  if ~isempty(GRAHANDLE)
    iY = round(get(wgts.CoronalSldr,'Value'));
    
    if get(wgts.StatMapCheck,'Value')
        alpha = str2double(get(wgts.AlphaEdt,'string'));
        statmat = (STATS.p<alpha).*STATS.dat;
        anaimg = squeeze(ANA.dat(:,iY,:))';
        statimg = squeeze(statmat(:,iY,:))';
        imgrgb = subFuseImage(anaimg,statimg,MINV,MAXV,STATMINV,STATMAXV,...
                            gammav,gamma_stat,cmapstr);
        set(GRAHANDLE.coronal,'Cdata',imgrgb);
    else
        set(GRAHANDLE.coronal,'cdata',squeeze(ANA.dat(:,iY,:))');
    end
    set(GRAHANDLE.sagitalV,   'xdata',[iY iY]);
    set(GRAHANDLE.transverseH,'ydata',[iY iY]);
    set(wgts.CoronalEdt,'String',sprintf('%d',iY));
    ydata = get(GRAHANDLE.triCoronal,'ydata');
    ydata(:) = iY;
    set(GRAHANDLE.triCoronal,'ydata',ydata,'cdata',squeeze(ANA.dat(:,iY,:)));
    axes(wgts.CoronalAxs); 
    hold on; 

    if get(wgts.ShowAtlasCheck,'Value')
        RoiMatrix = plotRoiMatrix(RoiMatrix, iY, 'coronal','on', mfilename);
    elseif get(wgts.StatMapCheck,'Value')
        RoiMatrix = plotRoiMatrix(RoiMatrix, iY, 'coronal','off', mfilename);
    end
    hold off;
    setappdata(wgts.main,'RoiMatrix',RoiMatrix);
    
    
  end
  
 case {'slider-transverse'}
  GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
  TXTHANDLE = getappdata(wgts.main,'TXTHANDLE');
  RoiMatrix = getappdata(wgts.main,'RoiMatrix');
  showTextRoiMatrix('transverse', 0, TXTHANDLE, 0);
  gammav = str2double(get(wgts.GammaEdt,'String'));
  gamma_stat = str2double(get(wgts.StatGammaEdt,'String'));
  cmapstr = get(wgts.StatColormapCmb,'String');
  cmapstr = cmapstr{get(wgts.StatColormapCmb,'value')};
  
  if ~isempty(GRAHANDLE) 
    iZ = round(get(wgts.TransverseSldr,'Value'));
    
    if get(wgts.StatMapCheck,'Value')
        alpha = str2double(get(wgts.AlphaEdt,'string'));
        statmat = (STATS.p<alpha).*STATS.dat;
        anaimg = (ANA.dat(:,:,iZ))';
        statimg = (statmat(:,:,iZ))';
        imgrgb = subFuseImage(anaimg,statimg,MINV,MAXV,STATMINV,STATMAXV,...
                            gammav,gamma_stat,cmapstr);        
        set(GRAHANDLE.transverse,'Cdata',imgrgb);
    else
        set(GRAHANDLE.transverse,'cdata',squeeze(ANA.dat(:,:,iZ))');
    end
    set(GRAHANDLE.sagitalH,   'ydata',[iZ iZ]);
    set(GRAHANDLE.coronalH,   'ydata',[iZ iZ]);
    set(wgts.TransverseEdt,'String',sprintf('%d',iZ));
    zdata = get(GRAHANDLE.triTransverse,'zdata');
    zdata(:) = iZ;
    set(GRAHANDLE.triTransverse,'zdata',zdata,'cdata',squeeze(ANA.dat(:,:,iZ))');
    axes(wgts.TransverseAxs); 
    hold on;
      
    if get(wgts.ShowAtlasCheck,'Value')
        RoiMatrix = plotRoiMatrix(RoiMatrix, iZ, 'transverse','on', mfilename);
    else
        RoiMatrix = plotRoiMatrix(RoiMatrix, iZ, 'transverse','off', mfilename);
    end
    hold off;
    setappdata(wgts.main,'RoiMatrix',RoiMatrix);
    
  end

    
 case {'edit-sagital'}
  iX = str2double(get(wgts.SagitalEdt,'String'));
  if isempty(iX),
    iX = round(get(wgts.SagitalSldr,'Value'));
    set(wgts.SagitalEdt,'String',sprintf('%d',iX));
  else
    if iX < 0,
      iX = 1; 
      set(wgts.SagitalEdt,'String',sprintf('%d',iX));
    elseif iX > size(ANA.dat,1),
      iX = size(ANA.dat,1);
      set(wgts.SagitalEdt,'String',sprintf('%d',iX));
    end
    set(wgts.SagitalSldr,'Value',iX);
    OrthoView_Callback(hObject,'slider-sagital',[]);
  end
  
 case {'edit-coronal'}
  iY = str2double(get(wgts.CoronalEdt,'String'));
  if isempty(iY),
    iY = round(get(wgts.CoronalSldr,'Value'));
    set(wgts.CoronalEdt,'String',sprintf('%d',iY));
  else
    if iY < 0,
      iY = 1; 
      set(wgts.CoronalEdt,'String',sprintf('%d',iY));
    elseif iY > size(ANA.dat,1),
      iY = size(ANA.dat,1);
      set(wgts.CoronalEdt,'String',sprintf('%d',iY));
    end
    set(wgts.CoronalSldr,'Value',iY);
    OrthoView_Callback(hObject,'slider-coronal',[]);
  end
 
 case {'edit-transverse'}
  iZ = str2double(get(wgts.TransverseEdt,'String'));
  if isempty(iZ),
    iZ = round(get(wgts.TransverseSldr,'Value'));
    set(wgts.TransverseEdt,'String',sprintf('%d',iZ));
  else
    if iZ < 0,
      iZ = 1; 
      set(wgts.TransverseEdt,'String',sprintf('%d',iZ));
    elseif iZ > size(ANA.dat,1),
      iZ = size(ANA.dat,1);
      set(wgts.TransverseEdt,'String',sprintf('%d',iZ));
    end
    set(wgts.TransverseSldr,'Value',iZ);
    OrthoView_Callback(hObject,'slider-transverse',[]);
  end

 case {'set-gamma'}
  gammav = str2double(get(wgts.GammaEdt,'String'));
  if ~isempty(gammav),
    cmap = gray(256).^(1/gammav);
    axes(wgts.SagitalAxs);     colormap(cmap);
    axes(wgts.CoronalAxs);     colormap(cmap);
    axes(wgts.TransverseAxs);  colormap(cmap);
    axes(wgts.ColorbarAxs);    colormap(cmap);
  end
  
 case {'dir-reverse'}
  % note that image(),imagesc() reverse Y axies
  Xrev = get(wgts.XReverseCheck,'Value');
  Yrev = get(wgts.YReverseCheck,'Value');
  Zrev = get(wgts.ZReverseCheck,'Value');
  if Xrev == 0,
    corX = 'normal';   traX = 'normal';
  else
    corX = 'reverse';  traX = 'reverse';
  end
  if Yrev == 0,
    sagX = 'normal';   traY = 'reverse';
  else
    sagX = 'reverse';  traY = 'normal';
  end
  if Zrev == 0,
    sagY = 'reverse';  corY = 'reverse';
  else
    sagY = 'normal';   corY = 'normal';
  end
  set(wgts.SagitalAxs,   'xdir',sagX,'ydir',sagY);
  set(wgts.CoronalAxs,   'xdir',corX,'ydir',corY);
  set(wgts.TransverseAxs,'xdir',traX,'ydir',traY);

 case {'crosshair'}
  GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
  if ~isempty(GRAHANDLE),
    if get(wgts.CrosshairCheck,'value') == 0,
      set(GRAHANDLE.sagitalV,   'visible','off');
      set(GRAHANDLE.sagitalH,   'visible','off');
      set(GRAHANDLE.coronalV,   'visible','off');
      set(GRAHANDLE.coronalH,   'visible','off');
      set(GRAHANDLE.transverseV,'visible','off');
      set(GRAHANDLE.transverseH,'visible','off');
    else
      set(GRAHANDLE.sagitalV,   'visible','on');
      set(GRAHANDLE.sagitalH,   'visible','on');
      set(GRAHANDLE.coronalV,   'visible','on');
      set(GRAHANDLE.coronalH,   'visible','on');
      set(GRAHANDLE.transverseV,'visible','on');
      set(GRAHANDLE.transverseH,'visible','on');
    end
  end

 case {'showmri'}
  GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
  if ~isempty(GRAHANDLE),
    if get(wgts.ShowMRICheck,'value') == 0,
      set(GRAHANDLE.triTransverse,   'visible','off');
      set(GRAHANDLE.triCoronal,      'visible','off');
      set(GRAHANDLE.triSagital,      'visible','off');
    else
      set(GRAHANDLE.triTransverse,   'visible','on');
      set(GRAHANDLE.triCoronal,      'visible','on');
      set(GRAHANDLE.triSagital,      'visible','on');
    end
  end
  
 case {'showatlas'}
  RoiMatrix = getappdata(wgts.main,'RoiMatrix');
  if get(wgts.ShowAtlasCheck,'value')
      if get(wgts.StatMapCheck,'value')
          set(wgts.StatMapCheck,'Value',0)
          OrthoView_Callback(hObject,'showstatmaps',[]);
      elseif ~isempty(RoiMatrix),
          activeRoi = RoiMatrix(:,cell2mat(RoiMatrix(6,:)));
          for k = 1:size(activeRoi,2)
            if activeRoi{6,k}
                set(activeRoi{8,k},   'visible','on');
                set(activeRoi{9,k},   'visible','on');
                set(activeRoi{10,k},   'visible','on');
                set(activeRoi{11,k},   'visible','on');
            end
          end
      end
  else
      if ~isempty(RoiMatrix),
          activeRoi = RoiMatrix(:,cell2mat(RoiMatrix(6,:)));
          for k = 1:size(activeRoi,2)
            if activeRoi{6,k}
                set(activeRoi{8,k},   'visible','off');
                set(activeRoi{9,k},   'visible','off');
                set(activeRoi{10,k},   'visible','off');
                set(activeRoi{11,k},   'visible','off');
            end
          end
      end
  end
  
 case{'showstatmaps'}
    GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');
    RoiMatrix = getappdata(wgts.main,'RoiMatrix');
    if get(wgts.StatMapCheck,'Value')
        set(wgts.ShowAtlasCheck,'Value',0);
        OrthoView_Callback(hObject,'showatlas',[]);

        alpha = str2double(get(wgts.AlphaEdt,'string'));
        statmat = (STATS.p<alpha).*STATS.dat;
        cmapstr = get(wgts.StatColormapCmb,'String');
        cmapstr = cmapstr{get(wgts.StatColormapCmb,'value')};

        axes(wgts.TriplotAxs); 
        if isfield(GRAHANDLE,'triStat')
            set(GRAHANDLE.triStat,'visible','on');
        else
        GRAHANDLE.triStat = statmap3d(statmat, 0, 'on', cmapstr, mfilename);
        setappdata(wgts.main,'GRAHANDLE',GRAHANDLE)
        end
    elseif isfield(GRAHANDLE,'triStat')
        set(GRAHANDLE.triStat,'visible','off');
        if ~isempty(RoiMatrix)&&get(wgts.ShowAtlasCheck,'Value'),
          activeRoi = RoiMatrix(:,cell2mat(RoiMatrix(6,:)));
          for k = 1:size(activeRoi,2)
            if activeRoi{6,k}
                set(activeRoi{8,k},   'visible','on');
            end
          end
        end
    end
    OrthoView_Callback(hObject,'slider-sagital',[]);
    OrthoView_Callback(hObject,'slider-coronal',[]);
    OrthoView_Callback(hObject,'slider-transverse',[]);
  
 case {'button-sagital'}
    click = get(wgts.main,'SelectionType');
    if strcmpi(click,'alt') && get(wgts.CrosshairCheck,'Value') == 1,
        pt = round(get(wgts.SagitalAxs,'CurrentPoint'));
        iY = pt(1,1);  iZ = pt(1,2);
        if iY > 0 && iY <= size(ANA.dat,2),
            set(wgts.CoronalEdt,'String',sprintf('%d',iY));
            set(wgts.CoronalSldr,'Value',iY);
            OrthoView_Callback(hObject,'slider-coronal',[]);
        end
        if iZ > 0 && iZ <= size(ANA.dat,3),
          set(wgts.TransverseEdt,'String',sprintf('%d',iZ));
          set(wgts.TransverseSldr,'Value',iZ);
          OrthoView_Callback(hObject,'slider-transverse',[]);
        end
    elseif strcmpi(click,'open')
        iX = round(get(wgts.SagitalSldr,'Value'));
        src = wgts.SagitalAxs;
        hfig = figure('Position',[322    65   652   655], 'Color', get(wgts.main,'color'));
        ha = axes;
        copyobj(allchild(src),ha);
        axis image;
        set(ha,'xlim',get(src,'xlim'),'xdir',get(src,'xdir'),'xcolor',get(src,'xcolor'),...
               'ylim',get(src,'ylim'),'ydir',get(src,'ydir'),'ycolor',get(src,'ycolor'),...
               'clim',get(src,'clim'),'color',get(src,'color'),...
               'fontname',get(src,'fontname'),'fontsize',get(src,'fontsize'),...
               'fontweight',get(src,'fontweight'));
        title(sprintf('Sagital(Y-Z): %d',iX),'fontweight','bold','fontsize',10);
        set(hfig,'colormap',get(wgts.main,'colormap'));
        hc = colorbar;
        set(hc,'xcolor',get(src,'xcolor'),'ycolor',get(src,'ycolor'),...
               'fontname',get(src,'fontname'),'fontsize',get(src,'fontsize'),...
               'fontweight',get(src,'fontweight'));
        set(findobj(ha),'ButtonDownFcn','');
    elseif strcmpi(click,'normal')
        TXTHANDLE = getappdata(wgts.main,'TXTHANDLE');
        pt = round(get(wgts.SagitalAxs,'CurrentPoint'));
        RoiMatrix = getappdata(wgts.main,'RoiMatrix');
        showTextRoiMatrix(RoiMatrix(:,cell2mat(RoiMatrix(6,:))), pt, TXTHANDLE, 1);
%     setappdata(wgts.main,'Roi',Roi);
    end
  
 case {'button-coronal'}
    click = get(wgts.main,'SelectionType');
    if strcmpi(click,'alt') && get(wgts.CrosshairCheck,'Value') == 1,
        pt = round(get(wgts.CoronalAxs,'CurrentPoint'));
        iX = pt(1,1);  iZ = pt(1,2);
        if iX > 0 && iX <= size(ANA.dat,1),
          set(wgts.SagitalEdt,'String',sprintf('%d',iX));
          set(wgts.SagitalSldr,'Value',iX);
          OrthoView_Callback(hObject,'slider-sagital',[]);
        end
        if iZ > 0 && iZ <= size(ANA.dat,3),
          set(wgts.TransverseEdt,'String',sprintf('%d',iZ));
          set(wgts.TransverseSldr,'Value',iZ);
          OrthoView_Callback(hObject,'slider-transverse',[]);
        end
    elseif strcmpi(click,'open')
        iY = round(get(wgts.CoronalSldr,'Value'));
        src = wgts.CoronalAxs;
        hfig = figure('Position',[322    65   652   655], 'Color', get(wgts.main,'color'));
        ha = axes;
        copyobj(allchild(src),ha);
        axis image;
        set(ha,'xlim',get(src,'xlim'),'xdir',get(src,'xdir'),'xcolor',get(src,'xcolor'),...
               'ylim',get(src,'ylim'),'ydir',get(src,'ydir'),'ycolor',get(src,'ycolor'),...
               'clim',get(src,'clim'),'color',get(src,'color'),...
               'fontname',get(src,'fontname'),'fontsize',get(src,'fontsize'),...
               'fontweight',get(src,'fontweight'));
        title(sprintf('Coronal(X-Z): %d',iY),'fontweight','bold','fontsize',10);
        set(hfig,'colormap',get(wgts.main,'colormap'));
        hc = colorbar;
        set(hc,'xcolor',get(src,'xcolor'),'ycolor',get(src,'ycolor'),...
               'fontname',get(src,'fontname'),'fontsize',get(src,'fontsize'),...
               'fontweight',get(src,'fontweight'));
        set(findobj(ha),'ButtonDownFcn','');
    elseif strcmpi(click,'normal')
        TXTHANDLE = getappdata(wgts.main,'TXTHANDLE');
        pt = round(get(wgts.CoronalAxs,'CurrentPoint'));
        RoiMatrix = getappdata(wgts.main,'RoiMatrix');
        showTextRoiMatrix(RoiMatrix(:,cell2mat(RoiMatrix(6,:))), pt, TXTHANDLE, 1);
    %     setappdata(wgts.main,'Roi',Roi);
    end

 case {'button-transverse'}
    click = get(wgts.main,'SelectionType');
    if strcmpi(click,'alt') && get(wgts.CrosshairCheck,'Value') == 1,
        pt = round(get(wgts.TransverseAxs,'CurrentPoint'));
        iX = pt(1,1);  iY = pt(1,2);
        if iX > 0 && iX <= size(ANA.dat,1),
          set(wgts.SagitalEdt,'String',sprintf('%d',iX));
          set(wgts.SagitalSldr,'Value',iX);
          OrthoView_Callback(hObject,'slider-sagital',[]);
        end
        if iY > 0 && iY <= size(ANA.dat,2),
          set(wgts.CoronalEdt,'String',sprintf('%d',iY));
          set(wgts.CoronalSldr,'Value',iY);
          OrthoView_Callback(hObject,'slider-coronal',[]);
        end
    elseif strcmpi(click,'open')
        iZ = round(get(wgts.TransverseSldr,'Value'));
        src = wgts.TransverseAxs;
        hfig = figure('Position',[322    65   652   655], 'Color', get(wgts.main,'color'));
        ha = axes;
        copyobj(allchild(src),ha);
        axis image;
        set(ha,'xlim',get(src,'xlim'),'xdir',get(src,'xdir'),'xcolor',get(src,'xcolor'),...
               'ylim',get(src,'ylim'),'ydir',get(src,'ydir'),'ycolor',get(src,'ycolor'),...
               'clim',get(src,'clim'),'color',get(src,'color'),...
               'fontname',get(src,'fontname'),'fontsize',get(src,'fontsize'),...
               'fontweight',get(src,'fontweight'));
        title(sprintf('Transverse(X-Y): %d',iZ),'fontweight','bold','fontsize',10);
        set(hfig,'colormap',get(wgts.main,'colormap'));
        hc = colorbar;
        set(hc,'xcolor',get(src,'xcolor'),'ycolor',get(src,'ycolor'),...
               'fontname',get(src,'fontname'),'fontsize',get(src,'fontsize'),...
               'fontweight',get(src,'fontweight'));
        set(findobj(ha),'ButtonDownFcn','');
    elseif strcmpi(click,'normal')
        TXTHANDLE = getappdata(wgts.main,'TXTHANDLE');
        pt = round(get(wgts.TransverseAxs,'CurrentPoint'));
        RoiMatrix = getappdata(wgts.main,'RoiMatrix');
        showTextRoiMatrix(RoiMatrix(:,cell2mat(RoiMatrix(6,:))), pt, TXTHANDLE, 1);
    %     setappdata(wgts.main,'Roi',Roi);
    end
  
 case {'button-triplot'}
  click = get(wgts.main,'SelectionType');
  if strcmpi(click,'open'),
    set(wgts.main,'Pointer','watch')
    iX = round(get(wgts.SagitalSldr,'Value'));
    iY = round(get(wgts.CoronalSldr,'Value'));
    iZ = round(get(wgts.TransverseSldr,'Value'));
    src = wgts.TriplotAxs;
    hfig = figure('Position',[337    64   611   657], 'Color', get(wgts.main,'color'));
    ha = axes; 
    
    if get(wgts.ShowAtlasCheck,'value')
        RoiMatrix = getappdata(wgts.main,'RoiMatrix');
        indRoi = cell2mat(RoiMatrix(6,:));
        roi3dMatrix(RoiMatrix(:,indRoi), 1, 'on','roi', mfilename);
    elseif get(wgts.StatMapCheck,'value')
        alpha = str2double(get(wgts.AlphaEdt,'string'));
        statmat = (STATS.p<alpha).*STATS.dat;
        cmapstr = get(wgts.StatColormapCmb,'String');
        cmapstr = cmapstr{get(wgts.StatColormapCmb,'value')};
        hpatch = statmap3d(statmat, 1, 'on', cmapstr,mfilename);
        set(hpatch,'ambientstrength',0.2);
    end
    set(ha,'CameraViewAngle',9.2046);
    tmpview = get(src,'view');
    view(tmpview);  
    
    hsurfs = findobj(allchild(src),'Type','surface'); %3 surfaces, the 3 planes
    newsurfs = copyobj(hsurfs,ha);
    set(newsurfs,'FaceLighting','none');
    set(ha,'xlim',get(src,'xlim'),'xdir',get(src,'xdir'),'xcolor',get(src,'xcolor'),...
           'ylim',get(src,'ylim'),'ydir',get(src,'ydir'),'ycolor',get(src,'ycolor'),...
           'clim',get(src,'clim'),'zdir',get(src,'zdir'),'zcolor',get(src,'zcolor'),...
           'color',get(src,'color'),'fontname',get(src,'fontname'),...
           'fontsize',get(src,'fontsize'),'fontweight',get(src,'fontweight'));
    
    xlabel(get(get(src,'xlabel'),'string'));
    ylabel(get(get(src,'ylabel'),'string'));
    zlabel(get(get(src,'zlabel'),'string'));
    title(sprintf('Triplot (X,Y,Z)=(%d,%d,%d)',iX,iY,iZ),'fontweight','bold','fontsize',10);
    set(hfig,'colormap',get(wgts.main,'colormap'));
    grid on; axis equal vis3d tight;
    set(wgts.main,'Pointer','arrow')
  end
  
 otherwise
     
end


return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to handle lightbox view
function LightboxView_Callback(hObject,eventdata,handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wgts = guihandles(get(hObject,'Parent'));
ANA  = getappdata(wgts.main,'ANA');
MINV    = getappdata(wgts.main,'MINV');
MAXV    = getappdata(wgts.main,'MAXV');
CMAP    = getappdata(wgts.main,'CMAP');
ViewMode = get(wgts.ViewModeCmb,'String');
ViewMode = ViewMode{get(wgts.ViewModeCmb,'Value')};
switch lower(ViewMode),
 case {'lightbox-coronal'}
  iDimension = 2;
 case {'lightbox-sagittal'}
  iDimension = 1;
 case {'lightbox-transverse'}
  iDimension = 3;
 otherwise
  iDimension = 3;
end
nmaximages = size(ANA.dat,iDimension);

NCol = 5;
NRow = 4;
% NCol = 3;
% NRow = 3;

switch lower(eventdata),
 case {'init'}
  NPages = floor((nmaximages-1)/NCol/NRow)+1;
  tmptxt = {};
  for iPage = 1:NPages,
    tmptxt{iPage} = sprintf('Page%d: %d-%d',iPage,...
                            (iPage-1)*NCol*NRow+1,min([nmaximages,iPage*NCol*NRow]));
  end
  set(wgts.ViewPageList,'String',tmptxt,'Value',1);
  ViewMode = get(wgts.ViewModeCmb,'String');
  ViewMode = ViewMode{get(wgts.ViewModeCmb,'Value')};
  if strcmpi(ViewMode,'lightbox'),
    LightboxView_Callback(hObject,'redraw',handles);
  end
  
 case {'redraw'}
  axes(wgts.LightboxAxs);  cla;
  pagestr = get(wgts.ViewPageList,'String');
  pagestr = pagestr{get(wgts.ViewPageList,'Value')};
  ipage = sscanf(pagestr,'Page%d:');
  SLICES = (ipage-1)*NCol*NRow+1:min([nmaximages,ipage*NCol*NRow]);
  RoiMatrix = getappdata(wgts.main,'RoiMatrix');
  if iDimension == 1,
    nX = size(ANA.dat,2);  nY = size(ANA.dat,3);
    INFSTR = 'Sag';
    Xrev = get(wgts.YReverseCheck,'Value');
    Yrev = get(wgts.ZReverseCheck,'Value');
  elseif iDimension == 2,
    nX = size(ANA.dat,1);  nY = size(ANA.dat,3);
    INFSTR = 'Cor';
    Xrev = get(wgts.XReverseCheck,'Value');
    Yrev = get(wgts.ZReverseCheck,'Value');
  else
    nX = size(ANA.dat,1);  nY = size(ANA.dat,2);
    INFSTR = 'Trans';
    Xrev = get(wgts.XReverseCheck,'Value');
    Yrev = get(wgts.YReverseCheck,'Value');
  end
  X = [0:nX-1];  Y = [nY-1:-1:0];
  if Xrev > 0,  X = fliplr(X);  end
  if Yrev > 0,  Y = fliplr(Y);  end
  for N = 1:length(SLICES),
    iSlice = SLICES(N);
    if iDimension == 1,
      tmpimg = squeeze(ANA.dat(iSlice,:,:));
    elseif iDimension == 2,
      tmpimg = squeeze(ANA.dat(:,iSlice,:));
    else
      tmpimg = squeeze(ANA.dat(:,:,iSlice));
    end
    iCol = floor((N-1)/NRow)+1;
    iRow = mod((N-1),NRow)+1;
    offsX = nX*(iRow-1);
    offsY = nY*NCol - iCol*nY;
    tmpx = X + offsX;  tmpy = Y + offsY;
    imagesc(tmpx,tmpy,tmpimg');  hold on;
    plotRoiLightboxM(RoiMatrix(:,cell2mat(RoiMatrix(6,:))), iSlice, INFSTR, offsX, offsY, Xrev, Yrev);
    text(min(tmpx)+1,min(tmpy)+1,sprintf('%s=%d',INFSTR,iSlice),...
         'color',[0.9 0.9 0.5],'VerticalAlignment','bottom',...
         'FontName','Comic Sans MS','FontSize',8,'Fontweight','bold');
  end
  axes(wgts.LightboxAxs);  
%   colormap(gray(256));
  set(gca,'Tag','LightboxAxs','color','black');
  set(gca,'XTickLabel',{},'YTickLabel',{},'XTick',[],'YTick',[]);
  set(gca,'xlim',[0 nX*NRow],'ylim',[0 nY*NCol],'clim',[MINV MAXV]);
  set(gca,'YDir','normal');
  set(allchild(gca),...
      'ButtonDownFcn','atlas3dgui(''LightboxView_Callback'',gcbo,''button-lightbox'',guidata(gcbo))');

  
 case {'button-lightbox'}
  click = get(wgts.main,'SelectionType');
  if strcmpi(click,'open'),
    % double click
    subZoomInLightBox(wgts,ANA);
  end
  
  
 otherwise
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to zoom-in plot
function subZoomInLightBox(wgts,ANA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ViewMode = get(wgts.ViewModeCmb,'String');
ViewMode = ViewMode{get(wgts.ViewModeCmb,'Value')};
ANA.ds = [0.25 0.25 0.25];
switch lower(ViewMode),
 case {'lightbox-coronal'}
  DX = ANA.ds(1);  DY = ANA.ds(3);
  tmpstr = sprintf('CORONAL %s %s',ANA.session,ANA.grpname);
  tmpxlabel = 'X (mm)';  tmpylabel = 'Z (mm)';
 case {'lightbox-sagittal'}
  DX = ANA.ds(2);  DY = ANA.ds(3);
  tmpstr = sprintf('SAGITAL %s %s',ANA.session,ANA.grpname);
  tmpxlabel = 'Y (mm)';  tmpylabel = 'Z (mm)';
 case {'lightbox-transverse'}
  DX = ANA.ds(1);  DY = ANA.ds(2);
  tmpstr = sprintf('TRANSVERSE %s %s',ANA.session,ANA.grpname);
  tmpxlabel = 'X (mm)';  tmpylabel = 'Y (mm)';
end


tmpstr = sprintf('%s P<%s',tmpstr,get(wgts.AlphaEdt,'String'));


hfig = wgts.main + 1005;
hsrc = wgts.LightboxAxs;


figure(hfig);  clf;
pos = get(hfig,'pos');
pos = [pos(1)-(680-pos(3)) pos(2)-(500-pos(4)) 680 500];

set(hfig,'Name',tmpstr,'pos',pos,'color',get(wgts.main,'color'));
haxs = copyobj(hsrc,hfig);
set(haxs,'ButtonDownFcn','');  % clear callback function
set(hfig,'Colormap',get(wgts.main,'Colormap'));
h = findobj(haxs,'type','image');
set(h,'ButtonDownFcn','');  % clear callback function

for N = 1:length(h),
  set(h(N),'xdata',get(h(N),'xdata')*DX,'ydata',get(h(N),'ydata')*DY);
  nx = length(get(h(N),'xdata'));  ny = length(get(h(N),'ydata'));
end
h = findobj(haxs,'type','text');
for N = 1:length(h),
  tmppos = get(h(N),'pos');
  tmppos(1) = tmppos(1)*DX;  tmppos(2) = tmppos(2)*DY;
  set(h(N),'pos',tmppos);
end
set(haxs,'Position',[0.08 0.1 0.75 0.75],'units','normalized');
h = findobj(haxs,'type','line');
for N =1:length(h),
  set(h(N),'xdata',get(h(N),'xdata')*DX,'ydata',get(h(N),'ydata')*DY);
end

h = findobj(haxs,'tag','ScaleBar');
if ~isempty(h),
  %if length(h) > 1,
  %  delete(h(1:end-1));
  %end
  %h = h(end);
  for N = 1:length(h),
    tmppos = get(h(N),'pos');
    tmppos([1 3]) = tmppos([1 3])*DX;
    tmppos([2 4]) = tmppos([2 4])*DY;
    set(h(N),'pos',tmppos);
  end
end


set(haxs,'xlim',get(haxs,'xlim')*DX,'ylim',get(haxs,'ylim')*DY);
set(haxs,'xtick',[0 10 20 30 40 50 60 70 80 90 100 110 120]);
set(haxs,'ytick',[0 10 20 30 40 50 60 70 80 90 100 110 120]);
xlabel(tmpxlabel);  ylabel(tmpylabel);
title(haxs,strrep(tmpstr,'_','\_'));
daspect(haxs,[1 1 1]);
pos = get(haxs,'pos');
hbar = copyobj(wgts.ColorbarAxs,hfig);
set(hbar,'pos',[0.85 pos(2) 0.045 pos(4)]);    
  

%clear callbacks
set(haxs,'ButtonDownFcn','');
set(allchild(haxs),'ButtonDownFcn','');


% make font size bigger
set(haxs,'FontSize',10);
set(get(haxs,'title'),'FontSize',10);
set(get(haxs,'xlabel'),'FontSize',10);
set(get(haxs,'ylabel'),'FontSize',10);
set(hbar,'FontSize',10);
set(get(hbar,'ylabel'),'FontSize',10);


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fuse anatomy and functional images
function imgrgb = subFuseImage(anaimg,statimg,MINV,MAXV,MIN_STAT,MAX_STAT,GAMMA_ana,GAMMA_stat,CMAP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GAMMA_ana = 1;
% GAMMA_stat = 1;
% colormap = 'autumn';

anaimg = (anaimg - MINV) / (MAXV - MINV);
anaimg = round(anaimg*255) + 1; % +1 for matlab indexing
anaimg((anaimg(:) <   0)) =   1;
anaimg((anaimg(:) > 256)) = 256;
ana_cmap = gray(256).^(1/GAMMA_ana);

ANARGB = ind2rgb(anaimg,ana_cmap);

mask = (statimg~=0);
mask = repmat(mask,[1 1 3]);

statimg = (statimg - MIN_STAT) / (MAX_STAT - MIN_STAT);
statimg = round(statimg*255) + 1; % +1 for matlab indexing
statimg((statimg(:) <   0)) =   1;
statimg((statimg(:) > 256)) = 256;
stat_cmap = eval(sprintf('%s(%g)',CMAP,256));
stat_cmap = stat_cmap.^(1/GAMMA_stat);

STATRGB = ind2rgb(statimg,stat_cmap);

imgrgb = ANARGB;
imgrgb(mask) = STATRGB(mask);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to redraw crosshair over rois
function subRedrawCrosshair(wgts,ANA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nX = size(ANA.dat,1);  nY = size(ANA.dat,2);  nZ = size(ANA.dat,3); 

iX = get(wgts.SagitalSldr,'Value');
iY = get(wgts.CoronalSldr,'Value');
iZ = get(wgts.TransverseSldr,'Value');
GRAHANDLE = getappdata(wgts.main,'GRAHANDLE');

delete([GRAHANDLE.transverseH, GRAHANDLE.transverseV, GRAHANDLE.coronalV,...
      GRAHANDLE.coronalH, GRAHANDLE.sagitalH, GRAHANDLE.sagitalV]);
  
% draw crosshair(s)
axes(wgts.SagitalAxs);
  hSagV = line([iY iY],[ 1 nZ],'color','b','tag','line');
  hSagH = line([ 1 nY],[iZ iZ],'color','b','tag','line');
  set([hSagV hSagH],...
      'ButtonDownFcn','atlas3dgui(''OrthoView_Callback'',gcbo,''button-sagital'',guidata(gcbo))');
  axes(wgts.CoronalAxs);
  hCorV = line([iX iX],[ 1 nZ],'color','b','tag','line');
  hCorH = line([ 1 nX],[iZ iZ],'color','b','tag','line');
  set([hCorV hCorH],...
      'ButtonDownFcn','atlas3dgui(''OrthoView_Callback'',gcbo,''button-coronal'',guidata(gcbo))');
  axes(wgts.TransverseAxs);
  hTraV = line([iX iX],[ 1 nY],'color','b','tag','line');
  hTraH = line([ 1 nX],[iY iY],'color','b','tag','line');
  set([hTraV hTraH],...
      'ButtonDownFcn','atlas3dgui(''OrthoView_Callback'',gcbo,''button-transverse'',guidata(gcbo))');
  if get(wgts.CrosshairCheck,'Value') == 0,
    set([hSagV hSagH hCorV hCorH hTraV hTraH],'visible','off');
  end
  GRAHANDLE.sagitalV    = hSagV;
  GRAHANDLE.sagitalH    = hSagH;
  GRAHANDLE.coronalV    = hCorV;
  GRAHANDLE.coronalH    = hCorH;
  GRAHANDLE.transverseV = hTraV;
  GRAHANDLE.transverseH = hTraH;
  
  setappdata(wgts.main,'GRAHANDLE',GRAHANDLE);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to get a color map in  rgb
function cmap_rgb = subGetColormap(wgts,depth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% depth = 20;
ydat = 1:depth;
ydat = ydat';

cmapstr = get(wgts.StatColormapCmb,'String');
cmapstr = cmapstr{get(wgts.StatColormapCmb,'value')};
switch lower(cmapstr),
 case {'autumn','spring','winter','cool','jet','hot','hsv','copper'}
  eval(sprintf('cmap = %s(%g);',cmapstr,depth));
 case {'red'}
  cmap = zeros(depth,3);  cmap(:,1) = 1;
 case {'green'}
  cmap = zeros(depth,3);  cmap(:,2) = 1;
 case {'blue'}
  cmap = zeros(depth,3);  cmap(:,3) = 1;
 case {'yellow'}
  cmap = zeros(depth,3);  cmap(:,1) = 1;  cmap(:,2) = 1;
 otherwise
  cmap = autumn(depth);
end
gammav = str2double(get(wgts.StatGammaEdt,'String'));
if ~isempty(gammav),
  cmap = cmap.^(1/gammav);
end
cmap = flipud(cmap);
cmap_rgb= ind2rgb(ydat,cmap);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION bregma coordinates
function subBregmaCoord(wgts, bregma)

CoordAxs = getappdata(wgts.main,'CoordAxs');
CorX = CoordAxs.CorX;
dif = (CorX(2)-CorX(1));
offset = rem(bregma(1),dif);
set(wgts.CoronalAxs,'XTick',([0, CorX]+offset)); %centro en bregma
CorXlabel = (CorX - bregma(1)+offset).*0.25;
CorXlabel = [(CorXlabel(1)-(CorXlabel(2)-CorXlabel(1))), CorXlabel];
set(wgts.CoronalAxs,'XTickLabel',CorXlabel);

set(wgts.TriplotAxs,'XTick',([0, CorX]+offset));%triplot
set(wgts.TriplotAxs,'XTickLabel',CorXlabel);

CorY = CoordAxs.CorY;
dif = (CorY(2)-CorY(1));
offset = rem(bregma(3),dif);
set(wgts.CoronalAxs,'YTick',([0, CorY]+offset)); %centro en bregma
CorYlabel = (CorY - bregma(3)+offset).*0.25;
CorYlabel = [(CorYlabel(1)-(CorYlabel(2)-CorYlabel(1))), CorYlabel];
set(wgts.CoronalAxs,'YTickLabel',CorYlabel)

SagX = CoordAxs.SagX;
dif = (SagX(2)-SagX(1));
offset = rem(bregma(2),dif);
set(wgts.SagitalAxs,'XTick',([0, SagX]+offset)); %centro en bregma
SagXlabel = (SagX - bregma(2)+offset).*0.25;
SagXlabel = [(SagXlabel(1)-(SagXlabel(2)-SagXlabel(1))), SagXlabel];
set(wgts.SagitalAxs,'XTickLabel',SagXlabel);

set(wgts.TriplotAxs,'YTick',([0, SagX]+offset));%triplot
set(wgts.TriplotAxs,'YTickLabel',SagXlabel);

SagY = CoordAxs.SagY;
dif = (SagY(2)-SagY(1));
offset = rem(bregma(3),dif);
set(wgts.SagitalAxs,'YTick',([0, SagY]+offset)); %centro en bregma
SagYlabel = (SagY - bregma(3)+offset).*0.25;
SagYlabel = [(SagYlabel(1)-(SagYlabel(2)-SagYlabel(1))), SagYlabel];
set(wgts.SagitalAxs,'YTickLabel',SagYlabel)

set(wgts.TriplotAxs,'ZTick',([0, SagY]+offset));%triplot
set(wgts.TriplotAxs,'ZTickLabel',SagYlabel);

TranX = CoordAxs.TranX;
dif = (TranX(2)-TranX(1));
offset = rem(bregma(1),dif);
set(wgts.TransverseAxs,'XTick',([0, TranX]+offset)); %centro en bregma
TranXlabel = (TranX - bregma(1)+offset).*0.25;
TranXlabel = [(TranXlabel(1)-(TranXlabel(2)-TranXlabel(1))), TranXlabel];
set(wgts.TransverseAxs,'XTickLabel',TranXlabel);

TranY = CoordAxs.TranY;
dif = (TranY(2)-TranY(1));
offset = rem(bregma(2),dif);
set(wgts.TransverseAxs,'YTick',([0, TranY]+offset)); %centro en bregma
TranYlabel = (TranY - bregma(2)+offset).*0.25;
TranYlabel = [(TranYlabel(1)-(TranYlabel(2)-TranYlabel(1))), TranYlabel];
set(wgts.TransverseAxs,'YTickLabel',TranYlabel)

return
