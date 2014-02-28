function select_popupmenu(fig,mainFig)


hdlEdt = findobj(get(fig,'children'),'style','edit');
list = findobj(get(fig,'children'),'tag','classeslist');

idx = get(list,'string');
n = str2double(idx(get(list,'value')));
regions = cell(1,n);
hdlEdt = flipud(hdlEdt);
for i = 1:n
    regions{i} = get(hdlEdt(i)','string');
end
npixhdl = findobj(get(fig,'children'),'tag','npix');

idx = get(npixhdl,'string');
npix = str2double(idx(get(npixhdl,'value')));
close(fig)
roi_seg_menu(mainFig,regions,npix);