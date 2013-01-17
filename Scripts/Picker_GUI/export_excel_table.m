function export_excel_table(hSelectFig)


[filename pathname] = uiputfile('*.xls','Save as');

if ~isnumeric(filename)

    fprintf('Exporting data to %s...\n',fullfile(pathname,filename))

    table = findobj(get(hSelectFig,'children'),'tag','statistics');
    columns = get(table,'ColumnName');
    stats = get(table,'Data');
    [row col] = size(stats);
    data = cell(row + 1, col);
    data(1,:) = columns;
    data(2:end,:) = stats;
    
    save_data( fullfile(pathname,filename), data);
        
    disp('done.')
end
