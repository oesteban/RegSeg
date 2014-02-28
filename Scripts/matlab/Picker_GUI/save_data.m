function [ output_args ] = save_data( path, data )
%SAVE_DATA Wrapper to save xls files if MATLAB is either in linux or
%windows
    try
        xlswrite( path, data);
    catch exception
        % http://www.mathworks.ch/support/solutions/en/data/1-1190ZB/index.html?solution=1-1190ZB
        ex2 = cellfun(@ex_func, data, 'UniformOutput', 0 );
        size_ex2 = cellfun(@length,ex2,'UniformOutput',0);
        str_length = max(max(cell2mat(size_ex2)));
        ex3 = cellfun(@(x) ex_func2(x,str_length),ex2,'uniformoutput',0);
        fid = fopen( path,'wt');
        for i = 1:size(ex3,1)
            for j = 1:size(ex3,2)
                termchar = ',';
                
                if j == size(ex3,2)
                    termchar = '\n';
                end
                    
                fprintf(fid,['%s',termchar], ex2{i,j} );
            end
        end
        fclose(fid);
    end

end

