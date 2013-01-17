function [ out ] = ex_func( in )

in_datatype = class(in);

switch in_datatype
    case 'char'
        out = in;
    case 'double'
        out = num2str(in);
end

