function cptsreg = transformControlPoints(transpath, cpts)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% we here assume origins are the same

paramsfin = elastix_parameter_read(transpath);

% adjust the file and rewrite it as a new temp file 
[savepath, fname, fext] = fileparts(transpath);

new_temp_path = fullfile(savepath, sprintf('%s_temp_cpts%s', fname, fext));
elastix_paramStruct2txt(new_temp_path, paramsfin);

% call final transformix
cptsreg  = transformix(cpts,new_temp_path);
cptsreg  = cptsreg.OutputPoint;
delete(new_temp_path);
end