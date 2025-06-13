function avreg = transformAnnotationVolume(transpath, av, avscale)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% we here assume origins are the same

paramsfin = elastix_parameter_read(transpath);

% adjust the file and rewrite it as a new temp file 
paramsfin.ResultImagePixelType = 'double';
paramsfin.FinalBSplineInterpolationOrder = 0;

[savepath, fname, fext] = fileparts(transpath);

new_temp_path = fullfile(savepath, sprintf('%s_temp_annotation%s', fname, fext));
elastix_paramStruct2txt(new_temp_path, paramsfin);

% call final transformix
avreg  = transformix(av,new_temp_path, 'movingscale', avscale*[1 1 1], 'verbose', 0);

delete(new_temp_path);
end