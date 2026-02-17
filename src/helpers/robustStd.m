function robuststd = robustStd(xx)
%NAN Summary of this function goes here
%   Detailed explanation goes here

madfac = 1.4826;
robuststd = madfac*mad(xx,1,1);

end

