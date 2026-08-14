% Urban Lab - NERF empowered by imec, KU Leuven and VIB
% Mace Lab  - Max Planck institute of Neurobiology
% Authors:  G. MONTALDO, E. MACE
% Review & test: C.BRUNNER, M. GRILLET
% September 2020
%
% Correlates functional scans with a square stimulus.
%
% c=mapCorrelation(scanfus,t0,t1)
%   scanfus,    fus-structure of type fusvolume 
%   t1,         start of the square correlation window
%   t2,         end of the square correlation window
%   c,          fus-structure of the same type fusvolume. 
%               It contains the correlation with the square window. 
%
% example: example03_correlation
%%
function c=mapCorrelation(scanfus,t0,t1)
[nz,nx,ny,nt]=size(scanfus.Data);

% creates a normalized square window between t0 and t1 
stim=zeros(nt,1);
stim(t0:t1,:)=1;                         % activity is a square window
hrf = hemodynamicResponse(0.1,[1.5 10 0.5 1 20 0 16]); 
stim=filter(hrf,1,stim);                 % filter the activity by the hrf

% normalize stim
stim=stim-mean(stim);
stim=stim./sqrt(sum(stim.^2));

c = struct();
c = copyStructBtoA(c, scanfus);
c.Data=zeros(nz,nx,ny);
% c.VoxelSize=scanfus.VoxelSize;
% c.Type=scanfus.Type;
% c.Direction=scanfus.Direction;
datatocorr = reshape(scanfus.Data, [], nt);
[corrvals, pvals]   = corr(datatocorr', stim);
corrvals   = reshape(corrvals, [nz, nx, ny]);
pvals      = reshape(pvals, [nz, nx, ny]);

for iplane=1:ny
    fprintf('correlation plane %d\n',iplane);
    c.Data(:,:,iplane)=medfilter(corrvals(:,:,iplane),5);
    c.pData(:,:,iplane)= pvals(:,:,iplane);
end

end

% median filter to eliminate outliers
% if a point is 1.5 times higher than the std in a square of n points
% it is set to the median in the square.
function af=medfilter(a,n)
[nz,nx]=size(a);
af=a;
for iz=1+n:nz-n
    for ix=1+n:nx-n
        tmp=a([-n:n]+iz,[-n:n]+ix);
        st=std(tmp(:));
        m=median(tmp(:));   
        if abs(af(iz,ix)-m)>1.5*st
        af(iz,ix)=m;
        end
    end
end

end


