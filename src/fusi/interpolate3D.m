% Internal function used by registrationccf.m and registerData.m
% that interpolates the data of 'scan' to the same voxel size than the 'atlas'
% it also manipulates the order of the axis to fit with the atlas

function scanInt=interpolate3D(atlas,scanNoperm)

scan=permuteScan(scanNoperm,atlas);

dz=scan.VoxelSize(1);
dx=scan.VoxelSize(2);
dy=scan.VoxelSize(3);

dzint=atlas.VoxelSize(1);
dxint=atlas.VoxelSize(2);
dyint=atlas.VoxelSize(3);

[nz,nx,ny]=size(scan.Data);

n1x=round((nx-1)*dx/dxint)+1;
n1y=round((ny-1)*dy/dyint)+1;
n1z=round((nz-1)*dz/dzint)+1;

% warning!! X and Y are permuted in matlab meshgrid axis1=Y axis2=X axis3=Z
[Xq,Yq,Zq] = meshgrid( (0:n1x-1)*dxint/dx+1,(0:n1z-1)*dzint/dz+1 ,(0:n1y-1)*dyint/dy+1);
ai=interp3(scan.Data,Xq,Yq,Zq,'linear',0);

scanInt.Data=ai;
scanInt.VoxelSize=atlas.VoxelSize;
scanInt.Direction=scan.Direction;
end
