% internal function to permute the orientation of the scan 

function dataperm=permuteScan(data,atlas)

o=[data.Direction(1:3:end)' data.Direction(2:3:end)'];
f=[atlas.Direction(1:3:end)' atlas.Direction(2:3:end)'];

ndim=size(o,1);

axfind=zeros(ndim,1);
Perm=zeros(ndim,1);
Flip=zeros(ndim,1);
for i=1:ndim
    for j=1:3
        if strcmp(o(i,:),f(j,:))
            Perm(j)=i;
            axfind(j)=1;
        elseif strcmp(flip(o(i,:)),f(j,:))
            Perm(j)=i;
            Flip(i)=1;
            axfind(j)=1;
        end
    end
end

% if exist find the time axis
if ndim==4
    for i=1:ndim
        if strcmp(o(i,:),'TM')
            Perm(ndim)=i;
            axfind(4)=1;
        end
    end
end

%chk
if  prod(axfind)==0
    fprintf(' origin axis: %s\n destin axis: %s\n',axorig,axfinal);
    error('error in name of the axis')
end


tmp=data.Data;
for i=1:3
     if Flip(i)==1
        tmp=flip(tmp,i);
     end
end
dataperm.Data=permute(tmp,Perm);

dataperm.VoxelSize=data.VoxelSize(Perm(1:3));
dataperm.Direction=atlas.Direction; 
end