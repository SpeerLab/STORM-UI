
%fid = fopen('Z:\Colenso\05_25_18_sample5\acquisition\bins\IRbead_1_00_647mlist.hdf5','r');
%c = hread(fid,7);
%A = fread(fid,'float32');
h5disp('Z:\Colenso\05_25_18_sample5\acquisition\bins\Visbead_1_00_488mlist.hdf5')
X_value = h5read('Z:\Colenso\05_25_18_sample5\acquisition\bins\IRbead_1_00_647mlist.hdf5','/fr_24/x')
Y_value = h5read('Z:\Colenso\05_25_18_sample5\acquisition\bins\IRbead_1_00_647mlist.hdf5','/fr_24/y')

%   1     2           3        4        5       6       7
% index  background	category   error    frame #	Height	iterations
%   8               9       10          11              12
% 	significance	sum     track_ID	track_length    X
%   13           14          15          
%   xsigma       y           Z
cat = floor(A(1:end)/1.401e-45);
x = A(12:15:end);
y = A(14:end);
z = A(15:end);
h = A(6:end);
frame = A(5:end);
length = A(11:end);
valid = floor(A(13:18:end)/1.401e-45);



ind = find(length>4);%==(mode(length(length>max(length)/2))));
%ind = find(length==(mode(length(length>max(length)/2))));
N = size(ind,1);
% mol = struct('cat',{cat(ind)},'x',{x(ind)},'y',{y(ind)},'z',{z(ind)},'h',{h(ind)},'area',{area(ind)},'width',{width(ind)},'phi',{phi(ind)},         'Ax',{Ax(ind)},'bg',{bg(ind)},'I',{I(ind)},'frame',{frame(ind)},'length',{length(ind)},'link',{-1*ones(N,1)},'valid',{valid(ind)});
mol.cat = cat(ind);
mol.x = x(ind);
mol.y = y(ind);
mol.z = z(ind);
mol.h = h(ind);
mol.frame = frame(ind);
mol.length = length(ind);

clear A;

f = mol;