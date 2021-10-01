function f = readbinfile_cmosbins(filename)

fid = fopen(filename,'r');
c = fread(fid,12);
A = fread(fid,'float32');
fclose(fid);

%   1       2   3   4   5   6       7       8       9   10  11  12
% Cas44178	X	Y	Xc	Yc	Height	Area	Width	Phi	Ax	BG	I	
%   13  14      15      16      17  18
% Frame	Length	Link	Valid	Z	Zc
cat = floor(A(13:18:end)/1.401e-45);
x = A(4:18:end);
y = A(5:18:end);
z = A(18:18:end);
h = A(6:18:end);
area = A(7:18:end);
width = A(8:18:end);
phi = A(9:18:end);
Ax = A(10:18:end);
bg = A(11:18:end);
I = A(12:18:end);
frame = floor(A(15:18:end)/1.401e-45);
length = floor(A(16:18:end)/1.401e-45);
valid = floor(A(13:18:end)/1.401e-45);


N = min([size(bg,1),size(area,1) ]);
ind = find(length>4);%==(mode(length(length>max(length)/2))));
%ind = find(length==(mode(length(length>max(length)/2))));
N = size(ind,1);
% mol = struct('cat',{cat(ind)},'x',{x(ind)},'y',{y(ind)},'z',{z(ind)},'h',{h(ind)},'area',{area(ind)},'width',{width(ind)},'phi',{phi(ind)},         'Ax',{Ax(ind)},'bg',{bg(ind)},'I',{I(ind)},'frame',{frame(ind)},'length',{length(ind)},'link',{-1*ones(N,1)},'valid',{valid(ind)});
mol.cat = cat(ind);
mol.x = x(ind);
mol.y = y(ind);
mol.z = z(ind);
mol.h = h(ind);
mol.area = area(ind);
mol.width = width(ind);
mol.phi = phi(ind);
mol.Ax = Ax(ind);
mol.bg = bg(ind);
mol.I = I(ind);
mol.frame = frame(ind);
mol.length = length(ind);
mol.link = -1*ones(N,1);
mol.valid = valid(ind);

clear A;

f = mol;