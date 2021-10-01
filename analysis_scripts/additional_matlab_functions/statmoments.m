function [v, unv] = statmoments(p, n)
%Statmoments Computes statistical central moments of image histogram.
%[v, unv] = statmoments(p,n) computes up to the nth statistical central
%moment of a histogram whose components are in vector P. The length of P
%must equal 256 or 65536.
%
% The program outputs a vector V with V(1) = mean, V(2) = variance, v(3) =
% 3rd moment, ... V(n) = nth central moment. The random variable values are
% normalized to the range [0,1], so all moments are in this range.
%
% The program also outputs a vector unv containing the same moments as v,
% but using un-normalizaed random variable values (eg, 0 to 255 if
% length(p) = 2^8. For example, if length(p) = 256 and v(1)=0.5, then
% unv(1) would have value 127.5 (half of the [0 255] range).

Lp = length(p);
if (Lp ~= 256) && (Lp ~=65536)
    error('P must be a 256- or 65536-element vector.');
end
G = Lp - 1;

%Make sure the histogram has unit area, and convert it to a column vector.
p = p/sum(p); p = p(:);
%form a vector of all the possible values of the random variable.
z = 0:G;
%now normalize the z's to the range [0,1].
z=z./G;
% The mean.
m=z*p;
%center variables around the mean.
z = z - m;
%Compute the central moments.
v = zeros(1,n);
v(1) = m;
for j=2:n
    v(j) = (z.^j)*p;
end
if nargout>1
    %compute the uncentralized moments.
    unv = zeros(1:n);
    unv(1)=m.*G;
    for j= 2:n
        unv(j) = ((z*G).^j)*p;
    end
end


