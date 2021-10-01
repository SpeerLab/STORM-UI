function t = statxture(f, scale)
%statxture computes statsitical measures of texture in an image.
%T=statxture(f,scale) computes six measures of texture from an
%image(region) f. Parameter scale is a 6-dim row vector whoses elements
%multiply the 6 corresponding elements of t for scaling purposes. if scale
%is not provided, it defaults to all 1's.  The output t is 6x1 vector with
%the following elements:
%t(1) =  average gray level
%t(2) = average contrast
%t(3) = measure of smoothness
%t(4) = third moment
%t(5) = measure of uniformity
%t(6) = entropy

if nargin ==1
    scale(1:6) = 1;
else %make sure it's a row vecotr.
    scale(:)';
end
%obtain histogram and normalize it.
p=imhist(f);
p=p./numel(f);
L = length(p);

%compute the three moments. we need the unnormalized ones from function
%statmoments, these are in vector mu.
[v, mu] = statmoments(p,3);
%compute the six texture measures:
%Average gray level.
t(1) = mu(1);
%standard deviation.
t(2) = mu(2).^0.5;
%smoothness.
%first normalize the variance to [0 1] by dividing it by (L-1)^2
varn = mu(2)/(L - 1)^2;
t(3) = 1 - 1/(1 + varn);
%third moment (normalized by  (L-1)^2 also).
t(4) =  mu(3)/(L - 1)^2;
%uniformity
t(5) = sum(p.^2);
%entropy.
t(6) = -sum(p.*(log2(p + eps)));

%scale the values.
t = t.*scale;