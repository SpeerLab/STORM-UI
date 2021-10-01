function [x, y] = intline(x1,x2, y1, y2)
%Intline Integer-coordinate line drawing algorithm.
%[X, Y] = Intline(X1, X2, Y1, Y2) computes an approximation to the line
%segment joining (x1,Y1) and (X2,Y2) with integer coordinates. X1,X2,Y1, Y2
%should be integers.  Intline is regersible, that is, inline(X1,X2,Y1,Y2)
%produces the same results as flipud(intline(X2,X1,Y2,Y1)).

%Copyright 1993-2002. The mathworks, Inc.

dx = abs(x2-x1);
dy = abs(y2-y1);

%check for degenerate case.
if ((dx ==0) && (dy == 0))
    x= x1;
    y=y1;
    return;
end
flip = 0;
if (dx>= dy)
    if (x1 > x2)
        %Always "draw" from left to right.
        t = x1; x1 = x2; x2 = t;
        t= y1; y1 = y2; y2 = t;
        flip = 1;
    end
    m = (y2 - y1)/(x2 - x1);
    x = (x1:x2).';
    y = round(y1 + m*(x-x1));
else
    if (y1> y2)
        %always "draw" from bottom to top.
        t = x1; x1 = x2; x2 = t;
        t = y1; y1 = y2; y2 = t;
        flip = 1;
    end
    m = (x2 - x1)/(y2 - y1);
    y = (y1:y2).';
    x = round(x1 + m*(y-y1));
end
if (flip)
    x = flipud(x);
    y = flipud(y);
end
