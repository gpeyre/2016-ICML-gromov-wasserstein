function C = extract_levelset(t,U, ls)

% extract_levelset - extract the longest connected levelset curve 
%
%   C = extract_levelset(t,U, ls);
%
%   t is a vector indexing the X/Y locations of pixels
%
%   Copyright (c) 2015 Gabriel Peyre

[c] = contourc(t,t,U, [ls ls]);
C = [];
L = -Inf;
while not(isempty(c))
    l = c(2,1); c(:,1) = [];
    if l>L
        L = l; 
        C = c(:,1:l);
    end
    c(:,1:l) = [];
end

end