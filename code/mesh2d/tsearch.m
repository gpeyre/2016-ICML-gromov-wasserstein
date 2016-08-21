function i = tsearch(x,y,th,px,py)

% ADDED BY GABRIEL PEYRE

i = tsearchn([x(:), y(:)],th,[px(:), py(:)]);  

end