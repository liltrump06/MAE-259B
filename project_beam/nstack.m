function n = nstack(ypos,stack)
ny = -(2*stack(1)*ypos+stack(2));
nx = 1;
n = [nx;ny]/norm([nx;ny]);
