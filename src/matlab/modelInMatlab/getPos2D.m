function pos = getPos2D(stencil)
pos = zeros(1,size(stencil,1));
for i = 1:size(stencil,2)
    s = stencil(i);
    pos(i) = ((s.j-1)*s.xn+s.i)*s.dof - (s.dof-1 - s.c); 
end