function pos = getPos3D(stencil)
pos = zeros(1,size(stencil,1));
for i = 1:size(stencil,2)
    s = stencil(i);
    pos(i) = ((s.k-1)*s.xn*s.yn + (s.j-1)*s.xn + s.i)*s.dof - (s.dof-1 - s.c);
end
