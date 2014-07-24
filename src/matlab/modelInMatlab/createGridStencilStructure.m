% create stencil structure akin to dmda grid in petsc.
% Initializes the fields:
% Component(c) = 0; i=j=k=0;
% Input fields: xn, yn, zn and number of degrees of freedom: dof
function stencil = createGridStencilStructure(xn, yn,zn,dof)
    stencil = struct('c',0,'i',0,'j',0,'k',0,'xn',xn,'yn',yn,'zn',zn,'dof',dof);
end