function savevtk(array, filename)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    [ny, nx, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, 'SCALARS scalars double\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    for z=1:nz
        for y=1:ny
            for x=1:nx
                fprintf(fid, '%f ', array(y,x,z));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return