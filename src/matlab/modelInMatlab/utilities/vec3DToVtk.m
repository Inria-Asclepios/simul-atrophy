% function from the following script source:
% Copyright (c) 2011, W Thielicke
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


function vec3DToVtk(u,v,w,x,y,z,filename)
%Output file name
% filename=('TestFile.vtk');

%load the MATLAB 3D flow example
% load wind
tic

nr_of_elements=numel(x);
fid = fopen(filename, 'w'); 

%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
fclose(fid);

%append binary x,y,z data
fid = fopen(filename, 'a'); 
fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');

%append another ASCII sub header
fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
fprintf(fid, 'VECTORS velocity_vectors float\n');

%append binary u,v,w data
fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b');

%append another binary u,v,w data set
fprintf(fid, '\nVECTORS another_vector_set float\n'); %ASCII header
fwrite(fid, [reshape(u*10,1,nr_of_elements);  reshape(v*2,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b'); %binary data

%append some scalar data
fprintf(fid, '\nSCALARS EinLustigerSkalar float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fwrite (fid, reshape(sqrt(u.^2+v.^2+w.^2),1,nr_of_elements),'float','b'); %binary data

fclose(fid);
toc