clear;clc;

verts = [0 0 1 1; 0 1 0 1; 0 0 0 0];faces = [1 3; 2 2; 3 4];



numVerts = size(verts, 2);
numFaces = size(faces, 2);

jj= repmat(1:numFaces,[3,1]);  % 3 vtx in each face
jj=jj(:)';

mem = sparse( faces(:),jj,1,numVerts,numFaces);

mem3 = sparse(numVerts, numFaces);

for i = 1:size(verts, 2)
    mem3(i,:) = any(faces == i);
end