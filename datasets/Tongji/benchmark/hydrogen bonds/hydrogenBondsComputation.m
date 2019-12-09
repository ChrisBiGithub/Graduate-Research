% compute hydrogen bonds using HBPLUS
trainingFile = '../training1.txt';
fid = fopen(trainingFile,'r');
if fid > 0
    complexes = textscan(fid,'%s');
    trainingComplexes = complexes{1};
    trainingSize = length(trainingComplexes);
    trainingComplexes = cell2mat(trainingComplexes);
    trainingComplexes = trainingComplexes(:,1:4);
    fclose(fid);
else
    disp('There is no trainig complex.');
    return;
end

testingFile = '../testing.txt';
fid = fopen(testingFile,'r');
if fid > 0
    complexes = textscan(fid,'%s');
    testingComplexes = complexes{1};
    testingSize = length(testingComplexes);
    testingComplexes = cell2mat(testingComplexes);
    testingComplexes = testingComplexes(:,1:4);
    fclose(fid);
else
    disp('There is no target complex.');
    return;
end

allComplexes = [trainingComplexes;testingComplexes];
allComplexes = unique(allComplexes,'rows');

complexSize = size(allComplexes,1);
for i = 1:complexSize
    complex = allComplexes(i,:);
    system(['hbplus pdb/' complex '.pdb']);
end