% compute protein surface roughness
% read in the file containing training pdb ids
trainingFile = '../../training1.txt';
fid = fopen(trainingFile,'r');
if fid > 0
    complexes = textscan(fid,'%s');
    complexes = complexes{1};
    trainingSize = length(complexes);
    complexes = cell2mat(complexes);
    fclose(fid);
else
    disp('There is no trainig complex.');
    return;
end

trainingComplexes = char();
trainingProteins = char();
for i = 1:trainingSize
    complex = complexes(i,:);
    trainingComplexes(i,:) = complex(1:4);
    trainingProteins(i) = complex(6);
    trainingComplex = trainingComplexes(i,:);
end

% read in the file containing testing pdb ids
testingFile = '../../testing.txt';
fid = fopen(testingFile,'r');
if fid > 0
    complexes = textscan(fid,'%s');
    complexes = complexes{1};
    testingSize = length(complexes);
    complexes = cell2mat(complexes);
    fclose(fid);
else
    disp('There is no target complex.');
    return;
end

targetComplexes = char();
targetProteins = char();
for i = 1:testingSize
    complex = complexes(i,:);
    targetComplexes(i,:) = complex(1:4);
    targetProteins(i) = complex(6);
    targetComplex = targetComplexes(i,:);
end

complexes = [trainingComplexes;targetComplexes];
proteinChains = [trainingProteins';targetProteins'];
complexNum = size(complexes,1);

radius = 0.2:0.1:4;
for i = 1:complexNum
    complex = complexes(i,:);
    proteinChain = proteinChains(i);
    sasa = [];
    for j = 1:39
        sasaFile = ['res_sasa/res_sasa_' complex '_' proteinChain '_' num2str(radius(j)) '.dat'];
        fid = fopen(sasaFile);
        content = textscan(fid,'%d %s %s %f');
        proteinNames = content{2};
        sasa(:,j) = content{4};
        fclose(fid);
    end
    load(['../../precomputed/' complex '.mat']);
    proteinSurface = File.proteinSurface;
    nts = File.NT;
    chainID = cat(1,nts.chainID);
    resType = cat(1,nts.type);
    ind1 = chainID == proteinChain;
    ind2 = resType(:,1)=='a';
    ind = find(sum(ind1+ind2,2)==2);
    interactions = File.interaction2;
    interactions = interactions(ind);
    proteinSurface = proteinSurface(ind);
    if isfield(File,'resiRoughness')
        resiRoughness = File.resiRoughness;
    else
        resiRoughness = nan(File.NumNT,1);
    end
    k = 1;
    for j = 1:length(ind)
        resName = proteinNames{j};
        if strcmp(resName,'ASP') || strcmp(resName,'GLU') || strcmp(resName,'LYS') || strcmp(resName,'ARG')...
                ||strcmp(resName,'HIS') || strcmp(resName,'TYR') || strcmp(resName,'TRP') || strcmp(resName,'PHE')...
                ||strcmp(resName,'CYS') || strcmp(resName,'MET') || strcmp(resName,'SER') || strcmp(resName,'THR')...
                ||strcmp(resName,'ASN') || strcmp(resName,'GLN') || strcmp(resName,'GLY') || strcmp(resName,'ALA')...
                ||strcmp(resName,'VAL') || strcmp(resName,'LEU') || strcmp(resName,'ILE') || strcmp(resName,'PRO')
            sasas = sasa(j,:);
            p = polyfit(log(radius)',log(sasas)',1);
            resiRoughness(ind(k)) = 2-p(1);
            temp = nts(ind(k));
            disp([temp.resName ' ' resName]);
            k = k+1;
        end
    end
    File.resiRoughness = resiRoughness;
    save(['../../precomputed/' complex '.mat'], 'File');
end
