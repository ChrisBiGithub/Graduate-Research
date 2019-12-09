% find protein residues whose ASA losts are larger than 1A^2

threshold = 1;

% read in the file containing training pdb ids
% file = '../../../training1.txt';
file = '../../../testing.txt';
fid = fopen(file,'r');
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

fid = fopen('nominal maximum area.txt');
temp = textscan(fid,'%s %f');
aaName = cell2mat(temp{1});
nma = temp{2};

% download training pdb files if they do not exist in the benchmark folder
trainingComplexes = char();
trainingProteins = char();
for i = 1:trainingSize
    ASA = struct;
    complex = complexes(i,:);
    trainingComplexes(i,:) = complex(1:4);
    trainingProteins(i) = complex(6);
    trainingComplex = trainingComplexes(i,:);
    trainingProtein = trainingProteins(i);
    fileName1 = ['res_sasa_' trainingComplex '_' trainingProtein '.dat'];
    fileName2 = ['res_sasa_' trainingComplex '_' trainingProtein '_RNA.dat'];
    fid1 = fopen(fileName1,'r');
    fid2 = fopen(fileName2,'r');
    temp1 = textscan(fid1, '%d %s %s %f');
    temp2 = textscan(fid2, '%d %s %s %f');
    proteinNames = temp1{2};
    ASA1 = temp1{4};
    ASA2 = temp2{4};
    ASAlost = ASA1 - ASA2;
    ind = ASAlost > threshold;
    interaction = struct;
    k = 1;
    for j = 1:length(ind)
        resName = proteinNames{j};
        if strcmp(resName,'ASP') || strcmp(resName,'GLU') || strcmp(resName,'LYS') || strcmp(resName,'ARG')...
                ||strcmp(resName,'HIS') || strcmp(resName,'TYR') || strcmp(resName,'TRP') || strcmp(resName,'PHE')...
                ||strcmp(resName,'CYS') || strcmp(resName,'MET') || strcmp(resName,'SER') || strcmp(resName,'THR')...
                ||strcmp(resName,'ASN') || strcmp(resName,'GLN') || strcmp(resName,'GLY') || strcmp(resName,'ALA')...
                ||strcmp(resName,'VAL') || strcmp(resName,'LEU') || strcmp(resName,'ILE') || strcmp(resName,'PRO')
            interaction.name = proteinNames{j};
            interaction.chain = trainingProtein;
            interaction.ASAlost = ind(j);
            temp = repmat(interaction.name,20,1);
            tempInd = sum(abs(temp-aaName),2)==0;
            if sum(tempInd)>0
                interaction.surface = (ASA1(j)/nma(tempInd))>0.1;
            end

            ASA.interaction(k) = interaction;
            k = k+1;
        end
    end
    save([trainingComplex '_' trainingProtein '_ASA.mat'], 'ASA');
    fclose(fid1);
    fclose(fid2);
end