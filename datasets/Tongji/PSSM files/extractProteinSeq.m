% extract protein sequence
trainingFile = '../benchmark/training1.txt';
testingFile = '../benchmark/testing.txt';

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

trainingComplexes = complexes;

fid = fopen(testingFile,'r');
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

testingComplexes = complexes;

allComplexes = [trainingComplexes; testingComplexes];

proteinSize = size(allComplexes,1);
load('aacode.mat');
proteinSeqsFile = 'proteinSequences.txt';
for i = 49:50
    complex = allComplexes(i,1:4);
    chainid = allComplexes(i,6);
    fileName = ['../benchmark/precomputed/' complex '.mat'];
    load(fileName);
    NTs = File.NT;
    resType = cat(1,NTs.type);
    resChain = cat(1,NTs.chainID);
    ind1 = resType(:,1)==repmat('a',length(resType(:,1)),1);
    ind2 = resChain==repmat(chainid,length(resChain),1);
    ind = (ind1+ind2)==2;
    targetAAs = NTs(ind); % AAs on the selected chain
    seq = cat(1,targetAAs.resName);
    seqFASTA = '';
    seqLen = sum(ind);
    for j = 1:seqLen
        aaName = seq(j,:);
        ind = sum(abs(repmat(aaName,20,1)-aaCodes(:,1:3)),2)==0;
        code = aaCodes(ind,4);
        seqFASTA = [seqFASTA code];
    end
    fid = fopen(proteinSeqsFile,'w');
    fwrite(fid,seqFASTA);
    fclose(fid);
    dos(['psiblast -query proteinSequences.txt -db nr -inclusion_ethresh 0.001 -num_iterations 3 -out_ascii_pssm pssm' num2str(i) '.txt -outfmt "7 score"']);
end