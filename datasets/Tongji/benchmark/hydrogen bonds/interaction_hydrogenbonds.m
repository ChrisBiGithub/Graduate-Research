% find interactions casued by hydrogen bonds
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
    fid = fopen([complex '.hb2']);
    content = textscan(fid, '%s', 'Delimiter', '\n');
    content = content{1};
    hbonds = struct;
    interactionCount = 1;
    for j = 9:length(content)
        interactionStr = content{j};
        interactionResi1Chain = interactionStr(1);
        interactionResi1SeqNum = interactionStr(2:5);
        interactionResi1Type = interactionStr(7:9);
        interactionResi2Chain = interactionStr(15);
        interactionResi2SeqNum = interactionStr(16:19);
        interactionResi2Type = interactionStr(21:23);
        
        % remove white space in residue type
        interactionResi1Type = strtrim(interactionResi1Type);
        interactionResi2Type = strtrim(interactionResi2Type);
        
        switch interactionResi1Type
            case {'ASP','GLU','LYS','ARG','HIS','TYR','TRP','PHE','CYS',...
                    'MET','SER','THR','ASN','GLN','GLY','ALA','VAL','LEU','ILE','PRO'}
                switch interactionResi2Type
                    case {'A','U','G','C'}
                        resi1 = struct;
                        resi2 = struct;
                        interaction = struct;
                        resi1.chain = interactionResi1Chain;
                        resi2.chain = interactionResi2Chain;
                        resi1.resSeqNum = str2num(interactionResi1SeqNum);
                        resi2.resSeqNum = str2num(interactionResi2SeqNum);
                        resi1.type = interactionResi1Type;
                        resi2.type = interactionResi2Type;
                        interaction.resi1 = resi1;
                        interaction.resi2 = resi2;
                        hbonds.interaction{interactionCount} = interaction;
                        interactionCount = interactionCount + 1;
                end
            case {'A','U','G','C'}
                switch interactionResi2Type
                    case {'ASP','GLU','LYS','ARG','HIS','TYR','TRP','PHE','CYS',...
                    'MET','SER','THR','ASN','GLN','GLY','ALA','VAL','LEU','ILE','PRO'}
                        resi1 = struct;
                        resi2 = struct;
                        interaction = struct;
                        resi1.chain = interactionResi1Chain;
                        resi2.chain = interactionResi2Chain;
                        resi1.resSeqNum = str2num(interactionResi1SeqNum);
                        resi2.resSeqNum = str2num(interactionResi2SeqNum);
                        resi1.type = interactionResi1Type;
                        resi2.type = interactionResi2Type;
                        interaction.resi2 = resi1;
                        interaction.resi1 = resi2;
                        hbonds.interaction{interactionCount} = interaction;
                        interactionCount = interactionCount + 1;
                end
        end
    end
    save([complex '_hbonds.mat'],'hbonds');
    fclose(fid);
end