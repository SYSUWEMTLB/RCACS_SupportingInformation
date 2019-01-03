function [ModelData,StopInfo] = read_ECOM_GCMTWB(FileToRead,SkillInfo,iStep)

ModelData.ELE =[];
ModelData.Z =[];
ModelData.S =[];
ModelData.T =[];
ModelData.U =[];
ModelData.V =[];

% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag5';
        'single' [SkillInfo.VPTS,1] 'Z';
        'uint32'  [1]  'tag6';
        'uint32'  [1]  'tag7';
        'single' [SkillInfo.VPTS,SkillInfo.KMB1] 'U';
        'uint32'  [1]  'tag8';
        'uint32'  [1]  'tag9';
        'single' [SkillInfo.VPTS,SkillInfo.KMB1] 'V';
        'uint32'  [1]  'tag10';
        'uint32'  [1]  'tag11';
        'single' [SkillInfo.VPTS,SkillInfo.KMB1] 'S';
        'uint32'  [1]  'tag12';
        'uint32'  [1]  'tag13';
        'single' [SkillInfo.VPTS,SkillInfo.KMB1] 'T';
        'uint32'  [1]  'tag14';
        }, ...
        'Repeat',1,  ...
        'offset',(iStep-1)*4*(3+SkillInfo.VPTS+2+4*(SkillInfo.VPTS*SkillInfo.KMB1+2)) ... %数据块的大小
        );
    mdata = m.data;
catch ME
    StopInfo.StopStep = iStep-1;
    StopInfo.EndOfFile = 1;
    StopInfo.ME = ME;
    return;
end
StopInfo.StopStep = iStep;
StopInfo.EndOfFile = 0;

ModelData.TIME = double(mdata(1).time);
ModelData.Z = double(mdata(1).Z);
ModelData.U = double(mdata(1).U);
ModelData.V = double(mdata(1).V);
ModelData.S = double(mdata(1).S);
ModelData.T = double(mdata(1).T);

return
end