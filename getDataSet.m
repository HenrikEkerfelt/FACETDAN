%% getDataSet(dataSetID)
function [data, header] = getDataSet(dataSetID)

    dataSetID = num2str(dataSetID);
    
    
    %% Check if already looked up:
    storagePath = 'pathStorage.mat';
    fieldName = sprintf('%s%s','field',dataSetID);
    paths = struct;
    
    if exist(storagePath,'file')
        load(storagePath);
        
        if isfield(paths,fieldName)
            dataSetInfo = getfield(paths,fieldName);
            if exist(dataSetInfo.path, 'file')
                load(dataSetInfo.path);
                data = data.raw;
                header = dataSetInfo.header;
                return
            end
        end
    end

    
    %% Search for dataset
    % find ~ -name "*21190" 2>/dev/null
    [foundDir, outp] = unix(sprintf( 'find ~ -name "*%s" 2>/dev/null', dataSetID ));
    if ~foundDir
        error('In getDataSet.m: Could not find the dataSet');
        
    end
    lines = splitlines(outp);
    if length(lines) == 2
        path = lines{1};
    else 
        error('In getDataSet.m: Data set not unique');
    end
    
    %% Find .mat file
    expr = 'E\d{3}_\d{5}$';
    startIdx = regexp(path,expr);
    dsName = sprintf('%s%s',path(startIdx:end),'.mat');
    
    matfile = sprintf('%s/%s',path,dsName);
    if exist(matfile)
        load(matfile);
    else
        error('In getDataSet.m: Could not find the .mat file')
    end
    
    %% Find header
    expr = '/nas/nas-li20-pm00/';
    idx = regexp(path,expr);
    if isempty(idx)
        error('In getDataSet.m: Something is wrong with the find header, idx is empty' )
    elseif length(idx) > 1
        error('In getDataSet.m: Something is wrong with the find header, idx gives multiple matches' )
    end
    if idx ~= 1 
        header = path(1:idx-1);
    else
        header = '';
    end
    
    %% Save to known data sets
    dataSetInfo.path = matfile;
    dataSetInfo.header = header;
    paths = setfield(paths,fieldName,dataSetInfo);
    
    save(storagePath, 'paths')
    
    data = data.raw;

    
        
end