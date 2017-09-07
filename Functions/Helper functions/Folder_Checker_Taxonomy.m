function Folder_Checker_Taxonomy(sourcepath, savepath)
    
if ~exist(savepath,'dir') 
    mkdir(savepath);
end

if ~exist(sourcepath,'dir') 
    mkdir(sourcepath);
    mkdir([sourcepath, filesep 'IV']);
    mkdir([sourcepath, filesep 'data']);
    mkdir([sourcepath, filesep 'metadata']);
    mkdir([sourcepath, filesep 'pictures']);
else
    if ~exist([sourcepath, filesep 'IV'],'dir')
        mkdir([sourcepath, filesep 'IV']);
    end
    if ~exist([sourcepath, filesep 'data'],'dir') 
        mkdir([sourcepath, filesep 'data']);
    end
    if ~exist([sourcepath, filesep 'metadata'],'dir') 
        mkdir([sourcepath, filesep 'metadata']);
    end
    if ~exist([sourcepath, filesep 'pictures'],'dir') 
        mkdir([sourcepath, filesep 'pictures']);
    end
end