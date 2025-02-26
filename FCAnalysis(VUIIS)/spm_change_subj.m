function spm_change_subj(orig_jobfile, orig_subj, new_jobfile, new_subj)
    % Load an SPM12 job structure from the specified file, and change
    % all occurrences of the string orig_subj to new_subj in paths,
    % filenames, etc.  Using this command, you can create a job file
    % manually in the GUI once, then convert it for each new subject
    % with a single matlab command (as long as directory structures and
    % file name conventions are the same for all subjects).  The output
    % job file new_jobfile will not be overwritten.  Example:
    % $Log: spm_change_subj.m, v $
    % Revision 1.42016/04/12 

    if ~exist(orig_jobfile, 'file')
        error(['Existing job file ' orig_jobfile ' not found.']);
    end
    
    load(orig_jobfile);
    matlabbatch = itemchange(matlabbatch, orig_subj, new_subj);
    save(new_jobfile, 'matlabbatch');
    
    function newitem = itemchange(item,origstr,newstr)

        iteminfo = whos('item');
        switch iteminfo.class

          case 'cell'
            for l = 1:length(item)
                g = item{l};
                g = itemchange(g,origstr,newstr);
                item{l} = g;
            end

          case 'struct'
            for l = 1:length(item)
                fnames = fieldnames(item(l));
                for f = 1:length(fnames)
                    s = struct('type','.','subs',fnames{f});
                    g = subsref(item(l),s);
                    g = itemchange(g,origstr,newstr);
                    item(l) = subsasgn(item(l),s,g);
                end
            end

          case 'char'
            item = strrep(item,origstr,newstr);

          otherwise

        end
        newitem = item;