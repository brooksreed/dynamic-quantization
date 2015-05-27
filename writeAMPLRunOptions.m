function writeAMPLRunOptions(fid,formulation,generalOptions,solveOptions)
% writes general AMPL .run file beginning (including solver options)
% (following this - model, data, solve, output options)
%
% writeAMPLRunOptions(fid,formulation,generalOptions,solveOptions)
%
% fid = file handle ID to write to (from previous fopen call)
%
% formulation: 'LP' or 'Int' (each req specialized options)
%
% generalOptions (all defaults = 0)
%   struct with
%   '.noAMPLPresolve'
%   '.noCPLEXPresolve'
%   '.AMPLtimes'
%   '.parallelmode'
%   '.threads'
%   '.timeoutCPLEX'
%
% NOTE: solveOptions can have any combo of fields
%   correct fields for formulation are checked
%
% if formulation=='LP': solveOptions (OPTIONAL)
%   struct with (any combination of the following)
%   '.netopt' (0,1,2)
%   '.alg' ('primal','dual','barrier','auto')
%   '.PD'  ('primal','dual')
%
% if formulation=='Int': solveOptions (OPTIONAL)
%   struct with MIP solver options
%   '.<CPLEX directive/' = value
% (for example)
%   '.solutionlim' = %i
%   '.timeoutBB' = %i
%   '.mipemphasis' = {0,1,2}
%   '.relobjdiff' (B&B branch throw-away relaxation)
%   '.varselect'
%   '.probe'
%   '.mipgap'
%   '.mipcuts'
%   '.mipsearch'
%   '.mipalgorithm'
%
% EXTENDED GENERAL OPTIONS DESCRIPTION (from AMPL/CPLEX user manual)
%
% timeoutCPLEX: in seconds (overridden if solutionlim set)
%
% parallelmode: deterministic or opportunistic
% parmode=0;      % 0 for default (auto - depends on threads)
% parmode=-1;     % -1 for opportunistic mode
% parmode = 1;    % 1 for deterministic mode
%
% threads
% value 0 tells CPLEX to use as many threads as possible, subj to:
%   - allowed by license when parallelmode = -1
%   - maintaining deterministic algorithm when parallelmode = 0
% positive value specifies that number of threads should be used

% BR 1/30/2012
% changelog: author,date,change
%{
-   added CPLEX memoryemphasis hardcorded in (applies to LP and MIP)
-   2/7 - added auxfiles from ampl to temp dir
-   2/14/2012 - added capability to set fields with no '='
            (no '=' if field is empty)
-   2/27/2012 - changed so CPLEX timeout enforced ALWAYS
        (not overridden by solutionlimit...)
- 3/4/2013: fixed bug with LP solver options (ignoring alg)
%}

if(nargin<4)
    solveOptions=[];    % all defaults
    generalOptions=[];
end
if(nargin<3)
    solveOptions=[];
end

% DEFAULTS (NO SPECIFICATION) USED FOR ANY UNSET PARAMETERS

% generalOptions have defaults...
if(~isfield(generalOptions,'noAMPLPresolve'))
    generalOptions.noAMPLPresolve=0;
end
if(~isfield(generalOptions,'noCPLEXPresolve'))
    generalOptions.noCPLEXPresolve=0;
end
if(~isfield(generalOptions,'AMPLtimes'))
    generalOptions.AMPLtimes=0;
end
if(~isfield(generalOptions,'parallelmode'))
    generalOptions.parallelmode=0;
end
if(~isfield(generalOptions,'threads'))
    generalOptions.threads=0;
end
if(~isfield(generalOptions,'timeoutCPLEX'))
    generalOptions.timeoutCPLEX=0;
end


% Write AMPL options

fprintf(fid, 'option solver cplexamp;\n'); % this calls academic (full)
fprintf(fid, 'option show_stats 1;\n');
if(generalOptions.AMPLtimes)
    fprintf(fid, 'option times 1;\n');
end
%fprintf(fid,' option gentimes 1;\n');
%fprintf(fid, 'option auxfiles ''acfrsu'';\n');

% disable presolve
if(generalOptions.noAMPLPresolve)
    fprintf(fid, 'option presolve 0;\n');
end

% CPLEX OPTIONS (settings from start of fcn):
fprintf(fid, '\noption cplex_options ''timing 1 ''\n');
fprintf(fid, '''parallelmode=%d ''\n',generalOptions.parallelmode);
fprintf(fid, '''threads=%d ''\n',generalOptions.threads);

if(isfield(solveOptions,'memoryemphasis'))
    fprintf(fid, '''memoryemphasis=%i ''\n',solveOptions.memoryemphasis);
end

% CPLEX Timeout (solutionlim overrides)
%if((generalOptions.timeoutCPLEX>0) ...
%        && (~isfield(solveOptions,'solutionlim')))
if(generalOptions.timeoutCPLEX>0)
    fprintf(fid, '''timelimit=%d ''\n',generalOptions.timeoutCPLEX);
    %fprintf('CPLEX TIMEOUT: %d \n',generalOptions.timeoutCPLEX);
end

% get fields from solveOptions
% print options...(generic template)
if(~isempty(solveOptions))
    names = fieldnames(solveOptions)
else
    names=[];
end


check={'netopt','alg','PD'};
for i = 1:length(check)
    ind = find(strcmp(names,check{i}));
    if(ind)
        names{ind}='ignore';
    end
end

switch formulation
    case 'int'
        % MIP solver options
        
        % WARNINGS WHEN PARAMETER NOT RECOGNIZED...
        % (could just not be hard-coded yet)
        
        % test that LP options not set for int...
        %         if(...
        %                 max([max(strcmp(names,'alg'));
        %                 max(strcmp(names,'PD'));
        %                 max(strcmp(names,'netopt'))]))
        %             fprintf('\nWarning - LP options set with Int formulation\n')
        %         end
        
        
        if(~isfield(solveOptions,'mipdisplay'))
            fprintf(fid, '''mipdisplay=1 ''\n'); % display LP iterations
        end
        
        for i = 1:length(names)
            if(~strcmp(names{i},'ignore'))
                field = getfield(solveOptions,names{i});
                if(~isempty(field))
                    fprintf(fid, '''%s=%i ''\n',names{i},field);
                else
                    fprintf(fid, '''%s ''\n',cmd,names{i});
                end
                %fprintf(fid, '''%s = %i ''\n',names{i},...
                %    getfield(solveOptions,names{i}));
                
            end
        end
        
    case 'LP'
        
        % test that only  LP options set
        % if none of names = an LP option & names isn't empty, warning
        testVec=[max(strcmp(names,'alg'));
            max(strcmp(names,'PD'));
            max(strcmp(names,'netopt'))];
        if((~max(testVec)) && (~isempty(names)))
            fprintf('\nWarning - no LP options set\n')
        end
        
        % LP solver choices here
        
        if(isfield(solveOptions,'netopt'))
            fprintf(fid, '''netopt=%i ''\n',getfield(solveOptions,'netopt'));
        end
        
        % tell CPLEX which algorithm to use
        if(isfield(solveOptions,'alg'))
            alg=solveOptions.alg;
            switch alg
                case 'primal'
                    fprintf(fid, '''primalopt ''\n');
                    disp('Forced Primal Opt')
                case 'dual'
                    fprintf(fid, '''dualopt ''\n');
                    disp('Forced Dual Opt')
                case 'barrier'
                    fprintf(fid, '''baropt ''\n');
                    disp('Forced Barrier Opt')
                case 'auto'
                    %fprintf(fid, '''autoopt ''\n');
                    disp('Auto Opt')
            end
        end
        
        % tell CPLEX which formulation to solve
        if(isfield(solveOptions,'PD'))
            PD=solveOptions.PD;
            switch PD
                case 'primal'
                    fprintf(fid, '''primal ''\n');
                case 'dual'
                    fprintf(fid, '''dual ''\n');
            end
        end
        
        % print other options as specified...
        % (I think CPLEX ignores MIP options for LP)
        for i = 1:length(names)
            if(~strcmp(names{i},'ignore'))
                
                field = getfield(solveOptions,names{i});
                
                if(~isempty(field))
                    fprintf(fid, '''%s=%i ''\n',names{i},field);
                else
                    fprintf(fid, '''%s ''\n',cmd,names{i});
                end
            end
        end
        
end
fprintf(fid,';\n'); % finish options line


