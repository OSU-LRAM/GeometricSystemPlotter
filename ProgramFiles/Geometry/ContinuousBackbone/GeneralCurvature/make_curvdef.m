function curv_fun = make_curvdef(curv_fun_string,paramlist,name,attempt_analytic,syspath,sysplotterpath)
% This function allows a user to specify the curvature of a backbone via a
% function contained in a string or symbolic expression. It pre-processes this curvature
% definition so that it can be integrated into the full backbone geometry
% via the function backbone_from_general_curvature. The results from this
% preprocessing are stored in an mfile curv_* located in the directory
% specified by the syspath input.
%
%
% Inputs:
%
% curv_fun_string: A string or symbolic expression describing the curvature
%   of the body as a function of position s along the body and any set of
%   shape variables. The position s is defined on a normalized unit-length
%   body centered at zero (in the range [-.5 .5]), and curvature is
%   interpreted to mean dtheta/ds of the tangent-line orientation when moving
%   in the positive-s direction.
%
% paramlist: A cell array of strings, each of which corresponds to one of
%   the shape variables used in defining curv_fun_string. The order of
%   these strings in the cell array determines the order in which they
%   appear in the curvature function created from curv_fun_string
%
% name: A string that will be used to name the curvature function created
%   from curv_fun_str.
%
% attempt_analytic (optional) : When creating a curv_ function,
%   make_curvdef will attempt to treat curv_fun_string as an analytic
%   function to the greatest extent possible, resulting in faster
%   performance when this function is evaluated (because it will require
%   fewer numerical integrations). Not all functions can be handled in this
%   way, and make_curvdef will fall back to a numerical function
%   specification if its attempt fails. For complicated curv_fun_str
%   inputs, this function can take significant time to determine that it
%   needs to fall back. To bypass the checks and directly use the numerical
%   methods, a 4-element boolean vector can be provided in this input.
%   Each element of this vector tells make_curvdef whether to attempt an
%   analytic form for the four elements of the output function. (see
%   description of these four elements in the "output" documentation below)
%
% syspath (optional unless make_curvdef is called
%   from within a sysf_ file) : Directory where output function should be placed.
%   If not specified, the directory is taken as the current user's Systems
%   directory, loaded from sysplotter_config. If make_curvdef  is called
%   from within a sysf_ file, this argument should be 'pathnames.syspath'
%
% sysplotterpath (optional unless make_curvdef is called
%   from within a sysf_ file) : Location of the sysplotter directory in
%   which the template for the curv_ file is located. If not specified, the
%   directory path is loaded from sysplotter_config. If make_curvdef  is called
%   from within a sysf_ file, this argument should be 'pathnames.sysplotterpath'
%
%
% Outputs:
%
% curv_fun: The output of this function is a handle to the curvature
%   function it creates, whose mfile was saved into the syspath directory
%   (loaded from sysplotter_confg or specified in the input).
%
%
%   This function takes in two arguments:
%
%       params: a vector of shape parameter values
%
%       mode: a string specifying what aspect of the curvature and related functions should be
%           returned. This can be:
%
%           'curvature': The curvature of points on the backbone, as specified in
%               curv_fun_str
%           'orientation': The orientation of tangent lines on the
%               backbone, relative to the tangent line at the midpoint of
%               the backbone. (This is the integral of curvature with
%               respect to the position s along the unit-normalized
%               backbone
%           'dcurvature': The derivative of the curvature with respect to
%               the shape parameters.
%           'dcurvature_int': The integral with respect to normalized
%               backbone position s of the derivative of curvature with respect to
%               the shape parameters
%
%   The output of curv_fun *is itself a function*. It takes arguments of s
%       in the range of [-.5 .5] (unit backbone centered at zero) and
%       returns the value of the quantity previous specified in "mode" as
%       an output.
%       


% Load the paths where files should be found

if  ~exist('syspath','var')
    
    load('sysplotter_config','syspath')
    
end


if ~exist('sysplotterpath','var')
    
    load('sysplotter_config','sysplotterpath')
    
end

% Initialize a cell structure to hold the list of functions to operate on
flist = {};

% Symbolically calculate all of the requested functions unless told
% otherwise
if ~exist('attempt_analytic','var')
    attempt_analytic = [1,1,1,1];
end

% Turn the curvature function from a string into a symbolic expression
if attempt_analytic(1)
    syms(paramlist{:},'s');
    curv_fun = eval(curv_fun_string);
    flist = [flist,{curv_fun}];
else
    flist = [flist,{[]}];
end

% Integrate the curvature function along the backbone
if attempt_analytic(2)
    int_curv_ds_fun = int(curv_fun,'s');
    flist = [flist,{int_curv_ds_fun}];
else
    flist = [flist,{[]}];
end

% Take the derivative of the curvature with respect to each parameter
if attempt_analytic(3)
    d_curv_dp_fun = jacobian(curv_fun,sym(paramlist));
    flist = [flist,{d_curv_dp_fun}];
else
    flist = [flist,{[]}];
end

% Integrate the parameter-derivative of the curvature function along the
% backbone
if attempt_analytic(4)
    int_d_curv_dp_ds_fun = int(d_curv_dp_fun,'s');
    flist = [flist,{int_d_curv_dp_ds_fun}];
else
    flist = [flist,{[]}];
end



% %%%%%%
% %%%%%%
%%%
%%%
% Print functions out to mfiles in a tempdir

% List of all the function names
fnamelist = {'curv_fun','int_curv_ds_fun','d_curv_dp_fun','int_d_curv_dp_ds_fun'};

% Make a temporary directory to hold mfiles for curvature functions
curvdef_tempdir = fullfile(sysplotterpath,'Utilities','curvature_mode_toolbox',...
'curvdef_tempdir');

mkdir(curvdef_tempdir)

% Test these functions to see if we can turn them into matlabFunctions
valid_mfile = ones(size(fnamelist));
for idx = 1:numel(fnamelist)
   
    if attempt_analytic(idx)
        
        %%%%%
        % See if we can turn this into a matlab function
        ftest = matlabFunction(flist{idx},'vars',['s';paramlist(:)]);
        testvars = num2cell(zeros(size(paramlist)));

        no_analytical = 0;
        try
            ftest(0,testvars{:}); % Attempt to call the function generated by matlabfunction
        catch ME
            % If matlab couldn't make an analytical form of the equation, mark
            % this
            if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
                no_analytical = 1;
            end 
        end
        
    else
       
        % Mark as no analytical form if none was attempted
        no_analytical = 1;
        
    end
     
    
    %%%%
    % If the matlabfunction passes the analytical-form test, save it to a
    % file and then concatenate it into the curvdef file. If it fails the
    % test, evaluate it at high density, save that to a file, and append a
    % file that loads this data and interpolates it.
            
    if ~no_analytical
        % Save the function to an mfile in the tempdir

        matlabFunction(flist{idx},'file',fullfile(curvdef_tempdir,fnamelist{idx}),'vars',['s';paramlist(:)]);
        
        switch fnamelist{idx}
            
            case 'curv_fun'
                
                curv_fun_call = {'output = @(s) curv_fun(s,params{:});'};

            case 'int_curv_ds_fun'
                
                int_curv_ds_call = {'output = @(s) int_curv_ds_fun(s,params{:});'}; 

            case 'd_curv_dp_fun'
                
                d_curv_dp_call = {'output = @(s) d_curv_dp_fun(s,params{:});'};
                             

            case 'int_d_curv_dp_ds_fun'
                
                int_d_curv_dp_ds_call = {'output = @(s) int_d_curv_dp_ds_fun(s,params{:});'}; 
                               
            otherwise
                
                error('fnamelist entry is not handled by switch statement')
                               
        end
        
    else
        
        % Mark that no mfile was created (so that later code does not
        % attempt to copy it in
        valid_mfile(idx) = 0;
        
        % Fill in the function calls
        switch fnamelist{idx} 
        
            
            case 'curv_fun'
                
                error('Curvature could not be turned into a matlab function. Either implement alternate behavior in make_curvdef, or correct curvature definition')
                             
            case 'int_curv_ds_fun'
                
                int_curv_ds_call = {'%% Padded length of unit backbone'...
                                   ,'all_limits = [-.51 0 .51];'...
                                   ,''...
                                   ,'%% Make dummy integration function'...
                                   ,['curv_fun_dummy = curv_' name '(cell2mat(params),''curvature'');']...
                                   ,'curvature = @(s,~) curv_fun_dummy(s);'...
                                   ,''...
                                   ,'%% Integral of the integrand function along s'...
                                   ,'output = ode_multistart(@ode45,curvature,all_limits,0,0);'}; 

            case 'd_curv_dp_fun'
                
                 d_curv_dp_call = {'%% Create a dummy function that takes in a vector of parameters'...
                                  ,'%% including s, and deals them into individual function parameters'...
                                  ,'curv_intermediate = @(all_params) vector_to_list_input(@curv_fun,all_params);'...
                                  ,''...
                                  ,'%% Create a function that takes the jacobian of the dummy function'...
                                  ,'fulljacobian = @(s) jacobianest(curv_intermediate,[s,params{:}]);'...
                                  ,''...
                                  ,'%% Create a function that truncates the s-derivative from the full jacobian'...
                                  ,'output = @(s) reshape_truncate_jacobian(fulljacobian(s));'};
                             

            case 'int_d_curv_dp_ds_fun'
                
                int_d_curv_dp_ds_call = {'%% Padded length of unit backbone'...
                                   ,'all_limits = [-.51 0 .51];'...
                                   ,''...
                                   ,'%% Make dummy integration function'...
                                   ,['d_curv_dp_fun_dummy = curv_' name '(cell2mat(params),''dcurvature'');']...
                                   ,'dcurvature = @(s,~) d_curv_dp_fun_dummy(s)'';'...
                                   ,''...
                                   ,'%% Integral of the integrand function along s'...
                                   ,'dummy_output = ode_multistart(@ode45,dcurvature,all_limits,0,zeros(size(params(:).'')));'...
                                   ,'output = @(t) transpose(dummy_output(t));'}; 
                               
            otherwise
                
                error('fnamelist entry is not handled by switch statement')
                               
        end              
        
    end
    

    
end



% Open the template file
fidi = fopen(fullfile(fileparts(which('make_curvdef')),...
    'make_curvdef_template.txt'));

% Create the output file
strrep(name,'curv_',''); % avoid writing curv_ twice if user already put it in
curv_name = ['curv_' name];
fido = fopen(fullfile(syspath,[curv_name '.m']),'w');

% Place the input function string as a comment in the top line of the file
fprintf(fido,'%% %s\n',char(curv_fun_string));

while ~feof(fidi)
    
	% Read the next line from the input file
	s = fgetl(fidi);
        
    
    % Replace the function calls
    switch strtrim(s)
        
        case 'AA_curv_fun_call'
            
            s = sprintf('\t\t%s\n',curv_fun_call{:});
            
        case 'AA_int_curv_ds_call'
            
             s = sprintf('\t\t%s\n',int_curv_ds_call{:});
            
        case 'AA_d_curv_dp_call'
            
             s = sprintf('\t\t%s\n',d_curv_dp_call{:});
            
        case 'AA_int_d_curv_dp_ds_call'
            
             s = sprintf('\t\t%s\n',int_d_curv_dp_ds_call{:});
            
        otherwise
            
            % Insert the filename
            s = strrep(s,'AA_CURVFILENAME',['curv_' name]);
            
    end


    % Copy put the line into the new file
    fprintf(fido,'%s\n',s);
    
end


% Close the template function
fclose(fidi);

% Concatenate any mfile functions onto the end of the template

for idx = 1:numel(flist)
    
    if valid_mfile(idx) %Only copy valid mfiles into the 

        % Insert a blank line
        fprintf(fido,'%s\n',[]);
        firstline = 1;

        % Load the matlab function
        fidi = fopen(fullfile(curvdef_tempdir,[fnamelist{idx} '.m']));



        % Copy all lines of the function into the curvdef file
        while ~feof(fidi)

            % Read the next line from the input file
            s = fgetl(fidi);

            % Print the line unless it is a comment or empty
            if ~strncmp(s,'%%',1) && ~strcmp(s,'')

                % If it is not the first line, indent
                if firstline
                    firstline = 0;
                else
                    fprintf(fido,'%s\t',[]);
                end

                fprintf(fido,'%s\n',s);           
            end

        end

        % Insert an end command and an extra blank line
        fprintf(fido,'%s\n\n','end');

        fclose(fidi);
    
    end

end



% cleanup
rmdir(curvdef_tempdir,'s')

fclose(fido);

% pass handle to newly created function as output
curv_fun = str2func(curv_name);

end