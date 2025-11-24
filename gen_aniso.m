    %===============================================================================
% gen_aniso.m is the main computation script
% no inputs for this script
%===============================================================================
clear all; 

% Dir_Name={'Mesaverde-30'};
% Dir_Name={'Mesaverde-60'};
% Dir_Name={'Mesaverde-HTI'};
% Dir_Name={'Mesaverde-0C'};
%Dir_Name={'Mesaverde-1'};
% Dir_Name={'Mesaverde-0-F'};
% Dir_Name={'Mesaverde-30-F'};
% Dir_Name={'Mesaverde-60-F'};
% Dir_Name={'Mesaverde-HTI-F'};
% Dir_Name={'Mesaverde-HTI-F1'};
% Dir_Name={'Mesaverde-HTI-F2'};
 Dir_Name={'Bakken-B'};
% Dir_Name={'Bakken-HTI'};
% Dir_Name={'Bakken-HTI-F'};
% Dir_Name={'Cotton-30'};
% Dir_Name={'Cotton-60'};
% Dir_Name={'Cotton-HTI'};
% Dir_Name={'Usari25-HTI'};
% Dir_Name={'Hornby1-HTI'};
% Dir_Name={'Minas-HTI'};
% Dir_Name={'Horne-HTI'};
% Dir_Name={'Barnett-HTI'};
% Dir_Name={'JH-HTI'};
% Dir_Name={'Deger-HTI'};
% Dir_Name={'JJ-HTI'};
% Dir_Name={'Austin-HTI'};    
% Dir_Name={'SlowForm-HTI'};
% Dir_Name={'FastForm-HTI'}; 
   
%===============================================================================
%===============================================================================
%===============================================================================
fprintf('\n\n============================================================================\n');
fprintf('Gen_Aniso Program has been started!\n');
fprintf(1,'\tThe used directory is %s\n',Dir_Name{1});
tStart_Prog=tic; % Setting the global timer for computations

    % Copy the input file in St1_3_SetModelUser_sp_SAFE.m work file located
    % in the directory contained gen_aniso.m
    run_dir=pwd;
    cd(Dir_Name{1});
    var=dir;
    for iname=1:size(var,1)
        if var(iname).isdir==0
           ivar=length(var(iname).name); 
           if strcmpi(var(iname).name(ivar-1:ivar),'.m'), break; end
        end
    end
    copyfile(var(iname).name,strcat(run_dir,'\St1_3_SetModelUser_sp_SAFE.m'));
    cd(run_dir); 

    %===============================================================================
    % Step 1. Initialization
    % E.g. set structure with parameters for computation of eigenfrequencies for
    % cylindrically layered structure with TTI layers
    %===============================================================================
    cd('private');
    InputParam = St1_SetModel();
    % Enter the PML layer parameters free on frequency
    InputParam.Model.AddDomain_Exist='none';
    InputParam.Model.AddDomainType=lower(InputParam.Model.AddDomainType); % pml, abc, abc+pml, same or none
    if strcmpi(InputParam.Model.AddDomainType,'abc+pml'), InputParam.Model.AddDomainType='pml+abc'; end
    if ~strcmpi(InputParam.Model.AddDomainType,'none')
       if strcmpi(InputParam.Model.AddDomainType,'pml') || strcmpi(InputParam.Model.AddDomainType,'abc') ...
          || strcmpi(InputParam.Model.AddDomainType,'pml+abc') || strcmpi(InputParam.Model.AddDomainType,'same')
          InputParam.Model.AddDomain_Exist='yes';
       else
          fprintf(1,'Error: AddDomainType is not set!!!\n');
       end
    end
    if strcmpi(InputParam.Model.AddDomain_Exist,'yes')
       if strcmpi(InputParam.Model.AddDomainLoc,'ext') % we have additional external layer
          InputParam.Model.DomainTheta(end+1)=InputParam.Model.DomainTheta(end);
          InputParam.Model.DomainEcc(end+1)=InputParam.Model.DomainEcc(end);
          InputParam.Model.DomainEccAngle(end+1)=InputParam.Model.DomainEccAngle(end);
          InputParam.Model.DomainParam(end+1)=InputParam.Model.DomainParam(end);
          InputParam.Model.DomainType{end+1}=InputParam.Model.DomainType{end};
          InputParam.Model.BCType{end}='SSstiff'; InputParam.Model.BCType{end+1}='rigid';
          InputParam.Model.DomainNth(end+1)=InputParam.Model.DomainNth(end);
       end
    end
    delete('output/*.*'); % delete old files in temp directory
    
    %===============================================================================
    % Step 2. Preparing data and structure for computations
    % E.g. assemble structure, which contains all parameters
    % necessary for the solution of the generalized eigenvalue problem
    %===============================================================================
    [CompStruct] = InputParam.Methods.St2_PrepareModel(InputParam);
    CompStruct.f_grid = CompStruct.Model.f_array; % Frequency grid
    
    % Cycle by frequency
    for if_grid=1:CompStruct.Model.N_disp
        CompStruct.if_grid=if_grid;
        fprintf('\tNow computing the %dth point on dispersion curve out of %d\n',if_grid,CompStruct.Model.N_disp);
        %========================================================================================
        %  Step 3. Running routines, which are necessary for the formulation of the problem
        %  E.g. assemble the spectral method matrices for eigenvalue problem
        %========================================================================================
        tStart_St3=tic; % Setting the local timer for computations
        [CompStruct,BasicMatrices,FEMatrices,FullMatrices] = CompStruct.Methods.St3_PrepareBasicMatrices(CompStruct,InputParam);
        fprintf('\t\tTime for St3 Program = %.1fs\n',toc(tStart_St3));
        
        % save Config, Methods, Model and etc. parameters in CompStruct.mat file (work file)
        fname = strcat(CompStruct.Config.root_path,'output\CompStruct.mat');
        CompStruct.ModelInitial=InputParam.Model;
        save(fname,'-struct','CompStruct','Config','Methods','ModelInitial','Model','Advanced','Misc','Data','Asymp');

        FEMatrices.DomainRx=CompStruct.Model.DomainRx;
        FEMatrices.DomainRy=CompStruct.Model.DomainRy;
        % Mesh and relevant different parameters are saved in FEMatrices-freq.mat file
        fname = strcat(CompStruct.Config.root_path,'output\FEMatrices-',num2str(CompStruct.f_grid(if_grid)),'.mat');
        %save(fname,'-struct','FEMatrices');
        save(fname,'-struct','FEMatrices', 'MMatrix_d', ...
             'MeshNodes','BoundaryEdges','MeshTri','MeshProps','PhysProp','DElements','DEMeshProps','DNodes','DNodesRem', ...
             'DNodesComp','DTakeFromVarPos','DPutToVarPos','DZeroVarPos','BNodes','BNodesFull','DomainRx','DomainRy');

        %========================================================================================
        % Step 4. Solving the problem 
        % E.g. computing the spectrum by the SAFE method
        %========================================================================================
        tStart_St4=tic; % Setting the local timer for computations
        [Results]=CompStruct.Methods.St4_ComputeSolution(CompStruct,BasicMatrices,FEMatrices,FullMatrices);
        cd(strcat(CompStruct.Config.root_path,'output\'));
        % Eigevalues and eigenvectors are saved in Results-freq.mat file
        dat_file = strcat('Results-',num2str(CompStruct.f_grid(if_grid)),'.mat'); 
        save(dat_file,'Results');
        cd(CompStruct.Config.root_path);
        
        fprintf('\t\tTime for St4 Program = %.1fs\n',toc(tStart_St4));
        
        %if isfield(CompStruct.Mesh,'Fig_handle')
        if strcmpi(CompStruct.Mesh.output,'yes')
           close(CompStruct.Mesh.Fig_handle);
        end
        
        
    end % end of cycle by frequency

    % move all results from "output" directory to "Dir_Name" directory
    movefile('output/*.*',Dir_Name{1});
    delete('St1_3_SetModelUser_sp_SAFE.m');

fprintf('Time for Gen_Aniso Program = %.1fs\n',toc(tStart_Prog));
fprintf('============================================================================\n');
