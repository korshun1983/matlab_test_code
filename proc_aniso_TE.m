% Proc_aniso_TE.m program runs after gen_aniso.m program and
% designed for processes the raw results of gen_aniso.m program.
% No inputs for this script

% 1. Calculate frequency dependencies of velocities and slownesses
% 2. Build polar grid (r,f)
% 3. On this polar grid calculated the dispacements ur and uf are calculated 
% 4. Calculate radial density of kinetic energy 
% 5. Calculate trigonometrical Fourier decomposition of ur,uf,uz (h=0,+-1, +-2, ..., +- NHarm-1)
% Results are saved in Slowness.mat and LayerTE-freq-layer.mat (freq in kHz, layer - its number)

clear all;

% Dir_Name={'Mesaverde-0C'};
Dir_Name={'Mesaverde-1'};
% Dir_Name={'Mesaverde-0-F'};
% Dir_Name={'Mesaverde-30-F'};
% Dir_Name={'Mesaverde-HTI-F'};
% Dir_Name={'Mesaverde-HTI-F1'};
% Dir_Name={'Mesaverde-HTI-F2'};
% Dir_Name={'Bakken-HTI-F'};

% Count for all layers, frequencies and eigenvalues
ProcTE.write_result='y'; % save (y) or not (n) results

ProcTE.DThetaInc=3; % Decrease azimuthal space parameter DTheta in DThetaInc times 
ProcTE.NHarm=4; % Number of modes of trigonometric Fourier transform (harmonics)

% Figure parameters
ProcTE.close_figs='n'; % close (y) or (n) previously opened figures
ProcTE.Polar_Grid='n'; % plot polar grid with mesh for each layer
ProcTE.num_fig=200;  % start figure number
ProcTE.font_size=12; % base fontname size

%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
fprintf('\n\n============================================================================\n');
fprintf('Proc_Aniso Program has been started!\n');
tStart_Prog=tic;

ProcTE.Dir_Name=Dir_Name;
curr_dir=pwd;
%run_prog=mfilename;
%full_path_run_prog=mfilename('fullpath');
%run_dir=full_path_run_prog(1:(length(full_path_run_prog)-length(run_prog)-1));

% Here main subprograms
CompStructA.Config.root_path=strcat(pwd,'\');

% Read Commin CompStruct data
fname_comp = char(strcat(ProcTE.Dir_Name,'\CompStruct.mat'));
fname_comp_add = char(strcat(ProcTE.Dir_Name,'\CompStruct_add.mat'));
if exist(fname_comp_add,'file')>0
    load(fname_comp_add);
else
    if exist(fname_comp,'file')>0
       load(fname_comp);
    else
       beep; beep; beep; fprintf(1,'File is not existed!\n'); return; 
    end
end
CompStructA.Config=Config;
CompStructA.Methods=Methods;
CompStructA.Model=Model;
CompStructA.Advanced=Advanced;
CompStructA.Misc=Misc;
CompStructA.Data=Data;
CompStructA.Asymp=Asymp;

Result_Freq=zeros(CompStructA.Model.N_disp,1);
ProcTE.freq_count=CompStructA.Model.f_array;
ProcTE.nfreqs=size(ProcTE.freq_count,2);

PlocTE.max_neigs=CompStructA.Advanced.num_eig_max;
Result_Velocity=zeros(CompStructA.Model.N_disp,PlocTE.max_neigs);
Result_Slowness=zeros(CompStructA.Model.N_disp,PlocTE.max_neigs);
Result_Attenuation=zeros(CompStructA.Model.N_disp,PlocTE.max_neigs);

ProcTE.max_neigs=CompStructA.Advanced.num_eig_max; % number of calculated eigs
ProcTE.neigs=1:ProcTE.max_neigs;

% Number of layers
ProcTE.nlayers_count=1:size(CompStructA.Model.DomainType,2);
ProcTE.nlayers=1:size(CompStructA.Model.DomainType,2);

if lower(ProcTE.close_figs)=='y', close all; end

% Give velocities in m/s
slow_coef=1.;
if strcmp(CompStructA.Config.SloUnits,'us/ft'), slow_coef=CompStructA.Misc.S_conv/1000.; end % 0.0254*12

%====================================================================================================================
iii_freq=0;
for ii_f=1:ProcTE.nfreqs % frequency cycle
%for ii_f=3
    %cd(run_dir);
    
    tStart=tic;
    fprintf(1,'Now fhe %d-th frequency (%.2f kHz) is processed out of %d\n', ...
            ii_f, ProcTE.freq_count(ii_f),ProcTE.nfreqs);
    ProcTE.ii_f=ii_f;
    iii_freq=iii_freq+1;
    
    % Read FEMatrices data (mesh, physical properties and etc.)
    fname = char(strcat(ProcTE.Dir_Name,'\FEMatrices-',num2str(ProcTE.freq_count(ii_f)),'.mat'));
    if exist(fname,'file')>0
        load(fname);
        FEMatricesA.MMatrix_d=MMatrix_d;
        FEMatricesA.MeshNodes=MeshNodes;
        FEMatricesA.BoundaryEdges=BoundaryEdges;
        FEMatricesA.MeshTri=MeshTri;
        FEMatricesA.MeshProps=MeshProps;
        FEMatricesA.PhysProp=PhysProp; 
        FEMatricesA.DElements=DElements;
        FEMatricesA.DEMeshProps=DEMeshProps; 
        FEMatricesA.DNodes=DNodes; 
        FEMatricesA.DNodesRem=DNodesRem;
        FEMatricesA.DNodesComp=DNodesComp; 
        FEMatricesA.DTakeFromVarPos=DTakeFromVarPos;
        FEMatricesA.DPutToVarPos=DPutToVarPos; 
        FEMatricesA.DZeroVarPos=DZeroVarPos; 
        FEMatricesA.BNodes=BNodes; 
        FEMatricesA.BNodesFull=BNodesFull;
        FEMatricesA.DomainRx=DomainRx;
        FEMatricesA.DomainRy=DomainRy;
    else
        beep; beep; beep; fprintf(1,'File is not existed!\n'); return;
    end

    %====================================================================================================================
    % Loading the results for the specific point of eigenvalue variable
    dat_file = char(strcat(ProcTE.Dir_Name,'\Results-',num2str(ProcTE.freq_count(ii_f)),'.mat'));
    load(dat_file);

    %====================================================================================================================
    % Calculate frequency dependencies of velocity and slowness
    var_mas=diag(Results.REig_vals);
    Result_Freq(ii_f)=Results.omega_val/(2.*pi);
    Result_Velocity(ii_f,1:PlocTE.max_neigs)=CompStructA.Misc.F_conv.*Results.omega_val./real(var_mas(1:PlocTE.max_neigs));
    Result_Slowness(ii_f,1:PlocTE.max_neigs)=slow_coef.*10.^6./Result_Velocity(ii_f,1:PlocTE.max_neigs);
    Result_Attenuation(ii_f,1:PlocTE.max_neigs)=-2.*imag(var_mas(1:PlocTE.max_neigs))./real(var_mas(1:PlocTE.max_neigs));

    %====================================================================================================================
    % Prepare index partion of layers
    ProcTE.stfn_layers=zeros(size(CompStructA.Model.DomainType,2),2);
    for ilayer=1:size(CompStructA.Model.DomainType,2)
        if strcmp(CompStructA.Model.DomainType(ilayer),'fluid')
           ivar=length(FEMatricesA.DNodesComp{ilayer});  % potential
        end
        if strcmp(CompStructA.Model.DomainType(ilayer),'HTTI')
           ivar=3*length(FEMatricesA.DNodesComp{ilayer}); % displacement
        end
        if ilayer==1
           ProcTE.stfn_layers(1,1)=1; ProcTE.stfn_layers(1,2)=ivar;
        else
           ProcTE.stfn_layers(ilayer,1)=ProcTE.stfn_layers(ilayer-1,2)+1;
           ProcTE.stfn_layers(ilayer,2)=ProcTE.stfn_layers(ilayer,1)+ivar-1;
        end
    end
    
    %====================================================================================================================
    % Calculate trigonometrical Fourier decomposition of ur,uf,uz and save in LayerTE-freq-layer.mat file
    % Prepare index partion of layers
    for ilayer=ProcTE.nlayers % layer cycle
    %for ilayer=1 % layer cycle
        ProcTE.ilayer=ilayer;
        fprintf(1,'\tNow fhe %d-th layer is processed out of %d (real %d) layers\n', ...
                ilayer, length(ProcTE.nlayers), length(ProcTE.nlayers_count));

        % Fluid layer
        if strcmp(CompStructA.Model.DomainType(ilayer),'fluid')
           cd('spectrum\');
           [ResFluid]=St61_Proc_Fluid_TE_sp_SAFE(CompStructA,FEMatricesA,ProcTE);
           cd('..'); %return;
           if lower(ProcTE.write_result)=='y'
              fname = char(strcat(ProcTE.Dir_Name,'\','LayerTE','-',num2str(ProcTE.freq_count(ii_f)),'-',num2str(ilayer)','.mat'));
              ResFluid.ProcTE=ProcTE;
              save(fname,'-struct','ResFluid','ProcTE','rr','TE_fer','cF_ur_rm_pm','cF_uf_rm_pm','cF_uz_rm_pm');
           end
        end

        % Solid layer
        if strcmp(CompStructA.Model.DomainType(ilayer),'HTTI')
           cd('spectrum\');
           [ResSolid]=St61_Proc_Solid_TE_sp_SAFE(CompStructA,FEMatricesA,ProcTE);
           cd('..'); %return;
           if lower(ProcTE.write_result)=='y'
              fname = char(strcat(ProcTE.Dir_Name,'\','LayerTE','-',num2str(ProcTE.freq_count(ii_f)),'-',num2str(ilayer)','.mat'));
              ResSolid.ProcTE=ProcTE;
              save(fname,'-struct','ResSolid','ProcTE','rr','TE_fer','cF_ur_rm_pm','cF_uf_rm_pm','cF_uz_rm_pm');
           end
        end

    end %  layer cycle

    fprintf('\tTime for this step = %.1fs\n',toc(tStart));
end % frequency cycle

%====================================================================================================================
% Save frequency dependencies of velocity and slowness in Slowness.mat file
fname_slow = char(strcat(ProcTE.Dir_Name,'\Slowness','.mat'));
save(fname_slow,'Result_Freq','Result_Velocity','Result_Slowness','Result_Attenuation');
cd(curr_dir);

fprintf('Time for Proc_Aniso TE = %.1fs\n',toc(tStart_Prog));
fprintf('============================================================================\n');
zzz=1; % convenient string for break point
