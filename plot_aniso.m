% plot_aniso.m only plots results specified below
% no inputs for this script
clear all;
Dir_Name={'Mesaverde-0-F'};
% Dir_Name={'Mesaverde-HTI-F'};
% Dir_Name={'Mesaverde-HTI-F1'};
% Dir_Name={'Mesaverde-HTI-F2'};
%Dir_Name={'TEST-BS'};

% Below 'y' corresponds to yes, 'n' corresponds to 'no'

% Slowness is plotted for all frequencies and eigenvalues in range (ms/f)
PlotAniso.vel_slo_attn_plot={'n','n','n'}; PlotAniso.vel_slo_attn_limits={[0 3000]; [50 300]; [0 1]};

% Plot input parameters of different maps, specified below: (freq) is frequency in kHz, (layer number), (neigs) is array of eigenvalues
PlotAniso.freq=9.5; PlotAniso.num_layer=2;
PlotAniso.neigs=[1:10 ];
%PlotAniso.neigs=[1:25];

% Plot maps of: 
% Absolute values of displacements in Cartesian coordinates (|ux|, |uy|, |uz|) or Cylindrical coordinates (|ur|, |uf|);
% Kinetic energy (te), potential energy (ue), kinetic energy minus potential energy (te), 
% Absolute values of components of Poyinting vector in Cartesian coordinates (|px|, |py|, |pz|) or  Cylindrical coordinates (|pr|, |pf})
% Absolute value of pressure in fluid or spherical part of full stress tensor in solid (|pp|); its polar angle (ap);
PlotAniso.fields={'ux', 'n'; 'uy', 'n'; 'uz', 'n'; 'ur', 'n'; 'uf', 'n'; 'te', 'y'; 'ue', 'n'; 'tu', 'n';
                  'px', 'n'; 'py', 'n'; 'pz', 'n'; 'pr', 'n'; 'pf', 'n'; 'pp', 'n'; 'ap', 'n' };

% Plot maps in 3D dimension (plot_3d) by setting the angle of the view (view_3D) and normalizing results(normalization)
PlotAniso.plot_3d='y'; PlotAniso.view_3D=[0, 90]; PlotAniso.normalization='n';

% Use contour plot to present results (plot_c), using (num_xy) points in x,y directions and (num_cont) contour lines
PlotAniso.plot_c='n'; PlotAniso.num_xy=[100; 100]; PlotAniso.num_cont=50;

% Use subplot possibility to plot maps (with exception of slowness results)
% n (rows) x m (columns) parameters are chosen automatically
PlotAniso.subplot='y'; 

% Figure parameters
PlotAniso.close_figs='n'; % close all figures
PlotAniso.num_fig=503; % start figure number
PlotAniso.font_size=8; % base fontname size
PlotAniso.linestyle=['o'; '-'; '-' ]; % styles of 3 lines
PlotAniso.linewidth=[1.; 1; 1];       % their widths
PlotAniso.color=['r','r','b'];        % their colors
PlotAniso.marker=['o'; 'o'; 'o']; % 3 markers
PlotAniso.markersize=[4; 10; 10]; % their size
PlotAniso.xlabel_show='n'; % show xlabel, text beside the x-axis 
   
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
%ColorOrder=[0 0 1; 0 0.5 0; 1 0 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25];
ColorOrder=[0 0 0];

PlotAniso.Dir_Name=Dir_Name;
% Current directory
curr_dir=pwd;
%run_prog=mfilename;
%full_path_run_prog=mfilename('fullpath');
%run_dir=full_path_run_prog(1:(length(full_path_run_prog)-length(run_prog)-1));
%cd(run_dir);

% Full input durectory name
text_var=char(PlotAniso.Dir_Name); key_full_dir='n';
if text_var(2)==':'
   PlotAniso.Full_Dir_Name=PlotAniso.Dir_Name;
else
   PlotAniso.Full_Dir_Name=strcat(pwd,'\',PlotAniso.Dir_Name);
end

% Read CompStruct and FEMatrices data 
fname = strcat(char(PlotAniso.Full_Dir_Name),'\CompStruct.mat');
if exist(fname,'file')>0
    load(fname);
    CompStructA.Config=Config;
    CompStructA.Methods=Methods;
    CompStructA.Model=Model;
    CompStructA.Advanced=Advanced;
    CompStructA.Misc=Misc;
    CompStructA.Data=Data;
    CompStructA.Asymp=Asymp;
else
    beep; beep; beep; fprintf(1,'File is not existed!\n'); return; 
end

PlotAniso.freq_count=CompStructA.Model.f_array;
find_freq='n';
for iff=1:CompStructA.Model.N_disp
    if abs(CompStructA.Model.f_array(iff)-PlotAniso.freq)<1.e-10
       PlotAniso.num_freq=iff; ifreq=iff; find_freq='y'; break;
    end
end
if find_freq=='n'
   beep; beep; beep; fprintf(1,'PlotAniso.freq did not find!\n'); return;
end

fname = strcat(char(PlotAniso.Dir_Name),'\FEMatrices-',num2str(CompStructA.Model.f_array(PlotAniso.num_freq)),'.mat');
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
%====================================================================================================================
%====================================================================================================================
% Build-in variables
% Number of layers
PlotAniso.nlayers=1:size(CompStructA.Model.DomainType,2);  
if ~ismember(PlotAniso.num_layer,PlotAniso.nlayers)
    beep; beep; beep; fprintf(1,'Error in PlotAniso.num_layer!\n'); return;
end

PlotAniso.max_neigs=CompStructA.Advanced.num_eig_max; % number of calculated eigs
if PlotAniso.neigs(length(PlotAniso.neigs))>PlotAniso.max_neigs
    beep; beep; beep; fprintf(1,'Error in PlotAniso.neigs!\n'); return;
end
    
if lower(PlotAniso.close_figs)=='y', close all; end

% Here main subprograms
CompStructA.Config.root_path=strcat(pwd,'\');
cd(strcat(CompStructA.Config.root_path,'spectrum\'));
PlotAniso.dir_input=PlotAniso.Dir_Name; % input durectory, may be anywhere

% Prepare index partion of layers
PlotAniso.stfn_layers=zeros(size(CompStructA.Model.DomainType,2),2);

for ilayer=1:size(CompStructA.Model.DomainType,2)
    if strcmp(CompStructA.Model.DomainType(ilayer),'fluid')
       ivar=length(FEMatricesA.DNodesComp{ilayer});  % potential
    end
    if strcmp(CompStructA.Model.DomainType(ilayer),'HTTI')
       ivar=3*length(FEMatricesA.DNodesComp{ilayer}); % displacement
    end
    if ilayer==1
       PlotAniso.stfn_layers(1,1)=1; PlotAniso.stfn_layers(1,2)=ivar;
    else
       PlotAniso.stfn_layers(ilayer,1)=PlotAniso.stfn_layers(ilayer-1,2)+1;
       PlotAniso.stfn_layers(ilayer,2)=PlotAniso.stfn_layers(ilayer,1)+ivar-1;
    end
end

%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
% Read eigs, calculate  velocities, slownesses and attenuations
[PlotAniso.Velocity,PlotAniso.Slowness,PlotAniso.Atten]=St51_Calc_VS_sp_SAFE(CompStructA,FEMatricesA,PlotAniso);

PlotAniso.SA=zeros(length(PlotAniso.max_neigs),2);
for im=1:PlotAniso.max_neigs
    PlotAniso.SA(im,:)=[PlotAniso.Slowness(im) PlotAniso.Atten(im)];
end
PlotAniso.SA=PlotAniso.SA';

% Plot phase velocities and slownesses for all calculated eigs!!!
if lower((PlotAniso.vel_slo_attn_plot{1}))=='y' % velocity
   PlotAniso.freq_count=CompStructA.Model.f_array;
   f_max=max(CompStructA.Model.f_array);

   % Plot velocities
   hf=figure(PlotAniso.num_fig+12); 
   set(gca,'FontSize',PlotAniso.font_size); %CompStructA.Model.f_max=10; PlotAniso.freq_count=1:10;
   set(gca,'ColorOrder',ColorOrder);
   hold on; box on; xlim([0  f_max+1]); 
   varv=0.0254*12e6./PlotAniso.vel_slo_attn_limits{1}; ylim(PlotAniso.vel_slo_attn_limits{1});
   plot(PlotAniso.freq_count,0.0254*12e6./PlotAniso.Slowness(1:length(PlotAniso.freq_count),1:PlotAniso.max_neigs), ...
       'Marker',PlotAniso.marker(1),'Markersize',PlotAniso.markersize(1), ...
       'LineStyle','none','Linewidth',PlotAniso.linewidth(1));

   freq_var=[0 PlotAniso.freq_count f_max+1];
   VP_val=ones(length(freq_var),1).*(CompStructA.Asymp.V_qP); 
   VSV_val=ones(length(freq_var),1).*(CompStructA.Asymp.V_qSV); 
   VSH_val =ones(length(freq_var),1).*(CompStructA.Asymp.V_SH);
   VSt_val =ones(length(freq_var),1).*(CompStructA.Asymp.V_St);
   
   plot(freq_var,VP_val*1e3,'LineStyle',PlotAniso.linestyle(2),'Linewidth',PlotAniso.linewidth(2),'Color',PlotAniso.color(2));
   plot(freq_var,VSV_val*1e3,'LineStyle',PlotAniso.linestyle(2),'Linewidth',PlotAniso.linewidth(2),'Color',PlotAniso.color(2));
   plot(freq_var,VSH_val*1e3,'LineStyle',PlotAniso.linestyle(3),'Linewidth',PlotAniso.linewidth(2),'Color',PlotAniso.color(3));
   xlabel('f (kHz)'); ylabel('Velocity m/s'); title('');
end
if lower((PlotAniso.vel_slo_attn_plot{2}))=='y' % slowness
   PlotAniso.freq_count=CompStructA.Model.f_array;
   f_max=max(CompStructA.Model.f_array);

   hf=figure(PlotAniso.num_fig+13); 
   set(gca,'FontSize',PlotAniso.font_size); %CompStructA.Model.f_max=10; PlotAniso.freq_count=1:10;
   set(gca,'ColorOrder',ColorOrder);
   hold on; box on; xlim([0  f_max+1]); ylim(PlotAniso.vel_slo_attn_limits{2});
   plot(PlotAniso.freq_count,PlotAniso.Slowness(1:length(PlotAniso.freq_count),1:PlotAniso.max_neigs), ...
       'Marker',PlotAniso.marker(1),'Markersize',PlotAniso.markersize(1), ...
       'LineStyle','none','Linewidth',PlotAniso.linewidth(1));

   freq_var=[0 PlotAniso.freq_count f_max+1];
   qP_val=ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_qP); 
   qSV_val=ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_qSV); 
   SH_val =ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_SH);
   St_val =ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_St);
   
   plot(freq_var,qP_val,'LineStyle',PlotAniso.linestyle(2),'Linewidth',PlotAniso.linewidth(2),'Color',PlotAniso.color(2));
   plot(freq_var,qSV_val,'LineStyle',PlotAniso.linestyle(2),'Linewidth',PlotAniso.linewidth(2),'Color',PlotAniso.color(2));
   plot(freq_var,SH_val ,'LineStyle',PlotAniso.linestyle(3),'Linewidth',PlotAniso.linewidth(2),'Color',PlotAniso.color(3));
   xlabel('f (kHz)'); ylabel('Slowness (\mus/f)'); title('');

end
if lower((PlotAniso.vel_slo_attn_plot{3}))=='y' % attenuation
   PlotAniso.freq_count=CompStructA.Model.f_array;
   f_max=max(CompStructA.Model.f_array);

   % Plot slownesses
   hf=figure(PlotAniso.num_fig+14); 
   set(gca,'FontSize',PlotAniso.font_size); 
   hold on; box on; xlim([0  f_max+1]); ylim(PlotAniso.vel_slo_attn_limits{3});
   PlotAniso.max_neigs_var=PlotAniso.max_neigs; %PlotAniso.max_neigs_var=5;
   plot(PlotAniso.freq_count,PlotAniso.Atten(1:length(PlotAniso.freq_count),1:PlotAniso.max_neigs_var), ...
       'Marker',PlotAniso.marker(1),'Markersize',PlotAniso.markersize(1), ...
       'LineStyle','none','Linewidth',PlotAniso.linewidth(1));

   xlabel('f (kHz)'); ylabel('1/Q'); title('');
end

%====================================================================================================================
%====================================================================================================================
if lower(PlotAniso.plot_3d)=='y' || lower(PlotAniso.plot_c)=='y'
   fprintf(1,'\n\n\n================================================\n');
   fprintf(1,'Now fhe %d-th frequency (%.2f kHz) is processed out of %d\n', ...
           PlotAniso.num_freq, PlotAniso.freq_count(PlotAniso.num_freq),length(PlotAniso.freq_count));

   % n (rows) x m (columns) parameters of subplot
   PlotAniso.subplot_nxm=[1 1];
   if lower(PlotAniso.subplot)=='y'
      vari=length(PlotAniso.neigs);
      
      if vari==1, n=1; m=1; end
      if vari==2, n=1; m=2; end
      if 3<=vari && vari<=4, n=2; m=2; end
      if 5<=vari && vari<=8, n=2; m=4; end
      if 9<=vari && vari<=12, n=3; m=4; end
      if 13<=vari && vari<=15, n=3; m=5; end
      if 16<=vari && vari<=18, n=3; m=6; end
      if 19<=vari && vari<=24, n=4; m=6; end
      if 25==vari, n=5; m=5; end
      if 25<vari && vari<=30, n=5; m=6; end
      if 31<=vari && vari<=35, n=5; m=7; end
      if 36<=vari && vari<=42, n=6; m=7; end
      if 43<=vari && vari<=48, n=6; m=8; end
      if 49<=vari && vari<=56, n=7; m=8; end
      if 57<=vari && vari<=63, n=7; m=9; end
      if 64<=vari && vari<=72, n=8; m=9; end
      if 73<=vari && vari<=80, n=8; m=10; end
      if 81<=vari && vari<=90, n=9; m=10; end
      if 91<=vari && vari<=99, n=9; m=11; end
      if 100<=vari && vari<=120, n=10; m=12; end
      if 121<=vari && vari<=132, n=11; m=12; end
      if 133<=vari && vari<=144, n=12; m=12; end
      if 145<=vari && vari<=156, n=12; m=13; end
      if 157<=vari && vari<=169, n=13; m=13; end
      if 170<=vari && vari<=182, n=13; m=14; end
      if 183<=vari && vari<=196, n=14; m=14; end
      if 197<=vari && vari<=210, n=14; m=15; end

      PlotAniso.subplot_nxm=[n m];
   end
   
   PlotAniso.iii_eig=0;
   % cycle on eigenvalues
   for ieig=PlotAniso.neigs
       PlotAniso.iii_eig=PlotAniso.iii_eig+1;
       PlotAniso.ieig=ieig;
           
       % Fuild: Calculation ux, uy, uz, ur, uf, T energy, U energy, T-P, 
       % poynting vector (px, py, pz, pr, pf), P and angle(P)
       if strcmp(CompStructA.Model.DomainType(PlotAniso.num_layer),'fluid')
          [ResFluid]=St52_Calc_Fluid_sp_SAFE(CompStructA,FEMatricesA,PlotAniso);
           St52_Plot_Fluid_sp_SAFE(CompStructA,FEMatricesA,PlotAniso,ResFluid);
       end

       % HTTI: Calculation ux, uy, uz, ur, ufi, T energy, U energy, T-P,
       % stress tensor, poynting vector (px, py, pz, pr, pf), spur and spur angle
       if strcmp(CompStructA.Model.DomainType(PlotAniso.num_layer),'HTTI')
          [ResSolid]=St52_Calc_Solid_sp_SAFE(CompStructA,FEMatricesA,PlotAniso);
          St52_Plot_Solid_sp_SAFE(CompStructA,FEMatricesA,PlotAniso,ResSolid);
       end
       
   end % cycle on eigenvalues

end % plot_3d or c

% Return in current directory
cd(curr_dir); colormap('jet');
zzz=1;
 

