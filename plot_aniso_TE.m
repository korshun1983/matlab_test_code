% plot_aniso_TE() plots various characteristics of the kinetic energy
% no inputs for this script
clear all;

% Dir_Name={'Mesaverde-HTI};
Dir_Name={'Mesaverde-30-F'};
% Dir_Name={'Mesaverde-30-F1'};
% Dir_Name={'Mesaverde-HTI-F'};
% Dir_Name={'Mesaverde-HTI-F1'};

% Below 'y' corresponds to yes, 'n' corresponds to 'no'

% What can be plotted for each specified eigenvalue (eig) and specified frequency (f):
% 1) Radial distributions (r) of kinetic energy, TE(f,eig,r,h), of all harmonic (h)
% 2) Radial distributions (r) of ITE(f,eig,r,h)=Integral_{r=r_min}^r r1*TE(f,eig,r1,h)dr1, 
%    normalized to max ITE(f,eig,r_max, h) in harmonics.
% 3) Radial distributions (r) of ITE(f,eig,r)=Sum_{h=-NHarm}^{h=NHarm}ITE(f,eig,r,h), 
%    normalized to max TE(f,eig,r_max).
% 4) Radial ddistributions (r) of kinetic energy, TE_f_e_Ir(f,eig,r), obtained by integrating 
%    kinetic energy density with respect to the angle, normalized to TE_f_e_Ir(f,eig,r_max)
% Number of harmonics is specified in proc_aniso.m, TE(f,eig,r,h) and TE_f_e_Ir(f,eig,r) 
% are also obtained by this program.
PlotTE.result=['n'; 'y'; 'y'; 'y';]; 

% Plot for these numbers of layers, frequency in kHz and array of eigenvalues for plotting 
PlotTE.nlayers=[1:2]; PlotTE.freq=8.25;  PlotTE.neigs=[1:5];

% Figure parameters
PlotTE.subplot='y'; % Use subplot possibility; n (rows) x m (columns) parameters are choosen automatically
PlotTE.close_figs='y'; % close all figures
PlotTE.num_fig=2300; % start figure number 
PlotTE.font_size=8; % base fontname size 
PlotTE.color=['r'; 'b'; 'k'; 'm'; ] ; % colors: r is 0th harmonic or 1st line, b is +- 1st harmonica or 2nd line, etc.
PlotTE.linewidth=[1.; 1; 1; 1; 1] ;   % line width: the same 
PlotTE.linestyle={'-'; ':'}; % line style: '-' is nonnegative harmonics, ':' is megative harmonics
PlotTE.xlabel_show='y'; % show xlabel, text beside the x-axis 
PlotTE.rl_ratio='y'; % x-axis is given in r/lambda, lambda is the wave length corresponded to eig

%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
PlotAniso.Dir_Name=Dir_Name;
curr_dir=pwd;
%run_prog=mfilename;
%full_path_run_prog=mfilename('fullpath');
%run_dir=full_path_run_prog(1:(length(full_path_run_prog)-length(run_prog)-1));
%cd(run_dir);

% Read CompStruct and FEMatrices data 
fname = char(strcat(PlotAniso.Dir_Name,'\CompStruct.mat'));
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

find_freq='n';
for iff=1:CompStructA.Model.N_disp
    if CompStructA.Model.f_array(iff)==PlotTE.freq
       PlotTE.nfreq=iff; ifreq=iff; find_freq='y'; break;
    end
end
if find_freq=='n'
   beep; beep; beep; fprintf(1,'PlotTE.freq did not find!\n'); return;
end

fname = strcat(char(PlotAniso.Dir_Name),'\FEMatrices-',num2str(CompStructA.Model.f_array(PlotTE.nfreq)),'.mat');
if exist(fname,'file')>0
    load(fname);
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
% Read results of proc_aniso.m program for each layer
PlotTE.nlayers_real=1:size(CompStructA.Model.DomainType,2);
ResLayer=cell(length(PlotTE.nlayers_real),1);
for ilayer=PlotTE.nlayers_real
    fname = strcat(char(PlotAniso.Dir_Name),'/','LayerTE','-',num2str(CompStructA.Model.f_array(PlotTE.nfreq)),'-',num2str(ilayer)','.mat');
    if exist(fname,'file')>0
        load(fname);
    else
        beep; beep; beep; fprintf(1,'File is not existed!\n'); return;
    end
    ResLayer{ilayer}.ProcTE=ProcTE;
    ResLayer{ilayer}.rr=rr;
    ResLayer{ilayer}.TE_f_e_r=TE_fer;
    ResLayer{ilayer}.cF_ur_rm_pm=cF_ur_rm_pm;
    ResLayer{ilayer}.cF_uf_rm_pm=cF_uf_rm_pm;
    ResLayer{ilayer}.cF_uz_rm_pm=cF_uz_rm_pm;
    zzz=1;
end

% Check input layers
for ilayer=PlotTE.nlayers
    if ~ismember(ilayer,PlotTE.nlayers_real)
       beep; beep; beep; fprintf(1,'Error in PlotTE.nlayers!\n'); return;
    end
end

%====================================================================================================================
%====================================================================================================================
% Read Velocities and Slownesses results
fname = strcat(char(PlotAniso.Dir_Name),'\Slowness.mat');
if exist(fname,'file')>0
    load(fname);
    PlotTE.freq_count=Result_Freq;
    PlotTE.Velocity=Result_Velocity;
    PlotTE.Slowness=Result_Slowness;
end

%====================================================================================================================
%====================================================================================================================
% Check input eigenvalues' array
PlotTE.max_neigs=CompStructA.Advanced.num_eig_max; % number of calculated eigs
PlotTE.neigs_count=1:PlotTE.max_neigs;
for ieig=PlotTE.neigs
    if ~ismember(ieig,PlotTE.neigs_count)
       beep; beep; beep; fprintf(1,'Error in PlotTE.neigs!\n'); return;
    end
end
    
if lower(PlotTE.close_figs)=='y', close all; end

%====================================================================================================================
%====================================================================================================================
% Set n (rows) x m (columns) parameters of subplot

PlotTE.subplot_nxm=[1 1];
if lower(PlotTE.subplot)=='y'
   vari=length(PlotTE.neigs);
   
   if vari==1, n=1; m=1; end
   if vari==2, n=1; m=2; end
   if 3<=vari && vari<=4, n=2; m=2; end
   if 5<=vari && vari<=8, n=2; m=4; end
   if 9<=vari && vari<=12, n=3; m=4; end
   if 13<=vari && vari<=15, n=3; m=5; end
   if 16<=vari && vari<=18, n=3; m=6; end
   if 19<=vari && vari<=24, n=4; m=6; end
   if 25<=vari && vari<=30, n=5; m=6; end
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

   PlotTE.subplot_nxm=[n m];
end

%====================================================================================================================
%====================================================================================================================
% Build r-axis 
NrrF=0; rrF=[];
for ilayer=PlotTE.nlayers
    NrrF=NrrF+size(ResLayer{ilayer}.rr,1);
    rrF=[rrF(1:end); ResLayer{ilayer}.rr];
end

%====================================================================================================================
%====================================================================================================================
% Prepare output arrays
NHarm=ResLayer{ilayer}.ProcTE.NHarm;
TE_f_eigs_r_h_p=zeros(length(PlotTE.neigs_count),NrrF,NHarm);
TE_f_eigs_r_h_m=zeros(length(PlotTE.neigs_count),NrrF,NHarm);
TE_f_eigs_Ir_h_p=zeros(length(PlotTE.neigs_count),NrrF-1,NHarm);
TE_f_eigs_Ir_h_m=zeros(length(PlotTE.neigs_count),NrrF-1,NHarm);
TE_f_eigs_Ir=zeros(length(PlotTE.neigs_count),NrrF-1);

%================================================================================================================
%================================================================================================================
%================================================================================================================
% 1) Calculate Kinetic Energy for specified f, eigenvalue, r and harmonic, TE(f,eig,r,h)
for ieig=PlotTE.neigs % cycle on eigenvalues
    for iharm=1:NHarm % cycle on harmonics
        TE_fer_h_p=[];
        TE_fer_h_m=[];
        for ilayer=PlotTE.nlayers
            varTE = 0.5.*FEMatricesA.PhysProp{ilayer}.rho*1000. ...
                   *(2.*pi.*PlotTE.freq.*CompStructA.Misc.F_conv).^2;

            var_p=  abs(ResLayer{ilayer}.cF_ur_rm_pm(ieig,:,iharm,1)).^2 ...
                   +abs(ResLayer{ilayer}.cF_uf_rm_pm(ieig,:,iharm,1)).^2 ...
                   +abs(ResLayer{ilayer}.cF_uz_rm_pm(ieig,:,iharm,1)).^2;
            var_p=varTE*var_p;

            var_m=  abs(ResLayer{ilayer}.cF_ur_rm_pm(ieig,:,iharm,2)).^2 ...
                   +abs(ResLayer{ilayer}.cF_uf_rm_pm(ieig,:,iharm,2)).^2 ...
                   +abs(ResLayer{ilayer}.cF_uz_rm_pm(ieig,:,iharm,2)).^2;
            var_m=varTE*var_m;

            TE_fer_h_p=[TE_fer_h_p(1:end); squeeze(var_p)'];
            TE_fer_h_m=[TE_fer_h_m(1:end); squeeze(var_m)'];
        end
        TE_f_eigs_r_h_p(ieig,:,iharm)=TE_fer_h_p;
        TE_f_eigs_r_h_m(ieig,:,iharm)=TE_fer_h_m;
    end % cycle on harmonics
end % cycle on eigenvalues

% plot the above calculated result
if lower(PlotTE.result(1))=='y'
   var_iii_eig=0;
   for ieig=PlotTE.neigs % cycle on eigenvalues
       % calculate figure number
       var_iii_eig=var_iii_eig+1; win_eig=1;
       num_fig=PlotTE.num_fig+100*(ilayer-1)+10*(ifreq-1)+(var_iii_eig-1);
       if lower(PlotTE.subplot)=='y' 
           num_fig=PlotTE.num_fig+100*(ilayer-1)+10*(ifreq-1);
           win_eig=var_iii_eig;
       end

       % build r-axis as r/lambda, lambda is the wave length
       lambda=(PlotTE.Velocity(ifreq,ieig)./1000)/PlotTE.freq;
       rrFL=rrF; if lower(PlotTE.rl_ratio)=='y', rrFL=rrFL/lambda; end

       hf=figure(num_fig);
       subplot(PlotTE.subplot_nxm(1),PlotTE.subplot_nxm(2),win_eig);
       set(gca,'FontSize',PlotTE.font_size);
       hold on; box on; 
       for iharm=2:NHarm % negative harmonics
           plot(rrFL,squeeze(TE_f_eigs_r_h_m(ieig,:,iharm)), ...
               'Linewidth',PlotTE.linewidth(iharm),'Color',PlotTE.color(iharm), ...
               'LineStyle',PlotTE.linestyle{2});
       end
       for iharm=1:NHarm % zero and positive harmonics
           plot(rrFL,squeeze(TE_f_eigs_r_h_p(ieig,:,iharm)), ...
               'Linewidth',PlotTE.linewidth(iharm),'Color',PlotTE.color(iharm), ...
               'LineStyle',PlotTE.linestyle{1});
       end
       if lower(PlotTE.xlabel_show)=='y'
          if lower(PlotTE.rl_ratio)=='y', xlabel('r/\lambda'); else xlabel('r (m)'); end
       end
       title(strcat('f=',num2str(PlotTE.freq),' kHz',', k_{',num2str(ieig),'}'));
   end % cycle on eigenvalues
end

%================================================================================================================
%================================================================================================================
%================================================================================================================
% 2) Calculate TE_f_eigs_Ir_h= Sum_{r=r_min}^{r} r1*TE(f,eig,r1,h) dr1
    
for ieig=PlotTE.neigs % cycle on eigenvalues
    for iharm=1:NHarm % cycle on harmonisc

        varl=squeeze(TE_f_eigs_r_h_p(ieig,1:NrrF-1,iharm))';
        varr=squeeze(TE_f_eigs_r_h_p(ieig,2:NrrF  ,iharm))';
        var_mas=0.5.*( rrF(1:NrrF-1).*varl+rrF(2:NrrF).*varr ).*(rrF(2:NrrF)-rrF(1:NrrF-1));
        TE_f_eigs_Ir_h_p(ieig,1,iharm)=var_mas(1);
        for ir=2:NrrF-1;
            TE_f_eigs_Ir_h_p(ieig,ir,iharm)= ...
               TE_f_eigs_Ir_h_p(ieig,ir-1,iharm)+var_mas(ir);
        end

        varl=squeeze(TE_f_eigs_r_h_m(ieig,1:NrrF-1,iharm))';
        varr=squeeze(TE_f_eigs_r_h_m(ieig,2:NrrF  ,iharm))';
        var_mas=0.5.*( rrF(1:NrrF-1).*varl+rrF(2:NrrF).*varr ).*(rrF(2:NrrF)-rrF(1:NrrF-1));
        TE_f_eigs_Ir_h_m(ieig,1,iharm)=var_mas(1);
        for ir=2:NrrF-1
            TE_f_eigs_Ir_h_m(ieig,ir,iharm)= ...
               TE_f_eigs_Ir_h_m(ieig,ir-1,iharm)+var_mas(ir);
        end

    end % cycle on harmonisc
end % cycle on eigenvalues

% plot normalized results for r=r_min, ..., r_max: 
% NTE_f_eigs_Ir_h= TE_f_eigs_Ir_h/(max of TE_f_eigs_Ir_h in harminics (h))
if lower(PlotTE.result(2))=='y'
    
   var_iii_eig=0;
   for ieig=PlotTE.neigs % cycle on eigenvalues
       % calculate figure number
       var_iii_eig=var_iii_eig+1; win_eig=1;
       num_fig=PlotTE.num_fig+1000+100*(ilayer-1)+10*(ifreq-1)+(var_iii_eig-1);
       if lower(PlotTE.subplot)=='y'
           num_fig=PlotTE.num_fig+1000+100*(ilayer-1)+10*(ifreq-1);
           win_eig=var_iii_eig;
       end

       % build r-axis as r/lambda, lambda is the wave length
       lambda=(PlotTE.Velocity(ifreq,ieig)./1000)/PlotTE.freq;
       rrFL=rrF; if lower(PlotTE.rl_ratio)=='y', rrFL=rrFL/lambda; end

       varmax=max(max(TE_f_eigs_Ir_h_p(ieig,NrrF-1,:)),max(TE_f_eigs_Ir_h_m(ieig,NrrF-1,:)));

       hf=figure(num_fig);
       subplot(PlotTE.subplot_nxm(1),PlotTE.subplot_nxm(2),win_eig);
       set(gca,'FontSize',PlotTE.font_size);
       hold on; box on;
       for iharm=2:NHarm % negative harmonics
           plot(rrFL(1:NrrF-1),squeeze(TE_f_eigs_Ir_h_m(ieig,:,iharm))./varmax, ...
               'Linewidth',PlotTE.linewidth(iharm),'Color',PlotTE.color(iharm), ...
               'LineStyle',PlotTE.linestyle{2});
       end
       for iharm=1:NHarm % zero and positive harmonics
           plot(rrFL(1:NrrF-1),squeeze(TE_f_eigs_Ir_h_p(ieig,:,iharm))./varmax, ...
               'Linewidth',PlotTE.linewidth(iharm),'Color',PlotTE.color(iharm), ...
               'LineStyle',PlotTE.linestyle{1});
       end
       if lower(PlotTE.xlabel_show)=='y'
          if lower(PlotTE.rl_ratio)=='y', xlabel('r/\lambda'); else  xlabel('r (m)'); end
       end
       title(strcat('f=',num2str(PlotTE.freq),' kHz',', k_{',num2str(ieig),'}'));
   end % cycle on eigenvalues

end 

%================================================================================================================
%================================================================================================================
%================================================================================================================
% 3) Calculate TE_f_eigs_Ir=Sum_{h=-NHarm}^{h=NHarm} (TE_f_eigs_Ir_h) and normalized results
if lower(PlotTE.result(3))=='y'

    for ieig=PlotTE.neigs % cycle on eigenvalues
       for ir=1:NrrF-1 % cycle on harmonics
           var_array=squeeze(TE_f_eigs_Ir_h_p(ieig,ir,:))+squeeze(TE_f_eigs_Ir_h_m(ieig,ir,:));
           var_array(1)=var_array(1)/2.;
           TE_f_eigs_Ir(ieig,ir)=sum(var_array);
       end
   end

   % Plot NTE_f_eigs_Ir=TE_f_eigs_Ir/TE_f_eigs_Ir(r_max)
   var_iii_eig=0;
   for ieig=PlotTE.neigs % cycle on eigenvalues
       % calculate figure number
       var_iii_eig=var_iii_eig+1; win_eig=1;
       num_fig=PlotTE.num_fig+3000+100*(ilayer-1)+10*(ifreq-1)+(var_iii_eig-1);
       if lower(PlotTE.subplot)=='y'
           num_fig=PlotTE.num_fig+3000+100*(ilayer-1)+10*(ifreq-1);
           win_eig=var_iii_eig;
       end

       % build r-axis as r/lambda, lambda is the wave length
       lambda=(PlotTE.Velocity(ifreq,ieig)./1000)/PlotTE.freq;
       rrFL=rrF; if lower(PlotTE.rl_ratio)=='y', rrFL=rrFL/lambda; end

       hf=figure(num_fig);
       subplot(PlotTE.subplot_nxm(1),PlotTE.subplot_nxm(2),win_eig);
       set(gca,'FontSize',PlotTE.font_size);
       hold on; box on; 
       plot(rrFL(1:NrrF-1),squeeze(TE_f_eigs_Ir(ieig,:))./TE_f_eigs_Ir(ieig,NrrF-1), ...
            'Linewidth',PlotTE.linewidth(iharm),'Color',PlotTE.color(2));
       if lower(PlotTE.xlabel_show)=='y'
          if lower(PlotTE.rl_ratio)=='y', xlabel('r/\lambda'); else   xlabel('r (m)'); end
       end
       title(strcat('f=',num2str(PlotTE.freq),' kHz',', k_{',num2str(ieig),'}'));
    end % cycle on eigenvalues

end % 3)

%================================================================================================================
%================================================================================================================
%================================================================================================================
% 4) Calculate TE_f_e_Ir(f,eig,r)=sum_{r=r_min}^{r} r1*TE_fer(f,eig,r1) dr1,
%    TE_fer(f,eig,r) was calculated proc_aniso.m by integrating kinetic enerdy density with respect to the angle.
%    NTE_f_e=TE_f_e_Ir(f,eig,r)/TE_f_e_Ir(f,eig,r_max)

if lower(PlotTE.result(4))=='y'

   % Prepare output arrays
   TE_f_e_Ir =zeros(length(PlotTE.neigs_count),NrrF-1);
   TE_f_e_NIr=zeros(length(PlotTE.neigs_count),NrrF-1);
   TE_ferA=zeros(length(PlotTE.neigs_count),NrrF);

   for ieig=PlotTE.neigs % cycle on eigenvalues

       TE_fer =[];
       for ilayer=PlotTE.nlayers
           TE_fer =[TE_fer(1:end); squeeze(ResLayer{ilayer}.TE_f_e_r(ieig,:))'];
       end
       TE_ferA(ieig,:)=TE_fer;

       varl=(TE_fer(1:NrrF-1));
       varr=(TE_fer(2:NrrF  ));
       var_mas=0.5.*( rrF(1:NrrF-1).*varl+rrF(2:NrrF).*varr ).*(rrF(2:NrrF)-rrF(1:NrrF-1));
       TE_f_e_Ir(ieig,1)=var_mas(1);
       for ir=2:NrrF-1
           TE_f_e_Ir(ieig,ir)= ...
               TE_f_e_Ir(ieig,ir-1)+var_mas(ir);
       end

       TE_f_e_NIr(ieig,:)= ...
           TE_f_e_Ir(ieig,:)./TE_f_e_Ir(ieig,NrrF-1);

   end % cycle on eigenvalues

   TE_test=zeros(length(PlotTE.neigs_count),1);
   var_iii_eig=0;
   for ieig=PlotTE.neigs % cycle on eigenvalues
        
       % calculate figure number
       var_iii_eig=var_iii_eig+1; win_eig=1;
       num_fig=PlotTE.num_fig+4000+100*(ilayer-1)+10*(ifreq-1)+(var_iii_eig-1);
       if lower(PlotTE.subplot)=='y'
           num_fig=PlotTE.num_fig+4000+100*(ilayer-1)+10*(ifreq-1);
           win_eig=var_iii_eig;
       end

       % build r-axis as r/lambda, lambda is the wave length
       lambda=(PlotTE.Velocity(ifreq,ieig)./1000)/PlotTE.freq;
       rrFL=rrF; if lower(PlotTE.rl_ratio)=='y', rrFL=rrFL/lambda; end

       hf=figure(num_fig);
       subplot(PlotTE.subplot_nxm(1),PlotTE.subplot_nxm(2),win_eig);
       set(gca,'FontSize',PlotTE.font_size);
       hold on; box on;
       plot(rrFL(1:NrrF-1),squeeze(TE_f_e_NIr(ieig,1:NrrF-1)), ...
           'Linewidth',PlotTE.linewidth(1),'Color',PlotTE.color(2));
       if lower(PlotTE.xlabel_show)=='y'
           if lower(PlotTE.rl_ratio)=='y', xlabel('r/\lambda'); else   xlabel('r (m)'); end
       end
       title(strcat('f=',num2str(PlotTE.freq),' kHz',', k_{',num2str(ieig),'}'));
        
   end % cycle on eigenvalues
   zzz=1;  

end % 4)
