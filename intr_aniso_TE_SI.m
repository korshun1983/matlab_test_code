clear all;

% Directory with processed results
% Dir_Name={'ISO-SAFE'};
% Dir_Name={'BS-45-SAFE'};
% Dir_Name={'BS-30-SAFE'};
% Dir_Name={'BS-60-SAFE'};
% Dir_Name={'BS-90-SAFE-PML'};
% Dir_Name={'BS-90-SAFE'};
% Dir_Name={'BS-45-SAFE'};
% Dir_Name={'AC-90-SAFE'};
Dir_Name={'BS-90-SAFE-D'};

% layer and frequency number for visualization, number of eigs for visualization
IntrTE3.freqs=[1:1:15]; 
IntrTE3.num_eigs_class_show=1;
IntrTE3.changer_12='y'; % if V_mode_1>qSV and V_mode_1 < V_mode_2, then change V_mode_1 and V_mode_2

IntrTE3.max_harm_limit=0.5; % For symmetry classification
% Classification #1 by TE_{PML or ABC}/TE_Total< TE_PML_Total_Limit (Nguyen, Treyssede)
IntrTE3.use_class_TE_NT='y';
% Classification #2 by max TE(PML or ABC)/max TE(Adjacent Region)
IntrTE3.use_class_PML_HTTI='n';

% plot parameters
IntrTE3.show_modes=['y';'y';'y';'y';'y';'n';'n']; % 0-mode,1-mode (p/m),2-mode (p/m),3-mode (p/m) 
IntrTE3.slowness='y'; 
IntrTE3.close_figs='y'; IntrTE3.num_fig=590; IntrTE3.font_size=15; IntrTE3.slowness_limits=[50 250];
IntrTE3.linestyle=['-'; '-'; '-'; '-']; IntrTE3.linewidth=[1; 1; 1; 1]; IntrTE3.color=['k'; 'b'; 'r'; 'm'];
IntrTE3.marker={'o'; 's'; 'd'; 'o'}; IntrTE3.markersize=[7; 7; 7; 7];

Title_Add=char(Dir_Name);
%====================================================================================================================
%====================================================================================================================
%====================================================================================================================
fprintf('\n\n============================================================================\n');
fprintf('Intr_Aniso_TE3 Program has been started!\n');
ColorOrder=[0 0 1; 0 0.5 0; 1 0 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25];
C_foor_n=1./(0.0254*12); y_label_text='Slowness (\mus/m)';

if lower(IntrTE3.close_figs)=='y', close all; end

if strcmp(IntrTE3.use_class_TE_NT,'y') && strcmp(IntrTE3.use_class_PML_HTTI,'y')
   IntrTE3.use_class_PML_HTTI='n';
   beep; beep; beep; 
   fprintf(1,'The first classification algorithm have been chosen!\n');
   fprintf(1,'Please, in future choose only one classification algorithm!\n'); return;
end

IntrTE3.Dir_Name=Dir_Name;
curr_dir=pwd;
% run_prog=mfilename;
% full_path_run_prog=mfilename('fullpath');
% run_dir=full_path_run_prog(1:(length(full_path_run_prog)-length(run_prog)-1));
% cd(run_dir);
% IntrTE3.dir_input=char(Dir_Name);

% Load Calculation Parameters: CompStruct and FEMatrices data 
fname = char(strcat(IntrTE3.Dir_Name,'\CompStruct.mat'));
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

% find input frerquencies in count frequency array and build freq number
% array
IntrTE3.nfreqs=[];
for iff3=1:length(IntrTE3.freqs)
    freq_check=IntrTE3.freqs(iff3);
    find_freq='n';
    for iff=1:CompStructA.Model.N_disp
        if CompStructA.Model.f_array(iff)==freq_check
           find_freq='y'; 
           IntrTE3.nfreqs=[IntrTE3.nfreqs; iff];
           break;
        end
    end
    if find_freq=='n'
       beep; beep; beep; fprintf(1,'At least One frequency did not find!\n'); return;
    end
end
IntrTE3.nfreqs=IntrTE3.nfreqs';
IntrTE3.freqs_count=CompStructA.Model.f_array;
IntrTE3.nfreqs_count=(1:CompStructA.Model.N_disp);

fname = char(strcat(IntrTE3.Dir_Name,'\FEMatrices-',num2str(IntrTE3.freqs(1)),'.mat'));
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

IntrTE3.nlayers_count=1:size(CompStructA.Model.DomainType,2); %=ProcTE.nlayers
ResLayer=cell(length(IntrTE3.nlayers_count),1);

% Load ProcTE struct
fname = char(strcat(IntrTE3.Dir_Name,'/','LayerTE','-',num2str(IntrTE3.freqs(1)),'-',num2str(1)','.mat'));
if exist(fname,'file')>0
    load(fname);
else
    beep; beep; beep; fprintf(1,'File is not existed!\n'); return;
end

%====================================================================================================================
%====================================================================================================================
% Load Slowness
fname = char(strcat(IntrTE3.Dir_Name,'/Slowness.mat'));
if exist(fname,'file')>0
    load(fname);
    IntrTE3.freqs_slow_count=Result_Freq;
    IntrTE3.Slowness=Result_Slowness;
    IntrTE3.Velocity=Result_Velocity;
    
    if strcmp(IntrTE3.slowness,'y')
       % Slowness
       f_max=max(CompStructA.Model.f_array);
       
       hf_full=figure(IntrTE3.num_fig+12); 
       set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
       hold on; box on; xlim([0  f_max+1]); ylim(IntrTE3.slowness_limits.*C_foor_n);
%        plot(IntrTE3.freqs_count,IntrTE3.Slowness(1:length(IntrTE3.freqs_count),1:CompStructA.Advanced.num_eig_max), ...
%            'Marker',char(IntrTE3.marker(1)),'Markersize',IntrTE3.markersize(1), ...
%            'LineStyle','none','Linewidth',IntrTE3.linewidth(1));
       plot(IntrTE3.freqs_count,...
           IntrTE3.Slowness(1:length(IntrTE3.freqs_count),1:CompStructA.Advanced.num_eig_max).*C_foor_n, ...
           'Marker',char(IntrTE3.marker(1)),'Markersize',IntrTE3.markersize(1), ...
           'LineStyle','none','Linewidth',IntrTE3.linewidth(1),'Color',IntrTE3.color(1));

       freq_var=[0 IntrTE3.freqs_count f_max+1];
       qP_val=ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_qP).*C_foor_n; 
       qSV_val=ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_qSV.*C_foor_n); 
       SH_val =ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_SH.*C_foor_n);
       St_val =ones(length(freq_var),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_St.*C_foor_n);
       plot(freq_var,qP_val,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(4));
       plot(freq_var,qSV_val,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(3));
       plot(freq_var,SH_val ,'LineStyle',IntrTE3.linestyle(3),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(2));
       %plot(freq_var,St_val,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(4));
       xlabel('f (kHz)'); ylabel(y_label_text); title('');
    end     
end

%====================================================================================================================
%====================================================================================================================
% Build input layers
IntrTE3.nlayers=ProcTE.nlayers;
IntrTE3.nlayers_count=ProcTE.nlayers_count; 
for ilayer=IntrTE3.nlayers
    if ~ismember(ilayer,IntrTE3.nlayers_count)
       beep; beep; beep; fprintf(1,'Error in IntrTE3.nlayers!\n'); return;
    end
end

% Build Eigs Array
IntrTE3.neigs=ProcTE.neigs;
IntrTE3.max_neigs=CompStructA.Advanced.num_eig_max; % number of calculated eigs
IntrTE3.neigs_count=ProcTE.neigs;
for ieig=IntrTE3.neigs
    if ~ismember(ieig,IntrTE3.neigs_count)
       beep; beep; beep; fprintf(1,'Error in IntrTE3.neigs!\n'); return;
    end
end


% Massive #3. 1-value of harmonic TE, 2-eig number ==> 
%             3-value of Classification #1, % 4-eig number, 5-slowness, 
%             6-value of Classification #2, % 7-eig number, 8-slowness, 
%             10 - frequency
Class_TE_0 =zeros(length(IntrTE3.nfreqs_count),length(IntrTE3.neigs_count),10);  % monople
Class_TE_1p=zeros(length(IntrTE3.nfreqs_count),length(IntrTE3.neigs_count),10);  % dypole +
Class_TE_1m=zeros(length(IntrTE3.nfreqs_count),length(IntrTE3.neigs_count),10);  % dypole -
Class_TE_2p=zeros(length(IntrTE3.nfreqs_count),length(IntrTE3.neigs_count),10);  % quadrupole + 
Class_TE_2m=zeros(length(IntrTE3.nfreqs_count),length(IntrTE3.neigs_count),10);  % quadrupole -
Class_TE_3p=zeros(length(IntrTE3.nfreqs_count),length(IntrTE3.neigs_count),10);  % 3pole + 
Class_TE_3m=zeros(length(IntrTE3.nfreqs_count),length(IntrTE3.neigs_count),10);  % 3pole -

%====================================================================================================================
% Freq Cycle
%====================================================================================================================
iff_1=0;
for iff=IntrTE3.nfreqs
    iff_1=iff_1+1;
    fprintf(1,'\tNow fhe %d-th frequency (%.2f kHz) is processed out of %d (real %d)\n', ...
            iff_1,IntrTE3.freqs_count(iff),length(IntrTE3.nfreqs),length(IntrTE3.nfreqs_count));

    % Load processing results
    for ilayer=IntrTE3.nlayers
        fname = char(strcat(IntrTE3.Dir_Name,'/','LayerTE','-',num2str(IntrTE3.freqs_count(iff)),'-',num2str(ilayer)','.mat'));
        if exist(fname,'file')>0
            load(fname);
        else
            beep; beep; beep; fprintf(1,'File is not existed!\n'); return;
        end
        ResLayer{ilayer}.ProcTE=ProcTE;
        ResLayer{ilayer}.rr=rr;
        ResLayer{ilayer}.TE_f_e_r=TE_fer(:,:);
        ResLayer{ilayer}.cF_ur_rm_pm=cF_ur_rm_pm(:,:,:,:);
        ResLayer{ilayer}.cF_uf_rm_pm=cF_uf_rm_pm(:,:,:,:);
        ResLayer{ilayer}.cF_uz_rm_pm=cF_uz_rm_pm(:,:,:,:);
        zzz=1;
    end

    %====================================================================================================================
    %====================================================================================================================
    % A) Calculate T Energy for specified f, eig, r and h (harmonic): TE(f,eig,r,h)
    NrrF=0; rrF=[];
    for ilayer=IntrTE3.nlayers
        NrrF=NrrF+size(ResLayer{ilayer}.rr,1);
        rrF=[rrF(1:end); ResLayer{ilayer}.rr];
    end
    NHarm=ResLayer{1}.ProcTE.NHarm;

    ifreq=iff; 
    
    TE_f_eigs_r_h_p=zeros(length(IntrTE3.neigs_count),NrrF,NHarm);
    TE_f_eigs_r_h_m=zeros(length(IntrTE3.neigs_count),NrrF,NHarm);
    TE_f_eigs_r_h=zeros(length(IntrTE3.neigs_count),NrrF,NHarm);

    for ieig=IntrTE3.neigs
        for iharm=1:NHarm

            TE_fer_h_p=[];
            TE_fer_h_m=[];
            for ilayer=IntrTE3.nlayers
                % rho*v^2/2
                varTE = 0.5.*FEMatricesA.PhysProp{ilayer}.rho*1000. ...
                       *(2.*pi.*IntrTE3.freqs_count(ifreq).*CompStructA.Misc.F_conv).^2;

                var_p=  abs(ResLayer{ilayer}.cF_ur_rm_pm(ieig,:,iharm,1)).^2 ...
                       +abs(ResLayer{ilayer}.cF_uf_rm_pm(ieig,:,iharm,1)).^2 ...
                       +abs(ResLayer{ilayer}.cF_uz_rm_pm(ieig,:,iharm,1)).^2;
                var_p=varTE*var_p;

                var_m=  abs(ResLayer{ilayer}.cF_ur_rm_pm(ieig,:,iharm,2)).^2 ...
                       +abs(ResLayer{ilayer}.cF_uf_rm_pm(ieig,:,iharm,2)).^2 ...
                       +abs(ResLayer{ilayer}.cF_uz_rm_pm(ieig,:,iharm,2)).^2;
                var_m=varTE*var_m;

                TE_fer_h_p=[TE_fer_h_p(1:end); squeeze(var_p')];
                TE_fer_h_m=[TE_fer_h_m(1:end); squeeze(var_m')];
            end
            TE_f_eigs_r_h_p(ieig,:,iharm)=TE_fer_h_p;
            TE_f_eigs_r_h_m(ieig,:,iharm)=TE_fer_h_m;

            var_const=1;
            %if iharm==1, var_const=0.5; end
            TE_f_eigs_r_h(ieig,:,iharm)=var_const.*(TE_fer_h_p+TE_fer_h_m);

        end % harms
    end % eigs
    
    %====================================================================================================================
    %====================================================================================================================
    % B)
    % B.1)Integrate TE(f,eig,r,h): TE(f,eig,Ir,h) = Sum_{r=r_min}^r r*TE(f,eig,r,h) dr, r_i=r_min, ..., r_max+r_pml
    % B.2) Normalize TE(f,eig,Ir,h) by max of h
    TE_f_eigs_Ir_h=zeros(length(IntrTE3.neigs_count),NrrF-1,NHarm); 
    TE_f_eigs_Ir_Nh=zeros(length(IntrTE3.neigs_count),NrrF-1,NHarm); 

    for ieig=IntrTE3.neigs
        for iharm=1:NHarm

            varl=squeeze(TE_f_eigs_r_h(ieig,1:NrrF-1,iharm));
            varr=squeeze(TE_f_eigs_r_h(ieig,2:NrrF  ,iharm));
            var_mas=0.5.*( rrF(1:NrrF-1).*varl'+rrF(2:NrrF).*varr' ).*(rrF(2:NrrF)-rrF(1:NrrF-1));
            TE_f_eigs_Ir_h(ieig,1,iharm)=var_mas(1);
            for ir=2:NrrF-1
                TE_f_eigs_Ir_h(ieig,ir,iharm)= ...
                    TE_f_eigs_Ir_h(ieig,ir-1,iharm)+var_mas(ir);
            end
            
        end % harms
        
        varmax=max(TE_f_eigs_Ir_h(ieig,NrrF-1,:));
        for iharm=1:NHarm
            TE_f_eigs_Ir_Nh(ieig,:,iharm)=(TE_f_eigs_Ir_h(ieig,:,iharm)./varmax);
        end
        
    end % eigs

    %====================================================================================================================
    %====================================================================================================================
    % C) For Polarization
    % C.1) Calculate TE_PZp(f,eig,h)=Sum_{i=1}^Nrr abs( f_(-n)+f_n )
    % C.2) Calculate TE_PZm(f,eig,h)=Sum_{i=1}^Nrr abs( f_(-n)-f_n )
    TE_PZp_f_e_h=zeros(length(IntrTE3.neigs_count),NHarm);
    TE_PZm_f_e_h=zeros(length(IntrTE3.neigs_count),NHarm);

    for ieig=IntrTE3.neigs
        for iharm=1:NHarm
        
            for ilayer=IntrTE3.nlayers
                var_p=squeeze(ResLayer{ilayer}.cF_ur_rm_pm(ieig,:,iharm,1));
                var_m=squeeze(ResLayer{ilayer}.cF_ur_rm_pm(ieig,:,iharm,2));
                
                var_sum_p_p_m=sum(abs(var_p+var_m));
                TE_PZp_f_e_h(ieig,iharm)=TE_PZp_f_e_h(ieig,iharm)+var_sum_p_p_m;
                
                var_sum_p_m_m=sum(abs(var_p-var_m));
                TE_PZm_f_e_h(ieig,iharm)=TE_PZm_f_e_h(ieig,iharm)+var_sum_p_m_m;
            end
            
        end % harms
    end % eigs

    %====================================================================================================================
    %====================================================================================================================
    % D) Full TE 
    % D.1) TE(f,eig,Ir)=sum_{r=r_min}^{r} r*TE_fer(f,eig,r) dr, r_i=r_min, ..., r_max+r_pml
    %      TE(f,eig,r) is calculated by integrating by angle in proc_aniso
    % D.2) TE_Total(f,eig)=TE(f,eig,r_max+r_pml)
    % D.3) TE_PML(f,eig)=TE(f,eig,r_max+r_pml)-TE(f,eig,r_max)
    % D.4) TE(f,eig,Nr_max)=Ir TE(f,eig,Ir)/TE_f_e_Ir(f,eig,Ir_max)
    TE_f_e_Ir   =zeros(length(IntrTE3.neigs_count),NrrF-1);
    TE_f_e_Total=zeros(length(IntrTE3.neigs_count),1);
    TE_f_e_PML=zeros(length(IntrTE3.neigs_count),1);
    TE_f_e_NIr  =zeros(length(IntrTE3.neigs_count),NrrF-1);
 
    for ieig=IntrTE3.neigs

        TE_fer =[];
        for ilayer=IntrTE3.nlayers
            TE_fer =[TE_fer(1:end); squeeze(ResLayer{ilayer}.TE_f_e_r(ieig,:)')];
        end

        varl=(TE_fer(1:NrrF-1));
        varr=(TE_fer(2:NrrF  ));
        var_mas=0.5.*( rrF(1:NrrF-1).*varl+rrF(2:NrrF).*varr ).*(rrF(2:NrrF)-rrF(1:NrrF-1));
        TE_f_e_Ir(ieig,1)=var_mas(1);
        for ir=2:NrrF-1
            TE_f_e_Ir(ieig,ir)= ...
                TE_f_e_Ir(ieig,ir-1)+var_mas(ir);
        end
        TE_f_e_Total(ieig)=TE_f_e_Ir(ieig,NrrF-1);
        
        var_ind=find(rrF==ResLayer{end}.rr(1));
        TE_f_e_PML(ieig)=TE_f_e_Ir(ieig,NrrF-1)-TE_f_e_Ir(ieig,var_ind-1);
        
        TE_f_e_NIr(ieig,:)= ...
            TE_f_e_Ir(ieig,:)./TE_f_e_Ir(ieig,NrrF-1);

    end % eigs

    %====================================================================================================================
    %====================================================================================================================
    % E) Calculate Max TE in each layer for specified f and eig
    TE_f_L_e=zeros(length(IntrTE3.neigs_count),length(IntrTE3.nlayers));
 
    for ieig=IntrTE3.neigs
        for ilayer=IntrTE3.nlayers
            TE_f_L_e(ieig,ilayer)=max(abs(ResLayer{ilayer}.TE_f_e_r(ieig,:)));
        end
    end % eigs

    %====================================================================================================================
    %====================================================================================================================
    % Classification by Symmetry and Polarization 

    for ieig=IntrTE3.neigs
        
        var_ind=(squeeze(TE_f_eigs_Ir_Nh(ieig,NrrF-1,:)>=IntrTE3.max_harm_limit));
        if isempty(var_ind), beep; beep; beep; fprintf(1,'Error!!!\n'); return; end
        %[var_res,var_ind]=max(TE_f_eigs_Ir_Nh(ieig,NrrF-1,:));
        if var_ind(1)==1 % monopole
           Class_TE_0(ifreq,ieig,1)=TE_f_eigs_Ir_Nh(ieig,NrrF-1,1);
           Class_TE_0(ifreq,ieig,2)=ieig;
        end
        
        if var_ind(2)==1 % dypole
           if TE_PZp_f_e_h(ieig,2)>=TE_PZm_f_e_h(ieig,2)
              Class_TE_1p(ifreq,ieig,1)=TE_f_eigs_Ir_Nh(ieig,NrrF-1,2);
              Class_TE_1p(ifreq,ieig,2)=ieig;
           else
              Class_TE_1m(ifreq,ieig,1)=TE_f_eigs_Ir_Nh(ieig,NrrF-1,2);
              Class_TE_1m(ifreq,ieig,2)=ieig;
           end
        end
        
        if var_ind(3)==1 % quadrupole
           if TE_PZp_f_e_h(ieig,3)>=TE_PZm_f_e_h(ieig,3)
              Class_TE_2p(ifreq,ieig,1)=TE_f_eigs_Ir_Nh(ieig,NrrF-1,3);
              Class_TE_2p(ifreq,ieig,2)=ieig;
           else
              Class_TE_2m(ifreq,ieig,1)=TE_f_eigs_Ir_Nh(ieig,NrrF-1,3);
              Class_TE_2m(ifreq,ieig,2)=ieig;
           end
        end

        if var_ind(4)==1 % quadrupole
           if TE_PZp_f_e_h(ieig,4)>=TE_PZm_f_e_h(ieig,4)
              Class_TE_3p(ifreq,ieig,1)=TE_f_eigs_Ir_Nh(ieig,NrrF-1,3);
              Class_TE_3p(ifreq,ieig,2)=ieig;
           else
              Class_TE_3m(ifreq,ieig,1)=TE_f_eigs_Ir_Nh(ieig,NrrF-1,3);
              Class_TE_3m(ifreq,ieig,2)=ieig;
           end
        end

    end % eigs

    eigs_array=1:IntrTE3.num_eigs_class_show;
    
    % Classification by TE_PML/TE_Total< TE_PML_Total_Limit (Nguyen, Treyssede)
    % For real eigs this parameter have small value
    Class_TE_NT=TE_f_e_PML(:)./TE_f_e_Total(:);
       
    % Classification by    max TE   /    max TE 
    %                    PML region   adjacent region
    % For real eigs this parameter have small value
    Class_TE_max_PML_HTTI=TE_f_L_e(:,end)./TE_f_L_e(:,end-1);
        
    %====================================================================================================================
    % monopole
    var_tmp_ind=Class_TE_0(ifreq,:,2);
    var_tmp_ind0=(Class_TE_0(ifreq,:,2)==0);
    var_tmp_ind(var_tmp_ind0)=[];
    [var_res1,var_ind1]=sort(Class_TE_NT(var_tmp_ind));
    Class_TE_0(ifreq,1:length(var_res1),3)=var_res1;
    Class_TE_0(ifreq,1:length(var_res1),4)=var_tmp_ind(var_ind1);
    Class_TE_0(ifreq,1:length(var_res1),5)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind1));

    [var_res2,var_ind2]=sort(Class_TE_max_PML_HTTI(var_tmp_ind));
    Class_TE_0(ifreq,1:length(var_res2),6)=var_res2;
    Class_TE_0(ifreq,1:length(var_res2),7)=var_tmp_ind(var_ind2);
    Class_TE_0(ifreq,1:length(var_res2),8)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind2));
    
    Class_TE_0(ifreq,1:length(var_ind),10)=IntrTE3.freqs_count(ifreq);

    mode_sort='ascend';
    % dypole +
    var_tmp_ind=Class_TE_1p(ifreq,:,2);
    var_tmp_ind0=(Class_TE_1p(ifreq,:,2)==0);
    var_tmp_ind(var_tmp_ind0)=[];
    [var_res1,var_ind1]=sort(Class_TE_NT(var_tmp_ind),mode_sort);
    Class_TE_1p(ifreq,1:length(var_res1),3)=var_res1;
    Class_TE_1p(ifreq,1:length(var_res1),4)=var_tmp_ind(var_ind1);
    Class_TE_1p(ifreq,1:length(var_res1),5)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind1));

    [var_res2,var_ind2]=sort(Class_TE_max_PML_HTTI(var_tmp_ind),mode_sort);
    Class_TE_1p(ifreq,1:length(var_res2),6)=var_res2;
    Class_TE_1p(ifreq,1:length(var_res2),7)=var_tmp_ind(var_ind2);
    Class_TE_1p(ifreq,1:length(var_res2),8)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind2));
    
    Class_TE_1p(ifreq,1:length(var_ind),10)=IntrTE3.freqs_count(ifreq);
    
    % dypole -
    var_tmp_ind=Class_TE_1m(ifreq,:,2);
    var_tmp_ind0=(Class_TE_1m(ifreq,:,2)==0);
    var_tmp_ind(var_tmp_ind0)=[];
    [var_res1,var_ind1]=sort(Class_TE_NT(var_tmp_ind),mode_sort);
    Class_TE_1m(ifreq,1:length(var_res1),3)=var_res1;
    Class_TE_1m(ifreq,1:length(var_res1),4)=var_tmp_ind(var_ind1);
    Class_TE_1m(ifreq,1:length(var_res1),5)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind1));
    
    [var_res2,var_ind2]=sort(Class_TE_max_PML_HTTI(var_tmp_ind),mode_sort);
    Class_TE_1m(ifreq,1:length(var_res2),6)=var_res2;
    Class_TE_1m(ifreq,1:length(var_res2),7)=var_tmp_ind(var_ind2);
    Class_TE_1m(ifreq,1:length(var_res2),8)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind2));

    Class_TE_1m(ifreq,1:length(var_ind),10)=IntrTE3.freqs_count(ifreq);

    % quadrupole +
    var_tmp_ind=Class_TE_2p(ifreq,:,2);
    var_tmp_ind0=(Class_TE_2p(ifreq,:,2)==0);
    var_tmp_ind(var_tmp_ind0)=[];
    [var_res1,var_ind1]=sort(Class_TE_NT(var_tmp_ind),mode_sort);
    Class_TE_2p(ifreq,1:length(var_res1),3)=var_res1;
    Class_TE_2p(ifreq,1:length(var_res1),4)=var_tmp_ind(var_ind1);
    Class_TE_2p(ifreq,1:length(var_res1),5)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind1));
    
    [var_res2,var_ind2]=sort(Class_TE_max_PML_HTTI(var_tmp_ind),mode_sort);
    Class_TE_2p(ifreq,1:length(var_res2),6)=var_res2;
    Class_TE_2p(ifreq,1:length(var_res2),7)=var_tmp_ind(var_ind2);
    Class_TE_2p(ifreq,1:length(var_res2),8)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind2));
    
    Class_TE_2p(ifreq,1:length(var_ind),10)=IntrTE3.freqs_count(ifreq);

    % quadrupole -
    var_tmp_ind=Class_TE_2m(ifreq,:,2);
    var_tmp_ind0=(Class_TE_2m(ifreq,:,2)==0);
    var_tmp_ind(var_tmp_ind0)=[];
    [var_res1,var_ind1]=sort(Class_TE_NT(var_tmp_ind),mode_sort);
    Class_TE_2m(ifreq,1:length(var_res1),3)=var_res1;
    Class_TE_2m(ifreq,1:length(var_res1),4)=var_tmp_ind(var_ind1);
    Class_TE_2m(ifreq,1:length(var_res1),5)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind1));
    
    [var_res2,var_ind2]=sort(Class_TE_max_PML_HTTI(var_tmp_ind),mode_sort);
    Class_TE_2m(ifreq,1:length(var_res2),6)=var_res2;
    Class_TE_2m(ifreq,1:length(var_res2),7)=var_tmp_ind(var_ind2);
    Class_TE_2m(ifreq,1:length(var_res2),8)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind2));
    
    Class_TE_2m(ifreq,1:length(var_ind),10)=IntrTE3.freqs_count(ifreq);

    % 3pole +
    var_tmp_ind=Class_TE_3p(ifreq,:,2);
    var_tmp_ind0=(Class_TE_3p(ifreq,:,2)==0);
    var_tmp_ind(var_tmp_ind0)=[];
    [var_res1,var_ind1]=sort(Class_TE_NT(var_tmp_ind),mode_sort);
    Class_TE_3p(ifreq,1:length(var_res1),3)=var_res1;
    Class_TE_3p(ifreq,1:length(var_res1),4)=var_tmp_ind(var_ind1);
    Class_TE_3p(ifreq,1:length(var_res1),5)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind1));

    [var_res2,var_ind2]=sort(Class_TE_max_PML_HTTI(var_tmp_ind),mode_sort);
    Class_TE_3p(ifreq,1:length(var_res2),6)=var_res2;
    Class_TE_3p(ifreq,1:length(var_res2),7)=var_tmp_ind(var_ind2);
    Class_TE_3p(ifreq,1:length(var_res2),8)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind2));
    
    Class_TE_3p(ifreq,1:length(var_ind),10)=IntrTE3.freqs_count(ifreq);

    % 3pole -
    var_tmp_ind=Class_TE_3m(ifreq,:,2);
    var_tmp_ind0=(Class_TE_3m(ifreq,:,2)==0);
    var_tmp_ind(var_tmp_ind0)=[];
    [var_res1,var_ind1]=sort(Class_TE_NT(var_tmp_ind),mode_sort);
    Class_TE_3m(ifreq,1:length(var_res1),3)=var_res1;
    Class_TE_3m(ifreq,1:length(var_res1),4)=var_tmp_ind(var_ind1);
    Class_TE_3m(ifreq,1:length(var_res1),5)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind1));

    [var_res2,var_ind2]=sort(Class_TE_max_PML_HTTI(var_tmp_ind),mode_sort);
    Class_TE_3m(ifreq,1:length(var_res2),6)=var_res2;
    Class_TE_3m(ifreq,1:length(var_res2),7)=var_tmp_ind(var_ind2);
    Class_TE_3m(ifreq,1:length(var_res2),8)=IntrTE3.Slowness(ifreq,var_tmp_ind(var_ind2));
    
    Class_TE_3m(ifreq,1:length(var_ind),10)=IntrTE3.freqs_count(ifreq);

end % Freqs Cycle

%====================================================================================================================
%====================================================================================================================
% Plot Results
freq_count=IntrTE3.freqs_count(IntrTE3.nfreqs);
if freq_count(1)-1>=0, freq_var=freq_count(1)-1; else freq_vsar=0; end
freq_countA=[freq_count(1)-1 freq_count freq_count(end)+1];
qSV_val=ones((length(freq_count)+2),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_qSV); 
SH_val =ones((length(freq_count)+2),1).*(CompStructA.Misc.S_conv./CompStructA.Asymp.V_SH);

if lower(IntrTE3.changer_12)=='y'
    iff_1=0;
    for iff=IntrTE3.nfreqs
        iff_1=iff_1+1;

        for iwave=1:7
            if iwave==1, Class_tmp=Class_TE_0; end
            if iwave==2, Class_tmp=Class_TE_1p; end
            if iwave==3, Class_tmp=Class_TE_1m; end
            if iwave==4, Class_tmp=Class_TE_2p; end
            if iwave==5, Class_tmp=Class_TE_2m; end
            if iwave==6, Class_tmp=Class_TE_3p; end
            if iwave==7, Class_tmp=Class_TE_3m; end

            if (Class_tmp(iff,2,5) > qSV_val(1)) && (Class_tmp(iff,1,5) < Class_tmp(iff,2,5))
               for iind=[3 4 5]
                   var=Class_tmp(iff,1,iind); 
                   Class_tmp(iff,1,iind)=Class_tmp(iff,2,iind); 
                   Class_tmp(iff,2,iind)=var;
               end
            end
            if (Class_tmp(iff,2,8) > qSV_val(1)) && (Class_tmp(iff,1,8) < Class_tmp(iff,2,8))
               for iind=[6 7 8]
                   var=Class_tmp(iff,1,iind); 
                   Class_tmp(iff,1,iind)=Class_tmp(iff,2,iind); 
                   Class_tmp(iff,2,iind)=var;
               end
            end

            if iwave==1, Class_TE_0=Class_tmp; end
            if iwave==2, Class_TE_1p=Class_tmp; end
            if iwave==3, Class_TE_1m=Class_tmp; end
            if iwave==4, Class_TE_2p=Class_tmp; end
            if iwave==5, Class_TE_2m=Class_tmp; end
            if iwave==6, Class_TE_3p=Class_tmp; end
            if iwave==7, Class_TE_3m=Class_tmp; end
        end
    end
end

if strcmp(IntrTE3.use_class_TE_NT,'y'), ikey=5; end
if strcmp(IntrTE3.use_class_PML_HTTI,'y'), ikey=8; end

% Stoneley
if lower(IntrTE3.show_modes(1)=='y')
    
    hf0=figure(IntrTE3.num_fig+10);
    set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
    hold on; box on; xlim([freq_countA(1) freq_countA(end)]); ylim(IntrTE3.slowness_limits.*C_foor_n);
    plot(IntrTE3.freqs_count,Class_TE_0(:,1:IntrTE3.num_eigs_class_show,ikey).*C_foor_n,'LineStyle','none', ...
        'Marker',char(IntrTE3.marker(1)),'MarkerSize',IntrTE3.markersize(1),'LineWidth',IntrTE3.linewidth(1));
    plot(freq_countA,SH_val.*C_foor_n ,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(2));
    plot(freq_countA,qSV_val.*C_foor_n,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(3));
    xlabel('f (kHz)'); ylabel(y_label_text); 
    title(strcat(Title_Add,' Monopole modes'));
end

% Dypole +
if lower(IntrTE3.show_modes(2)=='y')
    
    hf1p=figure(IntrTE3.num_fig+20);
    set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
    hold on; box on; xlim([freq_countA(1) freq_countA(end)]); ylim(IntrTE3.slowness_limits.*C_foor_n);
    plot(IntrTE3.freqs_count, Class_TE_1p(:,1:IntrTE3.num_eigs_class_show,ikey).*C_foor_n,'LineStyle','none', ...
        'Marker',char(IntrTE3.marker(2)),'MarkerSize',IntrTE3.markersize(2),'LineWidth',IntrTE3.linewidth(1));
    plot(freq_countA,SH_val.*C_foor_n ,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(2));
    plot(freq_countA,qSV_val.*C_foor_n,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(3));
    xlabel('f (kHz)'); ylabel(y_label_text); 
    title(strcat(Title_Add,' Flexural modes +'));
end

% Dypole -
if lower(IntrTE3.show_modes(3)=='y')
    
    hf1m=figure(IntrTE3.num_fig+21);
    set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
    hold on; box on; xlim([freq_countA(1) freq_countA(end)]); ylim(IntrTE3.slowness_limits.*C_foor_n);
    plot(IntrTE3.freqs_count, Class_TE_1m(:,1:IntrTE3.num_eigs_class_show,ikey).*C_foor_n, 'LineStyle','none', ...
        'Marker',char(IntrTE3.marker(2)),'MarkerSize',IntrTE3.markersize(2),'LineWidth',IntrTE3.linewidth(1));
    plot(freq_countA,SH_val.*C_foor_n ,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(1),'Color',IntrTE3.color(2));
    plot(freq_countA,qSV_val.*C_foor_n,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(1),'Color',IntrTE3.color(3));
    xlabel('f (kHz)'); ylabel(y_label_text);
    title(strcat(Title_Add,' Flexural modes -'));

end

% Quadupole +
if lower(IntrTE3.show_modes(4)=='y')
    
    hf2p=figure(IntrTE3.num_fig+30);
    set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
    hold on; box on; xlim([freq_countA(1) freq_countA(end)]); ylim(IntrTE3.slowness_limits.*C_foor_n);
    plot(IntrTE3.freqs_count, Class_TE_2p(:,1:IntrTE3.num_eigs_class_show,ikey).*C_foor_n,'LineStyle','none', ...
        'Marker',char(IntrTE3.marker(3)),'MarkerSize',IntrTE3.markersize(3),'LineWidth',IntrTE3.linewidth(2));
    plot(freq_countA,SH_val.*C_foor_n ,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(2));
    plot(freq_countA,qSV_val.*C_foor_n,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(2),'Color',IntrTE3.color(3));
    xlabel('f (kHz)'); ylabel(y_label_text);
    title(strcat(Title_Add,' Quadrupole modes +'));

end
     
% Quadupole -
if lower(IntrTE3.show_modes(5)=='y')
    
    hf2m=figure(IntrTE3.num_fig+31);
    set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
    hold on; box on; xlim([freq_countA(1) freq_countA(end)]); ylim(IntrTE3.slowness_limits.*C_foor_n);
    plot(IntrTE3.freqs_count, Class_TE_2m(:,1:IntrTE3.num_eigs_class_show,ikey).*C_foor_n,'LineStyle','none', ...
        'Marker',char(IntrTE3.marker(3)),'MarkerSize',IntrTE3.markersize(3),'LineWidth',IntrTE3.linewidth(3));
    plot(freq_countA,SH_val.*C_foor_n ,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(3),'Color',IntrTE3.color(2));
    plot(freq_countA,qSV_val.*C_foor_n,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(3),'Color',IntrTE3.color(3));
    xlabel('f (kHz)'); ylabel(y_label_text);
    title(strcat(Title_Add,' Quadrupole modes -'));

end

% 3pole +
if lower(IntrTE3.show_modes(6)=='y')
    
    hf3p=figure(IntrTE3.num_fig+40);
    set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
    hold on; box on; xlim([freq_countA(1) freq_countA(end)]); ylim(IntrTE3.slowness_limits.*C_foor_n);
    plot(IntrTE3.freqs_count, Class_TE_3p(:,1:IntrTE3.num_eigs_class_show,ikey).*C_foor_n,'LineStyle','none', ...
        'Marker',char(IntrTE3.marker(4)),'MarkerSize',IntrTE3.markersize(4),'LineWidth',IntrTE3.linewidth(4));
    plot(freq_countA,SH_val.*C_foor_n ,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(4),'Color',IntrTE3.color(2));
    plot(freq_countA,qSV_va.*C_foor_nl,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(4),'Color',IntrTE3.color(3));
    xlabel('f (kHz)'); ylabel(y_label_text');
    title(strcat(Title_Add,' 3pole +'));
    
end

% 3pole -
if lower(IntrTE3.show_modes(7)=='y')
    
    hf3p=figure(IntrTE3.num_fig+41);
    set(gca,'FontSize',IntrTE3.font_size); set(gca,'ColorOrder',ColorOrder);
    hold on; box on; xlim([freq_countA(1) freq_countA(end)]); ylim(IntrTE3.slowness_limits.*C_foor_n);
    plot(IntrTE3.freqs_count, Class_TE_3m(:,1:IntrTE3.num_eigs_class_show,ikey).*C_foor_n,'LineStyle','none', ...
        'Marker',char(IntrTE3.marker(4)),'MarkerSize',IntrTE3.markersize(4),'LineWidth',IntrTE3.linewidth(4));
    plot(freq_countA,SH_val.*C_foor_n,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(4),'Color',IntrTE3.color(2));
    plot(freq_countA,qSV_val.*C_foor_n,'LineStyle',IntrTE3.linestyle(2),'Linewidth',IntrTE3.linewidth(4),'Color',IntrTE3.color(3));
    xlabel('f (kHz)'); ylabel(y_label_text);
    title(strcat(Title_Add,' 3pole -'));

end
%====================================================================================================================
%====================================================================================================================

cd(curr_dir);
fprintf('============================================================================\n');
zzz=1;



