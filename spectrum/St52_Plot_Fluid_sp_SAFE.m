%===============================================================================
% Fuild: Plot ux, uy, uz, ur, uf, T energy, U energy, T-P, P and angle of P
% poynting vector (px,py,pz,pr,pf)
%===============================================================================

function St52_Plot_Fluid_sp_SAFE(CompStructA,FEMatricesA,PlotAniso,ResFluid)

    text1_add=num2str(PlotAniso.freq);
    text2_add=num2str(PlotAniso.ieig);
    ResFluid.title_add=strcat(text1_add,' kHz, k','_{',text2_add,'}');
    Result.title_add=ResFluid.title_add;
    Result.title_add_file=Result.title_add;
   
    %===============================================================================
    %===============================================================================
    % 3d Surface
    if lower(PlotAniso.plot_3d)=='y'
       key_plot_c_save=PlotAniso.plot_c;
       if lower(PlotAniso.plot_c)=='y', PlotAniso.plot_c='n'; end
   
       if lower(PlotAniso.subplot)=='y'
          text21_add=num2str(PlotAniso.neigs(1));
          text22_add=num2str(PlotAniso.neigs(length(PlotAniso.neigs)));
          ResFluid.title_add_file=strcat('3D,F=',text1_add,',Eigs=',text21_add,'-',text22_add);
          Result.title_add_file=ResFluid.title_add_file;
       end
        
       
       %============================================================================
       % plot at 10 nodes
  
       Result.xx=ResFluid.xx_10; Result.yy=ResFluid.yy_10; 
       Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
       Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);
       
       Result.field='ux';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+200;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |u_x|:';
              Result.function=abs(ResFluid.ux_10);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |u_x|:'; 
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='uy';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+300;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |u_y|:';
              Result.function=abs(ResFluid.uy_10);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |u_y|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='uz';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+400;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |u_z|:';
              Result.function=abs(ResFluid.uz_10);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |u_z|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='ur';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+500;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              %Result.title_show='Fluid, Surface |u_r|:';
              Result.function=abs(ResFluid.ur_10);
              Result.function=imag(ResFluid.ur_10);              
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |u_r|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='uf';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+600;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |u_f|:';
              Result.function=abs(ResFluid.uf_10);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |u_f|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='te';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+700;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface T Energy:';
              Result.function=ResFluid.TE_10;
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized T Energy:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='ue';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+800;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface U Energy:';
              Result.function=ResFluid.UE_10;
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized U Energy:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='tu';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+900;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface T-U Energy:';
              Result.function=ResFluid.TmU_10;
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized T-U Energy:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='px';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+1000;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |Poyn_x|:';
              Result.function=abs(ResFluid.PV_10(:,1));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |Poyn_x|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='py';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+1100;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |Poyn_y|:';
              Result.function=abs(ResFluid.PV_10(:,2));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |Poyn_y|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pz';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+1200;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |Poyn_z|:';
              Result.function=abs(ResFluid.PV_10(:,3));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |Poyn_z|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pr';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+1300;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |Poyn_r|:';
              Result.function=abs(ResFluid.PV_10(:,4));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |Poyn_r|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pf';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+1400;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface |Poyn_f|:';
              Result.function=abs(ResFluid.PV_10(:,5));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Fluid, Surface of Normalized |Poyn_f|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       %============================================================================
       % Plot pressure and its angle for all nodes
       Result.xx=ResFluid.xx; Result.yy=ResFluid.yy; 
       Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
       Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);
       
       Result.field='pp';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Fluid, Surface of |Pressure|:';
              Result.function=abs(ResFluid.PP);
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Fluid, Surface of Normalized |Pressure|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
                  
       end
       
       Result.field='ap';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+100;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
              Result.title_show='Fluid, Surface of Angle of Pressure:';
              Result.function=ResFluid.fiPP;
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Fluid, Surface of Normalized Angle of Pressure:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
                  
       end

       if lower(key_plot_c_save)=='y', PlotAniso.plot_c='y'; end
       
   end % plot_3D

   %===============================================================================
   %===============================================================================
   if lower(PlotAniso.plot_c)=='y'
      key_plot_3d_save=PlotAniso.plot_3d;
      if lower(PlotAniso.plot_3d)=='y', PlotAniso.plot_3d='n'; end
      
      if lower(PlotAniso.subplot)=='y'
         text21_add=num2str(PlotAniso.neigs(1));
         text22_add=num2str(PlotAniso.neigs(length(PlotAniso.neigs)));
         ResFluid.title_add_file=strcat('C,F=',text1_add,',Eigs=',text21_add,'-',text22_add);
         Result.title_add_file=ResFluid.title_add_file;
      end

      %============================================================================
      % plot at 10 nodes
      Result.xx=ResFluid.xx_10; Result.yy=ResFluid.yy_10; 
      Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
      Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);
      dx=(Result.xx_max-Result.xx_min)/PlotAniso.num_xy(1); dy=(Result.yy_max-Result.yy_min)/PlotAniso.num_xy(2);
      [Result.qxx,Result.qyy]=meshgrid(Result.xx_min:dx:Result.xx_max,Result.yy_min:dy:Result.yy_max);
      
      Result.field='ux';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+1700;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |u_x|:';
             Result.function=abs(ResFluid.ux_10);
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |u_x|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='uy';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+1800;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |u_y|:';
             Result.function=abs(ResFluid.uy_10);
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |u_y|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='uz';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+1900;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |u_z|:';
             Result.function=abs(ResFluid.uz_10);
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |u_z|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='ur';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2000;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |u_r|:';
             Result.function=abs(ResFluid.ur_10);
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |u_r|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='uf';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2100;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |u_f|:';
             Result.function=abs(ResFluid.uf_10);
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |u_f|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='te';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2200;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of T Energy:';
             Result.function=ResFluid.TE_10;
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized T Energy:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='ue';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2300;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of U Energy:';
             Result.function=ResFluid.UE_10;
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized U Energy:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='tu';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2400;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of T-U Energy:';
             Result.function=ResFluid.TmU_10;
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized T-U Energy:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end

      Result.field='px';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2500;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |Poyn_x|:';
             Result.function=abs(ResFluid.PV_10(:,1));
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |Poyn_x|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='py';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2600;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |Poyn_y|:';
             Result.function=abs(ResFluid.PV_10(:,2));
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |Poyn_y|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='pz';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2700;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |Poyn_z|:';
             Result.function=abs(ResFluid.PV_10(:,3));
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |Poyn_z|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='pr';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2800;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |Poyn_r|:';
             Result.function=abs(ResFluid.PV_10(:,4));
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |Poyn_r|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='pf';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+2900;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |Poyn_f|:';
             Result.function=abs(ResFluid.PV_10(:,5));
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |Poyn_f|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end

      %============================================================================ 
      % Plot pressure and its angle for all nodes
      Result.xx=ResFluid.xx; Result.yy=ResFluid.yy; 
      Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
      Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);
      dx=(Result.xx_max-Result.xx_min)/PlotAniso.num_xy(1); dy=(Result.yy_max-Result.yy_min)/PlotAniso.num_xy(2);
      [Result.qxx,Result.qyy]=meshgrid(Result.xx_min:dx:Result.xx_max,Result.yy_min:dy:Result.yy_max);
      
      Result.field='pp';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+1500;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of |Pressure|:';
             Result.function=abs(ResFluid.PP);
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized |Pressure|:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      Result.field='ap';
      for ifield=1:size(PlotAniso.fields,1)
          if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
             Result.num_fig=PlotAniso.num_fig+1600;
             if lower(PlotAniso.subplot)=='n'
                Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
             end
             
             Result.title_show='Fluid, Contours of Angle of Pressure:';
             Result.function=ResFluid.fiPP;
             if lower(PlotAniso.normalization)=='y'
                Result.title_show='Fluid, Contours of Normalized Angle of Pressure:';
                Result.function=Result.function./max(abs(Result.function));
             end
              
             St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
          end
      end
      
      if lower(key_plot_3d_save)=='y', PlotAniso.plot_3d='y'; end
   
   end % contour
       

end % program
