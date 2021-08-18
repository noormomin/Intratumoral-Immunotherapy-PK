function Momin2021
%Written by Noor Momin, October 22, 2019
%Updated by Noor Momin, April 9, 2020
%Updated by Noor Momin, Sept 11, 2020
%Updated by Noor Momin, May 25, 2021 (for paper submission)
%Updated by Noor Momin, August 18, 2021 (for paper revision)

clc; clear all; close all
warning('off')

%Description:
%INPUT PARAMETERS PAYLOAD AND ANTIGEN ()
%   %   Antigen/Payload Affinity - [M]
        Kd = 9.3E-9;
%   %   Payload Molecular Weight - [kDa]
        MW =  94;
%   %   Payload serum half-life  - [days]
        Th_L = log(2)/(0.019*24);%5.1;
%   %   Antigen turnover at QSS, intratumor half life - [days]
        Th_C = 30; %6;
%   %   Payload Uptake via Antigen Turnover - (0:no; 1:yes)
        A = 1;
%   %   Payload Intratumoral Injected Dose - [M]
        % 0.1 nmol dose in 20 uL--> 0.2 void fraction for a 50 mm^3 tumor, so 10 of 20 uL --> 0.05
        % nmol/50E-6 uL
        L = 1E-6;
%   %   Antigen Intratumoral Concentration - [M]
        C = 2E-7;
       
       
%INPUT PARAMETERS RECEPTOR (QSS)      
%   %   Receptor Presence  - (0:none; 1:yes)
        receptor = 1;             
%   %   Intrautmoral Receptor+ Cells - [cells/mL]
        cells = 150E3;
%   %   Intrautmoral Receptor Density on Cells - [receptors/cell]
        NR = 1000;
%   %   Receptor Engaged Cell Proliferation/Expansion - [1/day]
        %Kosmrlj et al Nature Letters, Supplement        
        expand = 1;
%   %   Receptor Rate Constants
        kon_R = 1.26E6;                     % [1/M/s] - on rate for IL-2 and IL-2RB from NKTR paper
        koff_R = 0.301;                      % [1/s] - off rate for IL-2 and IL-2RB NKTR paper
        kendo_R = 4E-2/(60);                % [1/s] - endocytic rate of IL-2R with ligand
        
% ODE solver options
options = odeset('RelTol',1e-14,'AbsTol',[1e-14]);
tspan = [0 259200];                        %72hours

%% Figures Generated in MATLAB Momin et al.
% % % % % % % % % % % % % % % READ ME % % % % % % % % % % % % % % % % % % %
% Figures were generated using the following code and colorized in Illustrator
% To generate figures, first compile the section above, then uncomment the
% figure you'd like to generate before running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Figure 2C 
        %Counter
        ticker = 0;
        ticker = ticker + 1;
        A = 0; 
        % % % % %Molecular Weight, Affinity and Fractional Activity Arrays
        elem =10;                           %number of elements in the array iteration
        MW_array = logspace(0, 3, elem);    %k iterations through for-loop
        Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
        Act_array = ones(elem);             
        Sys_Act_array = ones(elem);

        for k = 1:elem
                MW = MW_array(k);
                for g = 1:elem
                    Kd = Kd_array(g);
                    [p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
                    [t,y] = ode15s(@odefun,tspan,y0,options,p);
                    Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);
                    Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);
                end            
        end
        
        figure
        contourf(Kd_array,MW_array,Act_array)
        colormap('jet')
        caxis([0, 1.4])
        colorbar('eastoutside')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca, 'fontsize', 18)            
        xlabel('Collagen Affinity (M)')
        ylabel('Molecular Weight (kDa)')
        set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
        title('Figure 2b')
        hold on
        scatter([13.8E-9, 130E-9, 5E-6 221E-9, 5E-6], [137, 137, 137, 32, 32],100, 'ok','MarkerFaceColor', 'w')
        hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplementary Figure 8c
%         figure 
%         contourf(Kd_array,MW_array,Sys_Act_array)
%         colorbar('eastoutside')
%         set(gca,'XScale','log')
%         set(gca,'YScale','log')
%         set(gca, 'fontsize', 18)            
%         xlabel('Collagen Affinity (M)')
%         ylabel('Molecular Weight (kDa)')
%         title('Supplementary Figure 9')
%         set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplemental Figure 7a
%         [p, y0_LLM] = Inputs(14E-9,137,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLLM,yLLM] = ode15s(@odefun,tspan,y0_LLM,options,p);
% 
%         [p, y0_LLxM] = Inputs(130E-9,137,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLLxM,yLLxM] = ode15s(@odefun,tspan,y0_LLxM,options,p);
% 
%         [p, y0_LxLxM] = Inputs(5E-6,137,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLxLxM,yLxLxM] = ode15s(@odefun,tspan,y0_LxLxM,options,p);
% 
%         [p, y0_L] = Inputs(221E-9,32,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tL,yL] = ode15s(@odefun,tspan,y0_L,options,p);
% 
%         [p, y0_Lx] = Inputs(5E-6,32,[],Th_C,0,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [tLx,yLx] = ode15s(@odefun,tspan,y0_L,options,p);
% 
%         L_tumorLLM = yLLM(:,1) + yLLM(:,3) + yLLM(:,7) + yLLM(:,8);
%         L_tumorLLxM = yLLxM(:,1) + yLLxM(:,3) + yLLxM(:,7) + yLLxM(:,8);
%         L_tumorLxLxM = yLxLxM(:,1) + yLxLxM(:,3) + yLxLxM(:,7) + yLxLxM(:,8);
%         L_tumorL = yL(:,1) + yL(:,3) + yL(:,7) + yL(:,8);
%         L_tumorLx = yLx(:,1) + yLx(:,3) + yLx(:,7) + yLx(:,8);
% 
%         L0 = 1E-6;
%         ID_LLM = L_tumorLLM/(L0)*100;
%         ID_LLxM = L_tumorLLxM/(L0)*100;
%         ID_LxLxM = L_tumorLxLxM/(L0)*100;
%         ID_L = L_tumorL/(L0)*100;
%         ID_Lx = L_tumorLx/(L0)*100;
% 
%         figure
%         subplot(5,1,1), plot(tLLM/(60*60), ID_LLM)
%         ylim([0, 100])
%         title('Supplemental Figure 7a')
%         subplot(5,1,2), plot(tLLxM/(60*60), ID_LLxM)
%         ylim([0, 100])
%         subplot(5,1,3), plot(tLxLxM/(60*60), ID_LxLxM)
%         ylim([0, 100])
%         subplot(5,1,4), plot(tL/(60*60), ID_L)
%         ylim([0, 100])
%         subplot(5,1,5), plot(tLx/(60*60), ID_Lx)
%         ylim([0, 100])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplemental Figure 7b
%         EC50tumor = koff_R/kon_R;
%         A_LLM = 1./(1+EC50tumor./L_tumorLLM);
%         A_LLxM = 1./(1+EC50tumor./L_tumorLLxM);
%         A_LxLxM = 1./(1+EC50tumor./L_tumorLxLxM);
%         A_L = 1./(1+EC50tumor./L_tumorL);
%         A_Lx = 1./(1+EC50tumor./L_tumorLx);
% 
%         figure
%         subplot(5,1,1), plot(tLLM/(60*60), A_LLM)
%         ylim([0, 1])
%         title('Supplemental Figure 7b')
%         subplot(5,1,2), plot(tLLxM/(60*60), A_LLxM)
%         ylim([0, 1])
%         subplot(5,1,3), plot(tLxLxM/(60*60), A_LxLxM)
%         ylim([0, 1])
%         subplot(5,1,4), plot(tL/(60*60), A_L)
%         ylim([0, 1])
%         subplot(5,1,5), plot(tLx/(60*60), A_Lx)
%         ylim([0, 1])
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% %Supplementary Figure 9b
%         LIDcircV = 5E-8;
%         L_circLLM = 100*yLLM(:,4)/(LIDcircV);
%         L_circLLxM = 100*yLLxM(:,4)/(LIDcircV);
%         L_circLxLxM = 100*yLxLxM(:,4)/(LIDcircV);
%         L_circL = 100*yL(:,4)/(LIDcircV);
%         L_circLx = 100*yLx(:,4)/(LIDcircV);
%         
%         figure
%         subplot(5,1,1), plot(tLLM/(60*60), L_circLLM)
%         %ylim([0, 100])
%         title('Supplemental Figure 9b')
%         subplot(5,1,2), plot(tLLxM/(60*60), L_circLLxM)
%         %ylim([0, 100])
%         subplot(5,1,3), plot(tLxLxM/(60*60), L_circLxLxM)
%         %ylim([0, 100])
%         subplot(5,1,4), plot(tL/(60*60), L_circL)
%         %ylim([0, 100])
%         subplot(5,1,5), plot(tLx/(60*60),  L_circLx)
%         %ylim([0, 100])
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplemental Figure 8a
%         L0 = 1E-6;
%         L_RInterLLM = yLLM(:,9)/L0*100;
%         L_RInterLLxM = yLLxM(:,9)/L0*100;
%         L_RInterLxLxM = yLxLxM(:,9)/L0*100;
%         L_RInterL = yL(:,9)/L0*100;
%         L_RInterLx = yLx(:,9)/L0*100;
%         
%         figure
%         subplot(5,1,1), plot(tLLM/(60*60), L_RInterLLM)
%         %ylim([0, 100])
%         title('Supplemental Figure 8')
%         subplot(5,1,2), plot(tLLxM/(60*60), L_RInterLLxM)
%         %ylim([0, 100])
%         subplot(5,1,3), plot(tLxLxM/(60*60), L_RInterLxLxM)
%         %ylim([0, 100])
%         subplot(5,1,4), plot(tL/(60*60), L_RInterL)
%         %ylim([0, 100])
%         subplot(5,1,5), plot(tLx/(60*60),  L_RInterLx)
%         %ylim([0, 100])
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Supplementary Figure 8b
%         ticker = 0;
%         ticker = ticker + 1;
%         A = 0 
%         % % % % %Molecular Weight, Affinity and Fractional Activity Arrays
%         elem =10; %40                          %number of elements in the array iteration
%         MW_array = logspace(0, 3, elem);    %k iterations through for-loop
%         Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
%         Act_array = ones(elem);
% 
%         for k = 1:elem
%                 MW = MW_array(k);
%                 for g = 1:elem
%                     Kd = Kd_array(g);
%                     [p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,1,1,0,kon_R,koff_R,0,L,C);
%                     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%                     Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);
%                 end            
%         end
% 
%         figure
%         contourf(Kd_array,MW_array,Act_array)
%         colormap('jet')
%         %caxis([0,2.5])
%         colorbar('eastoutside')
%         set(gca,'XScale','log')
%         set(gca,'YScale','log')
%         set(gca, 'fontsize', 18)            
%         xlabel('Collagen Affinity (M)')
%         ylabel('Molecular Weight (kDa)')
%         set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
%         title('Supplementary Figure 8b')        
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Figure 10 Sensitivity Analysis
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(1E-9,500,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% figure
% subplot(3,3,1), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(100E-9,500,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,2), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(5E-6,500,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,3), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(1E-9,100,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,4), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(100E-9,100,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,5), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(5E-6,100,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,6), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(1E-9,10,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,7), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(100E-9,10,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,8), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])
% 
% elem = 6; %has to be an even number for integer logspacing
% Iog = [[],[],[],Th_C,[],[],cells,NR,expand,kon_R,koff_R,kendo_R,L,C]
% Act_array = ones(length(Iog),elem+1)
% for i = 1:length(Iog)
%     %subplot(3,5,i), title(num2str(i))
%     k = 0
% for n = logspace(-elem/2, elem/2,elem+1)
%     %for n = [0.0001, 0.001, 0.01, 0.1, 1]
%     k = k+1;
%     Imat = Iog;
%     Imat(i) = Imat(i)*n;
%     [p, y0] = Inputs(5E-6,10,[],Imat(1),1,1,Imat(2),Imat(3),Imat(4),Imat(5),Imat(6),Imat(7),Imat(8),Imat(9));
%     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%     Act_array(i,k) = tumorA(t,y,kon_R,koff_R,tspan);  
% end
% end
% subplot(3,3,9), heatmap(logspace(-elem/2, elem/2,elem+1), {'Collagen Half-Life','Tumor IL-2R+ cells/mL','IL-2R per cell', 'IL-2R+ Cell Expansion Rate','k_o_n_,_R','k_o_f_f_,_R','k_i_n_t_e_r_,_R','[Injected Protein]_0_,_t_u_m_o_r','[Collagen]_0_,_t_u_m_o_r'}, Act_array, 'CellLabelFormat','%.1f')
% colormap('jet')
% caxis([0,1.4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Figure 2C Exposure Instead
%         %Counter
%         ticker = 0;
%         ticker = ticker + 1;
%         A = 0; 
%         % % % % %Molecular Weight, Affinity and Fractional Activity Arrays
%         elem =100;                           %number of elements in the array iteration
%         MW_array = logspace(0, 3, elem);    %k iterations through for-loop
%         Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
%         Exp_array = ones(elem);             
%         %Sys_Act_array = ones(elem);
% 
%         for k = 1:elem
%                 MW = MW_array(k);
%                 for g = 1:elem
%                     Kd = Kd_array(g);
%                     [p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%                     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%                     Exp_array(k,g) = tumorExposure(t,y,kon_R,koff_R,tspan);
%          %           Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);
%                 end            
%         end
% % 
% %         [Xkd,Ymw] = meshgrid(1:elem,1:elem)
% %         [Xkd2,Ymw2] = meshgrid(1:0.01:elem, 1:0.1:elem);
% %         outData = interp2(Xkd, Ymw, Act_array, Xkd2, Ymw2, 'linear');
% % 
% % 
% %         
%         figure
%         contourf(Kd_array,MW_array,Exp_array,1000,'edgecolor','none')
%         colormap('jet')
%         %caxis([0,2.5])
%         colorbar('eastoutside')
%         set(gca,'XScale','log')
%         set(gca,'YScale','log')
%         set(gca, 'fontsize', 18)    
%         xlabel('Collagen Affinity (M)')
%         ylabel('Molecular Weight (kDa)')
%         set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
%         title('Half Life Hours')
%         hold on
%         scatter([20E-9, 158E-9, 5E-6 102E-9, 5E-6, 5E-6, 5E-6], [95, 95, 95, 15, 15,150,1],100, 'ok','MarkerFaceColor', 'w')
%         hold off
        
        
%         [p, y0] = Inputs(20E-9,95,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act1 = tumorExposure(t,y,kon_R,koff_R,tspan)
%         
%         [p, y0] = Inputs(158E-9,95,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act2 = tumorExposure(t,y,kon_R,koff_R,tspan)
%         
%         [p, y0] = Inputs(5E-6,95,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act3 = tumorExposure(t,y,kon_R,koff_R,tspan)
%         
%         [p, y0] = Inputs(102E-6,15,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act4 = tumorExposure(t,y,kon_R,koff_R,tspan)
% 
%         [p, y0] = Inputs(5E-6,15,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act5 = tumorExposure(t,y,kon_R,koff_R,tspan)
% 
%         [p, y0] = Inputs(5E-6,150,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act6 = tumorExposure(t,y,kon_R,koff_R,tspan)
% 
%         [p, y0] = Inputs(5E-6,1,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act7 = tumorExposure(t,y,kon_R,koff_R,tspan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %EIIIB Figure
%         %Counter
%         C = 2E-8;
%         ticker = 0;
%         ticker = ticker + 1;
%         A = 0; 
%         % % % % %Molecular Weight, Affinity and Fractional Activity Arrays
%         elem =20;                           %number of elements in the array iteration
%         MW_array = logspace(0, 3, elem);    %k iterations through for-loop
%         Kd_array = logspace(-10,-5,elem);   %g iterations through for-loop
%         Act_array = ones(elem);             
%         Sys_Act_array = ones(elem);
% 
%         for k = 1:elem
%                 MW = MW_array(k);
%                 for g = 1:elem
%                     Kd = Kd_array(g);
%                     [p, y0] = Inputs(Kd,MW,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%                     [t,y] = ode15s(@odefun,tspan,y0,options,p);
%                     Act_array(k,g) = tumorA(t,y,kon_R,koff_R,tspan);
%                     Sys_Act_array(k,g) = sysA(t,y,kon_R,koff_R,tspan);
%                 end            
%         end
%         
%         figure
%         contourf(Kd_array,MW_array,Act_array)
%         colormap('jet')
%         caxis([0, 1.4])
%         colorbar('eastoutside')
%         set(gca,'XScale','log')
%         set(gca,'YScale','log')
%         set(gca, 'fontsize', 18)            
%         xlabel('Collagen Affinity (M)')
%         ylabel('Molecular Weight (kDa)')
%         set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
%         title('Figure 2b')
%         hold on
%         scatter([1E-9, 5E-6, 1E-9, 5E-6], [98, 98, 30, 30],100, 'ok','MarkerFaceColor', 'w')
%         hold off
%         
%         [p, y0] = Inputs(1E-9,98,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act1 = tumorA(t,y,kon_R,koff_R,tspan)
%         
%         [p, y0] = Inputs(1E-9,30,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act2 = tumorA(t,y,kon_R,koff_R,tspan)
%         
%     [p, y0] = Inputs(5E-6,98,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act3 = tumorA(t,y,kon_R,koff_R,tspan)
%         
%         [p, y0] = Inputs(5E-6,30,[],Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C);
%         [t,y] = ode15s(@odefun,tspan,y0,options,p);
%         Act4 = tumorA(t,y,kon_R,koff_R,tspan)
%         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Additional Functions Used in this Code
%Fractional Activity in Tumor Function
    function cA = tumorA(t1,y1,kon_R, koff_R,tspan)
        EC50tumor = koff_R/kon_R;
        L_tumor = y1(:,1) + y1(:,3) + y1(:,7) + y1(:,8);
        Activity_L = 100./(1+EC50tumor./L_tumor);
        cA = trapz(t1, Activity_L)/(tspan(end)*100)*tspan(end)/(60*60*24);
    end

    function thalf = tumorExposure(t1,y1,kon_R, koff_R,tspan)
        EC50tumor = koff_R/kon_R;
        L_tumor = (y1(:,1) + y1(:,3) + y1(:,7) + y1(:,8))/(1E-6);
        f = fit(t1,L_tumor,'exp1');
        thalf = (log(2)/-f.b)/(60*60);
    end

%Seconds to Days Function
    function day = d(t)   %convert s to days
        day = t/(24*60*60);
    end

%Fractional Activity in Circulation Function
 function cA = sysA(t1,y1,kon_R,koff_R,tspan)
        EC50sys = koff_R/kon_R;
        L_sys = y1(:,4);
        Activity_sys = 100./(1+EC50sys./L_sys);
        cA = trapz(t1, Activity_sys)/(tspan(end)*100)*tspan(end)/(60*60*24);
    end




end

function ydot = odefun(~,y,p)
% Collect param values in cell array and redefine params with names
paramsCell = num2cell(p);
[k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17]=paramsCell{:};

% Collect y-vals in cell array and redefine y-vals with names
yCell = num2cell(y);
[L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng, R] = yCell{:};

dL = k2*LC - k1*C*L/k5 - k12*L*(R)/k5 + k13*LR + k3*k5*(Lesc-L)-k17*L;
dCint = k6*C;
dLC = k1*C*L/k5 -k2*LC - k6*A*LC - k12*LC*(R)/k5 + k13*LCR;
dLCint = k6*A*LC;
dLesc = k3*k5*(L/k5-Lesc)*(100E-6*k5/2E-3)-k4*Lesc- k17*Lesc; 
dLR = k12*L*(R)/k5  - k13*LR - k1*LR*C/k5 + k2*LCR - k14*(LR);
dLCR = k12*LC*(R)/k5 - k13*LCR - k2*LCR + k1*LR*C/k5;
dC = k2*LC - k1*C*L/k5  + dCint + dLCint - k6*C + k2*LCR - k1*LR*C/k5;
dLRint = k14*(LR);
dReng= LR+LCR;
dR = dLRint + k16*(dReng+R)*(1-(dReng+R)/(k15*10000)) - k12*L*R/k5 + k13*LR -k12*LC*R/k5 + k13*LCR;

ydot = [dL; dC; dLC; dLesc*k5; dCint; dLCint; dLR; dLCR; dLRint; dReng; dR].*(1/k5);
end

