function [f1]=Circies_Model()
    p = setup();
    %% Simulation parameters
    %Timespan
    tfin = 60*6;   %simulation final time
    step = 0.01;    %simulation step
    tspan = 0:step:tfin-step; 
    
    % options for ode function
    opti = odeset('AbsTol',1e-8,'RelTol',1e-6); %Sets Error Tolerance 
    
    %initial conditions
    Init = [1*10^14, 0, 0, 0, p.CN, 0, 0, 0];% AHL_ext, AHL_int, QscR(irrelevant), AHL.QscR-complex, Promotor Number, PAQ complex, mRNA, Protein
    
    [t0,x0] = ode23t(@(t,x) ODES(t,x,p),tspan, Init, opti);    

    
% %     % Figure 1 AHL Diffusion through Membrane
%     time = t0*60;
%     f1 = figure('Color',[1 1 1],'Name','AHL Diffusion through Membrane old');
%     
%     yyaxis left;
%     plot(time,x0(:,2),'LineWidth',2);
%     ylabel({'# of molecules inside'});
%     
%     yyaxis right;
%     plot(time,log10(x0(:,1)),'LineWidth',2);
%     ylabel({'# of molecules outside (log10)'});
%     
%     xlabel('Time (sec)')
%     legend('AHL_{int}','AHL_{ext}')
%     grid on
%     
% %     % Figure 2 AHL.QscR.complex formation
%     time = t0;
%     f2 = figure('Color',[1 1 1],'Name','AHL.QscR.complex formation');
%     
%     yyaxis left;
%     plot(time,x0(:,2),'LineWidth',2);
%     hold on
%     plot(time,x0(:,4),'LineWidth',2);
%     hold off
%     ylabel({'# of molecules inside'});
%     
%     yyaxis right;
%     plot(time,log10(x0(:,1)),'LineWidth',2);
%     ylabel({'# of molecules outside (log10)'});
%     
%     xlabel('Time (min)')
%     legend('AHL_{int}','AHL.QscR','AHL_{ext}')
%     grid on
% %     
% %     Figure 3 PAQ complex creation
%     time = t0;
%     f3 = figure('Color',[1 1 1],'Name','PAQ complex creation');
%     
%     plot(time,x0(:,2),'LineWidth',2);
%     hold on
%     plot(time,x0(:,4),'LineWidth',2);
%     plot(time,x0(:,5),'LineWidth',2);
%     plot(time,x0(:,6),'LineWidth',2);
%     hold off
%     
%     ylabel({'# of molecules inside'});
%     ylim([0,20])
%     xlabel('Time (min)')
%     legend('AHL_{int}','AHL.QscR', 'FreePromotor','PAQ')
%     grid on
%     
% %     Figure 4 DNA Transikription
%     time = t0;
%     f4 = figure('Color',[1 1 1],'Name','DNA Transkription -> mRNA');
%     
%     plot(time,x0(:,4),'LineWidth',2);
%     hold on
%     plot(time,x0(:,5),'LineWidth',2);
%     plot(time,x0(:,6),'LineWidth',2);
%     plot(time,x0(:,7),'LineWidth',2);
%     hold off
%     
%     ylabel({'# of molecules inside'});
%     ylim([0,20])
%     xlabel('Time (min)')
%     legend('AHL.QscR', 'FreePromotor','PAQ','mRNA')
%     grid on
%     
%     % Figure 5 mRNA Translation
%     time = t0/60;
%     f5 = figure('Color',[1 1 1],'Name','mRNA Translation -> Protein');
%     
%     plot(time,x0(:,5),'LineWidth',2);
%     hold on
%     plot(time,x0(:,6),'LineWidth',2);
%     plot(time,x0(:,7),'LineWidth',2);
%     hold off
%     ylabel({'# of molecules inside'});
%     ylim([0,20])
%     
%     yyaxis right
%     plot(time,x0(:,8),'LineWidth',2);
%     
%     ylabel({'Protein(# of molecules)'})
%     
%     xlabel('Time (hr)')
%     legend('FreePromotor','PAQ','mRNA', 'protein')
%     grid on
%     
% %     Figure 6 AHL_ext vs Protein
%     time = t0/60;
%     f6 = figure('Color',[1 1 1],'Name','AHL_{ext} vs Protein');
%     
%     plot(time,x0(:,1),'LineWidth',2);
%     hold on
%     plot(time,x0(:,8),'LineWidth',2);
%     hold off
%     
%     ylabel({'# of molecules)'});
%     
%     xlabel('Time (hr)')
%     legend('AHL_{ext}','protein')
%     grid on
    % Figure 7 Protein vs AHL Concentration ODE
%     tic;
%     start_conc = 14; %log10 of molecules per ml
%     fin_conc = 17; %log10 x molecules per ml
%     step = 0.05;
%     
%     conc_span = start_conc:step:fin_conc;
%     
%     result=zeros(1,numel(conc_span));
%     
%     for i = 1:numel(conc_span) % For each defined conc. step
%         test_conc = 10^conc_span(i); % Calculate the number of molecules per ml
%         Init(1,1) = test_conc; %Set the initial concentration to the initial value
%         [~,x0] = ode23t(@(t,x) ODES(t,x,p),tspan, Init, opti); %and run the whole model for this concentration
%         result(i) = x0(end,8); %Save the last value of the protein ammount (column 8)
%     end
%     
%     f7 = figure('Color',[1 1 1],'Name','# Protein vs #AHL_ext/ml');
%     plot(conc_span,result,'LineWidth',2);
%     
%     ylabel('# Protein')
%     xlabel('# AHL_{ext} (log10)')
%     grid on
%     disp("Elapsed time for Fig7 " + toc + " seconds.");
%     Figure 8 Protein vs AHL: ODE vs Hill
    tic;
    start_conc = 13; %log10 of molecules per ml
    fin_conc = 17; %log10 x molecules per ml
    step = 0.05;
    
    conc_span = start_conc:step:fin_conc;
    
    result_ode=zeros(1,numel(conc_span));
    result_hill_old = zeros(1,numel(conc_span));
    result_hill_new = zeros(1,numel(conc_span));
    
    
    for i = 1:numel(conc_span) % For each defined conc. step
        test_conc = 10^conc_span(i); % Calculate the number of molecules per ml
        Init(1,1) = test_conc; %Set the initial concentration to the initial value
        [~,x0] = ode23t(@(t,x) ODES(t,x,p),tspan, Init, opti); %and run the whole model for this concentration
        result_ode(i) = x0(end,8); %Save the last value of the protein ammount (column 8)
        
        test_conc = test_conc * 6.5*10^(-16);
        max_expr = ((p.k_transk * p.k_transl * p.CN) / (p.d_mRNA * p.d_protein));        
        affinity = (p.k_AHL_off / p.k_AHL_on) * (p.k_AQ_off / p.k_AQ_on) * (p.CN / test_conc^p.k_AQ_nh);
        result_hill_old(i) = max_expr * (p.basal + (1 - p.basal) * ((test_conc^p.k_AQ_nh) / (affinity^p.k_AQ_nh + test_conc^p.k_AQ_nh)));
        affinity = (p.k_AQ_off / p.k_AQ_on) * (p.k_AHL_off / p.k_AHL_on);
        result_hill_new(i) = max_expr * (p.basal + (1 - p.basal) * ((test_conc^p.k_AQ_nh) / (affinity^p.k_AQ_nh + test_conc^p.k_AQ_nh)));
    end
    deviation = zeros(1,numel(conc_span));
    for i = 1:numel(conc_span)
        deviation(i) = (result_hill_new(i) / result_ode(i));
    end
    
    f8 = figure('Color',[1 1 1],'Name','# Protein vs #AHL_ext/ml');
    plot(conc_span,result_ode,'LineWidth',2);
    hold on
    plot(conc_span,result_hill_old,'LineWidth',2);
    plot(conc_span,result_hill_new,'LineWidth',2);
    hold off
    
    legend('ODE','HILL_{old}','HILL_{new}')
    
    ylabel('# Protein')
    xlabel('# AHL_{ext} (log10)')
    grid on
    disp("Elapsed time for Fig8 " + toc + " seconds.");
    disp("std" + std(deviation));


end

function [dxdt]=ODES(t,x,p)
    %x(1) AHL_ext = changing
    dxdt(1,1) = -p.k_diff * p.rho * x(1) + p.k_diff * x(2);
    
    %x(1) AHL_ext = const.
    %dxdt(1,1) = 0; % Assuming ext. conc. constant
    
    %x(2) AHL_int
    dxdt(2,1) = +p.k_diff * p.rho * x(1) - p.k_diff * x(2); %Influence of the diffusion
%         -p.k_AHL_on * x(2) + p.k_AHL_off * x(4); %Influence of the AHL.QscR Complex
    
    %x(3) QscR conc = const.
    %QscR conc is currently not used, we assume there is enugh.
    dxdt(3,1) = 0;
    
    %x(4) AHL.QscR complex
    dxdt(4,1) = +p.k_AHL_on * x(2)^p.k_AHL_nh - p.k_AHL_off * x(4)...; %Influence of AHL_int concentration
        -p.k_AQ_on * x(5) * x(4)^p.k_AQ_nh +p.k_AQ_off * x(6); %Influence of the Binding to Promotor

    %x(5) Number of free Promotors
    dxdt(5,1) = -p.k_AQ_on * x(5) * x(4)^p.k_AQ_nh +p.k_AQ_off * x(6); %Influence of the PAQ complex
    
    %x(6) PAQ complex
    dxdt(6,1) = +p.k_AQ_on * x(5) * x(4)^p.k_AQ_nh -p.k_AQ_off * x(6);
    
    %x(7) mRNA
    dxdt(7,1) = +p.k_transk * x(6) + p.k_transk * p.basal * x(5) - p.d_mRNA * x(7);
    
    %x(8) protein
    dxdt(8,1) = +p.k_transl * x(7) - p.d_protein * x(8);
end

function [p] = setup()
    %%Parameters for E.coli
    
    %local density, (cells/ml)/OD600 E. coli
    %OD600 of ON Culture ~ 2
    %Data: http://book.bionumbers.org/what-is-the-concentration-of-bacterial-cells-in-a-saturated-culture/
    
    local_density = 1.3 * 10^9;
    OD600 = 2;
    
    %volume of a cell 0.65um^3/cell 
    v_cell = 6.5 * 10^(-16); % L
    
    %{
    Rho is the factor that describes the volume percent in a culture that
    is occupied by cells. This gives a conversion factor from external
    concentration to internal concentration.
    
    Dockery, J. D., Keener, J. P. 2001 A mathematical model for quorum
    sensing in Pseudomonas aeruginosa. Bull Math Biol. 63(1):95-116.
    https://pubmed.ncbi.nlm.nih.gov/11146885/

    rho = (cells / 1 ml) * (ml / cell)
    %}
    
%     p.rho = local_density * v_cell * OD600;
    p.rho = v_cell;
    
    %{
    D is the diffusion coefficient through the membrane [min^-1].
    D = S * P / V_cell, where S is the surface area of the cell,
    P is the Permeability of the cell, and V is the volume of the cell.
    
    https://www.biorxiv.org/content/10.1101/106229v1.full.pdf Page 7

    For E. coli with S = 22 µm^2, P = 3 * exp(-3) µm/min, V_cell = 0.65µl^3,
    D=5
    %}
    
    p.k_diff = 5;
    
    %apparent dissociation constant for AHL + QscR (AQ) [half saturation molecules]
    %~ 3.1nM https://journals.asm.org/doi/10.1128/JB.01041-10
    %Our volume is the vol of one cell: 6,5*10^-16 l/cell * 3.1*10^-9 mol/l -> 2.015*10^-24 mol/cell
    %Molecules: 2.015*10^-24 mol/cell * 6.022*10^23 molecules/mol = 1.213433 molecules / cell
    p.k_AHL_on = 0.1; % Fixme: arbitrary, has to be determined in the lab and ajusted for actual time
    p.k_AHL_off = 0.1213433; % k_AHL_on * 121.3433 NOTE: this one remains correct in relation too k_AHL_on
    
    %Hill Coeficient for AHL binding to QscR https://journals.asm.org/doi/10.1128/JB.01041-10
    p.k_AHL_nh = 1.0; % Expected Value from Paper for Hill coeficient
    
    %apparent dissociation constant for Promotor + AQ [half saturation molecules]
    %~ 2.2nM https://journals.asm.org/doi/10.1128/JB.01041-10
    %6,5*10^-16 * 2.2*10^-9 mol/l l/cell -> 1.43x10^-24 mol/cell
    %Molecules: 1.43x10^-24 mol/cell * 6.022*10^23 molecules/mol = 0.861146 molecules / cell
    p.k_AQ_on = 0.1; % Fixme: arbitrary, has to be determined in the lab and ajusted for actual time
    p.k_AQ_off = 0.0861146; % k_AQ_on * 0.861146
    
    %Hill Coeficient for AHL.QscR complex bonding to promotor https://journals.asm.org/doi/10.1128/JB.01041-10
    p.k_AQ_nh = 1.8;
    
   %Copy Number of Plasmid
    p.CN = 17;
    
    % Transkription rate [au] https://www.denovodna.com:4433/shared/dvXq9kDxF9aDKfUh7hxuPZoZldJ89EQz
    %p.k_transk = 397.39; % FIXME: Arbitrary Units, find ref value (nt/s maybe)
    
    % Translation rate [au] https://www.denovodna.com:4433/shared/3qEKCv9GJEMkSPjyRiyWOSrU6pw2nsIM
    % p.k_transl = 21.56; %FIXME: Arbitrary Units, find ref value (aa/s maybe)
    
    % Transkription rate in E. coli ~ 10-100 nt/s http://book.bionumbers.org/what-is-faster-transcription-or-translation/
    % 600 - 6000 nt/min
    % EGFP length ~ 714
    gene_length = 714;
    transk_rate = 2700; % FIXME: median
    p.k_transk = 1/(gene_length / transk_rate); % Einheit mRNA pro Minute
    
    % Translation rate in E. coli ~ 10-20 aa/s http://book.bionumbers.org/what-is-faster-transcription-or-translation/
    % 600 - 1200 aa/min
    prot_length = gene_length / 3;
    transl_rate = 900; % FIXME: median
    p.k_transl = 1/(prot_length/transl_rate); % Einheit Proteine pro Minute
    
    % mRNA decay constant [min^-1]
    % T_1/2 of mRNA in E. coli K10: 2.1 min https://www.ncbi.nlm.nih.gov/pmc/articles/PMC365694/table/tbl1/
    % d_mRNA = ln(2)/T_1/2
    
    p.d_mRNA = 0.3301;
    
    % protein decay constat [min^-1]
    % T_1/2 of an abnormal protein in E. coli ~ 45 min https://bionumbers.hms.harvard.edu/files/Protein%20half-lives%20in%20E.%20coli.pdf
    % d_protein = ln(2)/T_1/2
    
    p.d_protein = 0.015;
    
    % Basal expression [%]
    % Basal expression is the result of a leaky promotor
    basal_percent = 20;
    p.basal = basal_percent / 100;
   
    
end