%% ========================================================================
%  MAIN SCRIPT - SC PROJECT NELAHTP 2025/2026
%  Solenoide HTS CICC stacked-slotted a double-pancake @ 15 T / 50 kA
%  ------------------------------------------------------------------------
%  STRUTTURA
%    STEP 1 : Data acquisition and measurement analysis
%    STEP 2 : Geometry design del DP centrale -> matrici P ed E DEFINITIVE
%             (+ visualizzazione 3D macro/micro)
%    STEP 3 : Calcolo del campo magnetico (+ T_cs e margine di temperatura)
%    STEP 4 : Sistema matriciale DAE (+ dimostrazione dell'equiripartizione)
%    STEP 5 : CICC layout - sezione trasversale e centroide (R_core derivato)
%    STEP 6 : AC losses - sorgente di calore Q(s) per il termico
%    STEP 7 : Cooling needs - GANDALF ridotto slow-transient 1D
%             -> convergenza della portata minima, frontiera T_in-m_dot,
%                margine di 5 K e chiusura Loop B
%
%  NOMENCLATURA
%    TURN  = giro della spirale radiale DENTRO un double-pancake (r cresce, z fisso)
%    LAYER = un double-pancake impilato assialmente rispetto agli altri
%
%  GRADI DI LIBERTA'
%    N_layer : FISSATO dall'altezza del solenoide
%    N_turn  : GRADO DI LIBERTA' -> iterato per rispettare B_pk = 15 T
%
%  DIPENDENZE (stessa cartella o sul path MATLAB)
%    dataAcquisition.m, openItems.m                      (Step 1)
%    segmentMagneticField1d.m, Mu0.m,
%    inductanceMatrix1d.m, directSolve.m                 (SCfiles, Prof. Freschi)
%
%  PURE SCRIPT: nessuna definizione di funzione prima della fine del file.
% ========================================================================
clear; clc; close all;

% ------------------------------------------------------------------------
%  PROFILING. Alcuni blocchi costano minuti: la matrice delle induttanze
%  (Neumann, integral2 su ogni coppia di archi, O(Ne^2)), l'integrazione
%  temporale del DAE (25 passi impliciti x Picard, con L densa nel blocco
%  M/dt) e lo Step 7 (ricerca portata, convergenza e frontiera operativa;
%  termo-idraulici completi). Senza misurarli non si sa dove intervenire se
%  il run diventa troppo lento, e si finisce per ottimizzare la parte
%  sbagliata. Ogni blocco pesante viene cronometrato e il riepilogo finale
%  li ordina per costo.
%  I timer NON sono annidati: ogni sezione e' misurata una volta sola, cosi'
%  la somma delle voci si confronta direttamente col tempo totale e la
%  differenza (grafici, stampe, blocchi leggeri) resta visibile come residuo.
tmr        = struct('name',{},'sec',{});
tWall      = tic;                       % cronometro dell'intero script

%  CHIUSURA DEL LOOP B IN DUE RUN. T_op serve allo Step 3 ma e' il risultato
%  dello Step 7: e' un punto fisso. Il primo run parte dalla derivazione dello
%  Step 3 (il massimo compatibile con il margine anche a lift factor -20%),
%  la 7.6 dice quanto vale davvero, e il secondo run lo impone qui. Lasciare
%  [] per il primo giro; mettere il valore stampato dalla 7.6 per il secondo.
%  Due run bastano: T_op entra solo in Ic (perdite di isteresi) e in n, ed
%  entrambi dipendono da T molto piu' debolmente di quanto T dipenda da loro.
T_op_override = 16.46;    % [K] oppure [] per usare la derivazione dello Step 3

fprintf('=========================================================\n');
fprintf(' SC PROJECT NELAHTP 2025/2026 - main\n');
fprintf('=========================================================\n');


%% ========================================================================
%  STEP 1 - DATA ACQUISITION AND MEASUREMENT ANALYSIS
%  ========================================================================
fprintf('\n\n########## STEP 1 - DATA ACQUISITION ##########\n');

data = dataAcquisition(true);
openItems(data);

I_nom    = data.op.I_rated.val;      % [A]   50 kA
R_in     = data.op.ID_wind.val/2;    % [m]   0.5 m
H_max    = data.op.H_max.val;        % [m]   1.5 m
B_target = data.op.B_target.val;     % [T]   15 T
dIdt     = data.op.dIdt.val;         % [A/s] 200 A/s
dT_req   = data.op.dT_margin.val;    % [K]   5 K

%  DURATA DELLA RAMPA - definita QUI, dalla fisica dello scenario, perche' la
%  usano lo Step 4 (integrazione del DAE), lo Step 6 (dB/dt nelle perdite AC) e
%  lo Step 7 (finestra in cui la sorgente e' accesa). Definirla una volta sola
%  evita che tre Step la ricavino ciascuno per conto proprio e che una modifica
%  allo scenario ne aggiorni solo due su tre.
t_ramp = I_nom/dIdt;                 % [s] 250 s

fprintf('\n--- Parametri operativi estratti ---\n');
fprintf('  I_nom = %.0f kA | R_in = %.3f m | H_max = %.2f m | B_target = %.1f T\n', ...
    I_nom/1e3, R_in, H_max, B_target);
fprintf('  rampa: %.0f A/s -> t_ramp = %.0f s (usata da Step 4, 6 e 7)\n', dIdt, t_ramp);


%% ========================================================================
%  STEP 2 - GEOMETRY DESIGN DEL DP CENTRALE (matrici P ed E definitive)
%  ========================================================================
fprintf('\n\n########## STEP 2 - GEOMETRY DESIGN ##########\n');

% ------------------------------------------------------------------------
%  2.1  Parametri del cavo (PLACEHOLDER: da confermare nello Step 5 - layout)
% ------------------------------------------------------------------------
cab.width        = 0.025;   % [m] ingombro esterno CICC -> passo radiale tra turn
cab.N_slots      = 6;       % [-] slot per cross-section (stacked-slotted, ENEA-like)
cab.R_core       = 0.007276;% [m] centroide dello stack, coerente con la sezione
                            %     derivata nello Step 5 per N_tapes_slot = 80.
                            %     Se cambia la geometria dello stack, lo Step 5
                            %     ricalcola R_core e verifica la coerenza.
cab.N_twist      = 6;       % [-] giri di twist interi per turn.
                            %     La mesh DAE viene scelta in modo da mantenere
                            %     almeno 4 segmenti per periodo di twist.
cab.N_tapes_slot = 80;      % [-] nastri per slot.
                            %     La scelta riduce la corrente per nastro e
                            %     aumenta il margine elettromagnetico; la verifica
                            %     effettiva resta quella di Step 3 e Step 7.
cab.tapeNormal   = 'radial';% orientamento dello stack nello slot: 'radial'
                            % oppure 'azimuthal'.

% T_op_guess NON e' piu' un valore scelto a mano: viene CALCOLATO nello
% Step 3.3 come il MASSIMO valore che rispetta ancora il margine di 5 K
% richiesto dal brief, anche nello scenario peggiore di dispersione dei lift
% factor (-20%, dichiarata dal fornitore). Vedi Step 3.3 per la derivazione.

% ------------------------------------------------------------------------
%  GEOMETRIA DERIVATA DELLO STACK - tre grandezze DISTINTE, da non confondere:
%    R_core  = POSIZIONE: distanza del centroide dello stack dall'asse del cavo.
%              E' un dato di layout (PLACEHOLDER: lo fissa lo Step 5).
%    r0      = DIMENSIONE: raggio del conduttore equivalente che rappresenta lo
%              stack. Si DERIVA dalla sezione reale dello stack, non da R_core.
%    r_near  = bordo INTERNO dello stack = R_core - depth/2 (usato per il margine).
%  Confondere r0 con R_core rende i N_slots strand compenetrati: la mutua tra
%  adiacenti eguaglia l'autoinduttanza, L diventa singolare e le correnti del
%  DAE esplodono. r0 governa anche il self-field: B_self = mu0*I/(2*pi*r0).
depth_stack = cab.N_tapes_slot * data.tape.t_total.val;        % spessore stack
A_stack     = data.tape.width.val * depth_stack;               % sezione stack
cab.r0      = sqrt(A_stack/pi);                                % raggio equivalente
r_near      = cab.R_core - depth_stack/2;                      % bordo interno
sep_strand  = 2*cab.R_core*sin(pi/cab.N_slots);                % distanza tra strand

fprintf('\n--- 2.1b Geometria derivata dello stack ---\n');
fprintf('  sezione stack = %.2f x %.2f mm = %.2f mm^2\n', ...
    data.tape.width.val*1000, depth_stack*1000, A_stack*1e6);
fprintf('  r0 (dimensione) = %.3f mm | R_core (posizione) = %.2f mm | r_near = %.2f mm\n', ...
    cab.r0*1000, cab.R_core*1000, r_near*1000);
fprintf('  distanza tra strand adiacenti = %.3f mm -> r0/sep = %.3f\n', ...
    sep_strand*1000, cab.r0/sep_strand);
if r_near <= 0
    error('step2:stackTooThick', ...
        'Stack (%.2f mm) piu'' spesso di 2*R_core (%.2f mm): rivedere cab.R_core.', ...
        depth_stack*1000, 2*cab.R_core*1000);
end
if cab.r0 >= sep_strand/2
    error('step2:strandsOverlap', ...
        ['GEOMETRIA INCONSISTENTE: r0 = %.2f mm >= meta'' della distanza tra strand ' ...
         '(%.2f mm). Gli strand si compenetrano, L risulta singolare e le correnti ' ...
         'del DAE sono prive di senso. Aumentare cab.R_core o ridurre lo stack.'], ...
        cab.r0*1000, sep_strand/2*1000);
elseif cab.r0 > 0.4*sep_strand
    warning('step2:strandsClose', ...
        ['r0 = %.2f mm vicino a meta'' della distanza tra strand (%.2f mm): ' ...
         'approssimazione a filo sottile al limite. cab.R_core (%.1f mm) e'' ' ...
         'probabilmente troppo piccolo per %d nastri/slot -> rivedere nello Step 5.'], ...
        cab.r0*1000, sep_strand/2*1000, cab.R_core*1000, cab.N_tapes_slot);
end

disc.N_seg_field = 32;      % [-] lati/turn per la mappa di campo
%  Con N_twist = 6, N_seg_dae = 24 fornisce 4 segmenti per periodo di twist.
%  La matrice L costa circa O(Ne^2), quindi questo valore bilancia accuratezza
%  della discretizzazione del twist e costo di calcolo.
disc.N_seg_dae   = 24;      % [-] lati/turn per L e DAE

% ------------------------------------------------------------------------
%  2.1c  RISOLUZIONE DEL TWIST NELLA MESH
% ------------------------------------------------------------------------
%  La mesh deve avere abbastanza segmenti per periodo di twist, altrimenti
%  l'elica e' ALIASATA: gli strand non campionano uniformemente la sezione, non
%  risultano elettricamente equivalenti e l'equiripartizione fallisce.
%  Il criterio minimo adottato e' 4 segmenti per periodo di twist; valori
%  maggiori riducono l'aliasing ma aumentano rapidamente il costo di L.
segPerTwist_dae   = disc.N_seg_dae   / cab.N_twist;
segPerTwist_field = disc.N_seg_field / cab.N_twist;
fprintf('\n--- 2.1c Risoluzione del twist nella mesh ---\n');
fprintf('  mesh DAE  : %.1f segmenti per periodo di twist\n', segPerTwist_dae);
fprintf('  mesh campo: %.1f segmenti per periodo di twist\n', segPerTwist_field);
if segPerTwist_dae < 4
    warning('step2:twistAliased', ...
        ['Solo %.1f segmenti per twist nella mesh DAE: il twist e'' ALIASATO. Gli ' ...
         'strand non campionano uniformemente la sezione, non risultano equivalenti ' ...
         'e l''equiripartizione fallisce (scarto atteso >> 10%%). Aumentare ' ...
         'disc.N_seg_dae (costo di L ~ N^2) oppure ridurre cab.N_twist.'], segPerTwist_dae);
elseif segPerTwist_dae < 8
    fprintf('  NOTA: %.1f seg/twist sono sufficienti (scarto DAE atteso ~2%%).\n', segPerTwist_dae);
    fprintf('        Con disc.N_seg_dae = %d si scenderebbe sotto 1%% (costo di L x4).\n', ...
        2*disc.N_seg_dae);
end

% ------------------------------------------------------------------------
%  2.2  N_layer FISSATO dall'altezza del solenoide
% ------------------------------------------------------------------------
%  Ogni DP occupa 2*cable in z (due pancake affiancati). Si sceglie il massimo
%  numero DISPARI che sta in H_max: con un numero dispari esiste un layer
%  esattamente a z = 0, il piano equatoriale su cui il brief chiede di operare
%  il double-pancake.
%  N_layer E' UN VINCOLO, NON UN GRADO DI LIBERTA'.
%  Il brief dice "up to 1.5 m in height": l'altezza e' lo spazio disponibile e il
%  solenoide lo riempie. N_layer discende quindi INTERAMENTE dalla formula
%  (massimo dispari che sta in H_max) e l'UNICO grado di liberta' resta N_turn.

N_layer_max = floor(H_max / (2*cab.width));
N_layer_max = N_layer_max - (1 - mod(N_layer_max,2));   % massimo dispari ammesso
N_layer     = N_layer_max;      % DERIVATO da H_max: nessun valore scelto a mano
if mod(N_layer,2) == 0
    error('step2:evenLayers','N_layer deve essere DISPARI (serve un DP esattamente a z = 0).');
end
idxCentral  = (N_layer+1)/2;                            % indice del DP a z = 0
z_layers    = ((0:N_layer-1)' - (N_layer-1)/2) * 2*cab.width;

fprintf('\n--- 2.2 Impilamento assiale (N_layer fissato da H_max) ---\n');
fprintf('  H_max = %.2f m, spessore DP = %.0f mm -> N_layer = %d (dispari, DERIVATO)\n', ...
    H_max, 2*cab.width*1000, N_layer);
fprintf('  riempimento assiale = %.1f%% di H_max\n', ...
    N_layer*2*cab.width/H_max*100);
fprintf('  DP centrale: indice %d, z = %.4f m (deve essere 0)\n', ...
    idxCentral, z_layers(idxCentral));
fprintf('  estensione assiale occupata = %.3f m (<= %.2f m)\n', ...
    max(z_layers)-min(z_layers)+2*cab.width, H_max);

% ------------------------------------------------------------------------
%  2.3  N_turn: grado di liberta' iterato sul vincolo 15 T
% ------------------------------------------------------------------------
%  CRITERIO DI SCELTA. Il brief dice "sized to REACH a maximum magnetic field on
%  the conductor of 15 T": 15 T e' un OBIETTIVO DI PROGETTO, non un tetto da non
%  superare. N_turn e' intero e grossolano (un turn sposta B_pk di ~2 T), quindi
%  il target non e' centrabile esattamente: si sceglie il N_turn che minimizza
%  |B_pk - B_target|, non il massimo che sta sotto.
%  PERCHE' NON "il piu' grande sotto il target": con N_layer = 29 (formula) il
%  valore piu' vicino e' N_turn = 7 (B_pk ~ 15.1 T, +0.9%); il piu' grande
%  strettamente sotto sarebbe N_turn = 6 (~13.2 T, -12%). Scartare un progetto
%  che manca il target dello 0.9% per adottarne uno che lo manca del 12% e'
%  sbagliato: lo scarto residuo dello 0.9% e' peraltro dentro l'incertezza del
%  modello EM (filo sottile, B_self analitico, stack ridotto a un conduttore).
%  Il segno dello scarto viene comunque riportato, perche' un B_pk sopra target
%  consuma margine di temperatura e va tracciato fino allo Step 7.
fprintf('\n--- 2.3 Dimensionamento di N_turn sul target B_pk = %.1f T ---\n', B_target);

N_turn_list = 3:16;
Bpk_list = nan(size(N_turn_list));
t0 = tic;
for k = 1:numel(N_turn_list)
    tK = tic;
    Bpk_list(k) = peakFieldForNturn(N_turn_list(k), N_layer, idxCentral, z_layers, ...
                                    R_in, cab, disc.N_seg_field, I_nom);
    fprintf('   N_turn = %2d  ->  B_pk = %6.2f T  (scarto dal target %+6.2f %%)  [%.1f s]\n', ...
        N_turn_list(k), Bpk_list(k), (Bpk_list(k)-B_target)/B_target*100, toc(tK));
    if Bpk_list(k) > 1.15*B_target, break; end     % oltre: inutile continuare
end
tmr = addTimer(tmr, 'Step 2.3  scansione N_turn (Biot-Savart)', toc(t0));
valid = ~isnan(Bpk_list);
[~, kRel]  = min(abs(Bpk_list(valid) - B_target));
idxValid   = find(valid);
kSel       = idxValid(kRel);
N_turn = N_turn_list(kSel);
Bpk_design = Bpk_list(kSel);
errB   = (Bpk_design - B_target)/B_target*100;
fprintf('  --> SCELTO N_turn = %d  (B_pk = %.2f T, scarto %+.2f %% dal target)\n', ...
    N_turn, Bpk_design, errB);
if kSel > 1
    fprintf('      alternativa sotto: N_turn = %d -> %.2f T (%+.2f %%)\n', ...
        N_turn_list(kSel-1), Bpk_list(kSel-1), (Bpk_list(kSel-1)-B_target)/B_target*100);
end
if kSel < numel(N_turn_list) && ~isnan(Bpk_list(kSel+1))
    fprintf('      alternativa sopra: N_turn = %d -> %.2f T (%+.2f %%)\n', ...
        N_turn_list(kSel+1), Bpk_list(kSel+1), (Bpk_list(kSel+1)-B_target)/B_target*100);
end
if abs(errB) > 5
    warning('step2:targetMiss', ...
        ['Il miglior N_turn manca il target del %.1f %%: nessuna combinazione ' ...
         '(N_layer = %d fissato da H_max, N_turn intero) centra i %.1f T. ' ...
         'Serve agire su cab.width, che cambia sia il passo radiale sia N_layer.'], ...
        errB, N_layer, B_target);
elseif errB > 0
    fprintf('      NOTA: B_pk e'' %+.2f %% SOPRA il target -> consuma margine di\n', errB);
    fprintf('            temperatura. Tracciato fino allo Step 7 (portata minima).\n');
end

r_first = R_in;  r_last = R_in + cab.width*N_turn;
fprintf('  twist: %d giri interi/turn -> l_p = %.3f m (turn interno) ... %.3f m (esterno)\n', ...
    cab.N_twist, 2*pi*r_first/cab.N_twist, 2*pi*r_last/cab.N_twist);
fprintf('  ampere-turn totali = %.3e A\n', I_nom*N_layer*2*N_turn);

% ------------------------------------------------------------------------
%  2.4  MATRICI P ed E DEFINITIVE
% ------------------------------------------------------------------------
fprintf('\n--- 2.4 Matrici P ed E definitive ---\n');

% (a) MACRO: tutti i DP TRANNE il centrale (esclusione per INDICE di layer:
%     rimuove esattamente e solo il conduttore che il micro ricostruisce).
[P_macro, E_macro] = buildMacroStack(N_layer, idxCentral, z_layers, N_turn, ...
                                     disc.N_seg_field, R_in, cab.width);
fprintf('  P_macro: %5d nodi | E_macro: %5d segmenti  (%d layer attivi)\n', ...
    size(P_macro,1), size(E_macro,1), N_layer-1);

% (b) MICRO: il DP centrale, N_slots strand equivalenti twistati (1 per slot)
[P_center, th_center, s_center] = buildDPcenterline(N_turn, disc.N_seg_field, ...
                                                    R_in, cab.width, z_layers(idxCentral));
[Tv, Nv, Bv] = frenetFrame(P_center);
[P_micro, E_micro, strandID, faceInj, faceExt] = buildTwistedStrands( ...
    P_center, th_center, Tv, Nv, Bv, cab);
Nn_c = size(P_center,1);            % nodi della LINEA CENTRALE

fprintf('  P_micro: %5d nodi | E_micro: %5d segmenti  (%d strand)\n', ...
    size(P_micro,1), size(E_micro,1), cab.N_slots);
fprintf('  linea centrale: %d nodi | facce BC: %d iniezione, %d estrazione\n', ...
    Nn_c, numel(faceInj), numel(faceExt));

I_slot = I_nom / cab.N_slots;   % equiripartizione (dimostrata nello Step 4)
fprintf('  corrente per strand (equiripartizione) = %.1f A\n', I_slot);

% ------------------------------------------------------------------------
%  2.5  VISUALIZZAZIONE GEOMETRICA 3D (macro e micro)
% ------------------------------------------------------------------------
fprintf('\n--- 2.5 Visualizzazione geometrica 3D ---\n');

nodesPerLayer  = size(P_macro,1) / (N_layer - 1);
nodesPerStrand = size(P_micro,1) / cab.N_slots;

% ============ FIGURA 1: MACRO - distribuzione dei cavi nel solenoide =====
figure('Name','GEOMETRIA MACRO - solenoide completo','Color','w', ...
       'Position',[60 60 1400 480]);

% (a) vista 3D: tutti i DP impilati, colorati per quota
subplot(1,3,1); hold on; grid on; box on;
if exist('turbo','file'), cmapL = turbo(N_layer); else, cmapL = parula(N_layer); end
kL = 0;
for iz = 1:N_layer
    if iz == idxCentral, continue; end
    kL = kL + 1;
    idx = (kL-1)*nodesPerLayer + (1:nodesPerLayer);
    plot3(P_macro(idx,1), P_macro(idx,2), P_macro(idx,3), ...
        'Color',[cmapL(iz,:) 0.45], 'LineWidth',0.4);
end
plot3(P_center(:,1), P_center(:,2), P_center(:,3), 'r-', 'LineWidth',2.5);
axis equal; view(38,20);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('%d DP impilati (%d turn radiali ciascuno)\nrosso = DP centrale a z = 0', ...
    N_layer, N_turn));

% (b) spaccato: quarto di solenoide, per vedere il build interno
subplot(1,3,2); hold on; grid on; box on;
kL = 0;
for iz = 1:N_layer
    if iz == idxCentral, continue; end
    kL = kL + 1;
    idx = (kL-1)*nodesPerLayer + (1:nodesPerLayer);
    Pl = P_macro(idx,:);
    Pl(~(Pl(:,1)>=0 & Pl(:,2)>=0), :) = NaN;      % taglio a un quarto
    plot3(Pl(:,1), Pl(:,2), Pl(:,3), 'Color',[cmapL(iz,:) 0.7], 'LineWidth',0.6);
end
Pc = P_center;  Pc(~(Pc(:,1)>=0 & Pc(:,2)>=0), :) = NaN;
plot3(Pc(:,1), Pc(:,2), Pc(:,3), 'r-', 'LineWidth',2.5);
axis equal; view(42,24);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Spaccato (quarto): build radiale e assiale');

% (c) proiezione r-z: dove sta davvero il rame
subplot(1,3,3); hold on; grid on; box on;
r_macro = hypot(P_macro(:,1), P_macro(:,2));
plot(r_macro, P_macro(:,3), '.', 'MarkerSize',1, 'Color',[0.25 0.45 0.75]);
plot(hypot(P_center(:,1),P_center(:,2)), P_center(:,3), 'r.', 'MarkerSize',6);
xlabel('r [m]'); ylabel('z [m]'); yline(0,'k--');
title(sprintf('Proiezione r-z\nR_{in} = %.2f m, R_{out} = %.2f m', R_in, r_last));
legend('altri DP','DP centrale','Location','eastoutside'); axis tight;

% ============ FIGURA 2: MICRO - distribuzione degli strand nel cavo ======
figure('Name','GEOMETRIA MICRO - DP centrale twistato','Color','w', ...
       'Position',[80 80 1400 820]);
cmapS = lines(cab.N_slots);

% (a) vista 3D dell'intero DP centrale con i suoi strand
subplot(2,3,[1 2]); hold on; grid on; box on;
for j = 1:cab.N_slots
    idx = (j-1)*nodesPerStrand + (1:nodesPerStrand);
    plot3(P_micro(idx,1), P_micro(idx,2), P_micro(idx,3), ...
        'Color',cmapS(j,:), 'LineWidth',0.7);
end
plot3(P_center(:,1), P_center(:,2), P_center(:,3), 'k--', 'LineWidth',1.0);
axis equal; view(40,26);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('DP centrale: %d strand equivalenti twistati\n(nero tratteggiato = asse del cavo)', ...
    cab.N_slots));

% (b) zoom 3D su un passo di twist: si vede l'elica dei 6 strand
subplot(2,3,3); hold on; grid on; box on;
lp_inner = 2*pi*R_in / cab.N_twist;
mZoom = (s_center >= 0) & (s_center <= 1.05*lp_inner);
for j = 1:cab.N_slots
    idx = (j-1)*nodesPerStrand + find(mZoom);
    plot3(P_micro(idx,1), P_micro(idx,2), P_micro(idx,3), ...
        'Color',cmapS(j,:), 'LineWidth',2.0);
end
idxC = find(mZoom);
plot3(P_center(idxC,1), P_center(idxC,2), P_center(idxC,3), 'k--','LineWidth',1.2);
axis equal; view(44,28);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('Zoom: 1 passo di twist (l_p = %.3f m)', lp_inner));

% (c) sezioni trasversali: la rotazione degli slot lungo il cavo
nCut = 3;
sCut = linspace(0, lp_inner, nCut+1);  sCut(end) = [];
for c = 1:nCut
    subplot(2,3,3+c); hold on; grid on; box on; axis equal;
    [~, kN] = min(abs(s_center - sCut(c)));
    thPlot = linspace(0,2*pi,120);      % angolo di disegno (nome esplicito: 'th'
                                        % collideva con la struct dello Step 7)
    plot(cab.width/2*1000*cos(thPlot), cab.width/2*1000*sin(thPlot), 'k-','LineWidth',1.2);
    plot(cab.R_core*1000*cos(thPlot), cab.R_core*1000*sin(thPlot), 'k:','LineWidth',0.8);
    for j = 1:cab.N_slots
        idx = (j-1)*nodesPerStrand + kN;
        off3 = P_micro(idx,:) - P_center(kN,:);
        u = dot(off3, Nv(kN,:));  v = dot(off3, Bv(kN,:));
        plot(u*1000, v*1000, 'o','MarkerSize',11, ...
            'MarkerFaceColor',cmapS(j,:), 'MarkerEdgeColor','k');
        text(u*1000, v*1000, sprintf('%d',j), 'HorizontalAlignment','center', ...
            'FontSize',8,'Color','w','FontWeight','bold');
    end
    xlim([-1 1]*cab.width/2*1000*1.15); ylim([-1 1]*cab.width/2*1000*1.15);
    xlabel('N [mm]'); ylabel('B [mm]');
    title(sprintf('sezione a s = %.2f m', sCut(c)),'FontSize',9);
end

% (d) verifica visiva della commensurabilita': angolo di twist vs turn
subplot(2,3,6); hold on; grid on; box on;
alpha_check = cab.N_twist*th_center/(2*pi);       % giri di twist compiuti
turn_idx    = th_center/(2*pi);                   % turn percorsi
plot(turn_idx, alpha_check, 'b-','LineWidth',1.4);
plot(0:2*N_turn, cab.N_twist*(0:2*N_turn), 'ro','MarkerSize',5,'MarkerFaceColor','r');
xlabel('turn percorsi'); ylabel('giri di twist compiuti');
title(sprintf('Commensurabilita'': %d twist INTERI per turn', cab.N_twist));
legend('twist(\theta)','multipli interi','Location','northwest');

fprintf('  FIG 1: geometria MACRO (%d DP, vista 3D + spaccato + proiezione r-z)\n', N_layer);
fprintf('  FIG 2: geometria MICRO (%d strand twistati, 3D + zoom + sezioni)\n', cab.N_slots);


%% ========================================================================
%  STEP 3 - CALCOLO DEL CAMPO MAGNETICO
%  ========================================================================
fprintf('\n\n########## STEP 3 - CAMPO MAGNETICO ##########\n');

% ------------------------------------------------------------------------
%  3.1  Decomposizione in tre contributi fisicamente disgiunti
% ------------------------------------------------------------------------
%   B_ext  : dagli altri DP del solenoide (macro, DP centrale escluso)
%   B_mut  : dagli ALTRI strand dello stesso DP (esclude il proprio strand,
%            per non incappare nella singolarita' che segmentMagneticField1d
%            azzera silenziosamente quando Q coincide con un nodo sorgente)
%   B_self : contributo del PROPRIO strand, analitico alla superficie
fprintf('\n--- 3.1 Decomposizione del campo ---\n');

t0 = tic;
[B_tot, B_ext, B_mut, B_self_scalar] = computeDPField( ...
    P_macro, E_macro, I_nom, P_micro, E_micro, strandID, I_slot, cab, ...
    (1:size(P_micro,1))');
tmr = addTimer(tmr, 'Step 3.1  mappa di campo sul DP (macro+mutuo)', toc(t0));

Bmag_ext = vecnorm(B_ext,2,2);
Bmag_mut = vecnorm(B_mut,2,2);
Bmag_tot = vecnorm(B_tot,2,2) + B_self_scalar;

fprintf('  B_ext  (background dagli altri %2d DP) : max = %6.2f T\n', N_layer-1, max(Bmag_ext));
fprintf('  B_mut  (dagli altri %d strand del DP)  : max = %6.2f T\n', cab.N_slots-1, max(Bmag_mut));
fprintf('  B_self (proprio strand, analitico)    :       %6.2f T\n', B_self_scalar);
fprintf('  B_TOT  sul conduttore                 : max = %6.2f T  (target %.1f T)\n', ...
    max(Bmag_tot), B_target);

[Bpk, idxPk] = max(Bmag_tot);
rPk = hypot(P_micro(idxPk,1), P_micro(idxPk,2));
fprintf('  picco in r = %.4f m (R_in = %.3f m), z = %+.4f m\n', rPk, R_in, P_micro(idxPk,3));
if abs(rPk - R_in) > 1.5*cab.width
    warning('step3:peakLoc','Il picco NON e'' sul turn interno: verificare la geometria.');
end

% ------------------------------------------------------------------------
%  3.2  Distribuzione spaziale del campo
% ------------------------------------------------------------------------
fprintf('\n--- 3.2 Distribuzione spaziale ---\n');

s_nodes = repmat(s_center, cab.N_slots, 1);      % ascissa curvilinea per nodo
r_nodes = hypot(P_micro(:,1), P_micro(:,2));     % raggio per nodo

% th_center vive sulla LINEA CENTRALE (Nn_c x 1); va replicato su tutti i nodi
% degli N_slots strand per essere compatibile con gli indici di slot.
theta_nodes = repmat(th_center, cab.N_slots, 1);
slot_nodes  = strandNodeSlot(cab.N_slots, Nn_c);
alpha_nodes = twistAngle(theta_nodes, slot_nodes, cab);

Nrep = repmat(Nv, cab.N_slots, 1);
Brep = repmat(Bv, cab.N_slots, 1);

switch lower(cab.tapeNormal)
    case 'radial',    nFace =  Nrep.*cos(alpha_nodes) + Brep.*sin(alpha_nodes);
    case 'azimuthal', nFace = -Nrep.*sin(alpha_nodes) + Brep.*cos(alpha_nodes);
    otherwise, error('step3:tapeNormal','cab.tapeNormal: ''radial'' o ''azimuthal''.');
end
nFace = nFace ./ vecnorm(nFace,2,2);

B_perp    = abs(sum(B_tot .* nFace, 2));
B_par     = sqrt(max(vecnorm(B_tot,2,2).^2 - B_perp.^2, 0));
theta_deg = atan2d(B_perp, B_par);   % 0 = B nel piano | 90 = B perpendicolare

% Il twist modula localmente le componenti B_perp e B_par, soprattutto dove
% il campo di background del solenoide e' piccolo. B_self e' aggiunto in modulo
% al campo totale usato per T_cs e non entra nella decomposizione B_perp/B_par;
% per questo sqrt(B_perp^2+B_par^2) non coincide necessariamente con Bmag_tot.

fprintf('  orientamento stack: ''%s''\n', cab.tapeNormal);
fprintf('  B_perp: max = %.2f T | B_par: max = %.2f T | theta: %.1f ... %.1f deg\n', ...
    max(B_perp), max(B_par), min(theta_deg), max(theta_deg));

edgesR   = R_in + (0:N_turn)*cab.width;
Bprofile = nan(N_turn,1);  rProfile = nan(N_turn,1);
for k = 1:N_turn
    m = (r_nodes >= edgesR(k)) & (r_nodes < edgesR(k+1));
    if any(m), Bprofile(k) = max(Bmag_tot(m)); rProfile(k) = mean(r_nodes(m)); end
end
fprintf('  profilo radiale (B_max per turn):\n');
for k = 1:N_turn
    if ~isnan(Bprofile(k))
        fprintf('    turn %2d (r = %.3f m): %5.2f T\n', k, rProfile(k), Bprofile(k));
    end
end

% ------------------------------------------------------------------------
%  3.3  TEMPERATURA DI CURRENT SHARING E MARGINE
% ------------------------------------------------------------------------
%  T_cs e' definita implicitamente da:   Ic(T_cs, B) = I_tape
%  ed e' risolta per BISEZIONE da data.Ic.Tcs (Step 1). Non dipende da T_op:
%  T_op entra solo nel confronto finale  margine = T_cs - T_op >= 5 K.
%
%  DUE SCALE DI CORRENTE, da non confondere:
%    - I_slot  = I_nom/N_slots            -> corrente dell'INTERO stack (modello EM)
%    - I_tape  = I_slot/N_tapes_slot      -> corrente del SINGOLO nastro (per T_cs)
%  Chi va in current sharing per primo e' il singolo nastro, quindi T_cs va
%  calcolata con I_tape e con l'Ic del singolo nastro (NON moltiplicato).
%
%  CORREZIONE NEAR-EDGE: B_tot e' valutato al CENTROIDE dello stack (R_core),
%  ma il nastro che transisce per primo e' quello al bordo INTERNO, dove il
%  campo VICINO (B_self + B_mut) e' piu' intenso. Usare il centroide
%  sottostima il campo -> sovrastima Ic -> sovrastima T_cs e il margine.
%  Il campo LONTANO (B_ext) e' invariato: varia su scala di metri, non di mm.
fprintf('\n--- 3.3 Temperatura di current sharing e margine ---\n');

% depth_stack, r_near, cab.r0 e sep_strand sono gia' stati calcolati in 2.1b.
% Correzione near-edge: si applica SOLO ai termini prodotti da sorgenti esterne
% allo strand, perche' un nastro al bordo interno dello stack e' piu' vicino a
% quelle sorgenti. B_self NON va corretto: e' gia' valutato alla superficie del
% conduttore (raggio r0), che e' per definizione il massimo interno allo stack.
corr_mut  = cab.R_core / r_near;   % proxy conservativo per B_mut

I_tape_op = I_slot / cab.N_tapes_slot;                    % corrente per NASTRO
Bmag_near = Bmag_ext + corr_mut*Bmag_mut + B_self_scalar;

% T_cs su TUTTO il profilo (pronta per T_op(x) del modello termico, Step 7)
% T_cs con il FLAG di certificazione: sopra 8 T il Data Book misura solo 4.2 K
% e 20 K, quindi in molti punti del profilo la vera T_cs sta oltre i dati e il
% solver restituisce il bordo come LIMITE INFERIORE. I due casi vanno tenuti
% distinti: un margine calcolato da un limite inferiore e' esso stesso un
% limite inferiore, non una misura.
Tcs_profile = zeros(size(Bmag_near));
Tcs_cert    = false(size(Bmag_near));
t0 = tic;
for iq = 1:numel(Bmag_near)
    [Tcs_profile(iq), st] = data.Ic.Tcs(I_tape_op, Bmag_near(iq));
    Tcs_cert(iq) = strcmp(st,'certified');
end
tmr = addTimer(tmr, 'Step 3.3  T_cs per bisezione su ogni nodo', toc(t0));
nCert = nnz(Tcs_cert);
fprintf('  T_cs: %d nodi su %d CERTIFICATI (radice dentro i dati), %d limiti inferiori\n', ...
    nCert, numel(Tcs_cert), numel(Tcs_cert)-nCert);
[TrLo, TrHi] = data.Ic.Trange(Bpk);
fprintf('  a B_pk = %.2f T il Data Book misura solo T = %.1f ... %.1f K\n', Bpk, TrLo, TrHi);
% Il punto critico va cercato tra i nodi CERTIFICATI (T_cs realmente risolta,
% non un limite inferiore). Questa ricerca NON dipende da T_op: il minimo di
% Tcs_profile tra i nodi certificati e' lo stesso qualunque sia T_op, perche'
% T_op sposta solo un OFFSET costante sul margine.
Tcs_search = Tcs_profile;
Tcs_search(~Tcs_cert) = Inf;
if all(isinf(Tcs_search))
    warning('step3:noCertified', ...
        ['NESSUN nodo ha T_cs certificata: impossibile derivare T_op. Ridurre ' ...
         'cab.N_tapes_slot (piu'' corrente per nastro abbassa T_cs dentro la ' ...
         'finestra misurata dal Data Book).']);
    [~, idxCrit] = min(Tcs_profile);
else
    [~, idxCrit] = min(Tcs_search);
end
rCrit = hypot(P_micro(idxCrit,1), P_micro(idxCrit,2));
Bc    = Bmag_near(idxCrit);

fprintf('  nastri/slot = %d -> stack %.2f mm, r_near = %.2f mm (corr B_mut +%.1f%%)\n', ...
    cab.N_tapes_slot, depth_stack*1000, r_near*1000, (corr_mut-1)*100);
fprintf('  corrente per nastro I_tape = %.1f A  (I_slot = %.0f A / %d nastri)\n', ...
    I_tape_op, I_slot, cab.N_tapes_slot);
fprintf('  punto CRITICO (Tcs minima certificata): r = %.4f m, z = %+.4f m, B_near = %.2f T\n', ...
    rCrit, P_micro(idxCrit,3), Bc);

% ------------------------------------------------------------------------
%  T_op DERIVATO: il MASSIMO valore che rispetta il margine di 5 K richiesto
%  dal brief anche nello scenario PEGGIORE di dispersione dei lift factor
%  (-20%, dichiarata per iscritto dal fornitore).
%     T_op_guess = Tcs_caso_peggiore(punto critico) - dT_req
%  Qualunque T_op SOPRA questo valore violerebbe il margine nello scenario
%  peggiore; qualunque valore SOTTO e' conservativo ma spreca refrigerazione.
% ------------------------------------------------------------------------
if ~isfield(data,'lfx')
    error('step3:noSpreadData', ...
        ['data.lfx.spread non disponibile: impossibile derivare T_op senza la ' ...
         'dispersione dei lift factor dichiarata dal fornitore.']);
end
spread_worst = data.lfx.spread.val(2);
[Tcs_worst, st_worst] = data.Ic.Tcs(I_tape_op/(1-spread_worst), Bc);

if strcmp(st_worst,'coldCap')
    error('step3:designInfeasible', ...
        ['Anche alla T piu'' fredda misurata e nello scenario -%.0f%%, il nastro ' ...
         'sarebbe gia'' in current sharing: NESSUN T_op soddisfa il margine. ' ...
         'Aumentare cab.N_tapes_slot.'], spread_worst*100);
end
if strcmp(st_worst,'lowerBound')
    warning('step3:worstCaseUncertified', ...
        ['La T_cs nello scenario -%.0f%% eccede la T piu'' calda misurata a questo ' ...
         'campo: il valore usato per derivare T_op (%.2f K) e'' un LIMITE INFERIORE, ' ...
         'quindi T_op derivato e'' conservativo per difetto.'], spread_worst*100, Tcs_worst);
end

T_op_guess = Tcs_worst - dT_req;
T_op_was_overridden = ~isempty(T_op_override);
if T_op_was_overridden
    fprintf('\n  [LOOP B] T_op IMPOSTA dall''esterno a %.2f K (derivata sarebbe %.2f K)\n', ...
        T_op_override, T_op_guess);
    fprintf('  Le perdite di isteresi useranno Ic a questa temperatura: essendo\n');
    fprintf('  piu'' alta, Ic e'' piu'' bassa e la magnetizzazione SCENDE.\n');
    T_op_guess = T_op_override;
end

%  Se T_op e' imposto tramite T_op_override, il codice lo usa come valore
%  operativo assunto a monte; se non e' imposto, T_op viene derivato dal
%  margine richiesto nello scenario peggiore.
if T_op_was_overridden
    fprintf('\n  --- T_op IMPOSTO (Loop B), margine nello scenario peggiore (-%.0f%%) ---\n', ...
        spread_worst*100);
    fprintf('  T_cs nello scenario peggiore = %.2f K [%s]\n', Tcs_worst, st_worst);
    fprintf('  --> T_op_guess = %.2f K (IMPOSTO: NON e'' detto che rispetti %.0f K di margine)\n', ...
        T_op_guess, dT_req);
else
    fprintf('\n  --- T_op derivato dal caso peggiore (-%.0f%%) ---\n', spread_worst*100);
    fprintf('  T_cs nello scenario peggiore = %.2f K [%s]\n', Tcs_worst, st_worst);
    fprintf('  --> T_op_guess = %.2f K (MASSIMO che rispetta %.0f K di margine anche a -%.0f%%)\n', ...
        T_op_guess, dT_req, spread_worst*100);
end

TinMin = data.op.Tin_range.val(1);
if T_op_guess < TinMin
    warning('step3:TopBelowInlet', ...
        ['T_op derivato (%.2f K) e'' SOTTO la T di ingresso elio minima del brief ' ...
         '(%.1f K): margine non raggiungibile nemmeno col refrigerante piu'' freddo. ' ...
         'Aumentare cab.N_tapes_slot.'], T_op_guess, TinMin);
end

% --- margine con il T_op cosi' determinato ---------------------------------
margin_profile = Tcs_profile - T_op_guess;
margin_search  = margin_profile;
margin_search(~Tcs_cert) = Inf;
if all(isinf(margin_search))
    [margin_min, idxCrit2] = min(margin_profile);
else
    [margin_min, idxCrit2] = min(margin_search);
end
if idxCrit2 ~= idxCrit
    warning('step3:idxMismatch', ...
        'Nodo di margine minimo (%d) diverso dal nodo di Tcs minima (%d).', idxCrit2, idxCrit);
end

fprintf('  Ic del nastro a (T_op=%.2f K, B_pk) = %.1f A -> I_tape/Ic = %.3f\n', ...
    T_op_guess, data.Ic.eval(T_op_guess, Bpk), I_tape_op/data.Ic.eval(T_op_guess, Bpk));
fprintf('  T_cs = %.2f K | T_op = %.2f K -> MARGINE MINIMO = %.2f K (richiesto >= %.0f K)\n', ...
    Tcs_profile(idxCrit), T_op_guess, margin_min, dT_req);
Tcs_centroid = data.Ic.Tcs(I_tape_op, Bpk);
fprintf('  [confronto: al centroide T_cs = %.2f K, margine %.2f K -> ottimistico]\n', ...
    Tcs_centroid, Tcs_centroid - T_op_guess);

fprintf('\n  --- Robustezza alla dispersione dei lift factor ---\n');
for sc = [1, 1-data.lfx.spread.val(1), 1-spread_worst]
    if abs(sc-(1-spread_worst)) < 1e-12
        Tcs_sc = Tcs_worst;
    else
        Tcs_sc = data.Ic.Tcs(I_tape_op/sc, Bc);
    end
    fprintf('    LF x %.2f -> T_cs = %5.2f K, margine = %5.2f K  [%s]\n', ...
        sc, Tcs_sc, Tcs_sc-T_op_guess, ternary(Tcs_sc-T_op_guess>=dT_req-1e-9,'OK','SOTTO SOGLIA'));
end

if margin_min < dT_req
    warning('step3:margin', ...
        ['Margine %.2f K < %.0f K. Aumentare cab.N_tapes_slot o abbassare T_op.'], ...
        margin_min, dT_req);
else
    if Tcs_cert(idxCrit)
        fprintf('  --> MARGINE RISPETTATO e CERTIFICATO dai dati misurati\n');
    else
        fprintf('  --> margine rispettato ma da LIMITE INFERIORE (non certificato)\n');
    end
end

% ------------------------------------------------------------------------
%  3.3b  PROFILO DELLA CORRENTE CRITICA LUNGO IL DP CENTRALE
% ------------------------------------------------------------------------
%  Ic non e' un numero, e' un profilo: lungo il double-pancake il campo passa
%  da poco piu' di 1 T sul turn esterno a oltre 15 T su quello interno, e Ic
%  cambia di conseguenza di quasi un ordine di grandezza. Fin qui il codice
%  usava Ic solo dentro altri calcoli (T_cs in 3.3, la resistenza del DAE in
%  4.3, l'isteresi in 6.3) senza mai guardarla per quello che e'.
%
%  Questa sezione la tira fuori. Serve a tre cose concrete:
%    - vedere QUALE giro della spirale detta il limite, che non e'
%      necessariamente quello a campo massimo, perche' anche T_cs e la
%      corrente per nastro entrano nel bilancio;
%    - leggere la frazione di Ic effettivamente usata, I_tape/Ic, che e' il
%      parametro con cui si giudica se il conduttore e' sfruttato bene o
%      sovradimensionato, ed e' anche quello che governa la validita' della
%      power law nello Step 4;
%    - avere il numero da confrontare fra le due varianti del progetto
%      (tabella dei lift factor contro correlazione K_M2).
%
%  Ic e' valutata a T_op_guess, cioe' alla temperatura di esercizio derivata
%  in 3.3, e con il campo NEAR-EDGE: e' lo stesso campo che alimenta T_cs, per
%  non avere due profili di Ic incoerenti nello stesso programma.
fprintf('\n--- 3.3b Profilo della corrente critica lungo il DP ---\n');
t0 = tic;
Ic_prof = arrayfun(@(b) data.Ic.eval(T_op_guess, b), Bmag_near);   % [A] per NASTRO
tmr = addTimer(tmr, 'Step 3.3b profilo di Ic sul DP', toc(t0));

load_prof = I_tape_op ./ Ic_prof;              % frazione di Ic usata
[Ic_min, iIcMin] = min(Ic_prof);
[Ic_max, iIcMax] = max(Ic_prof);
fprintf('  Ic per NASTRO a T_op = %.2f K: %.1f ... %.1f A (escursione %.1fx)\n', ...
    T_op_guess, Ic_min, Ic_max, Ic_max/Ic_min);
fprintf('    minimo a s = %6.2f m, r = %.4f m, B_near = %5.2f T\n', ...
    s_nodes(iIcMin), r_nodes(iIcMin), Bmag_near(iIcMin));
fprintf('    massimo a s = %6.2f m, r = %.4f m, B_near = %5.2f T\n', ...
    s_nodes(iIcMax), r_nodes(iIcMax), Bmag_near(iIcMax));
fprintf('  Ic per STRAND (x%d nastri): %.0f ... %.0f A   contro I_slot = %.0f A\n', ...
    cab.N_tapes_slot, Ic_min*cab.N_tapes_slot, Ic_max*cab.N_tapes_slot, I_slot);
fprintf('  frazione di Ic usata, I_tape/Ic: %.3f ... %.3f\n', ...
    min(load_prof), max(load_prof));
%  Se in qualche nodo la corrente di esercizio superasse Ic, il conduttore
%  sarebbe gia' oltre il criterio a T_op e tutto il resto (margine, perdite,
%  power law) perderebbe senso. Meglio saperlo qui che allo Step 7.
if max(load_prof) >= 1
    warning('step3:overIc', ...
        ['In %d nodi su %d la corrente di esercizio supera Ic a T_op: il ' ...
         'conduttore e'' gia'' oltre il criterio. Aumentare cab.N_tapes_slot.'], ...
        nnz(load_prof >= 1), numel(load_prof));
elseif max(load_prof) > 0.9
    fprintf('    NOTA: si arriva al %.0f%% di Ic. Sopra il 90%% la power law e''\n', ...
        max(load_prof)*100);
    fprintf('    ripida e la resistenza dell''arco diventa sensibile a n (vedi 4.3b).\n');
end

%  PROFILO PER TURN. E' la lettura utile al progetto: dice quale giro della
%  spirale e' il piu' sollecitato. edgesR e rProfile sono gia' stati costruiti
%  in 3.2 per il profilo di campo, quindi le due tabelle sono confrontabili
%  riga per riga.
IcTurn = nan(N_turn,1);  loadTurn = nan(N_turn,1);  TcsTurn = nan(N_turn,1);
for k = 1:N_turn
    m = (r_nodes >= edgesR(k)) & (r_nodes < edgesR(k+1));
    if any(m)
        IcTurn(k)   = min(Ic_prof(m));
        loadTurn(k) = max(load_prof(m));
        TcsTurn(k)  = min(Tcs_profile(m));
    end
end
fprintf('  profilo per turn (valori piu'' sfavorevoli del giro):\n');
fprintf('    %4s %9s %9s %12s %11s %10s\n', ...
    'turn','r [m]','B_max [T]','Ic [A/nastro]','I_tape/Ic','T_cs [K]');
for k = 1:N_turn
    if ~isnan(IcTurn(k))
        fprintf('    %4d %9.3f %9.2f %12.1f %11.3f %10.2f\n', ...
            k, rProfile(k), Bprofile(k), IcTurn(k), loadTurn(k), TcsTurn(k));
    end
end
fprintf('  il nodo critico del margine (nodo %d, s = %.2f m) ha Ic = %.1f A/nastro\n', ...
    idxCrit, s_nodes(idxCrit), Ic_prof(idxCrit));
%  Il turn a Ic minima e quello a margine minimo NON devono per forza
%  coincidere: Ic minima segue il campo, il margine segue T_cs - T_op. Se
%  coincidono e' un caso, e vale la pena vederlo stampato.
[~, kIcW]  = min(IcTurn);
[~, kMarW] = min(TcsTurn);
if kIcW == kMarW
    fprintf('    turn con Ic minima e turn con T_cs minima coincidono (turn %d)\n', kIcW);
else
    fprintf('    turn con Ic minima: %d | turn con T_cs minima: %d (non coincidono:\n', ...
        kIcW, kMarW);
    fprintf('    Ic segue il campo, il margine segue T_cs meno T_op)\n');
end

% --- figura dedicata ------------------------------------------------------
figure('Name','STEP 3.3b - profilo della corrente critica','Color','w', ...
       'Position',[110 110 1320 430]);

subplot(1,3,1); hold on; grid on; box on;
plot(s_nodes, Ic_prof, '.', 'MarkerSize',3, 'Color',[0 0.45 0.85]);
plot(s_nodes(idxCrit), Ic_prof(idxCrit), 'kp','MarkerSize',13,'MarkerFaceColor','y');
yline(I_tape_op,'r--','I_{tape} di esercizio','LineWidth',1.2);
xlabel('s [m]'); ylabel('I_c per nastro [A]');
title(sprintf('I_c lungo il conduttore (T_{op} = %.2f K)', T_op_guess));
legend({'I_c(s)','nodo critico del margine','corrente di esercizio'}, ...
    'Location','best','FontSize',7);

subplot(1,3,2); hold on; grid on; box on;
yyaxis left;  plot(rProfile, IcTurn,'s-','LineWidth',1.4);  ylabel('I_c minima del turn [A]');
yyaxis right; plot(rProfile, Bprofile,'o-','LineWidth',1.2); ylabel('B_{max} del turn [T]');
xlabel('r [m]');
title('Profilo radiale: I_c contro campo');

subplot(1,3,3); hold on; grid on; box on;
plot(s_nodes, load_prof, '.', 'MarkerSize',3, 'Color',[0.85 0.33 0.10]);
yline(1,'r--','I_{tape} = I_c','LineWidth',1.2);
yline(max(load_prof),'k:',sprintf('max %.3f', max(load_prof)));
xlabel('s [m]'); ylabel('I_{tape} / I_c  [-]');
title('Frazione di I_c utilizzata');
ylim([0, max(1.05, 1.1*max(load_prof))]);

% ------------------------------------------------------------------------
%  3.4  Grafici del campo e del margine
% ------------------------------------------------------------------------
figure('Name','STEP 3 - campo e margine','Color','w','Position',[100 100 1300 760]);

subplot(2,3,1);
plot(N_turn_list(valid), Bpk_list(valid),'o-','LineWidth',1.4); hold on;
yline(B_target,'r--','15 T'); xline(N_turn,'k:','scelto');
xlabel('N_{turn} radiali'); ylabel('B_{pk} [T]');
title('Dimensionamento su 15 T'); grid on;

subplot(2,3,2);
plot(rProfile, Bprofile,'s-','LineWidth',1.4); hold on; yline(B_target,'r--');
xlabel('r [m]'); ylabel('B_{max} per turn [T]'); title('Profilo radiale'); grid on;

subplot(2,3,3);
plot(s_nodes, Bmag_tot,'.','MarkerSize',3);
xlabel('s [m]'); ylabel('|B| [T]'); title('B lungo il conduttore'); grid on;

subplot(2,3,4);
plot(s_nodes, B_perp,'.','MarkerSize',3); hold on; plot(s_nodes, B_par,'.','MarkerSize',3);
xlabel('s [m]'); ylabel('[T]'); legend('B_\perp','B_\parallel','Location','best');
title('Componenti nella terna del nastro'); grid on;

subplot(2,3,5);
% I nodi CERTIFICATI e quelli a LIMITE INFERIORE vanno distinti: il tratto
% piatto al bordo dei dati NON e' un valore calcolato, e' il limite della
% finestra misurata a quel campo. Confonderli rende il grafico ingannevole.
plot(s_nodes(Tcs_cert),  Tcs_profile(Tcs_cert), '.','MarkerSize',4,'Color',[0 0.45 0.85]); hold on;
plot(s_nodes(~Tcs_cert), Tcs_profile(~Tcs_cert),'.','MarkerSize',4,'Color',[0.85 0.5 0]);
yline(T_op_guess,'r--','T_{op}'); yline(T_op_guess+dT_req,'g--','T_{op}+5K');
plot(s_nodes(idxCrit), Tcs_profile(idxCrit),'kp','MarkerSize',13,'MarkerFaceColor','y');
xlabel('s [m]'); ylabel('T_{cs} [K]');
title('T_{cs} lungo il conduttore'); grid on;
legend({'certificata','limite inferiore','T_{op}','T_{op}+5K','punto critico'}, ...
    'Location','best','FontSize',7);

subplot(2,3,6);
plot(s_nodes(Tcs_cert),  margin_profile(Tcs_cert), '.','MarkerSize',4,'Color',[0 0.45 0.85]); hold on;
plot(s_nodes(~Tcs_cert), margin_profile(~Tcs_cert),'.','MarkerSize',4,'Color',[0.85 0.5 0]);
yline(dT_req,'r--','5 K richiesti');
plot(s_nodes(idxCrit), margin_min,'kp','MarkerSize',13,'MarkerFaceColor','y');
xlabel('s [m]'); ylabel('\DeltaT margine [K]');
title(sprintf('Margine (min certificato = %.2f K)', margin_min)); grid on;
legend({'certificato','limite inferiore','soglia','punto critico'}, ...
    'Location','best','FontSize',7);


%% ========================================================================
%  STEP 4 - SISTEMA MATRICIALE DAE
%  ========================================================================
%  Formulazione [Freschi s28-s33]:
%     [L 0; 0 0] d/dt [i;phi] + [R A; -A' G] [i;phi] = [0; i_s]
%  M singolare (DAE), R = R(i,B,T) non lineare dalla power law E-J.
%  BC [Freschi s34-s37]: faccia di iniezione a potenziale fisso (Dirichlet),
%  faccia di estrazione equipotenziale flottante -> directSolve('equivalue').
%
%  SCOPO: dimostrare che le correnti si ripartiscono UGUALMENTE tra gli strand.
%  Non e' una tautologia: la ripartizione uniforme discende dal fatto che, con
%  twist commensurato, ogni strand attraversa le stesse condizioni locali (B, e
%  quindi Ic e R) in ordine diverso ma con lo stesso INTEGRALE lungo il percorso.
%  Il DAE lo verifica risolvendo il circuito senza imporre nulla a priori.
fprintf('\n\n########## STEP 4 - SISTEMA DAE ##########\n');

% ------------------------------------------------------------------------
%  4.1  Mesh ridotta per L e DAE
% ------------------------------------------------------------------------
[P_c2, th_c2, ~] = buildDPcenterline(N_turn, disc.N_seg_dae, R_in, cab.width, ...
                                     z_layers(idxCentral));
[T2, N2, B2] = frenetFrame(P_c2);
[P_dae, E_dae, sid_dae, faceInj_dae, faceExt_dae] = buildTwistedStrands( ...
    P_c2, th_c2, T2, N2, B2, cab);

Ne = size(E_dae,1);  Nn = size(P_dae,1);
fprintf('\n--- 4.1 Mesh DAE ---\n');
fprintf('  %d nodi, %d archi (mesh campo: %d archi)\n', Nn, Ne, size(E_micro,1));

% ------------------------------------------------------------------------
%  4.2  Matrice delle induttanze L
% ------------------------------------------------------------------------
%  r0 (cab.r0) e' stato calcolato in 2.1b dalla sezione REALE dello stack, ed e'
%  lo stesso raggio usato per il self-field in 3.1: le due cose devono essere
%  coerenti, perche' descrivono lo stesso conduttore equivalente.
r0 = cab.r0;
fprintf('\n--- 4.2 Matrice induttanze L ---\n');
fprintf('  r0 = %.3f mm (sezione stack %.2f mm^2) | sep strand = %.3f mm\n', ...
    r0*1000, A_stack*1e6, sep_strand*1000);
fprintf('  calcolo L (Neumann, integral2)... puo'' richiedere minuti\n');
tic;
L = inductanceMatrix1d(P_dae, E_dae, r0);
L = (L + L.')/2;                   % simmetrizza (drift floating-point)
tL = toc;
tmr = addTimer(tmr, 'Step 4.2  matrice induttanze L (Neumann)', tL);
infoL = whos('L');
fprintf('  L: %dx%d, %.1f MB, %.1f s\n', size(L,1), size(L,2), infoL.bytes/1024^2, tL);

% --- diagnostica: EQUIVALENZA DEGLI STRAND (la causa fisica dell'equiripartizione)
% Il flusso concatenato da uno strand e' governato dalla somma della sua riga di
% L. Se gli strand sono equivalenti (twist ben risolto), le somme aggregate per
% strand coincidono e la corrente si divide in parti uguali. La dispersione fra
% queste somme e' una PREVISIONE quantitativa dello scarto atteso nel DAE.
rowSumL = sum(L, 2);
Lstrand = accumarray(sid_dae, rowSumL, [cab.N_slots 1], @sum);
fprintf('  flusso concatenato per strand [uWb/A]: ');
fprintf('%.3f ', Lstrand*1e6); fprintf('\n');
sprdL = (max(Lstrand)-min(Lstrand))/mean(Lstrand)*100;
fprintf('  dispersione induttiva tra strand: %.3f %%', sprdL);
if sprdL < 1
    fprintf('  -> strand equivalenti: attesa equiripartizione\n');
else
    fprintf('  -> strand NON equivalenti\n');
    fprintf('     Lo scarto del DAE sara'' di quest''ordine. Causa tipica: twist\n');
    fprintf('     aliasato nella mesh (vedi 2.1c).\n');
end

% --- diagnostica: L deve essere definita positiva e ben condizionata --------
% Se gli strand fossero geometricamente indistinguibili, L avrebbe autovalori
% nulli/negativi e il modo circolante non sarebbe vincolato -> correnti assurde.
evL = eig(full(L(1:min(200,size(L,1)), 1:min(200,size(L,1)))));  % campione
fprintf('  autovalori di L (campione): min = %.3e, max = %.3e\n', min(real(evL)), max(real(evL)));
if min(real(evL)) <= 0
    warning('step4:Lsingular', ...
        ['L non definita positiva (autovalore minimo %.2e): il modo di corrente ' ...
         'circolante non e'' vincolato e le correnti del DAE non sono affidabili. ' ...
         'Verificare r0 e la separazione tra strand.'], min(real(evL)));
end

% ------------------------------------------------------------------------
%  4.3  Incidenza, conduttanze trasversali, resistenze non lineari
% ------------------------------------------------------------------------
fprintf('\n--- 4.3 Assemblaggio del sistema ---\n');
A    = buildIncidence(E_dae, Nn);      % Ne x Nn: A(e,start)=-1, A(e,end)=+1
Gmat = sparse(Nn, Nn);                 % G = 0 (scelta di scope: niente conduttanze
                                       % trasversali, tutta la corrente nel SC)

segLen = vecnorm(P_dae(E_dae(:,2),:) - P_dae(E_dae(:,1),:), 2, 2);
Bedge  = interpFieldOnEdges(P_dae, E_dae, P_micro, Bmag_tot);   % campo medio per arco

E0  = data.basis.E0.val;

% ========================================================================
%  ESPONENTE n DELLA POWER LAW - tabella n(T,B) dello Step 1
% ========================================================================
%  n e' valutato arco per arco con data.nEJ.eval(T,B), sulla stessa condizione
%  locale usata per la legge E-J. Non si usa n_lab_Ltap.
nEJ = arrayfun(@(b) data.nEJ.eval(T_op_guess, b), Bedge);

fprintf('\n--- 4.3b Esponente n della power law ---\n');
fprintf('  n(T,B) dalla tabella del dataset del fornitore, valutato arco per arco\n');
fprintf('  a T_op = %.2f K: n = %.1f ... %.1f (media %.1f)\n', ...
    T_op_guess, min(nEJ), max(nEJ), mean(nEJ));
fprintf('  valore sintetico al punto di progetto (data.basis.n) = %g\n', data.basis.n.val);
if isfield(data.basis,'n_lab77')
    fprintf('  controllo di laboratorio: n = %.1f a 76.9 K campo proprio (vedi openItems)\n', ...
        data.basis.n_lab77.val);
end
if max(Bedge) > max(data.nEJ.Bgrid)
    fprintf('  NOTA: B arriva a %.1f T ma la tabella si ferma a %.0f T -> oltre il bordo\n', ...
        max(Bedge), max(data.nEJ.Bgrid));
    fprintf('        l''evaluatore CLAMPA. Alle basse T n si appiattisce ad alto campo,\n');
    fprintf('        quindi il clamp e'' stabile, ma resta un''estrapolazione dichiarata.\n');
end

% Ic dello STRAND = N_tapes * Ic del singolo nastro. Usare l'Ic del singolo
% nastro con la corrente dell'intero strand darebbe I/Ic ~ 60 e, elevato a
% n-1, farebbe esplodere R numericamente.
IcTape = arrayfun(@(b) data.Ic.eval(T_op_guess, b), Bedge);
IcEdge = cab.N_tapes_slot * IcTape;
fprintf('  E0 = %.1e V/m | n = %.1f-%.1f (tabella) | T_op = %.1f K\n', ...
    E0, min(nEJ), max(nEJ), T_op_guess);
fprintf('  Ic strand: min = %.0f A, max = %.0f A | I_slot/Ic_min = %.3f\n', ...
    min(IcEdge), max(IcEdge), I_slot/min(IcEdge));

%  VERIFICA CHE LA SCELTA DI n SIA IRRILEVANTE. Si calcola la tensione
%  resistiva sul DP agli estremi della banda plausibile di n. Se la differenza
%  resta molto sotto qualunque soglia di interesse, l'incertezza su n non
%  influenza il progetto e non serve raffinarla.
rIIc = I_slot/min(IcEdge);
L_dp_res = sum(segLen(sid_dae == 1));
fprintf('  sensibilita'' a n (I/Ic = %.3f, lunghezza strand %.1f m):\n', rIIc, L_dp_res);
for nTest = [data.basis.n_sweep.val(1), mean(nEJ), data.basis.n_sweep.val(2)]
    fprintf('    n = %5.1f -> tensione resistiva sul DP = %.2e V\n', ...
        nTest, E0*rIIc^nTest*L_dp_res);
end
if E0*rIIc^data.basis.n_sweep.val(1)*L_dp_res < 1e-6
    fprintf('    -> sotto il microvolt anche col n piu'' basso: la scelta di n NON\n');
    fprintf('       influenza il progetto. La ripartizione e'' governata da L.\n');
else
    warning('step4:nMatters', ...
        ['La tensione resistiva non e'' trascurabile (%.2e V) al n piu'' basso ' ...
         'della banda: qui la scelta di n conta e l''estrapolazione a 15 T va ' ...
         'sostituita con una misura a campo pieno.'], ...
         E0*rIIc^data.basis.n_sweep.val(1)*L_dp_res);
end

% --- diagnostica: perche' ci aspettiamo l'equiripartizione -----------------
% Ogni strand attraversa le stesse condizioni locali in ordine diverso: se il
% twist e' commensurato, l'integrale di R lungo ciascuno strand e' identico.
Rlin_edge = (E0 .* segLen ./ IcEdge);        % resistenza a I = Ic (riferimento)
Rsum = accumarray(sid_dae, Rlin_edge, [cab.N_slots 1], @sum);
fprintf('\n  Somma di R lungo ciascuno strand [nOhm]: ');
fprintf('%.4f ', Rsum*1e9); fprintf('\n');
sprd = (max(Rsum)-min(Rsum))/mean(Rsum)*100;
fprintf('  dispersione tra strand: %.4f %%  ', sprd);
if sprd < 1
    fprintf('-> gli strand sono ELETTRICAMENTE EQUIVALENTI\n');
    fprintf('     (e'' questa la ragione fisica dell''equiripartizione)\n');
else
    fprintf('-> ATTENZIONE: strand non equivalenti, attendersi ripartizione disuniforme\n');
end

% ------------------------------------------------------------------------
%  4.4  Integrazione implicita sulla rampa
% ------------------------------------------------------------------------
%  CONVENZIONE DEI SEGNI (verificata): la seconda riga a blocchi e'
%      -A'*i + G*phi = i_s
%  e (-A'*i)_nodo = corrente netta USCENTE dal nodo attraverso gli archi.
%  Quindi i_s > 0 = corrente INIETTATA nel nodo dall'esterno.
%  Qui la faccia di iniezione e' vincolata a phi = 0 (Dirichlet: la sua
%  equazione di KCL viene sostituita, la corrente vi entra implicitamente) e
%  sulla faccia di estrazione si impone i_s < 0 (corrente ESTRATTA).
%  Con i_s > 0 sulla faccia di estrazione le correnti d'arco uscirebbero
%  NEGATIVE (verso opposto agli archi), facendo fallire la verifica.
%  t_ramp e' gia' definito nello Step 1 come I_nom/dIdt: la durata della rampa
%  e' un dato dello scenario, non un parametro dell'integrazione. Qui lo Step 4
%  integra ESATTAMENTE sulla rampa, quindi i due coincidono - ma per scelta,
%  non per definizione (se un giorno si aggiungesse un flat-top al DAE, t_end
%  crescerebbe e t_ramp dovrebbe restare 250 s).
t_end  = t_ramp;                       % 250 s: il DAE copre la sola rampa
%  PASSO TEMPORALE DEL DAE ELETTROMAGNETICO.
%  Si mantiene dt = 10 s, cioe' 25 passi sulla rampa 0 -> 50 kA.
%  Questo valore riguarda SOLO lo Step 4.4 e non viene modificato dallo
%  studio di convergenza termoidraulico dello Step 7.
dt     = 10;                           % [s]
nSteps = round(t_end/dt);
x = zeros(Ne+Nn,1);                    % [i; phi]  (parte da zero: nessun bias)
M = blkdiag(sparse(L), sparse(Nn,Nn));

fprintf('\n--- 4.4 Integrazione temporale ---\n');
fprintf('  rampa 0 -> %.0f kA in %.0f s, dt = %.1f s (%d passi)\n', ...
    I_nom/1e3, t_end, dt, nSteps);

%  CONDIZIONAMENTO E RAFFINAMENTO ITERATIVO
%  Il sistema e' un problema a PUNTO SELLA (blocco (2,2) nullo perche' G = 0) e
%  mescola scale molto diverse: L/dt ~ 1e-8, A ~ 1, R dalla power law che
%  attraversa decadi durante la rampa. Il numero di condizione stimato risulta
%  O(1e17), al limite della doppia precisione -> viene monitorato con cond1est
%  (terzo output di directSolve) e il residuo ||Ax-b||/||b|| e' riportato in 4.5.
%
%  directSolve puo' applicare iterativeRefinement passando 'iterMAX' > 0.
%  ATTENZIONE: iterativeRefinement rispetta il proprio flag 'verbose', ma al suo
%  interno richiama directSolve SENZA propagarlo (riga 50 del file fornito),
%  quindi ogni passo di raffinamento stampa comunque. Con 25 passi temporali x
%  fino a 20 iterazioni di Picard sarebbero centinaia di righe: il raffinamento
%  e' quindi DISATTIVATO di default. Abilitarlo (nIterRef > 0) solo se il
%  residuo riportato in 4.5 risulta elevato, accettando le stampe.
nIterRef = 0;      % passi di raffinamento iterativo (0 = disabilitato, vedi sopra)

hist_t = zeros(nSteps,1);
hist_i = zeros(nSteps, cab.N_slots);
condEst = NaN;
%  AVANZAMENTO. Con dt = 10 s i passi sono 25 e ognuno risolve un sistema che
%  contiene il blocco DENSO L. Si stampa una riga circa ogni dieci per cento,
%  con la stima del tempo residuo estrapolata dal ritmo medio finora.
t0      = tic;
nPicTot = 0;
nextMark = 0.10;
for it = 1:nSteps
    tNow = it*dt;
    Inow = min(dIdt*tNow, I_nom);
    i_s = zeros(Nn,1);
    i_s(faceExt_dae) = -Inow/numel(faceExt_dae);   % ESTRAZIONE (segno corretto)
    for pic = 1:20                                 % Picard sulla non linearita' R(i)
        iEdge = x(1:Ne);
        Rvec  = resistanceEJ(iEdge, IcEdge, segLen, E0, nEJ);
        K     = [spdiags(Rvec,0,Ne,Ne), sparse(A); -sparse(A).', Gmat];
        Asys  = M/dt + K;
        b     = M/dt*x + [zeros(Ne,1); i_s];
        if it == 1 && pic == 1
            % DIAGNOSTICA (solo al primo solve): chiedendo il terzo output
            % directSolve invoca cond1est e restituisce la stima del numero
            % di condizione in norma 1. Serve a sapere se fidarsi del solve.
            [xNew, ~, condEst] = directSolve(Asys, b, ...
                            Ne + faceInj_dae, zeros(numel(faceInj_dae),1), ...
                            'equivalue', {Ne + faceExt_dae}, 'method','lu', ...
                            'iterMAX', nIterRef, 'verbose','n');
            fprintf('  numero di condizione stimato (cond1est) = %.2e\n', condEst);
            if condEst > 1e14
                fprintf('  -> sistema mal condizionato (atteso: punto sella con G = 0):\n');
                fprintf('     raffinamento iterativo attivo (%d passi).\n', nIterRef);
            end
        else
            xNew  = directSolve(Asys, b, ...
                            Ne + faceInj_dae, zeros(numel(faceInj_dae),1), ...
                            'equivalue', {Ne + faceExt_dae}, 'method','lu', ...
                            'iterMAX', nIterRef, 'verbose','n');
        end
        if norm(xNew - x)/max(norm(xNew),eps) < 1e-8, x = xNew; break; end
        x = xNew;
    end
    nPicTot = nPicTot + pic;
    hist_t(it) = tNow;
    hist_i(it,:) = accumarray(sid_dae, x(1:Ne), [cab.N_slots 1], @mean).';
    if it/nSteps >= nextMark
        el = toc(t0);
        fprintf('    avanzamento %3.0f%%  (%4d/%d passi, %6.1f s trascorsi, ~%5.1f s residui)\n', ...
            it/nSteps*100, it, nSteps, el, el*(nSteps-it)/it);
        nextMark = nextMark + 0.10;
    end
end
tDAE = toc(t0);
tmr = addTimer(tmr, 'Step 4.4  integrazione temporale DAE', tDAE);
fprintf('  %d passi in %.1f s (%.3f s/passo) | %.1f iterazioni di Picard per passo\n', ...
    nSteps, tDAE, tDAE/nSteps, nPicTot/nSteps);

iFinal  = x(1:Ne);
iStrand = accumarray(sid_dae, iFinal, [cab.N_slots 1], @mean);

% ------------------------------------------------------------------------
%  4.5  DIMOSTRAZIONE DELL'EQUIRIPARTIZIONE
% ------------------------------------------------------------------------
fprintf('\n--- 4.5 Dimostrazione dell''equiripartizione ---\n');
fprintf('  correnti finali per strand [A]: '); fprintf('%.1f ', iStrand); fprintf('\n');
fprintf('  atteso I_nom/N_slots           : %.1f A\n', I_slot);
devMax = max(abs(iStrand - I_slot))/I_slot*100;
fprintf('  scarto massimo                 : %.4f %%\n', devMax);

% (a) conservazione della corrente totale
fprintf('  somma delle correnti = %.1f A (imposta: %.1f A, errore %.3e)\n', ...
    sum(iStrand), I_nom, abs(sum(iStrand)-I_nom)/I_nom);

% (a-bis) residuo del sistema lineare risolto (verifica del solve, non della fisica).
%  ATTENZIONE: directSolve RISCRIVE le righe dei vincoli (potenziale imposto
%  sulle facce di iniezione, equipotenzialita' su quelle di estrazione), quindi
%  il sistema risolto NON e' Asys*x = b su quelle righe. Misurare il residuo
%  sulle matrici grezze dava ||Ax-b||/||b|| ~ 1, un numero privo di significato:
%  la soluzione non deve soddisfare quelle righe. Il residuo si valuta sulle
%  righe LIBERE, e i vincoli si verificano a parte.
bcRows = [Ne + faceInj_dae(:); Ne + faceExt_dae(:)];
freeR  = true(size(b));  freeR(bcRows) = false;
rvec   = Asys*x - b;
resid  = norm(rvec(freeR))/max(norm(b(freeR)), eps);
vInj   = x(Ne + faceInj_dae);
vExt   = x(Ne + faceExt_dae);
fprintf('  residuo sulle righe libere ||Ax-b||/||b|| = %.3e  (cond stimato %.2e)\n', ...
    resid, condEst);
fprintf('  vincoli: |V| max sulle facce di iniezione = %.2e V | dispersione su\n', ...
    max(abs(vInj)));
fprintf('    quelle di estrazione = %.2e V (devono essere ~0)\n', max(vExt)-min(vExt));

% (b) uniformita' lungo ciascuno strand (continuita' in serie)
spreadInStrand = zeros(cab.N_slots,1);
for j = 1:cab.N_slots
    ij = iFinal(sid_dae == j);
    spreadInStrand(j) = (max(ij)-min(ij))/mean(ij)*100;
end
fprintf('  variazione della corrente LUNGO ciascuno strand: max %.2e %%\n', ...
    max(spreadInStrand));

fprintf('  dispersione induttiva prevista in 4.2: %.3f %% (confronto con lo scarto)\n', sprdL);
if devMax < 2
    fprintf('\n  ==> EQUIRIPARTIZIONE DIMOSTRATA\n');
    fprintf('      Il DAE e'' stato risolto senza imporre la ripartizione: le correnti\n');
    fprintf('      risultano uguali perche'' il twist commensurato rende gli strand\n');
    fprintf('      elettricamente equivalenti (stessa somma di R, vedi 4.3).\n');
else
    warning('step4:noEqui', ...
        ['Scarto %.2f%% dall''equiripartizione (atteso ~%.2f%% dalla dispersione ' ...
         'induttiva). Se molto maggiore di 2%%, verificare la risoluzione del ' ...
         'twist in 2.1c.'], devMax, sprdL);
end

% ------------------------------------------------------------------------
%  4.6  CONTROPROVA: il test non e' una tautologia
% ------------------------------------------------------------------------
%  Si ripete l'ultimo passo con gli strand deliberatamente NON equivalenti
%  (Ic perturbata strand per strand, come se il twist non fosse commensurato).
%  Se il solutore restituisse comunque correnti uguali, la verifica di 4.5
%  sarebbe priva di significato.
runControl = true;
if runControl
    fprintf('\n--- 4.6 Controprova (strand resi NON equivalenti) ---\n');
    t0 = tic;
    pert = 1 + 0.30*(sid_dae - (cab.N_slots+1)/2)/cab.N_slots;   % +/- 15% circa
    IcPert = IcEdge .* pert;
    xC = zeros(Ne+Nn,1);
    i_sC = zeros(Nn,1);  i_sC(faceExt_dae) = -I_nom/numel(faceExt_dae);
    for pic = 1:40
        RvecC = resistanceEJ(xC(1:Ne), IcPert, segLen, E0, nEJ);
        KC    = [spdiags(RvecC,0,Ne,Ne), sparse(A); -sparse(A).', Gmat];
        bC    = [zeros(Ne,1); i_sC];
        xNewC = directSolve(KC, bC, Ne + faceInj_dae, zeros(numel(faceInj_dae),1), ...
                            'equivalue', {Ne + faceExt_dae}, 'method','lu', ...
                            'verbose','n');
        if norm(xNewC - xC)/max(norm(xNewC),eps) < 1e-8, xC = xNewC; break; end
        xC = xNewC;
    end
    iStrandC = accumarray(sid_dae, xC(1:Ne), [cab.N_slots 1], @mean);
    devC = max(abs(iStrandC - I_slot))/I_slot*100;
    fprintf('  correnti con strand perturbati [A]: '); fprintf('%.1f ', iStrandC); fprintf('\n');
    fprintf('  scarto: %.2f %%  (nel caso commensurato era %.4f %%)\n', devC, devMax);
    if devC > 5*max(devMax, 1e-3)
        fprintf('  ==> Il modello DISTINGUE i due casi: la verifica di 4.5 e'' significativa.\n');
    else
        warning('step4:control', ...
            'La controprova non distingue i casi: il test di equiripartizione non e'' conclusivo.');
    end
    tmr = addTimer(tmr, 'Step 4.6  controprova (strand perturbati)', toc(t0));
end

% ------------------------------------------------------------------------
%  4.7  Grafici del DAE
% ------------------------------------------------------------------------
figure('Name','STEP 4 - DAE ed equiripartizione','Color','w','Position',[120 120 1300 420]);

subplot(1,3,1); hold on; grid on; box on;
cmapS4 = lines(cab.N_slots);
for j = 1:cab.N_slots
    plot(hist_t, hist_i(:,j), '-', 'Color', cmapS4(j,:), 'LineWidth',1.3);
end
plot(hist_t, dIdt*hist_t/cab.N_slots, 'k--','LineWidth',1.0);
xlabel('t [s]'); ylabel('corrente per strand [A]');
title('Correnti durante la rampa');
lg = arrayfun(@(j) sprintf('strand %d',j), 1:cab.N_slots, 'UniformOutput', false);
legend([lg, {'I(t)/N'}], 'Location','northwest','FontSize',7);

subplot(1,3,2); hold on; grid on; box on;
bar(1:cab.N_slots, iStrand); yline(I_slot,'r--','I_{nom}/N_{slots}','LineWidth',1.2);
xlabel('strand'); ylabel('corrente finale [A]');
title(sprintf('Ripartizione finale (scarto %.1e %%)', devMax));

subplot(1,3,3); hold on; grid on; box on;
if runControl
    bar(1:cab.N_slots, [iStrand, iStrandC]);
    legend('twist commensurato','strand perturbati','Location','best','FontSize',8);
else
    bar(1:cab.N_slots, iStrand);
end
yline(I_slot,'r--','LineWidth',1.2);
xlabel('strand'); ylabel('corrente [A]'); title('Controprova');


%% ========================================================================
%  SALVATAGGIO
%  ========================================================================
save('sc_project_results.mat', ...
     'data','I_nom','R_in','H_max','B_target','dIdt','dT_req','T_op_guess', ...   % Step 1
     'cab','disc','N_layer','idxCentral','z_layers','N_turn', ...                 % Step 2
     'P_macro','E_macro','P_micro','E_micro','strandID','faceInj','faceExt', ...
     'P_center','th_center','s_center','Tv','Nv','Bv','I_slot', ...
     'B_tot','B_ext','B_mut','B_self_scalar','Bmag_tot','Bpk','idxPk', ...        % Step 3
     'B_perp','B_par','theta_deg','s_nodes','r_nodes','Bprofile','rProfile', ...
     'Bmag_near','Tcs_profile','Tcs_cert','nCert','margin_profile','margin_min','idxCrit', ...
     'r_near','depth_stack','I_tape_op', ...
     'L','A','iStrand','iFinal','IcEdge','Bedge','hist_t','hist_i','Rsum','condEst', ...
     'Lstrand','sprdL','segPerTwist_dae', ...
     'r0','A_stack','sep_strand','r_near','depth_stack', ...    % Step 4
     '-v7.3');
fprintf('\n\n--> Salvato sc_project_results.mat\n');

%% ------------------------- RIEPILOGO ------------------------------------
fprintf('\n=========================================================\n');
fprintf(' RIEPILOGO DEL DESIGN\n');
fprintf('=========================================================\n');
fprintf('  Solenoide : %d DP x %d turn radiali, R_in = %.2f m, H = %.2f m\n', ...
    N_layer, N_turn, R_in, max(z_layers)-min(z_layers)+2*cab.width);
fprintf('  Cavo      : %d slot, %d nastri/slot, twist %d giri/turn\n', ...
    cab.N_slots, cab.N_tapes_slot, cab.N_twist);
fprintf('  Correnti  : %.0f kA totali -> %.0f A/slot -> %.1f A/nastro\n', ...
    I_nom/1e3, I_slot, I_tape_op);
fprintf('  Campo     : B_pk = %.2f T (target %.1f T)\n', Bpk, B_target);
if Tcs_cert(idxCrit), stt = 'CERTIFICATA dai dati'; else, stt = 'solo limite inferiore'; end
fprintf('  Margine   : T_cs = %.2f K (%s), T_op = %.1f K -> %.2f K (richiesto %.0f K)\n', ...
    Tcs_profile(idxCrit), stt, T_op_guess, margin_min, dT_req);
fprintf('              nodi con T_cs certificata: %d/%d\n', nCert, numel(Tcs_cert));
fprintf('  DAE       : scarto dall''equiripartizione = %.2e %%\n', devMax);
fprintf('---------------------------------------------------------\n');
fprintf('  n dalla tabella n(T,B) dello Step 1, valutata arco per arco (4.3b).\n');
fprintf('=========================================================\n');


%% ========================================================================
%  STEP 5 - CICC LAYOUT: sezione trasversale e RISOLUZIONE DEL CENTROIDE
%  ========================================================================
%  R_core, r0 e r_near sono grandezze distinte:
%    R_core = posizione del centroide dello stack rispetto all'asse del cavo;
%    r0     = raggio equivalente associato alla sezione dello stack;
%    r_near = bordo interno dello stack.
%  La sezione viene costruita dall'esterno verso l'interno e R_core e'
%  ricavato geometricamente dallo slot.

fprintf('\n\n########## STEP 5 - CICC LAYOUT ##########\n');

% ------------------------------------------------------------------------
%  5.1  Parametri costruttivi
% ------------------------------------------------------------------------
% Questi tre spessori NON derivano dal Data Book o dal brief:
% sono ASSUNZIONI GEOMETRICHE di progetto della sezione CICC.
lay.t_jacket  = 2.0e-3;   % [m] ASSUNTO: jacket in acciaio
lay.t_web_out = 1.0e-3;   % [m] ASSUNTO: parete Cu tra slot e jacket
lay.t_web_in  = 1.0e-3;   % [m] ASSUNTO: parete Cu tra slot e canale centrale

% Anche i giochi nulli sono una SCELTA DI TOPOLOGIA, non grandezze derivate.
% La profondita' e la larghezza dello slot sono invece ricavate da:
%   slot_d = N_tapes_slot * t_tape + clr_r
%   slot_w = tape_width + clr_t
% Con clr_r = clr_t = 0 lo stack riempie esattamente lo slot e l'He scorre
% soltanto nel canale centrale.
lay.clr_r     = 0;        % [m] ASSUNTO: gioco radiale nullo
lay.clr_t     = 0;        % [m] ASSUNTO: gioco tangenziale nullo

% ------------------------------------------------------------------------
%  5.2  Sezione derivata -> R_core (RISOLUZIONE DEL CENTROIDE)
% ------------------------------------------------------------------------
r_cab          = cab.width/2;
lay.slot_d     = depth_stack + lay.clr_r;            % profondita' slot
lay.slot_w     = data.tape.width.val + lay.clr_t;    % larghezza slot
r_former_out   = r_cab        - lay.t_jacket;
r_slot_out     = r_former_out - lay.t_web_out;
R_core_derived = r_slot_out   - lay.slot_d/2;        % <<< CENTROIDE DERIVATO
r_slot_in      = r_slot_out   - lay.slot_d;
r_channel      = r_slot_in    - lay.t_web_in;

fprintf('\n--- 5.2 Sezione trasversale e centroide ---\n');
fprintf('  r_cavo            = %5.2f mm\n', r_cab*1e3);
fprintf('  jacket SS   %5.2f -> r_former_out = %5.2f mm\n', lay.t_jacket*1e3, r_former_out*1e3);
fprintf('  web esterno %5.2f -> slot esterno = %5.2f mm\n', lay.t_web_out*1e3, r_slot_out*1e3);
fprintf('  slot (%d nastri) %5.2f mm\n', cab.N_tapes_slot, lay.slot_d*1e3);
fprintf('     -> R_core (CENTROIDE) = %5.3f mm   <<< DERIVATO dalla sezione\n', ...
    R_core_derived*1e3);
fprintf('        r0 (dimensione)    = %5.3f mm   (dalla sezione dello stack)\n', cab.r0*1e3);
fprintf('        r_near (bordo int) = %5.3f mm   (per il margine)\n', ...
    (R_core_derived-depth_stack/2)*1e3);
fprintf('  web interno %5.2f -> canale He    = %5.2f mm (diam %.2f mm)\n', ...
    lay.t_web_in*1e3, r_channel*1e3, 2*r_channel*1e3);

if r_channel <= 0
    error('step5:noChannel', ...
        'La sezione non chiude: nessuno spazio per il canale. Ridurre stack/jacket/web.');
end

% ------------------------------------------------------------------------
%  5.3  Verifiche geometriche
% ------------------------------------------------------------------------
sep_der    = 2*R_core_derived*sin(pi/cab.N_slots);
arc_slot   = 2*pi*R_core_derived/cab.N_slots;
web_tang   = arc_slot - lay.slot_w;
r_near_der = R_core_derived - depth_stack/2;

fprintf('\n--- 5.3 Verifiche geometriche ---\n');
fprintf('  fit angolare: arco/slot %.2f mm - slot %.2f mm = web %.2f mm  [%s]\n', ...
    arc_slot*1e3, lay.slot_w*1e3, web_tang*1e3, ternary(web_tang>0,'OK','NON CI STA'));
fprintf('  r0/sep = %.3f (errore >=0.5, warning >0.4)                    [%s]\n', ...
    cab.r0/sep_der, ternary(cab.r0/sep_der<0.4,'OK','LIMITE'));
fprintf('  r_near = %.2f mm > 0                                          [%s]\n', ...
    r_near_der*1e3, ternary(r_near_der>0,'OK','NEGATIVO'));
if web_tang <= 0
    error('step5:slotsDontFit','Gli slot non entrano nella circonferenza a R_core = %.2f mm.', ...
        R_core_derived*1e3);
end

% ------------------------------------------------------------------------
%  5.4  Aree e perimetri (input per il modello termico dello Step 7)
% ------------------------------------------------------------------------
A_cable   = pi*r_cab^2;
A_channel = pi*r_channel^2;
%  A_void e' NULLA per costruzione (slot = stack, vedi 5.1). Si calcola lo
%  stesso e si verifica, invece di darlo per scontato: se qualcuno rimettesse
%  un gioco in 5.1 senza accorgersi delle conseguenze, il controllo lo becca.
A_void    = cab.N_slots*(lay.slot_w*lay.slot_d - A_stack);
if A_void > 1e-12
    error('step5:voidNotZero', ...
        ['A_void = %.3f mm^2 invece di 0: gli slot non sono riempiti ' ...
         'esattamente dallo stack. Con giochi non nulli il modello a UN solo ' ...
         'canale dello Step 7 non e'' piu'' valido, perche'' ci sarebbe elio ' ...
         'anche negli slot. Azzerare lay.clr_r e lay.clr_t, oppure estendere ' ...
         'lo Step 7 a due regioni di fluido.'], A_void*1e6);
end
A_He      = A_channel;         % l'elio scorre SOLO nel canale centrale
A_SC      = cab.N_slots*A_stack;
A_SS      = pi*(r_cab^2 - r_former_out^2);
A_Cu      = pi*(r_former_out^2 - r_channel^2) - cab.N_slots*lay.slot_w*lay.slot_d;
P_He_channel = 2*pi*r_channel;
P_stack_He   = cab.N_slots*2*(lay.slot_w + lay.slot_d);
P_Cu_SS      = 2*pi*r_former_out;
P_He_wet     = P_He_channel;                   % bagnato SOLO il canale
%  DIAMETRO IDRAULICO. Con l'elio confinato nel canale centrale, che e' un
%  condotto circolare, 4*A/P e 2*r coincidono identicamente: non c'e' nessuna
%  scelta da fare. Si calcola nella forma generale e si verifica l'identita',
%  cosi' il controllo resta valido anche se la sezione del canale cambiasse.
D_h          = 4*A_He/P_He_wet;                % [m]
if abs(D_h - 2*r_channel) > 1e-9
    error('step5:DhMismatch', ...
        'D_h = %.4f mm diverso da 2*r_canale = %.4f mm: la regione di elio non e'' il solo canale.', ...
        D_h*1e3, 2*r_channel*1e3);
end

fprintf('\n--- 5.4 Aree e perimetri (per lo Step 7) ---\n');
fprintf('  %-13s %8s %8s\n','regione','[mm^2]','[%sez]');
for kk = 1:6
    switch kk
        case 1, nm='canale He'; vv=A_channel;
        case 2, nm='void slot'; vv=A_void;
        case 3, nm='He TOTALE'; vv=A_He;
        case 4, nm='stack SC';  vv=A_SC;
        case 5, nm='former Cu'; vv=A_Cu;
        case 6, nm='jacket SS'; vv=A_SS;
    end
    fprintf('  %-13s %8.2f %8.1f\n', nm, vv*1e6, vv/A_cable*100);
end
fprintf('  D_h = %.3f mm (= diametro del canale, unica regione di elio)\n', D_h*1e3);
fprintf('  P_bagnato = %.2f mm (solo parete del canale)\n', P_He_wet*1e3);
fprintf('  P_stack-former = %.2f mm (contatto SOLIDO-SOLIDO, non bagnato)\n', P_stack_He*1e3);

% ------------------------------------------------------------------------
%  5.5  Verifica idraulica preliminare
% ------------------------------------------------------------------------
he.rho = 65.0;   % [kg/m3] He supercritico ~10 K, 5 bar [letteratura NIST/CoolProp]
he.mu  = 1.5e-6; % [Pa.s]
mdot   = data.op.mdot_max.val;
L_dp   = s_center(end);
v_He   = mdot/(he.rho*A_He);
Re_He  = he.rho*v_He*D_h/he.mu;
f_D    = 0.316*Re_He^-0.25;
dp_He  = f_D*(L_dp/D_h)*he.rho*v_He^2/2;
t_transit = L_dp/v_He;

fprintf('\n--- 5.5 Idraulica preliminare ---\n');
fprintf('  mdot = %.1f g/s | L_DP = %.2f m\n', mdot*1e3, L_dp);
fprintf('  v = %.2f m/s | Re = %.2e (%s) | f = %.4f\n', v_He, Re_He, ...
    ternary(Re_He>4000,'turbolento','laminare'), f_D);
fprintf('  dp = %.4f bar  [%s]\n', dp_He/1e5, ternary(dp_He/1e5<10,'OK','ECCESSIVO'));
fprintf('  transito He = %.1f s (rampa %.0f s -> transitorio lento)\n', t_transit, I_nom/dIdt);
fprintf('  NOTA: stima con rho e mu COSTANTI. Lo Step 7 rifa'' il calcolo con le\n');
fprintf('        proprieta'' reali dell''He lungo il canale (rho varia di ~4x).\n');

% ------------------------------------------------------------------------
%  5.6  CHIUSURA DEL LOOP A
% ------------------------------------------------------------------------
fprintf('\n--- 5.6 Chiusura del Loop A ---\n');
dR = abs(R_core_derived - cab.R_core);
fprintf('  R_core usato in Step 2 : %.3f mm\n', cab.R_core*1e3);
fprintf('  R_core derivato qui    : %.3f mm  (differenza %.3f mm)\n', R_core_derived*1e3, dR*1e3);
if dR < 0.2e-3
    fprintf('  ==> LOOP A CHIUSO: geometria autoconsistente.\n');
else
    warning('step5:loopA', ...
        'R_core incoerente (%.2f mm). Impostare cab.R_core = %.5f nello Step 2.1 e rilanciare.', ...
        dR*1e3, R_core_derived);
end

% ------------------------------------------------------------------------
%  5.7  Disegno della sezione
% ------------------------------------------------------------------------
figure('Name','STEP 5 - sezione CICC','Color','w','Position',[140 140 660 640]);
hold on; axis equal; grid on; box on;
thc = linspace(0,2*pi,200);
fill(r_cab*1e3*cos(thc),        r_cab*1e3*sin(thc),        [0.75 0.75 0.78]);
fill(r_former_out*1e3*cos(thc), r_former_out*1e3*sin(thc), [0.85 0.55 0.30]);
fill(r_channel*1e3*cos(thc),    r_channel*1e3*sin(thc),    [0.60 0.80 0.95]);
for j = 1:cab.N_slots
    a  = (j-1)*2*pi/cab.N_slots;
    xc = R_core_derived*cos(a)*1e3;  yc = R_core_derived*sin(a)*1e3;
    Rm = [cos(a) -sin(a); sin(a) cos(a)];
    cor = Rm*[-lay.slot_d/2 lay.slot_d/2 lay.slot_d/2 -lay.slot_d/2;
              -lay.slot_w/2 -lay.slot_w/2 lay.slot_w/2 lay.slot_w/2]*1e3;
    fill(xc+cor(1,:), yc+cor(2,:), [0.95 0.95 0.55]);
    cs = Rm*[-depth_stack/2 depth_stack/2 depth_stack/2 -depth_stack/2;
             -data.tape.width.val/2 -data.tape.width.val/2 ...
              data.tape.width.val/2  data.tape.width.val/2]*1e3;
    fill(xc+cs(1,:), yc+cs(2,:), [0.15 0.15 0.15]);
    plot(xc, yc, 'r+', 'MarkerSize',9, 'LineWidth',1.5);   % centroide
end
xlabel('[mm]'); ylabel('[mm]');
title(sprintf(['Sezione CICC: %d slot x %d nastri\n' ...
    'R_{core} = %.2f mm (centroide, +), canale %.1f mm, D = %.0f mm'], ...
    cab.N_slots, cab.N_tapes_slot, R_core_derived*1e3, 2*r_channel*1e3, cab.width*1e3));
legend({'jacket SS','former Cu','canale He','slot','stack SC','centroide'}, ...
    'Location','eastoutside','FontSize',8);

save('step5_layout.mat','lay','R_core_derived','r_channel','r_former_out', ...
     'A_He','A_channel','A_void','A_SC','A_Cu','A_SS', ...
     'P_He_channel','P_stack_He','P_Cu_SS','D_h','v_He','Re_He','dp_He','t_transit');
fprintf('\n--> Salvato step5_layout.mat (aree e perimetri per lo Step 7)\n');

%% ========================================================================
%  STEP 6 - AC LOSSES durante la rampa 0 -> 50 kA
%  ========================================================================
%  Calcola la sorgente di calore Q(s) [W/m] lungo il conduttore, che alimenta
%  il modello termo-idraulico dello Step 7.
%
%  PERCHE' FORMULE CHIUSE E NON IL DAE DELLO STEP 4:
%  nello Step 4 la scelta di scope e' G = 0 (nessuna conduttanza trasversale).
%  Con G = 0 il DAE NON PUO' dare le perdite di accoppiamento, che sono per
%  definizione le correnti indotte che circolano ATTRAVERSO quelle conduttanze.
%  Analogamente non vede l'isteresi intra-nastro, perche' un arco lumped non ha
%  profilo di penetrazione del flusso. Si usano quindi formule di letteratura.
%
%  TERMINI
%    (a) isteresi di MAGNETIZZAZIONE (stato critico)
%    (b) isteresi di TRASPORTO (Norris)
%    (c) accoppiamento inter-strand (modello tau)
%    (d) eddy nei metalli (former Cu, jacket SS)
%
%  ENERGIA PER CICLO vs POTENZA - la distinzione che spiega l'apparente
%  contraddizione con la regola di scala P_hyst ~ |dB/dt|, P_coup ~ (dB/dt)^2:
%    - Le formule (a) e (b) danno una ENERGIA per ciclo [J/m]. Nel modello di
%      stato critico le correnti di schermo scorrono sempre a +/-Jc quale che sia
%      la velocita' di rampa: l'area del ciclo di isteresi dipende SOLO
%      dall'escursione di campo/corrente, non da quanto tempo ci si mette.
%      L'energia per ciclo e' quindi indipendente da dB/dt.
%    - La POTENZA e' quell'energia divisa per la durata: P = Q/t_rampa, e siccome
%      t_rampa = escursione/(dB/dt) si ottiene P_hyst proporzionale a |dB/dt|.
%      La regola di scala vale quindi per la POTENZA, ed e' esattamente cio' che
%      il codice produce (Q [J/m] calcolata dall'escursione, poi / t_ramp).
%    - Per l'accoppiamento e' l'opposto: le correnti indotte sono proporzionali a
%      dB/dt e dissipano I^2*R, quindi P ~ (dB/dt)^2 e l'energia per ciclo ~ dB/dt,
%      cioe' SVANISCE a rampa lenta. E' questa differenza che permette di separare
%      i due contributi sperimentalmente rampando a ~0.01 T/s (isteresi pura).
%  A 0.2 kA/s (dB/dt ~ 0.06 T/s) ci si aspetta quindi isteresi dominante: il
%  codice lo VERIFICA numericamente invece di assumerlo.
fprintf('\n\n########## STEP 6 - AC LOSSES ##########\n');
t0 = tic;

% ------------------------------------------------------------------------
%  6.1  Scenario: UNA rampa monotona
% ------------------------------------------------------------------------
%  Rampa singola 0 -> 50 kA in 250 s, nessun ciclo, nessun flat-top.
%
%  LE FORMULE ADOTTATE DANNO DIRETTAMENTE LA POTENZA, non l'energia per ciclo:
%     P_hyst = (2/(3*pi)) * Jc * d_eff * |dB/dt| * A_sc
%     P_coup = (n*tau/mu0) * (dB/dt)^2 * A
%  Entrambe contengono dB/dt, quindi si applicano a una rampa monotona senza
%  dover convertire da "energia per ciclo" a potenza. Sparisce cosi' il
%  La formula usata e' direttamente quella definita per la rampa considerata;
%  non e' necessario introdurre fattori correttivi di ciclo.
%
%  IPOTESI DI STATO VERGINE: il magnete parte smagnetizzato (raffreddato a
%  campo nullo). Una risalita successiva senza superare Tc costerebbe di piu'.
%  6.2  Regime di penetrazione: l'asintoto usato in 6.3 e' valido?
% ------------------------------------------------------------------------
%  Campo caratteristico di una strip sottile:  Bc = mu0*Ic/(pi*w)
%                          [Brandt & Indenbom, Phys. Rev. B 48, 12893 (1993)]
%  Se l'escursione >> Bc la strip e' COMPLETAMENTE PENETRATA e la perdita per
%  ciclo tende all'asintoto lineare. Va verificato, non assunto.
Ic_at_op = data.Ic.eval(T_op_guess, Bpk);
Bc_strip = Mu0()*Ic_at_op/(pi*data.tape.width.val);
fprintf('\n--- 6.2 Regime di penetrazione ---\n');
fprintf('  Ic(T_op,B_pk) = %.1f A -> Bc = mu0*Ic/(pi*w) = %.4f T\n', Ic_at_op, Bc_strip);
fprintf('  escursione = %.2f T -> rapporto = %.0f  [%s]\n', Bpk, Bpk/Bc_strip, ...
    ternary(Bpk/Bc_strip>10,'PENETRAZIONE COMPLETA: asintoto valido', ...
                            'regime intermedio: asintoto SOVRASTIMA'));
if Bpk/Bc_strip <= 10
    warning('step6:partialPenetration', ...
        'Rapporto %.1f < 10: strip non completamente penetrata, asintoto sovrastima.', ...
        Bpk/Bc_strip);
end

%  QUANTO DELLA RAMPA E' IN REGIME DI PENETRAZIONE PARZIALE? La verifica
%  sopra e' al campo FINALE, ma la 6.3c integra la potenza anche sui primi
%  istanti della rampa, dove il campo locale e' piu' basso e l'asintoto puo'
%  non valere ancora. Si CALCOLA la frazione di rampa (in lambda, che qui
%  coincide con la frazione di tempo, rampa lineare) in cui il nodo di picco
%  ha campo sotto 10*Bc, invece di assumerla piccola.
lam_soglia_pen = min(10*Bc_strip/Bpk, 1);
fprintf('  regime di penetrazione parziale (B < 10*Bc) per lambda < %.4f,\n', lam_soglia_pen);
fprintf('    cioe'' il primo %.1f%% della rampa (in campo/tempo, rampa lineare)\n', ...
    lam_soglia_pen*100);
if lam_soglia_pen > 0.05
    warning('step6:partialPenetrationRamp', ...
        ['La finestra a penetrazione parziale copre il %.1f%% della rampa: non e'' ' ...
         'piu'' ovviamente trascurabile, andrebbe pesata nell''integrale invece di ' ...
         'usare solo l''asintoto.'], lam_soglia_pen*100);
else
    fprintf('    -> sotto il 5%%: trascurabile, l''asintoto di piena penetrazione usato\n');
    fprintf('       in 6.3/6.3c resta valido per la quasi totalita'' della rampa\n');
end

% ------------------------------------------------------------------------
%  6.3  Isteresi (termine dominante)
% ------------------------------------------------------------------------
%     P_hyst = (2/(3*pi)) * Jc * d_eff * |dB/dt| * A_sc          [W/m]
%
%  DUE LETTURE DEI SIMBOLI, e conviene esplicitarle perche' l'una senza l'altra
%  porta a sbagliare di ordini di grandezza:
%    Jc * A_sc = Ic  -- Jc e' la densita' di corrente critica SULLA SEZIONE
%                       SUPERCONDUTTRICE, quindi il prodotto con A_sc e'
%                       semplicemente la corrente critica. Si usa direttamente
%                       Ic, evitando di dover conoscere lo spessore dello
%                       strato REBCO (2.5 um) e la sua incertezza.
%    d_eff          -- dimensione efficace di penetrazione del flusso. Per un
%                       nastro in campo PERPENDICOLARE le correnti di schermo
%                       circolano sulla LARGHEZZA, quindi d_eff = w = 4 mm.
%                       Non lo spessore: quello vale per campo parallelo, e
%                       darebbe un risultato 1000 volte piu' piccolo.
%
%  Si usa B_perp, la componente normale alla faccia calcolata nello Step 3
%  nella terna del nastro: e' quella che penetra la lamina. La componente
%  parallela offre al flusso solo lo spessore del film.
%
%  Il termine e' proporzionale a Ic, quindi cresce dove il campo e' MINORE.
%  Il prodotto Ic*B_perp ha percio' un massimo intermedio lungo il DP, non
%  agli estremi: e' il motivo per cui il profilo va calcolato nodo per nodo
%  invece che al solo punto di picco.
fprintf('\n--- 6.3 Isteresi ---\n');
Ic_local = arrayfun(@(b) data.Ic.eval(T_op_guess, b), Bmag_near);   % [A] per nastro
acl.d_eff = data.tape.width.val;    % [m] campo perpendicolare -> larghezza del nastro
dBdt_perp = B_perp / t_ramp;        % [T/s] locale (B proporzionale a I lungo la rampa)

P_hyst_tape = (2/(3*pi)) * Ic_local .* acl.d_eff .* dBdt_perp;   % [W/m] per nastro
P_mag_nodo  = P_hyst_tape * cab.N_tapes_slot;                    % [W/m] per strand
fprintf('  Ic locale nastro: %.1f ... %.1f A | dB_perp/dt: %.4f ... %.4f T/s\n', ...
    min(Ic_local), max(Ic_local), min(dBdt_perp), max(dBdt_perp));
fprintf('  d_eff = %.1f mm (larghezza del nastro, campo perpendicolare)\n', acl.d_eff*1e3);
fprintf('  P_isteresi per strand: max %.4f W/m\n', max(P_mag_nodo));

% ------------------------------------------------------------------------
%  6.3c  ISTERESI: DIPENDENZA TEMPORALE Q(s,t)
% ------------------------------------------------------------------------
%  Le righe sopra usano Ic al campo FINALE della rampa: e' il valore che
%  serve a tutti gli usi puntuali di Ic_local nel resto del file (grafici,
%  Bc_strip, la stampa di sintesi). Qui si costruisce IN PIU' la versione
%  che dipende dal tempo, perche' durante la rampa il campo locale non e'
%  quello finale mai raggiunto: e' lambda(t)*B_near(s), con
%
%     lambda(t) = I(t)/I_nom,   0 -> 1 durante la rampa lineare
%
%  dB_perp/dt resta costante (la rampa e' lineare in I, quindi in B), ma
%  Ic(T_op, campo) NO: cresce quando lambda scende, perche' il campo scende.
%  Usare Ic(B_finale) per tutta la rampa vuol dire usare il valore piu'
%  BASSO di Ic in ogni istante tranne l'ultimo, quindi SOTTOSTIMARE la
%  potenza media. Non e' una correzione di dettaglio: il fattore fra la
%  potenza mediata sulla rampa e quella al campo finale, stampato sotto, e'
%  vicino a 2, non a qualche punto percentuale.
%
%  cool.lam_grid e' definita qui (non in Step 7) perche' e' il campo, non il
%  raffreddamento, a impostare la fisica del problema: lo Step 7 la eredita
%  e la usa solo per interpolare nel tempo.
cool.n_lam    = 21;                          % livelli di rampa, 0 e 1 inclusi
cool.lam_grid = linspace(0, 1, cool.n_lam).';
lam_row       = cool.lam_grid.';             % riga, per il prodotto esterno sul campo

%  data.Ic.eval accetta solo B SCALARE (e' per questo che 6.3 la chiama con
%  arrayfun su Bmag_near): la stessa idiom si applica qui sulla matrice
%  campo, che arrayfun tratta elemento per elemento preservando la forma.
Bmat_hyst  = Bmag_near * lam_row;                              % [nodi x n_lam]
Ic_local_t = arrayfun(@(b) data.Ic.eval(T_op_guess, b), Bmat_hyst);
P_hyst_tape_t = (2/(3*pi)) * Ic_local_t .* acl.d_eff .* dBdt_perp;  % dBdt_perp costante
P_mag_nodo_t  = P_hyst_tape_t * cab.N_tapes_slot;               % [nodi x n_lam]

%  Fattore effettivo al nodo di picco, e fattore integrale pesato sulle
%  potenze VERE del run (non un peso assunto a priori): serve a confrontare
%  col conto analitico fatto a mano prima di scrivere questo codice.
%  trapz, non mean: sulla griglia lambda uniforme le due differiscono di
%  qualche punto percentuale, e trapz e' l'integrale coerente con
%  l'interpolazione lineare usata poi in gandalfSolve.
fatt_hyst_pk = trapz(cool.lam_grid, P_mag_nodo_t(idxPk,:)) / P_mag_nodo_t(idxPk,end);
fprintf('  fattore temporale isteresi (media rampa / campo finale), nodo di picco: %.2f\n', ...
    fatt_hyst_pk);
if fatt_hyst_pk < 1.2
    warning('step6:hystTimeFlat', ...
        ['Il fattore temporale sull''isteresi (%.2f) e'' insolitamente vicino a 1: ' ...
         'verificare che lam_row stia davvero campionando 0..1 e che Bmat_hyst sia ' ...
         'il prodotto esterno atteso, non una copia di Bmag_near.'], fatt_hyst_pk);
end

% ------------------------------------------------------------------------
%  6.4  Isteresi di TRASPORTO (Norris, strip)
% ------------------------------------------------------------------------
%   Q = (mu0*Ic^2/pi)*[(1-F)ln(1-F) + (1+F)ln(1+F) - F^2] ,  F = I_tape/Ic
%                                   [W.T. Norris, J. Phys. D 3, 489 (1970)]
%
%  NON E' UNO DEI DUE TERMINI DELLE FORMULE ADOTTATE: quelle coprono l'isteresi
%  di MAGNETIZZAZIONE (correnti di schermo indotte dal campo) e l'accoppiamento.
%  L'isteresi di TRASPORTO nasce invece dalla rampa della corrente che il cavo
%  porta, ed e' un meccanismo distinto che si somma. Si tiene, ma con due
%  precauzioni: e' definita per CICLO, quindi va divisa per la durata della
%  rampa; e il codice verifica sotto che resti trascurabile, perche' se non lo
%  fosse la sua conversione approssimata da ciclo a rampa peserebbe.
fprintf('\n--- 6.4 Isteresi di trasporto (Norris) ---\n');
F_norris  = min(I_tape_op ./ Ic_local, 0.999);     % clamp: F->1 diverge
Q_tr_tape = (Mu0()*Ic_local.^2/pi) .* ...
            ((1-F_norris).*log(1-F_norris) + (1+F_norris).*log(1+F_norris) - F_norris.^2);
P_tr_nodo = Q_tr_tape * cab.N_tapes_slot / t_ramp;
fprintf('  F = I_tape/Ic: %.3f ... %.3f\n', min(F_norris), max(F_norris));
fprintf('  P_trasporto per strand: max %.6f W/m (%.2f %% dell''isteresi)\n', ...
    max(P_tr_nodo), max(P_tr_nodo)/max(P_mag_nodo)*100);
if max(P_tr_nodo) > 0.05*max(P_mag_nodo)
    warning('step6:norris', ...
        ['Il termine di trasporto vale il %.1f%% dell''isteresi: non e'' piu'' ' ...
         'trascurabile, e la sua conversione da ciclo a rampa (divisione per ' ...
         't_ramp) diventa un''approssimazione che pesa sul totale.'], ...
        max(P_tr_nodo)/max(P_mag_nodo)*100);
end

% ------------------------------------------------------------------------
%  6.5  Accoppiamento inter-strand
% ------------------------------------------------------------------------
%     P_coup = (n * tau / mu0) * (dB/dt)^2 * A ,   tau = (mu0/2)*(l_p/2pi)^2/rho_eff
%
%  n e' il fattore di forma della geometria: n = 2 per campo trasverso, che e'
%  il caso qui (il campo del solenoide e' perpendicolare all'asse del cavo).
%  Con n = 2 la formula coincide con quella classica gia' usata prima: la
%  differenza e' che ora il fattore di forma e' esplicito invece che assorbito
%  nella costante.
%
% ========================================================================
%  IL CONTATTO ELETTRICO STACK-FORMER: UNA TENSIONE DI PROGETTO
% ========================================================================
%  Lo slot e' a contatto col rame del former PER PROTEZIONE AL QUENCH: se il
%  nastro transisce, la corrente deve poter passare nel rame, che ha sezione e
%  capacita' termica molto maggiori, invece di restare confinata nel
%  superconduttore e bruciarlo.
%
%  Quel contatto pero' e' la STESSA interfaccia che determina rho_eff, e le
%  due esigenze tirano in direzioni opposte:
%    - protezione al quench  -> vuole contatto elettrico BUONO (rho_eff bassa)
%    - perdite di accoppiamento -> vogliono contatto CATTIVO (rho_eff alta)
%  Non e' un dettaglio da nascondere in un'assunzione: e' un compromesso di
%  progetto, e va deciso guardando i numeri.
%
%  SCELTA ADOTTATA: si privilegia la protezione al quench, quindi rho_eff e'
%  posta al valore METALLICO, calcolato dai dati del Data Book (rho_Cu a 10 K
%  con magnetoresistenza, per il fattore geometrico del percorso attorno agli
%  slot). E' il valore piu' basso fisicamente possibile e quindi il caso
%  PEGGIORE per le perdite: se il progetto chiude cosi', chiude a maggior
%  ragione con un contatto meno buono.
%
%  CONSEGUENZA, ed e' pesante: con questa rho_eff l'accoppiamento diventa il
%  termine DOMINANTE, alcune volte l'isteresi. La leva per riportarlo sotto
%  controllo e' il PASSO DI TWIST, perche' P_coup varia con l_p^2: il codice
%  calcola sotto quale l_p servirebbe e lo confronta con quello attuale.
%  RAPPORTO CON G = 0 DELLO STEP 4 - non c'e' contraddizione, e vale la pena
%  dirlo bene perche' e' facile confondersi.
%  G = 0 non dipende dall'esistenza di un percorso trasversale: dipende dal
%  RAPPORTO FRA LE RESISTIVITA'. In esercizio nominale il nastro e' nello stato
%  superconduttivo, con resistivita' effettivamente nulla; il percorso
%  alternativo attraverso il rame, per quanto ben collegato, ha resistivita'
%  finita. La corrente di TRASPORTO non ha quindi nessuna ragione di lasciare
%  il superconduttore, e porla tutta nel SC non e' un'approssimazione ma la
%  descrizione corretta del caso nominale, che e' quello che stiamo
%  caratterizzando. Il percorso nel rame conta quando il nastro sviluppa
%  resistenza - current sharing e quench - che e' uno scenario diverso e non
%  oggetto di questo studio.
%
%  E le correnti di ACCOPPIAMENTO, che invece attraversano il rame? Sono un
%  problema distinto: non sono corrente di trasporto imposta ai terminali ma
%  correnti CIRCOLANTI indotte da dB/dt. Il DAE con G = 0 non puo'
%  rappresentarle - non ha maglie trasverse su cui possano chiudersi - ed e'
%  esattamente per questo che l'accoppiamento si calcola qui con una formula
%  analitica invece di leggerlo dal DAE. I due modelli descrivono due modi di
%  corrente diversi e non si sovrappongono.
fprintf('\n--- 6.5 Accoppiamento inter-strand ---\n');
%  RRR = 100: DATO DI PROGETTO fornito dal docente insieme al brief, non una
%  nostra scelta. Sta qui e non in dataAcquisition solo perche' e' arrivato
%  dopo; la sua collocazione naturale sarebbe data.op, col tag 'brief'.
%  Entra in tre punti: la resistivita' residua del rame (rho_Cu_0T), la
%  magnetoresistenza di Kohler (che dipende da B*RRR) e la scelta del set di
%  coefficienti NIST per k(T,RRR). Essendo un dato, i tre usi devono restare
%  agganciati a questa singola variabile: e' gia' cosi'.
acl.RRR       = 100;        % [-] RRR del rame del former [dato del docente]
acl.rho_Cu_0T = 1.55e-10;   % [Ohm*m] Cu a RRR = 100, ~10 K, campo nullo [NIST Monograph 177]
%  MAGNETORESISTENZA: FUNZIONE DEL CAMPO, NON PIU' UN NUMERO SOLO.
%  Regola di Kohler: l'aumento RELATIVO di resistivita' dipende dal solo
%  prodotto B*RRR. Si usa la forma approssimata data in
%    E. Metral et al., "Beam screen issues", arXiv:1108.1643, eq. (5):
%       delta_rho/rho_0 = 10^(-2.69) * (B*RRR)^1.055
%  scelta perche' e' ad accesso aperto ED e' tarata sullo stesso rame che
%  usiamo noi: la nota riporta rho_0 = 1.55e-10 Ohm*m a 20 K con RRR = 100,
%  esattamente acl.rho_Cu_0T. L'esponente 1.055 ~ 1 e' il ramo lineare di
%  alto campo del Kohler plot; a campo molto basso la curva vera e'
%  quadratica, quindi qui delta_rho e' leggermente sovrastimata (rho piu'
%  alta -> perdite piu' basse), ma quei nodi pesano poco perche' P_coup
%  va come B^2.
%  Piu' accurata sarebbe la forma polinomiale di 4o grado in log(B*S(T,RRR))
%  usata nei codici di quench, ma i suoi coefficienti stanno in fonti chiuse.
acl.MR_src    = 'Kohler approx., arXiv:1108.1643 eq.(5)';
acl.MR_of_B   = @(B) 1 + 10^(-2.69) * (abs(B)*acl.RRR).^1.055;
acl.MR_factor = acl.MR_of_B(max(Bmag_near(:)));   % valore al campo di picco
acl.rho_Cu    = acl.rho_Cu_0T * acl.MR_factor;
%  VERIFICA DELLA FORMULA contro i valori calcolati nella fonte stessa
%  (rho_0 = 1.55e-10, RRR = 100): 1.8e-10 a 0.535 T, 5.5e-10 a 8.33 T,
%  11.2e-10 a 20 T. Se la formula fosse trascritta male questo scatta.
mr_chk_B   = [0.535 8.33 20];
mr_chk_ref = [1.8 5.5 11.2]*1e-10;
mr_chk_our = 1.55e-10*(1 + 10^(-2.69)*(mr_chk_B*100).^1.055);
mr_chk_err = max(abs(mr_chk_our - mr_chk_ref)./mr_chk_ref);
if mr_chk_err > 0.05
    error('step6:kohler', ...
        'La Kohler non riproduce i valori della fonte (scarto %.1f%%).', mr_chk_err*100);
end
fprintf('  MR(B) da Kohler [%s]\n', acl.MR_src);
fprintf('    verifica contro i valori della fonte: scarto massimo %.1f%% OK\n', mr_chk_err*100);
%  geomFac CALCOLATO, non assunto: omogeneizzazione 2D della sezione (funzione
%  locale transverseGeomFactor). Vedi la sua intestazione per definizione e
%  validazione. La geometria arriva tutta dallo Step 5.
geoSec = struct('r_form_o', r_former_out, 'r_chan', r_channel, ...
                'r_slot_in', r_slot_in, 'r_slot_out', r_slot_out, ...
                'slot_w', lay.slot_w, 'slot_d', lay.slot_d, ...
                'R_core', R_core_derived, 'N_slots', cab.N_slots);
tSec = tic;
[acl.geomFac, gfInfo] = transverseGeomFactor(geoSec, 201);
fprintf('  geomFac = %.3f CALCOLATO sulla sezione (griglia %dx%d, %.1f s)\n', ...
    acl.geomFac, gfInfo.n, gfInfo.n, toc(tSec));
fprintf('    frazione di rame nel disco del former %.2f -> il resto e'' allungamento\n', ...
    gfInfo.area_frac);
fprintf('    del percorso attorno agli slot e al canale\n');
%  DIAGNOSTICA RICHIESTA: resistenza fra due strand ADIACENTI. Non e' la
%  grandezza che entra in rho_eff (li' serve una resistivita' omogeneizzata,
%  non una resistenza fra due elettrodi) ma dice quanto il rame fra due slot
%  vicini sia un percorso quasi diretto.
gf_pair = transverseGeomFactor(geoSec, 201, 'pair');
fprintf('    [diagnostica] percorso fra due strand ADIACENTI: fattore %.2f\n', gf_pair);
fprintf('    fra slot vicini il rame passa quasi dritto; il valore usato in\n');
fprintf('    rho_eff resta quello omogeneizzato\n');
if gfInfo.ref_err > 0.01
    warning('step6:geomFacRef', ...
        ['La validazione su disco di rame pieno da'' %.4f invece di 1: la griglia ' ...
         'sta introducendo un bias, infittire n.'], gfInfo.sigma_ref);
else
    fprintf('    validazione su disco di rame pieno: %.4f (atteso 1) OK\n', gfInfo.sigma_ref);
end
acl.rho_eff_metallic = acl.rho_Cu * acl.geomFac;   % percorso tutto in rame
acl.n_shape   = 2.0;        % [-] fattore di forma, campo trasverso

acl.rho_eff      = acl.rho_eff_metallic;   % contatto elettrico voluto (quench protection)
acl.rho_eff_band = [acl.rho_eff_metallic, 1e-6];   % da contatto metallico a resistivo

lp_inner  = 2*pi*R_in/cab.N_twist;    % passo di twist sul turn interno [m]
%  IL PASSO DI TWIST NON E' COSTANTE. twistAngle (Step 2.4) impone un numero
%  FISSO di giri di twist per turn a QUALUNQUE raggio: alpha = theta_slot +
%  N_twist*theta. Segue che il passo vero e' l_p(r) = 2*pi*r/N_twist, non il
%  valore calcolato al turn interno e riusato ovunque. Con P_coup
%  proporzionale a l_p^2, usare lp_inner al turn esterno (dove r e' maggiore,
%  circa 1.35 volte su questa geometria: R_out/R_in = (R_in + cab.width*
%  N_turn)/R_in) SOTTOSTIMA l'accoppiamento di circa il quadrato di quel
%  rapporto, cioe' un fattore ~1.8: e' il limite gia' segnalato in precedenza
%  come questione aperta, qui chiuso. Il valore esatto per QUESTA geometria
%  (non un numero fisso) e' stampato a runtime subito sotto.
%  r_nodes e' il raggio del nodo sullo STRAND (include il piccolo offset del
%  twist), non della linea centrale del cavo: la differenza e' dell'ordine
%  di R_core, cioe' millimetrica contro raggi di centinaia di millimetri, e
%  la si trascura dichiaratamente invece di ricostruire una seconda
%  geometria solo per questo termine.
lp_nodo = 2*pi*r_nodes/cab.N_twist;   % [m] passo di twist locale, per nodo
fprintf('  l_p locale: %.0f mm (turn interno) ... %.0f mm (turn esterno), contro\n', ...
    min(lp_nodo)*1e3, max(lp_nodo)*1e3);
%  rho_eff e tau sono LOCALI: il campo cambia di un fattore 15 lungo il
%  conduttore e con esso la magnetoresistenza; ora anche l_p e' locale. tau_c
%  resta stampato al picco per il controllo di quasi-stazionarieta'.
acl.rho_eff_nodo = acl.rho_Cu_0T * acl.MR_of_B(Bmag_near) * acl.geomFac;
tau_c_nodo = (Mu0()/2)*(lp_nodo/(2*pi)).^2 ./ acl.rho_eff_nodo;
tau_c     = max(tau_c_nodo(:));
dBdt_loc  = Bmag_near / t_ramp;       % [T/s] locale (B proporzionale a I)
%  Nella formula di accoppiamento si usa A_stack, cioe' la sezione conduttiva
%  associata allo stack per slot, non l'intero ingombro geometrico del CICC.
P_coup_nodo = (acl.n_shape*tau_c_nodo/Mu0()) .* dBdt_loc.^2 * A_stack;
fprintf('  effetto di l_p(s) sull''accoppiamento, al turn esterno: fattore %.2f\n', ...
    (max(lp_nodo)/lp_inner)^2);

fprintf('  rho_eff = %.2e Ohm*m (contatto METALLICO col former, per quench protection)\n', ...
    acl.rho_eff);
fprintf('    = rho_Cu(RRR %d) x MR(B) x %.3f (percorso attorno agli slot)\n', ...
    acl.RRR, acl.geomFac);
fprintf('    MR(B) da %.2f al turn esterno a %.2f al picco: NON e'' piu'' un numero\n', ...
    acl.MR_of_B(min(Bmag_near(:))), acl.MR_factor);
fprintf('    unico, e'' valutata sul campo di ogni nodo (Kohler, ramo lineare)\n');
fprintf('    -> rho_eff da %.2e a %.2e Ohm*m lungo il conduttore\n', ...
    min(acl.rho_eff_nodo(:)), max(acl.rho_eff_nodo(:)));
fprintf('  l_p = %.3f m | n = %.0f -> tau = %.3f s (MASSIMO lungo il conduttore,\n', ...
    lp_inner, acl.n_shape, tau_c);
fprintf('    che cade dove rho_eff e'' minima, cioe'' al campo MINIMO; al picco di\n');
fprintf('    campo tau vale %.3f s. Il quasi-stazionario si verifica sul massimo.\n', ...
    min(tau_c_nodo(:)));
fprintf('  P_accoppiamento per strand: max %.4f W/m\n', max(P_coup_nodo));

%  VALIDITA' DEL MODELLO QUASI-STAZIONARIO. La formula presuppone che le
%  correnti di accoppiamento siano in equilibrio con dB/dt, cioe' tau << durata
%  del transitorio. Con contatto metallico tau cresce di tre ordini di
%  grandezza rispetto a un contatto resistivo, quindi la verifica non e' piu'
%  automatica e va fatta.
fprintf('  tau/t_rampa = %.4f  [%s]\n', tau_c/t_ramp, ...
    ternary(tau_c < 0.1*t_ramp, 'quasi-stazionario OK', 'AL LIMITE'));
if tau_c > 0.1*t_ramp
    warning('step6:tauSlow', ...
        ['tau = %.2f s contro una rampa di %.0f s: le correnti di accoppiamento ' ...
         'non seguono piu'' quasi-staticamente il campo e la formula sovrastima. ' ...
         'Servirebbe un modello transitorio del circuito di accoppiamento.'], ...
        tau_c, t_ramp);
end

%  SENSIBILITA' e LEVA DI PROGETTO
fprintf('  sensibilita'' a rho_eff (per strand, sul turn a campo massimo):\n');
for rr = [acl.rho_eff_band(1), 1e-8, 1e-7, acl.rho_eff_band(2)]
    tt = (Mu0()/2)*(lp_inner/(2*pi))^2 / rr;
    pp = (acl.n_shape*tt/Mu0()) * max(dBdt_loc)^2 * A_stack;
    fprintf('    rho_eff = %.1e -> tau = %8.3f s -> P = %.4f W/m\n', rr, tt, pp);
end
P_coup_max_strand = max(P_coup_nodo);   % serve al controllo in 6.7

% ------------------------------------------------------------------------
%  6.5c  ACCOPPIAMENTO: DIPENDENZA TEMPORALE Q(s,t)
% ------------------------------------------------------------------------
%  Qui l'andamento nel tempo e' l'OPPOSTO di quello dell'isteresi, e vale la
%  pena dirlo esplicitamente perche' non e' intuitivo. La magnetoresistenza
%  MR(B) CRESCE con B, quindi rho_eff e' MINIMA a inizio rampa (B basso):
%  tau = (mu0/2)(l_p/2pi)^2/rho_eff e' quindi MASSIMO a inizio rampa, non a
%  fine. Le perdite di accoppiamento sono piu' alte all'INIZIO della rampa,
%  quando il conduttore e' piu' "pulito" elettricamente, non alla fine.
%  acl.MR_of_B e' gia' vettorizzata per costruzione (solo potenze e prodotti
%  elementwise), quindi il prodotto esterno lambda*B non richiede arrayfun.
Bmat_coup     = Bmag_near * lam_row;                              % [nodi x n_lam]
rho_eff_nodo_t = acl.rho_Cu_0T * acl.MR_of_B(Bmat_coup) * acl.geomFac;
tau_c_nodo_t   = (Mu0()/2)*(lp_nodo/(2*pi)).^2 ./ rho_eff_nodo_t;  % lp_nodo gia' locale
%  dB/dt e' costante nel tempo (rampa lineare): P_coup_nodo_t varia solo
%  attraverso tau_c_nodo_t.
P_coup_nodo_t = (acl.n_shape*tau_c_nodo_t/Mu0()) .* dBdt_loc.^2 * A_stack;

%  trapz, non mean: stessa ragione della 6.3c.
fatt_coup_pk = trapz(cool.lam_grid, P_coup_nodo_t(idxPk,:)) / P_coup_nodo_t(idxPk,end);
fprintf('  fattore temporale accoppiamento (media rampa / campo finale), nodo di picco: %.2f\n', ...
    fatt_coup_pk);
if fatt_coup_pk < 1.2
    warning('step6:coupTimeFlat', ...
        ['Il fattore temporale sull''accoppiamento (%.2f) e'' insolitamente vicino a ' ...
         '1: verificare Bmat_coup e acl.MR_of_B su input matrice.'], fatt_coup_pk);
end

%  VERIFICA QUASI-STAZIONARIA VERA. La verifica in 6.5 confrontava tau_c/
%  t_ramp usando tau AL CAMPO FINALE. Ma tau e' massimo dove rho_eff e'
%  minima, cioe' a CAMPO BASSO: se il transitorio delle correnti di
%  accoppiamento non e' trascurabile rispetto a t_ramp proprio nei primi
%  istanti (dove la 6.5c mostra che la potenza e' piu' alta, non piu' bassa),
%  l'ipotesi quasi-stazionaria su cui poggia P_coup = f(tau) va verificata
%  sul massimo VERO, non su quello al campo finale.
tau_c_max_t = max(tau_c_nodo_t(:));
fprintf('  tau massimo su TUTTA la rampa = %.3f s\n', tau_c_max_t);
fprintf('  tau_max/t_rampa = %.4f  [%s]\n', tau_c_max_t/t_ramp, ...
    ternary(tau_c_max_t < 0.1*t_ramp, 'quasi-stazionario OK su tutta la rampa', ...
                                       'AL LIMITE: la verifica al campo finale non bastava'));
if tau_c_max_t > 0.1*t_ramp
    warning('step6:tauSlowRamp', ...
        ['tau massimo sulla rampa (%.2f s, a campo basso) supera il 10%% di t_ramp ' ...
         '(%.0f s): le correnti di accoppiamento non seguono quasi-staticamente il ' ...
         'campo nella parte iniziale, dove pero'' la potenza di accoppiamento e'' ' ...
         'maggiore. Il fattore temporale stampato sopra e'' un limite superiore, non ' ...
         'un valore su cui contare.'], tau_c_max_t, t_ramp);
end

%  Fattore sull'ENERGIA della rampa al SOLO nodo di picco: e' una diagnostica
%  LOCALE, utile a vedere dove la correzione pesa di piu', ma NON il numero
%  da citare in relazione. Il fattore rappresentativo dell'intero DP,
%  integrato anche in s oltre che in t, e' stampato in 7.3 (usa q_th_t,
%  costruita a valle di P_ac_cable_t qui sopra): e' quello il numero
%  fisicamente corretto, perche' pesa ogni tratto di conduttore per la
%  potenza che davvero deposita, non solo il punto piu' caldo.
Etot_finale_pk = (P_mag_nodo(idxPk) + P_coup_nodo(idxPk)) * t_ramp;
Etot_vera_pk   = trapz(cool.lam_grid, P_mag_nodo_t(idxPk,:) + P_coup_nodo_t(idxPk,:)) * t_ramp;
fprintf('  fattore sull''energia AC al SOLO nodo di picco (diagnostica locale): %.2f\n', ...
    Etot_vera_pk/Etot_finale_pk);
fprintf('      (energia con Q(s) congelato al campo finale: %.2f kJ; con Q(s,t): %.2f kJ)\n', ...
    Etot_finale_pk/1e3, Etot_vera_pk/1e3);
fprintf('      il fattore sull''INTERO DP, quello da citare, e'' stampato in 7.3\n');

% ------------------------------------------------------------------------
%  6.6  Eddy nei metalli (former Cu, jacket SS)
% ------------------------------------------------------------------------
%   P/V = (dB/dt)^2 * w^2 / (12*rho)   (lastra piana)   [Wilson 1983]
%  Il rho del rame a 15 T include la MAGNETORESISTENZA: usare il valore a campo
%  nullo sovrastimerebbe le perdite (rho piu' basso -> piu' eddy).
fprintf('\n--- 6.6 Eddy nei metalli ---\n');
% acl.rho_Cu_0T, acl.MR_factor e acl.rho_Cu sono definite in 6.5, dove servono
% gia' per il limite inferiore di rho_eff.
acl.rho_SS    = 5.0e-7;     % [Ohm*m] acciaio inox a bassa T [NIST]
w_Cu = 2*(r_former_out - r_channel);
w_SS = lay.t_jacket;
P_Cu_nodo = dBdt_loc.^2 * w_Cu^2 / (12*acl.rho_Cu) * A_Cu / cab.N_slots;
P_SS_nodo = dBdt_loc.^2 * w_SS^2 / (12*acl.rho_SS) * A_SS / cab.N_slots;
fprintf('  rho_Cu(15T) = %.2e Ohm*m (x%.1f magnetoresistenza) | rho_SS = %.1e\n', ...
    acl.rho_Cu, acl.MR_factor, acl.rho_SS);
fprintf('  P_eddy Cu per strand: max %.6f W/m | SS: max %.6f W/m\n', ...
    max(P_Cu_nodo), max(P_SS_nodo));
if max(P_SS_nodo) < 0.01*max(P_mag_nodo)
    fprintf('  -> eddy SS < 1%% del termine dominante: trascurabile (VERIFICATO, non assunto)\n');
end

% ------------------------------------------------------------------------
%  6.7  Totale Q(s): sorgente per lo Step 7
% ------------------------------------------------------------------------
%  IL TOTALE USA SOLO ISTERESI E ACCOPPIAMENTO. Gli altri tre termini sono
%  calcolati sopra e si sono rivelati irrilevanti su questa macchina:
%  trasporto (Norris) e eddy nell'acciaio sotto lo 0.1%, eddy nel rame 0.1%.
%  Escluderli dal totale toglie rumore senza togliere informazione, ed e' la
%  stessa struttura delle due formule adottate: isteresi nel superconduttore
%  piu' accoppiamento fra gli stack.
%  Restano pero' CALCOLATI e VERIFICATI, non cancellati: la verifica sotto
%  scatta a ogni run e avvisa se uno di loro tornasse a contare. Serve perche'
%  i tre termini scalano diversamente dagli altri due - le eddy vanno con
%  (dB/dt)^2 e con la resistivita' dei metalli, il trasporto con I/Ic - quindi
%  la loro irrilevanza vale per QUESTI parametri, non per costruzione.
P_ac_nodo = P_mag_nodo + P_coup_nodo;                     % [W/m] per strand
P_minor   = P_tr_nodo + P_Cu_nodo + P_SS_nodo;            % termini esclusi

% Somma sugli N_slots strand che condividono la stessa ascissa s -> W/m di CAVO
nPerStrand = size(P_micro,1)/cab.N_slots;
P_ac_cable = zeros(nPerStrand,1);
for j = 1:cab.N_slots
    P_ac_cable = P_ac_cable + P_ac_nodo((j-1)*nPerStrand + (1:nPerStrand));
end
s_cable  = s_center;
P_ac_tot = trapz(s_cable, P_ac_cable);          % [W] potenza totale sul DP

%  VERSIONE TEMPO-DIPENDENTE, stessa somma sugli strand ma sulla matrice
%  Q(s,t) della 6.3c/6.5c: serve alla 7.3 per costruire q_th_t, il confronto
%  fra l'energia integrata sulla rampa e quella al campo finale.
P_ac_nodo_t  = P_mag_nodo_t + P_coup_nodo_t;               % [nodi x n_lam]

%  Somma separata dei due contributi sui N_slots strand -> [W/m] di cavo.
P_mag_cable_t  = zeros(nPerStrand, cool.n_lam);
P_coup_cable_t = zeros(nPerStrand, cool.n_lam);
P_ac_cable_t   = zeros(nPerStrand, cool.n_lam);
for j = 1:cab.N_slots
    idxj = (j-1)*nPerStrand + (1:nPerStrand);
    P_mag_cable_t  = P_mag_cable_t  + P_mag_nodo_t(idxj,:);
    P_coup_cable_t = P_coup_cable_t + P_coup_nodo_t(idxj,:);
    P_ac_cable_t   = P_ac_cable_t   + P_ac_nodo_t(idxj,:);
end

P_mag_tot_t  = trapz(s_cable, P_mag_cable_t , 1);   % [W]
P_coup_tot_t = trapz(s_cable, P_coup_cable_t, 1);   % [W]
P_ac_tot_t   = trapz(s_cable, P_ac_cable_t  , 1);   % [W] 1 x n_lam, potenza istante per istante
P_ac_avg     = trapz(cool.lam_grid, P_ac_tot_t);    % [W] potenza VERA, mediata sulla rampa

%  Per il grafico temporale si usa la potenza lineare media sul DP [W/m].
L_cable = s_cable(end) - s_cable(1);
if L_cable <= 0
    L_cable = 1;
end
q_mag_time  = P_mag_tot_t  / L_cable;               % [W/m]
q_coup_time = P_coup_tot_t / L_cable;               % [W/m]
t_loss      = cool.lam_grid * t_ramp;               % [s]
t_loss_pts  = 0:10:t_ramp;
if isempty(t_loss_pts) || t_loss_pts(end) < t_ramp
    t_loss_pts = [t_loss_pts t_ramp];
end
q_mag_pts  = interp1(t_loss, q_mag_time , t_loss_pts, 'linear');
q_coup_pts = interp1(t_loss, q_coup_time, t_loss_pts, 'linear');
%  Controllo di coerenza: l'ultima colonna di P_ac_tot_t (campo finale) deve
%  ricadere su P_ac_tot entro l'errore di interpolazione della mesh di campo.
errPtot = abs(P_ac_tot_t(end) - P_ac_tot)/max(P_ac_tot, eps);
if errPtot > 1e-6
    warning('step6:ptotMismatch', ...
        ['P_ac_tot_t al campo finale (%.3f W) non coincide con P_ac_tot (%.3f W), ' ...
         'scarto %.2e: la matrice tempo-dipendente non sta usando la stessa mesh ' ...
         'di P_ac_cable.'], P_ac_tot_t(end), P_ac_tot, errPtot);
end

%  CONTROLLO SULL'ACCOPPIAMENTO, l'analogo di quello che la 6.6 fa per l'eddy
%  nell'acciaio: verificare invece di assumere. Se il termine di accoppiamento
%  supera il 10% del totale, la specifica su rho_eff non e' piu' un dettaglio
%  e va portata nei requisiti di fabbricazione con una tolleranza.
%  La ripartizione va letta NELLO STESSO nodo, altrimenti si sommano i massimi
%  di termini che piccano in punti diversi e il totale supera il 100%.
[~, iPk]  = max(P_ac_cable);
frac_coup = P_coup_nodo(iPk)/(P_mag_nodo(iPk) + P_coup_nodo(iPk));
frac_hyst = 1 - frac_coup;
fprintf('  ripartizione NEL NODO di picco (s = %.1f m): isteresi %.1f %%, accoppiamento %.1f %%\n', ...
    s_cable(iPk), frac_hyst*100, frac_coup*100);

%  LEVA DI PROGETTO SUL PASSO DI TWIST. P_coup varia con l_p^2, quindi se
%  l'accoppiamento e' troppo alto la via non e' rinunciare al contatto col rame
%  (serve alla protezione al quench) ma TWISTARE PIU' STRETTO. Si calcola il
%  passo che riporterebbe l'accoppiamento sotto una frazione data dell'isteresi
%  e il numero di giri per turn corrispondente.
fprintf('\n  --- Come gestire l''accoppiamento: due leve, in ordine di costo ---\n');
P_hyst_pk = max(P_mag_nodo);
fprintf('  LEVA 1, LA PORTATA. L''accoppiamento e'' potenza da smaltire, e il budget\n');
fprintf('  di elio ha margine: lo Step 7 cerca la portata minima e la confronta con\n');
fprintf('  i %.0f g/s disponibili. Se ci sta, non serve toccare il cavo, ed e'' la\n', ...
    data.op.mdot_max.val*1e3);
fprintf('  soluzione piu'' economica: nessuna modifica di geometria, nessun impatto\n');
fprintf('  sulla protezione al quench, nessuna cascata sulla mesh del DAE.\n');
fprintf('  Il verdetto lo da'' lo Step 7.5, non questa sezione.\n');
fprintf('\n  LEVA 2, IL PASSO DI TWIST (P_coup varia con l_p^2). Da usare solo se la\n');
fprintf('  portata non basta: costa una revisione di disc.N_seg_dae e del tempo di\n');
fprintf('  calcolo, e un twist piu stretto aumenta la deformazione dei nastri.\n');
fprintf('  attuale: l_p = %.0f mm, N_twist = %d -> accoppiamento %.0f %% dell''isteresi\n', ...
    lp_inner*1e3, cab.N_twist, P_coup_max_strand/P_hyst_pk*100);
for fr = [0.10 0.20 0.30]
    lp_req = lp_inner*sqrt(fr*P_hyst_pk/max(P_coup_max_strand,eps));
    fprintf('    per accoppiamento <= %2.0f %% dell''isteresi: l_p <= %3.0f mm -> N_twist >= %2.0f\n', ...
        fr*100, lp_req*1e3, ceil(2*pi*R_in/lp_req));
end
fprintf('  (riferimento ENEA dalle slide: passo di twist 0.5 m)\n');
fprintf('  ATTENZIONE: alzare N_twist richiede di rivedere disc.N_seg_dae, perche''\n');
fprintf('  la 2.1c chiede almeno 4 segmenti per periodo di twist (ora %.1f).\n', segPerTwist_dae);

fprintf('\n  NOTA: che l''accoppiamento domini (%.0f %% del totale) NON e'' di per se''\n', frac_coup*100);
fprintf('  un problema. E'' la conseguenza attesa del contatto metallico voluto per\n');
fprintf('  la protezione al quench, ed e'' un problema solo se la potenza risultante\n');
fprintf('  non e'' smaltibile. Il criterio di accettazione e'' quindi la portata\n');
fprintf('  minima dello Step 7 contro i %.0f g/s disponibili, non la quota di\n', ...
    data.op.mdot_max.val*1e3);
fprintf('  accoppiamento sul totale. Il verdetto e'' rimandato li''.\n');
%  L'energia AC viene integrata sulla potenza tempo-dipendente della rampa,
%  non stimata a partire dal solo valore al campo finale.
E_ac_tot = P_ac_avg * t_ramp;                   % [J] energia VERA sulla rampa
fprintf('  potenza media sulla rampa (Q(s,t)) = %.1f W, contro %.1f W al campo finale\n', ...
    P_ac_avg, P_ac_tot);
fprintf('    fattore %.2f | energia vera sulla rampa = %.1f kJ\n', ...
    P_ac_avg/max(P_ac_tot,eps), E_ac_tot/1e3);

fprintf('\n--- 6.7 Totale ---\n');
tot = max(P_mag_nodo)+max(P_tr_nodo)+max(P_coup_nodo)+max(P_Cu_nodo)+max(P_SS_nodo);
fprintf('  %-22s %12s %11s\n','termine','max [W/m]','quota');
fprintf('  %-22s %12.6f %10.1f%%\n','magnetizzazione', max(P_mag_nodo), max(P_mag_nodo)/tot*100);
fprintf('  %-22s %12.6f %10.1f%%\n','trasporto (Norris)', max(P_tr_nodo), max(P_tr_nodo)/tot*100);
fprintf('  %-22s %12.6f %10.1f%%\n','accoppiamento', max(P_coup_nodo), max(P_coup_nodo)/tot*100);
fprintf('  %-22s %12.6f %10.1f%%\n','eddy Cu', max(P_Cu_nodo), max(P_Cu_nodo)/tot*100);
fprintf('  %-22s %12.6f %10.1f%%\n','eddy SS', max(P_SS_nodo), max(P_SS_nodo)/tot*100);
fprintf('  %s\n', repmat('-',1,48));
fprintf('  %-22s %12.6f %10.1f%%   <- ESCLUSI dal totale\n', ...
    'somma dei minori', max(P_minor), max(P_minor)/tot*100);
if max(P_minor) > 0.01*max(P_ac_nodo)
    warning('step6:minorTerms', ...
        ['I termini esclusi (trasporto + eddy) valgono il %.1f%% del totale: non ' ...
         'sono piu'' trascurabili e vanno rimessi in P_ac_nodo. Causa tipica: e'' ' ...
         'cambiato dB/dt, la resistivita'' dei metalli, o il rapporto I/Ic.'], ...
        max(P_minor)/max(P_ac_nodo)*100);
else
    fprintf('  -> minori sotto l''1%%: esclusione VERIFICATA a questo run, non assunta\n');
end
%  La media va pesata sulla LUNGHEZZA: i segmenti hanno passi diversi al variare
%  del raggio, quindi mean() sui nodi non e' l'integrale diviso L.
fprintf('  Q(s) di CAVO: max %.3f W/m, medio %.3f W/m (= P_tot/L)\n', ...
    max(P_ac_cable), P_ac_tot/(s_cable(end)-s_cable(1)));
fprintf('  potenza totale sul DP AL CAMPO FINALE = %.1f W (vedi sopra per la media\n', ...
    P_ac_tot);
fprintf('    vera sulla rampa e per l''energia, che usano P_ac_avg)\n');

% --- Bilancio termico preliminare: e' raffreddabile? ----------------------
%  Stima grossolana (lo Step 7 la fara' col modello vero e le proprieta' reali).
cp_He_approx = 5000;                 % [J/kg/K] He supercritico ~7 K, 5 bar [NIST]
%  P_ac_avg, non P_ac_tot: e' la potenza che la portata deve smaltire IN
%  MEDIA sulla rampa, coerente con quanto ora fa lo Step 7.
dT_He = P_ac_avg / (data.op.mdot_max.val * cp_He_approx);
fprintf('\n  --- Bilancio termico preliminare (indicativo) ---\n');
fprintf('  riscaldamento He stimato = P/(mdot*cp) = %.2f K\n', dT_He);
fprintf('  T_in min (brief) = %.1f K -> T uscita stimata = %.2f K (T_op ammessa %.2f K)\n', ...
    data.op.Tin_range.val(1), data.op.Tin_range.val(1)+dT_He, T_op_guess);
if data.op.Tin_range.val(1)+dT_He > T_op_guess
    warning('step6:coolingTight', ...
        ['Il riscaldamento stimato dell''elio (%.2f K) porterebbe l''uscita a %.2f K, ' ...
         'SOPRA la T_op ammessa (%.2f K). Lo Step 7 dovra'' verificarlo col modello ' ...
         'vero. Rimedi possibili: piu'' nastri (piu'' margine), percorso idraulico ' ...
         'piu'' corto, o revisione dello scenario di rampa.'], ...
        dT_He, data.op.Tin_range.val(1)+dT_He, T_op_guess);
end

% ------------------------------------------------------------------------
%  6.8  RIEPILOGO DELLE ASSUNZIONI DELLO STEP 6
% ------------------------------------------------------------------------
fprintf('\n--- 6.8 Assunzioni dello Step 6 ---\n');
fprintf('  A6.1 Scenario CHIUSO: una sola rampa monotona 0 -> %.0f kA da stato vergine,\n', I_nom/1e3);
fprintf('       nessun ciclo. Le formule adottate danno direttamente la POTENZA in\n');
fprintf('       funzione di dB/dt, quindi non serve nessuna conversione da energia\n');
fprintf('       per ciclo (vedi 6.1).\n');
fprintf('  A6.2 Formule chiuse e non DAE: con G = 0 il DAE non puo'' dare accoppiamento\n');
fprintf('       ne'' isteresi intra-nastro (arco lumped, nessun profilo di penetrazione).\n');
fprintf('  A6.3 Stato critico (Bean) = limite n->infinito della power law: le formule di\n');
fprintf('       isteresi NON usano n. n(T,B) resta usato solo nella R(i,B,T) dello Step 4.\n');
fprintf('  A6.4 Nastri dello stack trattati come INDIPENDENTI (perdita singola x %d).\n', cab.N_tapes_slot);
fprintf('       In uno stack reale la schermatura reciproca RIDUCE la magnetizzazione:\n');
fprintf('       questa stima e'' quindi CONSERVATIVA.\n');
fprintf('  A6.5 Solo B_perp genera magnetizzazione; B_par trascurata (offre al flusso lo\n');
fprintf('       spessore del film, %.1f um, non la larghezza).\n', data.tape.t_REBCO.val*1e6);
fprintf('  A6.6 rho_eff = %.1e Ohm*m, valore METALLICO: conseguenza della scelta di\n', acl.rho_eff);
fprintf('       mettere lo stack a contatto col rame per la protezione al quench.\n');
fprintf('       Non e'' un''assunzione libera ma il caso peggiore per le perdite,\n');
fprintf('       calcolato dai dati del Data Book. Il fattore geometrico (x%.3f) e''\n', acl.geomFac);
fprintf('       CALCOLATO sulla sezione (6.5) e la magnetoresistenza e'' MR(B) da\n');
fprintf('       Kohler, quindi nessuno dei due e'' piu'' assunto.\n');
fprintf('  A6.7 n = %.0f, fattore di forma per campo trasverso. Con n = 2 la formula\n', acl.n_shape);
fprintf('       coincide con quella classica dell''accoppiamento.\n');
fprintf('  A6.8 d_eff = %.1f mm = larghezza del nastro (campo perpendicolare). Con\n', acl.d_eff*1e3);
fprintf('       campo parallelo si userebbe lo spessore e il risultato cadrebbe di\n');
fprintf('       tre ordini di grandezza: la scelta va dichiarata, non sottintesa.\n');
fprintf('  A6.7b Magnetoresistenza Cu: MR(B) da Kohler [arXiv:1108.1643 eq.(5)],\n');
fprintf('       x%.2f al picco di campo, x%.2f al turn esterno.\n', ...
    acl.MR_factor, acl.MR_of_B(min(Bmag_near(:))));
fprintf('  A6.8 Equiripartizione I/%d per slot: dimostrata nello Step 4 (twist commensurato).\n', cab.N_slots);
fprintf('  A6.9 Asintoto di penetrazione completa verificato in 6.2 (rapporto %.0f >> 10).\n', Bpk/Bc_strip);

tmr = addTimer(tmr, 'Step 6    perdite AC (Ic locale su ogni nodo)', toc(t0));

% ------------------------------------------------------------------------
%  6.9  Grafici
% ------------------------------------------------------------------------
figure('Name','STEP 6 - AC losses','Color','w','Position',[160 160 1300 420]);

subplot(1,3,1); hold on; grid on; box on;
plot(s_nodes, max(P_mag_nodo,1e-12),'.','MarkerSize',3);
plot(s_nodes, max(P_coup_nodo,1e-12),'.','MarkerSize',3);
plot(s_nodes, max(P_Cu_nodo,1e-12),'.','MarkerSize',3);
set(gca,'YScale','log');
xlabel('s [m]'); ylabel('P per strand [W/m]');
title('Termini di perdita');
legend({'magnetizzazione','accoppiamento','eddy Cu'},'Location','best','FontSize',7);

subplot(1,3,2); hold on; grid on; box on;
plot(s_cable, P_ac_cable,'LineWidth',1.4);
xlabel('s [m]'); ylabel('Q [W/m] di cavo');
title(sprintf('Sorgente per lo Step 7 (max %.2f W/m)', max(P_ac_cable)));

subplot(1,3,3); hold on; grid on; box on;
hMag = plot(t_loss, q_mag_time, '-', 'LineWidth',1.6);
scatter(t_loss_pts, q_mag_pts, 28, 'o', ...
    'MarkerEdgeColor', hMag.Color, 'MarkerFaceColor', hMag.Color, ...
    'DisplayName','isteresi: punti ogni 10 s');
hCoup = plot(t_loss, q_coup_time, '-', 'LineWidth',1.6);
scatter(t_loss_pts, q_coup_pts, 28, 's', ...
    'MarkerEdgeColor', hCoup.Color, 'MarkerFaceColor', hCoup.Color, ...
    'DisplayName','accoppiamento: punti ogni 10 s');
xlabel('t [s]');
ylabel('potenza media lineare sul DP [W/m]');
title('AC losses nel tempo');
legend({'isteresi','isteresi: punti ogni 10 s', ...
        'accoppiamento','accoppiamento: punti ogni 10 s'}, ...
       'Location','best','FontSize',8);

save('step6_aclosses.mat','P_ac_cable','s_cable','P_ac_nodo','P_mag_nodo','P_tr_nodo', ...
     'P_coup_nodo','P_Cu_nodo','P_SS_nodo','P_ac_tot','E_ac_tot','acl','tau_c','Bc_strip', ...
     'P_ac_cable_t','P_ac_nodo_t','P_mag_nodo_t','P_coup_nodo_t','P_ac_tot_t','P_ac_avg', ...
     'tau_c_nodo_t','tau_c_max_t','cool');
fprintf('\n--> Salvato step6_aclosses.mat (Q(s) e Q(s,t) sorgente per lo Step 7)\n');

%% ========================================================================
%  STEP 7 - COOLING NEEDS: GANDALF ridotto, slow-transient 1D
%  ========================================================================
%  OBIETTIVI:
%    1) trovare la portata minima di He che garantisce
%         min_{s,t}[T_cs(s)-T_sc(s,t)] >= 5 K
%       a T_in = 4.5 K;
%    2) dimostrare la convergenza numerica DEL RISULTATO m_dot,min;
%    3) costruire la frontiera operativa T_in,max(m_dot), fino a 15 g/s.
%
%  MODELLO FISICO DI RIFERIMENTO.
%  Il framework resta quello GANDALF 1-region per CICC. Si adotta pero' la
%  riduzione slow-transient: l'energia di He e solidi e' transitoria, mentre
%  pressione e velocita' sono trattate quasi-stazionariamente.
%
%  TEMPO ACUSTICO.
%      t_ac = L_ramo/c
%  e' il tempo impiegato da una perturbazione di pressione a percorrere un
%  ramo alla velocita' del suono. Se t_ac << t_ramp, la dinamica acustica si
%  assesta molto prima dell'evoluzione termica. Questo giustifica la riduzione
%  di p e v a variabili quasi-stazionarie. Il tempo di transito del fluido,
%  L_ramo/v, resta invece nel problema energetico attraverso l'avvezione.
%
%  INCOGNITE TERMICHE TRANSITORIE:
%    T_sc(s,t)  stack HTS
%    T_Cu(s,t)  former di rame
%    T_j(s,t)   jacket SS304
%    T_He(s,t)  elio
%
%  CONTATTO STACK-FORMER:
%    R_contact,th = 0. Rimangono finite le resistenze di conduzione interne
%    allo stack e al former, quindi T_sc e T_Cu restano incognite distinte.
%
%  PROPRIETA':
%    He: tabella real-fluid embedded nel file;
%    Cu e SS304: correlazioni NIST gia' presenti;
%    Hastelloy C-276: cpHast(T) e kHast(T) gia' presenti nel main.
fprintf('\n\n########## STEP 7 - COOLING NEEDS (GANDALF RIDOTTO / SLOW TRANSIENT) ##########\n');
tStep7 = tic;

% ------------------------------------------------------------------------
%  7.1  Geometria, condizioni operative e parametri numerici
% ------------------------------------------------------------------------
cool.A_He     = A_channel;
cool.P_wet    = P_He_channel;
cool.D_h      = 4*cool.A_He/cool.P_wet;
cool.L        = s_center(end);
cool.feed     = 'center';
cool.T_in     = data.op.Tin_range.val(1);       % [K] punto base del brief
cool.p_out    = 10e5;                           % [Pa] riferimento, centro banda
cool.mdot_max = data.op.mdot_max.val;           % [kg/s] massimo totale sul DP
cool.RRR      = acl.RRR;
cool.MR_Cu    = acl.MR_factor;
cool.t_pre    = 0;
cool.t_end    = t_ramp + 100;                   % [s] include la coda post-rampa

if isfield(data.op,'pin_range')
    cool.p_range = toPascal(data.op.pin_range.val);
else
    cool.p_range = [5 20]*1e5;
end
if cool.p_out < cool.p_range(1) || cool.p_out > cool.p_range(2)
    cool.p_out = mean(cool.p_range);
end

% Tolleranze delle ricerche scalari.
cool.flow_tol      = 0.001e-3;   % [kg/s] 0.001 g/s sulla portata
cool.Tin_tol       = 0.02;       % [K] sulla frontiera T_in,max
cool.max_bisect    = 14;
cool.mdot_floor    = 0.25e-3;    % [kg/s] solo limite numerico inferiore

% Proprieta' equivalenti dello stack.
lyr = data.tape.layers.val;
ttot = data.tape.t_total.val;
cool.f_Cu_tape = (lyr.Cu_top + lyr.Cu_back)/ttot;
cool.f_Ag_tape = (lyr.Ag_top + lyr.Ag_back)/ttot;
cool.f_base_tape = max(1 - cool.f_Cu_tape - cool.f_Ag_tape, 0);
tbase = lyr.Hastelloy + lyr.REBCO + lyr.buffer;
if tbase <= 0
    error('step7:stackLayers','Spessore del substrato/REBCO/buffer non valido.');
end
cool.rho_base_tape = (8890*lyr.Hastelloy + 6300*lyr.REBCO + 5000*lyr.buffer)/tbase;
cool.rho_stack = cool.f_Cu_tape*8960 + cool.f_Ag_tape*10490 + ...
                 cool.f_base_tape*cool.rho_base_tape;

% Percorsi termici trasversi: non compare alcun termine 1/h_contact.
cool.P_sc       = cab.N_slots*2*lay.slot_d;
cool.ell_sc     = data.tape.width.val/2;
cool.ell_cu_in  = lay.t_web_in/2;
cool.ell_cu_out = lay.t_web_out/2;
cool.ell_j      = lay.t_jacket/2;
cool.contact    = 'perfect';

[~, kst_ref] = stackVolProps7(T_op_guess,cool);
kcu_ref = kCu(T_op_guess,cool.MR_Cu,cool.RRR);
kss_ref = kSS(T_op_guess);
cool.h_sc = 1/(cool.ell_sc/max(kst_ref,eps) + cool.ell_cu_in/max(kcu_ref,eps));
cool.h_cj = 1/(cool.ell_cu_out/max(kcu_ref,eps) + cool.ell_j/max(kss_ref,eps));

[rho0,~,cp0,c0,mu0,kHe0] = heProps(cool.T_in,cool.p_out); %#ok<ASGLU>
v0    = (cool.mdot_max/2)/(rho0*cool.A_He);
t_ac  = (cool.L/2)/c0;
t_adv = (cool.L/2)/max(v0,eps);

fprintf('\n--- 7.1 Modello e scale temporali ---\n');
fprintf('  L_DP = %.2f m | due rami da %.2f m | ingresso He a s = L/2\n',cool.L,cool.L/2);
fprintf('  A_He = %.2f mm^2 | D_h = %.3f mm | P_wet = %.2f mm\n', ...
    cool.A_He*1e6,cool.D_h*1e3,cool.P_wet*1e3);
fprintf('  T_in base = %.2f K | p_out riferimento = %.1f bar | mdot_max = %.1f g/s\n', ...
    cool.T_in,cool.p_out/1e5,cool.mdot_max*1e3);
fprintf('  t_ac = L_ramo/c = %.4f s | t_ramp = %.0f s | rapporto = %.3e\n', ...
    t_ac,t_ramp,t_ac/t_ramp);
fprintf('  tempo di transito a mdot_max = %.1f s\n',t_adv);
fprintf('  -> pressione/velocita'' quasi-stazionarie; energia pienamente transitoria.\n');
fprintf('  contatto stack-former: R_contact = 0 | h_eq rappresentativo = %.0f W/m^2/K\n', ...
    cool.h_sc);

% ------------------------------------------------------------------------
%  7.2  Sorgenti AC sulla griglia nativa dello Step 6
% ------------------------------------------------------------------------
if ~isfield(cool,'lam_grid') || ~exist('P_mag_nodo_t','var') || ~exist('P_coup_nodo_t','var')
    error('step7:missingLossHistory', ...
        'Mancano P_mag_nodo_t/P_coup_nodo_t o cool.lam_grid dallo Step 6.');
end

nPerStrand7 = size(P_mag_nodo_t,1)/cab.N_slots;
if abs(nPerStrand7-round(nPerStrand7)) > 1e-12
    error('step7:lossSize','Numero di nodi delle perdite non divisibile per N_slots.');
end
nPerStrand7 = round(nPerStrand7);
if numel(s_cable) ~= nPerStrand7
    error('step7:lossGrid','s_cable e matrice delle perdite hanno dimensioni incompatibili.');
end
nLam7 = numel(cool.lam_grid);
if size(P_mag_nodo_t,2) ~= nLam7 || size(P_coup_nodo_t,2) ~= nLam7
    error('step7:lossTime','Numero di livelli temporali incompatibile con cool.lam_grid.');
end

q_sc_cab_t = zeros(nPerStrand7,nLam7);
q_cu_cab_t = zeros(nPerStrand7,nLam7);
for j = 1:cab.N_slots
    idxj = (j-1)*nPerStrand7 + (1:nPerStrand7);
    q_sc_cab_t = q_sc_cab_t + P_mag_nodo_t(idxj,:);
    q_cu_cab_t = q_cu_cab_t + P_coup_nodo_t(idxj,:);
end
q_tot_cab_t = q_sc_cab_t + q_cu_cab_t;

if size(P_ac_cable_t,1) ~= nPerStrand7 || size(P_ac_cable_t,2) ~= nLam7
    error('step7:lossReference','P_ac_cable_t ha dimensioni incompatibili.');
end
errQ7 = max(abs(q_tot_cab_t(:)-P_ac_cable_t(:)))/max(max(P_ac_cable_t(:)),eps);
if errQ7 > 1e-8
    error('step7:sourceBalance','Le sorgenti Step 7 non ricostruiscono Q_AC: errore %.2e.',errQ7);
end

E7 = trapz(cool.lam_grid,trapz(s_cable,q_tot_cab_t,1))*t_ramp;
fprintf('\n--- 7.2 Sorgenti termiche ---\n');
fprintf('  al campo finale: isteresi %.1f W | coupling %.1f W | totale %.1f W\n', ...
    trapz(s_cable,q_sc_cab_t(:,end)),trapz(s_cable,q_cu_cab_t(:,end)), ...
    trapz(s_cable,q_tot_cab_t(:,end)));
fprintf('  energia integrata sulla rampa = %.2f kJ (Step 6: %.2f kJ)\n',E7/1e3,E_ac_tot/1e3);
fprintf('  errore di ricostruzione Q(s,t) = %.2e\n',errQ7);

% ------------------------------------------------------------------------
%  7.3  T_cs sulla griglia nativa
% ------------------------------------------------------------------------
nPerStrandT = numel(Tcs_profile)/cab.N_slots;
if abs(nPerStrandT-round(nPerStrandT)) > 1e-12 || round(nPerStrandT) ~= numel(s_center)
    error('step7:TcsSize','Tcs_profile non compatibile con s_center e N_slots.');
end
nPerStrandT = round(nPerStrandT);
Tcs_mat7  = reshape(Tcs_profile,nPerStrandT,cab.N_slots);
cert_mat7 = reshape(Tcs_cert,nPerStrandT,cab.N_slots);
[Tcs_s7,iWorst7] = min(Tcs_mat7,[],2);
cert_s7 = false(nPerStrandT,1);
for k = 1:nPerStrandT
    cert_s7(k) = cert_mat7(k,iWorst7(k));
end

fprintf('\n--- 7.3 Current sharing per il termico ---\n');
fprintf('  T_cs nativa: %.2f ... %.2f K | punti certificati %d/%d\n', ...
    min(Tcs_s7),max(Tcs_s7),nnz(cert_s7),numel(cert_s7));
fprintf('  Il profilo di pieno campo e'' usato in ogni istante: scelta conservativa.\n');

% ------------------------------------------------------------------------
%  7.4  CONVERGENZA SPAZIALE DELLA SOLUZIONE TERMICA
% ------------------------------------------------------------------------
% Lo studio di convergenza viene fatto a dt FISSO, variando solo dx.
% La soluzione con il passo spaziale piu' piccolo e' usata come riferimento.
%
% Si usa una griglia MOLTO piu' fitta per rendere il diagramma log-log piu'
% informativo. I passi sono distribuiti in scala logaritmica e addensati
% soprattutto verso il lato fine della discretizzazione. Il riferimento e'
% il dx piu' piccolo della serie.

% Passo temporale mantenuto FISSO nello studio di convergenza spaziale.
% In questo modo varia solo dx e l'errore misurato resta puramente spaziale.
conv.dt_fixed = 0.625;                 % [s]

dxTarget = sort(unique([ ...
    logspace(log10(0.005),log10(1.6),12), ...
    2*logspace(log10(0.005),log10(0.8),10), ...
    5*logspace(log10(0.005),log10(0.10),8)]),'descend');

% La griglia globale deve avere un numero DISPARI di nodi, cosi' il centro
% del DP e' un nodo esatto. Il numero di intervalli viene quindi reso pari.
nInt = 2*round(cool.L./(2*dxTarget));
nInt = max(nInt,2);
NxTmp = nInt + 1;

% Elimina eventuali duplicati prodotti dall'arrotondamento.
[conv.Nx,iaNx] = unique(NxTmp,'stable');
conv.dx_target = dxTarget(iaNx);
conv.dx = cool.L ./ (conv.Nx-1);
conv.nLevel = numel(conv.Nx);

fprintf('\n--- 7.4 Studio di convergenza spaziale della soluzione ---\n');
fprintf('  dt mantenuto fisso = %.3f s\n',conv.dt_fixed);
fprintf('  dx testati [m]: ');
fprintf('%.4g ',conv.dx);
fprintf('\n');

% ------------------------------------------------------------------------
%  7.4a  Soluzione di riferimento sul dx piu' piccolo
% ------------------------------------------------------------------------
coolRef = cool;
coolRef.Nx = conv.Nx(end);
coolRef.dt = conv.dt_fixed;

[s_ref,qsc_ref,qcu_ref,Tcs_ref,cert_ref] = makeThermalGrid7(coolRef.Nx,cool.L, ...
    s_cable,q_sc_cab_t,q_cu_cab_t,s_center,Tcs_s7,cert_s7); %#ok<ASGLU>

optRef = struct;
optRef.wantHist = true;
optRef.verbose = false;
optRef.guess = NaN;
optRef.halfWidth = NaN;

Fref = findMinFlow7(coolRef,s_ref,qsc_ref,qcu_ref,Tcs_ref,t_ramp, ...
    A_Cu,A_SC,A_SS,P_Cu_SS,dT_req,optRef);

if ~Fref.solverOK
    error('step7:referenceSolverFailure', ...
        'Il solver termico e'' fallito sulla griglia di riferimento.');
end
if ~Fref.feasible
    error('step7:referenceNoFeasibleFlow', ...
        'Nessuna portata <= %.1f g/s soddisfa il margine sulla griglia di riferimento.', ...
        cool.mdot_max*1e3);
end

mdot_min = Fref.mdot;
cool.mdot_design = mdot_min;

Rd = Fref.R;
Tcs_th = Tcs_ref;
cert_th = cert_ref;
s_th = s_ref;
q_sc_th_t = qsc_ref;
q_cu_th_t = qcu_ref;

margin_th = Tcs_th - Rd.Tc_max_t;
[margin7,iWorstNode] = min(margin_th);
Tc_max_global = max(Rd.Tc_max_t);
dp_th = Fref.dp;
p_in_th = Fref.p_in;

conv.dx_ref = conv.dx(end);
conv.Nx_ref = conv.Nx(end);
conv.mdot_ref = mdot_min;

% Soluzione termica di riferimento: massimi spazio-temporali dei 4 campi.
Uref = [Rd.Tj_max_t; Rd.Tcu_max_t; Rd.Tc_max_t; Rd.THe_max_t];
normUref = max(norm(Uref,2),eps);

fprintf('  riferimento: Nx=%d | dx_ref=%.5f m | mdot=%.4f g/s\n', ...
    conv.Nx_ref,conv.dx_ref,mdot_min*1e3);

% ------------------------------------------------------------------------
%  7.4b  Errori delle griglie piu' grossolane rispetto al riferimento
% ------------------------------------------------------------------------
conv.err_rel = nan(conv.nLevel,1);
conv.margin = nan(conv.nLevel,1);
conv.Tscmax = nan(conv.nLevel,1);
conv.dp = nan(conv.nLevel,1);
conv.results = cell(conv.nLevel,1);

for lev = 1:conv.nLevel
    if lev == conv.nLevel
        % riferimento
        conv.err_rel(lev) = 0;
        conv.margin(lev) = margin7;
        conv.Tscmax(lev) = Tc_max_global;
        conv.dp(lev) = dp_th;
        conv.results{lev} = Fref;
        continue
    end

    [sL,qscL,qcuL,TcsL,certL] = makeThermalGrid7(conv.Nx(lev),cool.L, ...
        s_cable,q_sc_cab_t,q_cu_cab_t,s_center,Tcs_s7,cert_s7); %#ok<ASGLU>

    coolL = cool;
    coolL.Nx = conv.Nx(lev);
    coolL.dt = conv.dt_fixed;

    % STESSO punto operativo della soluzione di riferimento.
    PL = evaluateThermalPoint7(mdot_min,cool.T_in,coolL,sL,qscL,qcuL,TcsL,t_ramp, ...
        A_Cu,A_SC,A_SS,P_Cu_SS,dT_req,false);

    if ~PL.solverOK
        error('step7:convergenceSolverFailure', ...
            'Il solver termico e'' fallito a Nx=%d, dx=%.5f m.', ...
            conv.Nx(lev),conv.dx(lev));
    end

    % Interpola i quattro profili sulla griglia piu' fine.
    TjL  = interp1(sL,PL.R.Tj_max_t, s_ref,'linear','extrap');
    TcuL = interp1(sL,PL.R.Tcu_max_t,s_ref,'linear','extrap');
    TscL = interp1(sL,PL.R.Tc_max_t, s_ref,'linear','extrap');
    THeL = interp1(sL,PL.R.THe_max_t,s_ref,'linear','extrap');

    UL = [TjL; TcuL; TscL; THeL];

    conv.err_rel(lev) = norm(UL-Uref,2)/normUref;
    conv.margin(lev) = PL.margin;
    conv.Tscmax(lev) = PL.Tscmax;
    conv.dp(lev) = PL.dp;
    conv.results{lev} = PL;

    fprintf('  Nx=%5d | dx=%8.5f m -> errore rel.=%.4e | ', ...
        conv.Nx(lev),conv.dx(lev),conv.err_rel(lev));
    fprintf('margine=%.4f K | T_sc,max=%.3f K\n', ...
        PL.margin,PL.Tscmax);
end

% Ordine osservato: fit sui soli punti con errore positivo e finito.
maskConv = isfinite(conv.err_rel) & conv.err_rel > 0;
if nnz(maskConv) >= 2
    pp = polyfit(log(conv.dx(maskConv)),log(conv.err_rel(maskConv)),1);
    conv.order = pp(1);
    conv.fitC = exp(pp(2));
else
    conv.order = NaN;
    conv.fitC = NaN;
end

% Errore del penultimo passo rispetto al riferimento.
conv.err_penultimate = conv.err_rel(end-1);
conv.pass = conv.err_penultimate < 0.01;       % 1 %

fprintf('  errore penultimo livello = %.4e (%.3f %%)\n', ...
    conv.err_penultimate,100*conv.err_penultimate);
fprintf('  ordine osservato log-log = %.3f\n',conv.order);
if conv.pass
    fprintf('  -> CONVERGENZA SPAZIALE ACCETTATA: errore < 1 %% rispetto al riferimento.\n');
else
    warning('step7:spatialConvergence', ...
        'Errore del penultimo livello = %.3f %%: la griglia non e'' ancora entro 1%%.', ...
        100*conv.err_penultimate);
end

% Campioni della curva margine-portata sulla griglia di riferimento.
md_plot = unique(max(cool.mdot_floor,min(cool.mdot_max, ...
          [0.75 0.90 1.00 1.10 1.25 1.50]*mdot_min)));
md_scan = md_plot(:);
mar_scan = nan(size(md_scan));
for k = 1:numel(md_scan)
    Pk = evaluateThermalPoint7(md_scan(k),cool.T_in,coolRef,s_ref,qsc_ref,qcu_ref, ...
        Tcs_ref,t_ramp,A_Cu,A_SC,A_SS,P_Cu_SS,dT_req,false);
    mar_scan(k) = Pk.margin;
end

% Il resto dello Step 7 usa la griglia di riferimento.
cool.Nx = coolRef.Nx;
cool.dt = coolRef.dt;

fprintf('\n  RISULTATO DI PROGETTO SULLA GRIGLIA DI RIFERIMENTO\n');
fprintf('  Nx=%d | dx=%.5f m | dt=%.3f s\n',cool.Nx,conv.dx_ref,cool.dt);
fprintf('  mdot minima = %.3f g/s su %.1f g/s disponibili\n',mdot_min*1e3,cool.mdot_max*1e3);
fprintf('  margine minimo = %.3f K a s = %.2f m\n',margin7,s_th(iWorstNode));
fprintf('  massimi spazio-temporali: T_sc=%.3f K | T_Cu=%.3f K | T_He=%.3f K\n', ...
    max(Rd.Tc_max_t),max(Rd.Tcu_max_t),max(Rd.THe_max_t));
fprintf('  dp ramo A/B = %.4f / %.4f bar | p_in,max = %.3f bar\n', ...
    Rd.dpA/1e5,Rd.dpB/1e5,p_in_th/1e5);

% ------------------------------------------------------------------------
%  7.5  FRONTIERA OPERATIVA T_in,max(m_dot)
% ------------------------------------------------------------------------
% Il brief consente T_in fra 4.5 e 50 K e m_dot <= 15 g/s. Non esiste un
% unico "ottimo" senza una funzione costo: si costruisce quindi la frontiera
% delle coppie che soddisfano esattamente il margine di 5 K.
%
% Per contenere il tempo di calcolo, la frontiera completa usa il penultimo
% livello spaziale della convergenza, mantenendo lo stesso dt fisso. Il punto piu' importante, T_in,max a 15 g/s,
% viene poi ricalcolato anche sulla griglia finale piu' fine.
env.Nx = conv.Nx(end-1);
env.dt = conv.dt_fixed;
cool_env = cool;
cool_env.Nx = env.Nx;
cool_env.dt = env.dt;

[s_env,qsc_env,qcu_env,Tcs_env,cert_env] = makeThermalGrid7(env.Nx,cool.L, ...
    s_cable,q_sc_cab_t,q_cu_cab_t,s_center,Tcs_s7,cert_s7); %#ok<ASGLU>

% Portate scelte per descrivere la frontiera.
md_env = [mdot_min 5e-3 7.5e-3 10e-3 12.5e-3 cool.mdot_max];
md_env = md_env(md_env >= mdot_min-1e-12 & md_env <= cool.mdot_max+1e-12);
md_env = unique(md_env,'sorted');

env.mdot = md_env(:);
env.Tin_max = nan(size(env.mdot));
env.margin = nan(size(env.mdot));
env.p_in = nan(size(env.mdot));
env.Tscmax = nan(size(env.mdot));

Tin_lo_brief = data.op.Tin_range.val(1);
Tin_hi_brief = data.op.Tin_range.val(2);
Tin_hard = min(Tin_hi_brief,min(Tcs_env)-dT_req);
Tin_hard = max(Tin_hard,Tin_lo_brief);

fprintf('\n--- 7.5 Frontiera operativa T_in,max(mdot) ---\n');
fprintf('  griglia frontiera: Nx=%d, dt=%.3f s\n',env.Nx,env.dt);
fprintf('  limite termico superiore preliminare: T_in <= min(T_cs)-5 K = %.2f K\n',Tin_hard);

for k = 1:numel(env.mdot)
    ET = findMaxTin7(env.mdot(k),cool_env,s_env,qsc_env,qcu_env,Tcs_env,t_ramp, ...
        A_Cu,A_SC,A_SS,P_Cu_SS,dT_req,Tin_lo_brief,Tin_hard);
    env.Tin_max(k) = ET.Tin;
    env.margin(k) = ET.margin;
    env.p_in(k) = ET.p_in;
    env.Tscmax(k) = ET.Tscmax;
    if ET.feasibleAtMinTin
        fprintf('  mdot=%6.3f g/s -> T_in,max=%6.3f K | margine=%6.3f K | p_in=%.3f bar\n', ...
            env.mdot(k)*1e3,env.Tin_max(k),env.margin(k),env.p_in(k)/1e5);
    else
        fprintf('  mdot=%6.3f g/s -> non soddisfa 5 K nemmeno a T_in=%.2f K\n', ...
            env.mdot(k)*1e3,Tin_lo_brief);
    end
end

% Punto estremo a 15 g/s verificato sulla griglia finale.
env15 = findMaxTin7(cool.mdot_max,cool,s_th,q_sc_th_t,q_cu_th_t,Tcs_th,t_ramp, ...
    A_Cu,A_SC,A_SS,P_Cu_SS,dT_req,Tin_lo_brief,min(Tin_hi_brief,min(Tcs_th)-dT_req));
env.Tin_max_15_fine = env15.Tin;
env.margin_15_fine = env15.margin;
env.p_in_15_fine = env15.p_in;
env.Tscmax_15_fine = env15.Tscmax;

fprintf('  verifica FINE a 15 g/s: T_in,max = %.3f K | margine = %.3f K | p_in = %.3f bar\n', ...
    env.Tin_max_15_fine,env.margin_15_fine,env.p_in_15_fine/1e5);
fprintf('  -> la curva rappresenta una frontiera di progetto, non un unico optimum economico.\n');

% ------------------------------------------------------------------------
%  7.6  Chiusura del Loop B
% ------------------------------------------------------------------------
fprintf('\n--- 7.6 Chiusura del Loop B ---\n');
fprintf('  T_op usata a monte per le perdite = %.3f K\n',T_op_guess);
fprintf('  T_sc massima al punto base dello Step 7 = %.3f K\n',Tc_max_global);
if Tc_max_global <= T_op_guess
    fprintf('  -> LOOP B CHIUSO: T_sc,max <= T_op. L''assunzione a monte era conservativa.\n');
else
    warning('step7:loopB', ...
        ['LOOP B NON CHIUSO: T_sc,max = %.3f K > T_op = %.3f K. ' ...
         'Rilanciare con T_op_override = %.2f K.'], ...
         Tc_max_global,T_op_guess,Tc_max_global);
end
if env.Tscmax_15_fine <= T_op_guess
    fprintf('  anche il punto caldo della frontiera a 15 g/s resta sotto T_op (%.3f <= %.3f K).\n', ...
        env.Tscmax_15_fine,T_op_guess);
else
    warning('step7:loopBEnvelope', ...
        'La frontiera a 15 g/s porta T_sc,max sopra T_op: le perdite vanno ricalcolate con T_op piu'' alto.');
end

% ------------------------------------------------------------------------
%  7.7  Verifiche idrauliche e validita' della riduzione slow-transient
% ------------------------------------------------------------------------
rel_dp = Rd.dp_imbal/max(mean([Rd.dpA Rd.dpB]),eps);
rel_flow_est = rel_dp/1.8;
fprintf('\n--- 7.7 Verifiche idrauliche ---\n');
fprintf('  punto base: p_out=%.3f bar | p_in,max=%.3f bar | dp=%.4f bar\n', ...
    cool.p_out/1e5,p_in_th/1e5,dp_th/1e5);
fprintf('  Re = %.2e ... %.2e | h_He = %.0f ... %.0f W/m^2/K\n', ...
    min(Rd.Re),max(Rd.Re),min(Rd.hHe),max(Rd.hHe));
fprintf('  v = %.3f ... %.3f m/s | tempo di transito = %.1f s\n', ...
    min(Rd.v),max(Rd.v),Rd.t_transit);
fprintf('  squilibrio dp fra i rami = %.3f %% -> correzione stimata al 50/50 = %.3f %%\n', ...
    rel_dp*100,rel_flow_est*100);
fprintf('  dp/p_ref = %.4f %% al punto base\n',dp_th/cool.p_out*100);
fprintf('  punto frontiera a 15 g/s: p_in = %.3f bar\n',env.p_in_15_fine/1e5);

if t_ac/t_ramp > 0.01
    warning('step7:acousticScale', ...
        't_ac/t_ramp = %.3e: la separazione di scala non e'' forte.',t_ac/t_ramp);
end
if dp_th/cool.p_out > 0.05 || (env.p_in_15_fine-cool.p_out)/cool.p_out > 0.05
    warning('step7:pressureApprox', ...
        'La perdita di pressione supera il 5%%: usare p locale anche nelle proprieta'' dell''He.');
end
if rel_flow_est > 0.02
    warning('step7:flowSplit', ...
        'Il 50/50 richiederebbe una correzione stimata superiore al 2%%.');
end

% ------------------------------------------------------------------------
%  7.8  Figure dello Step 7
% ------------------------------------------------------------------------
figure('Name','STEP 7 - GANDALF ridotto: soluzione finale','Color','w', ...
       'Position',[180 180 1180 760]);

subplot(2,2,1); hold on; grid on; box on;
plot(s_th,Rd.T_He,'LineWidth',1.4);
plot(s_th,Rd.T_sc,'LineWidth',1.5);
plot(s_th,Rd.T_c ,'LineWidth',1.4);
plot(s_th,Rd.T_j ,'LineWidth',1.2);
xlabel('s [m]'); ylabel('T a t_{end} [K]');
title(sprintf('Profili finali, t = %.0f s',cool.t_end));
legend({'He','stack','former','jacket'},'Location','best','FontSize',8);

subplot(2,2,2); hold on; grid on; box on;
plot(Rd.t_hist,Rd.THe_hist,'LineWidth',1.4);
plot(Rd.t_hist,Rd.Tc_hist ,'LineWidth',1.5);
xline(t_ramp,'k--','fine rampa');
xlabel('t [s]'); ylabel('T massima [K]');
title('Storia temporale');
legend({'He','stack'},'Location','best','FontSize',8);

subplot(2,2,3); hold on; grid on; box on;
hfill = fill([s_th;flipud(s_th)],[Rd.Tc_max_t;flipud(Tcs_th)],[0.8 0.95 0.8], ...
    'EdgeColor','none','FaceAlpha',0.5);
uistack(hfill,'bottom');
plot(s_th,Tcs_th,'LineWidth',1.5);
plot(s_th,Rd.Tc_max_t,'LineWidth',1.5);
xlabel('s [m]'); ylabel('T [K]');
title('T_{cs} contro T_{sc,max}');
legend({'T_{cs}-T_{sc,max}','T_{cs}','T_{sc,max}'},'Location','best','FontSize',8);

subplot(2,2,4); hold on; grid on; box on;
plot(s_th,margin_th,'LineWidth',1.6);
yline(dT_req,'r--','5 K richiesti');
plot(s_th(iWorstNode),margin7,'kp','MarkerSize',13,'MarkerFaceColor','y');
xlabel('s [m]'); ylabel('margine [K]');
title(sprintf('Margine termico minimo %.3f K',margin7));

figure('Name','STEP 7 - convergenza spaziale log-log','Color','w', ...
       'Position',[220 220 780 560]);
hold on; grid on; box on;

maskPlot = isfinite(conv.err_rel) & conv.err_rel > 0;
dxPlot = conv.dx(maskPlot);
errPlot = conv.err_rel(maskPlot);

% Ordina per dx crescente per una lettura naturale sul log-log.
[dxPlot,ordPlot] = sort(dxPlot);
errPlot = errPlot(ordPlot);

loglog(dxPlot,errPlot,'o-','LineWidth',1.7,'MarkerSize',7);

% Retta di best fit sullo stesso piano log-log.
if isfinite(conv.order)
    dxFit = logspace(log10(min(dxPlot)),log10(max(dxPlot)),100);
    errFit = conv.fitC*dxFit.^conv.order;
    loglog(dxFit,errFit,'--','LineWidth',1.3);
    legend({'errore numerico',sprintf('fit: O(dx^{%.2f})',conv.order)}, ...
           'Location','best');
else
    legend({'errore numerico'},'Location','best');
end

set(gca,'XScale','log','YScale','log');
xlabel('\Deltax [m]');
ylabel('||U-U_{ref}||_2 / ||U_{ref}||_2');
title(sprintf('Convergenza spaziale, dx_{ref} = %.4g m',conv.dx_ref));

figure('Name','STEP 7 - frontiera operativa Tin-mdot','Color','w','Position',[260 220 760 520]);
hold on; grid on; box on;
goodEnv = isfinite(env.Tin_max);
plot(env.mdot(goodEnv)*1e3,env.Tin_max(goodEnv),'o-','LineWidth',1.8);
plot(cool.mdot_max*1e3,env.Tin_max_15_fine,'kp','MarkerSize',12,'MarkerFaceColor','y');
xlabel('mdot [g/s]'); ylabel('T_{in,max} [K]');
title('Frontiera operativa: margine minimo = 5 K');
legend({'frontiera, griglia intermedia','15 g/s, griglia finale'},'Location','best');

% Sezione 2D: post-processing sulla soluzione base finale.
Tnodes7 = [Rd.THe_max_t,Rd.Tc_max_t,Rd.Tcu_max_t,Rd.Tj_max_t];
Tspan7 = max(Tnodes7,[],2)-min(Tnodes7,[],2);
[span7,iGrad7] = max(Tspan7);
[~,kst7] = stackVolProps7(Rd.Tc_max_t(iGrad7),cool);
geoSecT = struct('r_cab',r_cab,'r_form_o',r_former_out,'r_chan',r_channel, ...
                 'r_slot_in',r_slot_in,'r_slot_out',r_slot_out, ...
                 'slot_w',lay.slot_w,'slot_d',lay.slot_d, ...
                 'R_core',R_core_derived,'N_slots',cab.N_slots);
kmapT = struct('k_Cu',kCu(Rd.Tcu_max_t(iGrad7),cool.MR_Cu,cool.RRR), ...
               'k_SS',kSS(Rd.Tj_max_t(iGrad7)), ...
               'k_st',kst7,'A_Cu',A_Cu,'A_SC',A_SC);
srcT = struct('q_cu',max(q_cu_th_t(iGrad7,:)), ...
              'q_sc',max(q_sc_th_t(iGrad7,:)));
[Tmap7,X7,Y7,info7] = sectionThermal2D(geoSecT,kmapT,srcT, ...
                         Rd.hHe(iGrad7),Rd.THe_max_t(iGrad7),181);

figure('Name','STEP 7 - sezione critica 2D','Color','w','Position',[200 200 620 560]);
pcolor(X7*1e3,Y7*1e3,Tmap7); shading interp; axis equal tight; hold on;
cb7=colorbar; cb7.Label.String='T [K]';
angSec7=linspace(0,2*pi,200);
plot(r_cab*1e3*cos(angSec7),r_cab*1e3*sin(angSec7),'k-','LineWidth',1.2);
plot(r_former_out*1e3*cos(angSec7),r_former_out*1e3*sin(angSec7),'k--','LineWidth',0.8);
plot(r_channel*1e3*cos(angSec7),r_channel*1e3*sin(angSec7),'k-','LineWidth',1.2);
for j=1:cab.N_slots
    aa=(j-1)*2*pi/cab.N_slots;
    xc=R_core_derived*cos(aa)*1e3; yc=R_core_derived*sin(aa)*1e3;
    Rm=[cos(aa) -sin(aa); sin(aa) cos(aa)];
    cor=Rm*[-lay.slot_d/2 lay.slot_d/2 lay.slot_d/2 -lay.slot_d/2; ...
            -lay.slot_w/2 -lay.slot_w/2 lay.slot_w/2 lay.slot_w/2]*1e3;
    plot(xc+cor(1,[1:end 1]),yc+cor(2,[1:end 1]),'k-','LineWidth',1.0);
end
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Sezione a escursione massima: s=%.1f m, span=%.3f K',s_th(iGrad7),span7));
fprintf('  mappa 2D: T = %.3f ... %.3f K nella sezione scelta\n',info7.T_min,info7.T_max);

% Spirale del pancake colorata con la massima T_sc nel tempo.
[Pspir7,~,sSpir7] = buildDPcenterline(N_turn,disc.N_seg_field,R_in,cab.width,0);
nHalf7=(size(Pspir7,1)+1)/2;
xs7=Pspir7(1:nHalf7,1)*1e3; ys7=Pspir7(1:nHalf7,2)*1e3; ss7=sSpir7(1:nHalf7);
Tsp7=interp1(s_th,Rd.Tc_max_t,ss7,'linear','extrap');
THesp7=interp1(s_th,Rd.THe_max_t,ss7,'linear','extrap');
Msp7=interp1(s_th,margin_th,ss7,'linear','extrap');

figure('Name','STEP 7 - pancake: spirale in temperatura','Color','w','Position',[150 150 1250 560]);
subplot(1,2,1); hold on; axis equal; box on;
patch([xs7;nan],[ys7;nan],[Tsp7;nan],'EdgeColor','interp','FaceColor','none','LineWidth',7);
cb8=colorbar; cb8.Label.String='T_{sc,max} [K]';
plot(xs7(1),ys7(1),'ks','MarkerFaceColor','w','MarkerSize',8);
plot(xs7(end),ys7(end),'ko','MarkerFaceColor','k','MarkerSize',7);
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Pancake, mdot = %.3f g/s',mdot_min*1e3));

subplot(1,2,2); hold on; grid on; box on;
plot(ss7,THesp7,'LineWidth',1.4);
plot(ss7,Tsp7,'LineWidth',1.7);
plot(ss7,Msp7,'LineWidth',1.4);
yline(dT_req,'r--','margine richiesto');
xlabel('s [m]'); ylabel('T / margine [K]'); title('Profilo lungo il pancake');
legend({'He','stack max','margine'},'Location','best','FontSize',8);

% Profili finali, distinti dai massimi spazio-temporali.
q_sc_th = q_sc_th_t(:,end);
q_cu_th = q_cu_th_t(:,end);
q_j_th  = zeros(size(q_sc_th));
q_th    = q_sc_th + q_cu_th;

save('step7_cooling.mat','cool','conv','env','s_th','q_th','q_sc_th','q_cu_th','q_j_th', ...
     'q_sc_th_t','q_cu_th_t','P_ac_avg','E_ac_tot','Tcs_th','cert_th', ...
     'margin_th','margin7','mdot_min','Rd','md_scan','mar_scan', ...
     'dp_th','p_in_th','Tc_max_global','t_ac','t_adv');

tmr = addTimer(tmr,'Step 7    convergenza + frontiera operativa',toc(tStep7));
fprintf('\n--> Salvato step7_cooling.mat\n');

%% ------------------------- RIEPILOGO FINALE -----------------------------
fprintf('\n\n=========================================================\n');
fprintf(' RIEPILOGO DEL DESIGN (Step 1-7)\n');
fprintf('=========================================================\n');
fprintf('  Solenoide : %d DP x %d turn radiali, R_in = %.2f m, H = %.3f m / %.2f m\n', ...
    N_layer, N_turn, R_in, N_layer*2*cab.width, H_max);
fprintf('  Cavo      : D = %.0f mm, %d slot, %d nastri/slot, twist %d giri/turn\n', ...
    cab.width*1e3, cab.N_slots, cab.N_tapes_slot, cab.N_twist);
fprintf('  Sezione   : canale He %.2f mm, D_h %.3f mm, jacket SS %.1f mm\n', ...
    2*r_channel*1e3, D_h*1e3, lay.t_jacket*1e3);
fprintf('  Correnti  : %.0f kA -> %.0f A/slot -> %.1f A/nastro\n', ...
    I_nom/1e3, I_slot, I_tape_op);
fprintf('  Campo     : B_pk = %.2f T (target %.1f T, scarto %+.2f %%)\n', ...
    Bpk, B_target, (Bpk-B_target)/B_target*100);
fprintf('  DAE       : scarto dall''equiripartizione = %.2e %%\n', devMax);
fprintf('  AC losses : %.1f W (media rampa) / %.1f W (campo finale) sul DP, %.2f W/m\n', ...
    P_ac_avg, P_ac_tot, max(P_ac_cable));
fprintf('              di picco, %.1f kJ sulla rampa (Q(s,t), non congelata)\n', E_ac_tot/1e3);
fprintf('  RAFFREDD. : mdot minima = %.3f g/s su %.1f disponibili (%.0f %%)\n', ...
    cool.mdot_design*1e3,cool.mdot_max*1e3,cool.mdot_design/cool.mdot_max*100);
fprintf('              griglia finale Nx=%d, dt=%.3f s | convergenza: %s\n', ...
    cool.Nx,cool.dt,ternary(conv.pass,'OK','DA RAFFINARE'));
fprintf('              T_in = %.2f K, p_in = %.3f bar, dp = %.4f bar\n', ...
    cool.T_in,p_in_th/1e5,dp_th/1e5);
fprintf('              massimi spazio-temporali: T_sc %.2f K | T_Cu %.2f K | T_He %.2f K\n', ...
    max(Rd.Tc_max_t),max(Rd.Tcu_max_t),max(Rd.THe_max_t));
fprintf('              MARGINE MINIMO %.3f K (richiesto %.0f K)\n',margin7,dT_req);
fprintf('              profili a t_end: stack %.2f...%.2f K | former %.2f...%.2f K | He %.2f...%.2f K\n', ...
    min(Rd.T_sc),max(Rd.T_sc),min(Rd.T_c),max(Rd.T_c),min(Rd.T_He),max(Rd.T_He));
fprintf('  FRONTIERA : a 15 g/s si puo'' arrivare a T_in = %.2f K mantenendo 5 K di margine\n', ...
    env.Tin_max_15_fine);
fprintf('---------------------------------------------------------\n');
fprintf('  FONTI E ASSUNZIONI PRINCIPALI\n');
fprintf('  Materiali solidi:\n');
fprintf('    Cu OFHC: cp e k da correlazioni NIST; RRR = %d.\n',acl.RRR);
fprintf('    SS304: cp e k da correlazioni NIST.\n');
fprintf('    Hastelloy C-276: cp da Lu, Choi, Zhou (2008); sotto 10 K\n');
fprintf('      si usa il tratto misurato approssimato linearmente.\n');
fprintf('    k dell''Hastelloy: NIST Inconel 718 come surrogato Ni-base.\n');
fprintf('  Elio:\n');
fprintf('    proprieta'' real-fluid da tabella embedded in funzione di T e p.\n');
fprintf('  Step 7:\n');
fprintf('    modello GANDALF ridotto slow-transient: energia transitoria di\n');
fprintf('    stack/former/jacket/He; pressione e velocita'' quasi-stazionarie.\n');
fprintf('    contatto termico stack-former ideale: R_contact = 0.\n');
fprintf('    portata divisa nominalmente 50/50 fra i due rami e verificata tramite dp.\n');
fprintf('    criterio base: portata MINIMA con margine >= %.0f K e mdot <= %.1f g/s.\n', ...
    dT_req,cool.mdot_max*1e3);
fprintf('    la portata minima e'' verificata con convergenza spaziale log-log su dx, a dt fisso.\n');
fprintf('    e'' inoltre costruita la frontiera T_in,max(mdot) nel range del brief.\n');
fprintf('    nessun requisito su h_contact viene ricavato dal surplus di margine.\n');
fprintf('  Limiti dichiarati:\n');
fprintf('    proprieta'' He valutate alla pressione di riferimento se dp/p e'' piccolo;\n');
fprintf('    Ag approssimato come Cu e REBCO/buffer come substrato nello stack;\n');
fprintf('    T_cs puo'' essere un limite inferiore dove i dati non consentono una radice.\n');
fprintf('---------------------------------------------------------\n');
fprintf('  IL COMPROMESSO CENTRALE DEL PROGETTO:\n');
fprintf('   Lo stack e'' assunto elettricamente ben accoppiato al rame nel worst case di coupling.\n');
fprintf('   Questo rende rho_eff metallica e l''accoppiamento il termine dominante\n');
fprintf('   (%.0f %% delle perdite). NON e'' un difetto: e'' potenza da smaltire, e il\n', frac_coup*100);
fprintf('   criterio di accettazione e'' la PORTATA, non la quota di accoppiamento.\n');
fprintf('   Verdetto: %.2f g/s minimi su %.0f disponibili (%.0f %% del budget).\n', ...
    cool.mdot_design*1e3, cool.mdot_max*1e3, cool.mdot_design/cool.mdot_max*100);
if cool.mdot_design < 0.8*cool.mdot_max
    fprintf('   -> ci sta con margine: NON serve twistare piu'' stretto, il cavo resta\n');
    fprintf('      com''e'' e la protezione al quench e'' preservata.\n');
else
    fprintf('   -> il budget di portata e'' quasi esaurito. Valutare la seconda leva:\n');
    fprintf('      twistare piu'' stretto (numeri in 6.7), al costo di rivedere\n');
    fprintf('      disc.N_seg_dae e di piu'' deformazione sui nastri.\n');
end
fprintf('=========================================================\n');

%% ------------------------- PROFILING ------------------------------------
%  Dove se ne va il tempo. Le voci sono ordinate per costo decrescente.
%  Il RESIDUO (totale meno somma delle voci) raccoglie tutto cio' che non e'
%  cronometrato: grafici, stampe, Step 1, Step 5, assemblaggi leggeri. Se il
%  residuo diventasse la voce dominante vorrebbe dire che il collo di bottiglia
%  sta in un blocco non ancora strumentato, e li' andrebbe messo un timer.
tTot = toc(tWall);
fprintf('\n=========================================================\n');
fprintf(' PROFILING - dove se ne va il tempo\n');
fprintf('=========================================================\n');
if isempty(tmr)
    fprintf('  nessun timer registrato\n');
else
    [~, ord] = sort([tmr.sec], 'descend');
    fprintf('  %-46s %9s %7s\n', 'blocco', 'tempo [s]', 'quota');
    for k = ord
        fprintf('  %-46s %9.1f %6.1f%%\n', tmr(k).name, tmr(k).sec, tmr(k).sec/tTot*100);
    end
    resid = tTot - sum([tmr.sec]);
    fprintf('  %-46s %9.1f %6.1f%%\n', ...
        '(residuo: grafici, stampe, blocchi leggeri)', resid, resid/tTot*100);
end
fprintf('  ---------------------------------------------------------\n');
fprintf('  %-46s %9.1f\n', 'TOTALE', tTot);
fprintf('=========================================================\n');

%% ========================================================================
%  FUNZIONI LOCALI
%% ========================================================================
function tmr = addTimer(tmr, name, sec)
% Accumula una voce di profiling. Se il nome esiste gia' i tempi si SOMMANO,
% cosi' un blocco eseguito piu' volte non genera righe duplicate nel riepilogo.
k = find(strcmp({tmr.name}, name), 1);
if isempty(k)
    tmr(end+1).name = name;    %#ok<AGROW>
    tmr(end).sec    = sec;
else
    tmr(k).sec = tmr(k).sec + sec;
end
end

% ------------------------------------------------------------------------
function out = ternary(cond, a, b)
% Helper per stampare esiti di verifica in linea.
if cond, out = a; else, out = b; end
end

% ------------------------------------------------------------------------
function [P_center, theta, s] = buildDPcenterline(N_turn, N_seg, R_in, pitch, z0)
% Linea centrale del cavo di UN double-pancake: pancake inferiore (spirale in
% uscita, z0-pitch/2) + superiore (rientro, z0+pitch/2), conduttore continuo
% con transizione al raggio esterno. Ritorna anche l'angolo poloidale theta
% (serve al twist commensurato) e l'ascissa curvilinea s.
% TOPOLOGIA: il conduttore ENTRA e RIESCE dal raggio ESTERNO, e la transizione
% fra i due pancake avviene al raggio INTERNO.
%   - i lap joint (inter-pancake e terminali) stanno all'esterno, dove il campo
%     e' minimo e dove sono fisicamente accessibili;
%   - il raggio interno, dove il campo e' massimo, resta conduttore CONTINUO.
rmax = R_in + pitch*N_turn;

% pancake inferiore: dall'ESTERNO verso l'INTERNO
th1 = linspace(0, 2*pi*N_turn, N_turn*N_seg + 1)';
r1  = rmax - (pitch/(2*pi))*th1;                     % rmax -> R_in
Pc1 = [r1.*cos(th1), r1.*sin(th1), (z0 - pitch/2)*ones(size(th1))];

% pancake superiore: dall'INTERNO verso l'ESTERNO
th2 = linspace(2*pi*N_turn, 4*pi*N_turn, N_turn*N_seg + 1)';
r2  = R_in + (pitch/(2*pi))*(th2 - 2*pi*N_turn);     % R_in -> rmax
Pc2 = [r2.*cos(th2), r2.*sin(th2), (z0 + pitch/2)*ones(size(th2))];

P_center = [Pc1; Pc2(2:end,:)];       % rimuove il nodo doppio alla transizione
theta    = [th1; th2(2:end)];
d = diff(P_center,1,1);
s = [0; cumsum(sqrt(sum(d.^2,2)))];
end

% ------------------------------------------------------------------------
function [P, E] = buildMacroStack(N_layer, idxCentral, z_layers, N_turn, N_seg, R_in, pitch)
% Solenoide = N_layer double-pancake impilati assialmente, tutti identici a
% quello centrale. Il layer idxCentral e' ESCLUSO (lo modella il micro):
% l'esclusione per INDICE rimuove esattamente e solo il conduttore che il
% modello micro ricostruisce -> nessun doppio conteggio, nessun buco.
P = []; E = []; off = 0;
for iz = 1:N_layer
    if iz == idxCentral, continue; end
    Pc = buildDPcenterline(N_turn, N_seg, R_in, pitch, z_layers(iz));
    P = [P; Pc]; %#ok<AGROW>
    n = (1:size(Pc,1)-1)' + off;
    E = [E; n, n+1]; %#ok<AGROW>
    off = off + size(Pc,1);
end
end

% ------------------------------------------------------------------------
function [That, Nhat, Bhat] = frenetFrame(P_center)
% Terna locale lungo il cavo: T tangente, N radiale-esterno (rispetto a z),
% B binormale. Il piano (N,B) e' la sezione trasversale in cui stanno gli slot.
Nn = size(P_center,1);
d  = diff(P_center,1,1);
ds = sqrt(sum(d.^2,2));
That = zeros(Nn,3);
That(1:end-1,:) = d ./ ds;
That(end,:)     = That(end-1,:);
Zc   = repmat([0 0 1], Nn, 1);
Nhat = cross(Zc, That, 2);
Nhat = Nhat ./ vecnorm(Nhat,2,2);
Bhat = cross(That, Nhat, 2);
end

% ------------------------------------------------------------------------
function alpha = twistAngle(theta, slotIdx, cab)
% Angolo di twist COMMENSURATO, parametrizzato sull'angolo poloidale theta:
%   alpha = theta_slot + N_twist*theta,   N_twist INTERO
% Su ogni turn theta avanza di 2*pi -> il twist compie N_twist giri interi a
% QUALUNQUE raggio, quindi l_p(turn k) = 2*pi*r_k/N_twist e' per costruzione
% un divisore esatto del perimetro di quel turn.
%
% DIMENSIONI: theta e slotIdx devono avere la stessa lunghezza (o uno dei due
% essere scalare). theta e' definito sulla linea centrale (Nn_c x 1); per
% valutarlo su tutti i nodi degli strand va replicato: repmat(theta,N_slots,1).
if ~isscalar(theta) && ~isscalar(slotIdx) && numel(theta) ~= numel(slotIdx)
    error('twistAngle:sizeMismatch', ...
        ['theta (%d elementi) e slotIdx (%d elementi) devono avere la stessa ' ...
         'lunghezza. Replicare theta sui %d strand: repmat(theta, N_slots, 1).'], ...
        numel(theta), numel(slotIdx), cab.N_slots);
end
alpha = (slotIdx-1)*(2*pi/cab.N_slots) + cab.N_twist*theta;
end

% ------------------------------------------------------------------------
function slotIdx = strandNodeSlot(N_slots, Nn_c)
% Vettore (N_slots*Nn_c x 1) con l'indice di slot di ciascun nodo di P_micro,
% nello stesso ordine in cui buildTwistedStrands impila gli strand.
slotIdx = reshape(repmat(1:N_slots, Nn_c, 1), [], 1);
end

% ------------------------------------------------------------------------
function [P, E, strandID, faceInj, faceExt] = buildTwistedStrands(P_center, theta, That, Nhat, Bhat, cab)
% Uno strand equivalente per slot (lo stack di nastri e' ridotto a un solo
% conduttore: spessore ignorato, si conserva solo la posizione), twistato
% secondo il vincolo di commensurabilita'. Costruisce anche la topologia:
%   strandID        : per ogni arco, a quale strand appartiene
%   faceInj/faceExt : nodi della prima/ultima sezione (BC del DAE, Freschi s34-s37)
Nc = size(P_center,1);
P = []; E = []; strandID = []; off = 0;
faceInj = zeros(cab.N_slots,1);  faceExt = zeros(cab.N_slots,1);
for j = 1:cab.N_slots
    alpha = twistAngle(theta, j, cab);          % theta e' Nc x 1, j scalare -> ok
    o3 = Nhat.*(cab.R_core.*cos(alpha)) + Bhat.*(cab.R_core.*sin(alpha));
    P = [P; P_center + o3]; %#ok<AGROW>
    n = (1:Nc-1)' + off;
    E = [E; n, n+1]; %#ok<AGROW>
    strandID = [strandID; j*ones(Nc-1,1)]; %#ok<AGROW>
    faceInj(j) = 1  + off;
    faceExt(j) = Nc + off;
    off = off + Nc;
end
end

% ------------------------------------------------------------------------
function [B_tot, B_ext, B_mut, B_self] = computeDPField(P_macro, E_macro, I_macro, ...
                                    P_micro, E_micro, strandID, I_slot, cab, probeIdx)
% Campo sul DP centrale, in tre contributi fisicamente disgiunti.
%   B_ext : dagli altri DP (macro). Il DP centrale e' escluso dal macro.
%   B_mut : dagli ALTRI strand, calcolato strand per strand ESCLUDENDO gli
%           archi del proprio strand. Necessario: quando il punto di campo
%           coincide con un nodo sorgente, segmentMagneticField1d incontra una
%           forma 0/0 sui due segmenti adiacenti e la AZZERA silenziosamente
%           (fix isfinite) -> mancherebbe un pezzo di campo senza alcun avviso.
%   B_self: contributo del proprio strand, analitico alla superficie del
%           conduttore equivalente, mu0*I_slot/(2*pi*R_core). Sommato in
%           MODULO (scelta conservativa: la direzione ruota col twist).
% NB: segmentMagneticField1d restituisce H [A/m] -> conversione B = Mu0()*H.
Q = P_micro(probeIdx,:);
B_ext = Mu0() * segmentMagneticField1d(P_macro, E_macro, I_macro, Q);

B_mut = zeros(size(Q));
slotAll   = strandNodeSlot(cab.N_slots, size(P_micro,1)/cab.N_slots);
slotProbe = slotAll(probeIdx);
for j = 1:cab.N_slots
    rows = find(slotProbe == j);
    if isempty(rows), continue; end
    keep = (strandID ~= j);                     % esclude il proprio strand
    B_mut(rows,:) = Mu0() * segmentMagneticField1d(P_micro, E_micro(keep,:), ...
                                                   I_slot, Q(rows,:));
end
% Il campo che lo strand genera sulla PROPRIA superficie dipende dal raggio
% del conduttore (cab.r0), NON dalla sua distanza dall'asse del cavo (R_core).
% Questo valore e' gia' il massimo interno allo stack (il self-field e' nullo
% al centro del conduttore e massimo alla superficie): non va quindi corretto
% ulteriormente con il fattore near-edge, che riguarda solo i termini generati
% da sorgenti ESTERNE allo strand.
B_self = Mu0() * I_slot / (2*pi*cab.r0);
B_tot  = B_ext + B_mut;
end

% ------------------------------------------------------------------------
function Bpk = peakFieldForNturn(N_turn, N_layer, idxCentral, z_layers, R_in, cab, N_seg, I_nom)
% Picco di campo totale per un dato N_turn. Il campo e' campionato SOLO sui
% nodi del primo turn (raggio minore): per un solenoide il picco sul conduttore
% sta sempre li'. Riduce di ~N_turn volte il costo della ricerca.
[Pc, cool, ~]   = buildDPcenterline(N_turn, N_seg, R_in, cab.width, z_layers(idxCentral));
[T,N,B]       = frenetFrame(Pc);
[Pm, Em, sid] = buildTwistedStrands(Pc, cool, T, N, B, cab);
[Pma, Ema]    = buildMacroStack(N_layer, idxCentral, z_layers, N_turn, N_seg, R_in, cab.width);
% SONDA sul turn a RAGGIO MINORE, individuato PER RAGGIO e non per indice.
% Con la topologia (esterno -> interno -> esterno) il turn INIZIALE sta al raggio
% MASSIMO: campionare i primi N_seg nodi darebbe il campo MINIMO (sottostima ~5x)
% e il dimensionamento di N_turn sceglierebbe un valore enormemente eccessivo.
% La selezione per raggio e' robusta a qualunque topologia.
Nc    = size(Pc,1);
rr    = hypot(Pm(:,1), Pm(:,2));
probe = find(rr <= R_in + cab.width);
if isempty(probe)
    error('peakField:noProbe','Nessun nodo entro il turn interno: verificare la geometria.');
end
[Bt, ~, ~, Bs] = computeDPField(Pma, Ema, I_nom, Pm, Em, sid, I_nom/cab.N_slots, cab, probe);
Bpk = max(vecnorm(Bt,2,2) + Bs);
end

% ------------------------------------------------------------------------
function A = buildIncidence(E, Nn)
% Matrice di incidenza arco->nodo [Freschi s27]: A(e,start) = -1, A(e,end) = +1,
% cosi' che -A*phi restituisca la tensione lungo l'arco (phi_start - phi_end).
Ne = size(E,1);
A  = sparse([1:Ne, 1:Ne]', [E(:,1); E(:,2)], [-ones(Ne,1); ones(Ne,1)], Ne, Nn);
end

% ------------------------------------------------------------------------
function Rvec = resistanceEJ(iEdge, IcEdge, segLen, E0, nEJ)
% Resistenza non lineare dalla power law E-J [Freschi s5]:
%   E = E0 (I/Ic)^n  ->  R = (E0*l/Ic) * (|I|/Ic)^(n-1)
% nEJ puo' essere SCALARE o VETTORIALE (uno per arco): l'elevamento a potenza
% e' elementwise, quindi passare n(T,B) valutato arco per arco non richiede
% nessuna modifica qui. E' cio' che fa lo Step 4.3b, perche' il campo varia
% di un fattore 3 fra turn esterno e interno e con esso l'esponente.
% Floor su |I| e su R: a corrente esattamente nulla R sarebbe 0 e la matrice
% del sistema singolare.
x    = max(abs(iEdge), 1e-6) ./ IcEdge;
Rvec = (E0 .* segLen ./ IcEdge) .* x.^(nEJ-1);
Rvec = max(Rvec, 1e-15);
end

% ------------------------------------------------------------------------
function Bedge = interpFieldOnEdges(P_dae, E_dae, P_field, Bmag_field)
% Campo medio su ogni arco della mesh DAE, preso dal nodo piu' vicino della
% mesh fine su cui e' stata calcolata la mappa di campo.
mid   = 0.5*(P_dae(E_dae(:,1),:) + P_dae(E_dae(:,2),:));
idx   = nearestNode(P_field, mid);
Bedge = Bmag_field(idx);
end

% ------------------------------------------------------------------------
function idx = nearestNode(X, Y)
% Nearest neighbour senza toolbox: per ogni riga di Y l'indice della riga di X
% piu' vicina. A blocchi, per contenere la memoria.
idx = zeros(size(Y,1),1);
blk = 500;
nx2 = sum(X.^2,2)';
for i = 1:blk:size(Y,1)
    jj = i:min(i+blk-1, size(Y,1));
    d2 = nx2 - 2*Y(jj,:)*X.' + sum(Y(jj,:).^2,2);
    [~, idx(jj)] = min(d2, [], 2);
end
end

% ------------------------------------------------------------------------
function [s,qsc,qcu,Tcs,cert] = makeThermalGrid7(Nx,L,sLoss,qscLoss,qcuLoss,sTcs,Tcs0,cert0)
% Costruisce una griglia termica dispari e interpola SEMPRE dai dati nativi,
% evitando interpolazioni successive fra livelli di raffinamento.
if mod(Nx,2)==0
    error('makeThermalGrid7:mesh','Nx deve essere dispari per avere il nodo centrale.');
end
s = linspace(0,L,Nx).';
qsc = max(interp1(sLoss,qscLoss,s,'linear','extrap'),0);
qcu = max(interp1(sLoss,qcuLoss,s,'linear','extrap'),0);
Tcs = interp1(sTcs,Tcs0,s,'linear','extrap');
cert = interp1(sTcs,double(cert0),s,'nearest','extrap') > 0.5;
end

% ------------------------------------------------------------------------
function P = evaluateThermalPoint7(mdot,Tin,cool,s,qsc,qcu,Tcs,t_ramp, ...
                                   A_Cu,A_SC,A_SS,P_CuSS,dT_req,wantHist)
% Valuta un singolo punto operativo (mdot,Tin).
c = cool;
c.T_in = Tin;
R = slowCoolingSolve7(mdot,c,s,qsc,qcu,t_ramp,A_Cu,A_SC,A_SS,P_CuSS,wantHist);

P.R = R;
P.mdot = mdot;
P.Tin = Tin;
P.solverOK = R.ok;

if R.ok
    mprof = Tcs - R.Tc_max_t;
    P.margin = min(mprof);
    P.dp = max(R.dpA,R.dpB);
    P.p_in = c.p_out + P.dp;
    P.Tscmax = max(R.Tc_max_t);
    P.Tcumax = max(R.Tcu_max_t);
    P.Themax = max(R.THe_max_t);
    P.feasible = P.margin >= dT_req && ...
                 P.p_in <= c.p_range(2) && c.p_out >= c.p_range(1);
else
    P.margin = -Inf;
    P.dp = Inf;
    P.p_in = Inf;
    P.Tscmax = Inf;
    P.Tcumax = Inf;
    P.Themax = Inf;
    P.feasible = false;
end
end

% ------------------------------------------------------------------------
function F = findMinFlow7(cool,s,qsc,qcu,Tcs,t_ramp, ...
                          A_Cu,A_SC,A_SS,P_CuSS,dT_req,opt)
% Bisezione sulla portata. Cerca una parentesi [fail,pass] e la restringe
% fino alla tolleranza cool.flow_tol. Se e' disponibile una stima precedente
% la usa per limitare il numero di transitori alle griglie piu' fini.
if nargin < 12 || isempty(opt), opt=struct; end
if ~isfield(opt,'wantHist'), opt.wantHist=false; end
if ~isfield(opt,'verbose'), opt.verbose=false; end
if ~isfield(opt,'guess'), opt.guess=NaN; end
if ~isfield(opt,'halfWidth'), opt.halfWidth=NaN; end

mdFloor = cool.mdot_floor;
mdMax = cool.mdot_max;

if isfinite(opt.guess)
    if isfinite(opt.halfWidth)
        hw = opt.halfWidth;
    else
        hw = max(0.5e-3,0.2*opt.guess);
    end
    mdLo = max(mdFloor,opt.guess-hw);
    mdHi = min(mdMax,opt.guess+hw);
else
    mdLo = mdFloor;
    mdHi = mdMax;
end

PLo = evaluateThermalPoint7(mdLo,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
    A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);
PHi = evaluateThermalPoint7(mdHi,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
    A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);

% Espansione della parentesi se la stima iniziale non separa fail/pass.
step = max(mdHi-mdLo,0.5e-3);
while PLo.feasible && mdLo > mdFloor+1e-12
    mdHi = mdLo;
    PHi = PLo;
    mdLo = max(mdFloor,mdLo-step);
    step = 1.5*step;
    PLo = evaluateThermalPoint7(mdLo,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
        A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);
end
step = max(mdHi-mdLo,0.5e-3);
while ~PHi.feasible && mdHi < mdMax-1e-12
    mdLo = mdHi;
    PLo = PHi;
    mdHi = min(mdMax,mdHi+step);
    step = 1.5*step;
    PHi = evaluateThermalPoint7(mdHi,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
        A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);
end

testedM = [mdLo;mdHi];
testedMar = [PLo.margin;PHi.margin];

if ~PHi.feasible
    % Nemmeno la portata massima passa.
    Pfinal = evaluateThermalPoint7(mdMax,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
        A_Cu,A_SC,A_SS,P_CuSS,dT_req,opt.wantHist);
    F = packFlowResult7(Pfinal,testedM,testedMar,false);
    return
end

if PLo.feasible
    % Anche il floor numerico passa: il vero minimo e' sotto il dominio cercato.
    Pfinal = evaluateThermalPoint7(mdLo,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
        A_Cu,A_SC,A_SS,P_CuSS,dT_req,opt.wantHist);
    F = packFlowResult7(Pfinal,testedM,testedMar,true);
    F.lowerBoundActive = true;
    return
end

for ib = 1:cool.max_bisect
    if mdHi-mdLo <= cool.flow_tol
        break
    end
    md = 0.5*(mdLo+mdHi);
    Pm = evaluateThermalPoint7(md,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
        A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);
    testedM(end+1,1)=md; %#ok<AGROW>
    testedMar(end+1,1)=Pm.margin; %#ok<AGROW>
    if opt.verbose
        fprintf('    bisezione %2d: mdot=%7.4f g/s -> margine=%7.4f K\n', ...
            ib,md*1e3,Pm.margin);
    end
    if Pm.feasible
        mdHi=md; PHi=Pm;
    else
        mdLo=md; PLo=Pm;
    end
end

Pfinal = evaluateThermalPoint7(mdHi,cool.T_in,cool,s,qsc,qcu,Tcs,t_ramp, ...
    A_Cu,A_SC,A_SS,P_CuSS,dT_req,opt.wantHist);
F = packFlowResult7(Pfinal,testedM,testedMar,true);
F.lowerBoundActive = false;
end

% ------------------------------------------------------------------------
function F = packFlowResult7(P,testedM,testedMar,feasible)
F.mdot = P.mdot;
F.margin = P.margin;
F.dp = P.dp;
F.p_in = P.p_in;
F.Tscmax = P.Tscmax;
F.Tcumax = P.Tcumax;
F.Themax = P.Themax;
F.R = P.R;
F.solverOK = P.solverOK;
F.feasible = feasible && P.feasible;
F.mdot_tested = testedM;
F.margin_tested = testedMar;
F.lowerBoundActive = false;
end

% ------------------------------------------------------------------------
function E = findMaxTin7(mdot,cool,s,qsc,qcu,Tcs,t_ramp, ...
                         A_Cu,A_SC,A_SS,P_CuSS,dT_req,Tlo,Thi)
% Per una portata fissata trova la massima T_in che mantiene 5 K di margine.
% La funzione e' monotona nel range operativo: aumentando T_in il margine
% diminuisce. La bisezione restituisce quindi la frontiera termica.
Plo = evaluateThermalPoint7(mdot,Tlo,cool,s,qsc,qcu,Tcs,t_ramp, ...
    A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);

E.feasibleAtMinTin = Plo.feasible;
if ~Plo.feasible
    E.Tin = NaN;
    E.margin = Plo.margin;
    E.p_in = Plo.p_in;
    E.Tscmax = Plo.Tscmax;
    E.R = [];
    return
end

Thi = max(Thi,Tlo);
Phi = evaluateThermalPoint7(mdot,Thi,cool,s,qsc,qcu,Tcs,t_ramp, ...
    A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);

if Phi.feasible
    E.Tin = Thi;
    E.margin = Phi.margin;
    E.p_in = Phi.p_in;
    E.Tscmax = Phi.Tscmax;
    E.R = [];
    return
end

lo=Tlo; hi=Thi; Pbest=Plo;
for ib=1:cool.max_bisect
    if hi-lo <= cool.Tin_tol
        break
    end
    Tm=0.5*(lo+hi);
    Pm=evaluateThermalPoint7(mdot,Tm,cool,s,qsc,qcu,Tcs,t_ramp, ...
        A_Cu,A_SC,A_SS,P_CuSS,dT_req,false);
    if Pm.feasible
        lo=Tm; Pbest=Pm;
    else
        hi=Tm;
    end
end

E.Tin = lo;
E.margin = Pbest.margin;
E.p_in = Pbest.p_in;
E.Tscmax = Pbest.Tscmax;
E.R = [];
end

% ------------------------------------------------------------------------
function R = slowCoolingSolve7(mdot, cool, s, q_sc_t, q_cu_t, t_ramp, ...
                               A_Cu, A_SC, A_SS, P_CuSS, wantHist)
% Due rami in parallelo con ingresso comune al nodo centrale.
% Ogni ramo porta mdot/2. Il solver termico di ramo e' identico; cambiano
% solamente il profilo di sorgente e, di poco, le proprieta' locali.
N = numel(s);
if mod(N,2) == 0
    error('slowCoolingSolve7:mesh','La griglia globale deve avere un numero dispari di nodi.');
end
iC = (N+1)/2;
idxA = (iC:-1:1).';
idxB = (iC:N).';
xA = s(iC)-s(idxA);
xB = s(idxB)-s(iC);

RA = slowBranchSolve7(mdot/2,cool,xA,q_sc_t(idxA,:),q_cu_t(idxA,:), ...
                      t_ramp,A_Cu,A_SC,A_SS,P_CuSS,wantHist);
RB = slowBranchSolve7(mdot/2,cool,xB,q_sc_t(idxB,:),q_cu_t(idxB,:), ...
                      t_ramp,A_Cu,A_SC,A_SS,P_CuSS,wantHist);

R.ok = RA.ok && RB.ok;
fieldsFinal = {'T_sc','T_c','T_j','T_He','p','v','Re','hHe','rho'};
for k=1:numel(fieldsFinal)
    fn=fieldsFinal{k};
    R.(fn)=nan(N,1);
    R.(fn)(idxA)=RA.(fn);
    R.(fn)(idxB)=RB.(fn);
    R.(fn)(iC)=0.5*(RA.(fn)(1)+RB.(fn)(1));
end

fieldsMax = {'Tc_max_t','Tcu_max_t','Tj_max_t','THe_max_t'};
for k=1:numel(fieldsMax)
    fn=fieldsMax{k};
    R.(fn)=nan(N,1);
    R.(fn)(idxA)=RA.(fn);
    R.(fn)(idxB)=RB.(fn);
    R.(fn)(iC)=max(RA.(fn)(1),RB.(fn)(1));
end
R.T_Cu = R.T_c;
R.dpA = RA.dp_max;
R.dpB = RB.dp_max;
R.dp = max(R.dpA,R.dpB);
R.dp_imbal = abs(R.dpA-R.dpB);
R.t_transit = max(RA.t_transit,RB.t_transit);

if wantHist
    nt=min(numel(RA.t_hist),numel(RB.t_hist));
    R.t_hist=RA.t_hist(1:nt);
    R.THe_hist=max(RA.THe_hist(1:nt),RB.THe_hist(1:nt));
    R.Tc_hist=max(RA.Tc_hist(1:nt),RB.Tc_hist(1:nt));
else
    R.t_hist=[]; R.THe_hist=[]; R.Tc_hist=[];
end
end

% ------------------------------------------------------------------------
function R = slowBranchSolve7(mdot, cool, x, q_sc_t, q_cu_t, t_ramp, ...
                              A_Cu, A_SC, A_SS, P_CuSS, wantHist)
% Modello energetico 1D di un ramo, scritto con la stessa struttura
% numerica dell'esame.
%
% GRIGLIA:
%   xx = discretizzazione spaziale del ramo
%   nn = numero di nodi spaziali
%   nvar = 4 variabili per nodo
%
% Il vettore delle incognite e' INTERLACCIATO nello spazio:
%
%   Tm = [ Tj_1  TCu_1  Tsc_1  THe_1 ...
%          Tj_2  TCu_2  Tsc_2  THe_2 ...
%          ...
%          Tj_nn TCu_nn Tsc_nn THe_nn ]'
%
% quindi, al nodo i:
%   Ti     = jacket
%   Ti + 1 = rame
%   Ti + 2 = stack superconduttore
%   Ti + 3 = elio
% e il nodo spaziale successivo inizia quattro posizioni dopo.
%
% SPAZIO:
%   - jacket, rame, SC: differenze centrate per la conduzione assiale;
%   - He: differenza all'indietro (upwind) per il termine convettivo,
%         perche' la portata e' positiva dall'ingresso all'uscita.
%
% TEMPO:
%   Eulero implicito (Backward Euler), primo ordine:
%
%       Tmpiu1 = AA1 \ (AA2*Tm + bb)
%
% Le proprieta' termofisiche sono valutate al passo precedente e congelate
% durante il singolo solve lineare (linearizzazione semi-implicita).
%
% Rispetto all'esame cambia solo la fisica: qui ci sono quattro temperature,
% sorgenti AC tempo-dipendenti e proprieta' criogeniche dipendenti da T.

xx = x(:);
nn = length(xx);
nvar = 4;
ntot = nvar*nn;

if nn < 2
    error('slowBranchSolve7:grid','Servono almeno due nodi per ramo.');
end

dx = xx(2)-xx(1);
if max(abs(diff(xx)-dx)) > 1e-10*max(xx(end),1)
    error('slowBranchSolve7:grid','La griglia di ramo deve essere uniforme.');
end

if size(q_sc_t,1) ~= nn || size(q_cu_t,1) ~= nn
    error('slowBranchSolve7:source','Dimensione delle sorgenti incompatibile con la griglia.');
end

dt = cool.dt;
nStep = ceil(cool.t_end/dt);
tt = (0:nStep)'*dt;

% Indici nel vettore interlacciato.
ij  = (1:4:ntot)';      % jacket
icu = (2:4:ntot)';      % rame
isc = (3:4:ntot)';      % stack SC
ihe = (4:4:ntot)';      % elio

% Stato iniziale. Prima della rampa non ci sono sorgenti: T = T_in uniforme
% e' gia' la soluzione iniziale coerente.
Tm = cool.T_in*ones(ntot,1);

Tj  = Tm(ij);
Tcu = Tm(icu);
Tsc = Tm(isc);
The = Tm(ihe);

Tj_max  = Tj;
Tcu_max = Tcu;
Tsc_max = Tsc;
The_max = The;

dp_max = 0;
R.ok = true;

if wantHist
    histHe = nan(nStep+1,1);
    histSc = nan(nStep+1,1);
    histHe(1) = max(The);
    histSc(1) = max(Tsc);
else
    histHe = [];
    histSc = [];
end

for iter = 1:nStep
    t = tt(iter+1);

    if max(Tm) >= 49.5
        R.ok = false;
        break
    end

    % -------------------------------------------------------------
    % Sorgenti AC al tempo corrente
    % -------------------------------------------------------------
    if t <= t_ramp
        lam = max(0,min(1,t/t_ramp));
        qsc = interp1(cool.lam_grid,q_sc_t.',lam,'linear','extrap').';
        qcu = interp1(cool.lam_grid,q_cu_t.',lam,'linear','extrap').';
    else
        qsc = zeros(nn,1);
        qcu = zeros(nn,1);
    end

    % -------------------------------------------------------------
    % Proprieta' al passo precedente
    % -------------------------------------------------------------
    cp_j = cpSS(Tj);
    k_j  = kSS(Tj);
    Cj   = 7900*cp_j*A_SS;

    cp_cu = cpCu(Tcu);
    k_cu  = kCu(Tcu,cool.MR_Cu,cool.RRR);
    Ccu   = 8960*cp_cu*A_Cu;

    [rhoCp_sc,k_sc] = stackVolProps7(Tsc,cool);
    Csc = rhoCp_sc*A_SC;

    [rhoHe,~,cpHe,~,muHe,kHe] = heProps(The,cool.p_out*ones(nn,1));
    [fD,Nu] = tubeCorrelation7(mdot,cool.A_He,cool.D_h,cpHe,muHe,kHe);
    hHe = Nu.*kHe/cool.D_h;
    H = hHe*cool.P_wet;
    Che = rhoHe.*cpHe*cool.A_He;
    v = mdot./(rhoHe*cool.A_He);

    % -------------------------------------------------------------
    % Conduttanze trasverse alle interfacce fra materiali
    % -------------------------------------------------------------
    % Per ogni interfaccia si usa una resistenza termica in serie per unita'
    % di lunghezza. Con R_contact = 0:
    %
    %   q'_(SC-Cu) = Gsc*(Tsc - Tcu)
    %   Gsc = P_sc / (ell_sc/k_sc + ell_cu_in/k_cu)
    %
    %   q'_(Cu-j) = Gcj*(Tcu - Tj)
    %   Gcj = P_CuSS / (ell_cu_out/k_cu + ell_j/k_j)
    %
    % Per il fluido:
    %
    %   q'_(Cu-He) = H*(Tcu - The),    H = hHe*P_wet
    %
    % Lo stesso flusso entra con segno opposto nelle due equazioni accoppiate,
    % quindi lo scambio interno conserva energia.
    Gsc = cool.P_sc ./ ...
        (cool.ell_sc./max(k_sc,eps) + cool.ell_cu_in./max(k_cu,eps));

    Gcj = P_CuSS ./ ...
        (cool.ell_cu_out./max(k_cu,eps) + cool.ell_j./max(k_j,eps));

    % -------------------------------------------------------------
    % Coefficienti normalizzati, analoghi agli aa/beta/cc dell'esame
    % -------------------------------------------------------------
    [aPj,aWj,aEj]    = axialCoeffExam7(k_j*A_SS, Cj,  dt,dx);
    [aPcu,aWcu,aEcu] = axialCoeffExam7(k_cu*A_Cu,Ccu,dt,dx);
    [aPsc,aWsc,aEsc] = axialCoeffExam7(k_sc*A_SC,Csc,dt,dx);

    % scambio jacket <-> rame
    beta_j_cu  = dt*Gcj./Cj;
    beta_cu_j  = dt*Gcj./Ccu;

    % scambio rame <-> SC
    beta_cu_sc = dt*Gsc./Ccu;
    beta_sc_cu = dt*Gsc./Csc;

    % scambio rame <-> He
    beta_cu_he = dt*H./Ccu;
    beta_he_cu = dt*H./Che;

    % coefficiente convettivo upwind dell'elio
    cc = dt*(mdot*cpHe/dx)./Che;

    % -------------------------------------------------------------
    % Matrici del transitorio: stessa logica dell'esame
    % -------------------------------------------------------------
    AA1 = spalloc(ntot,ntot,18*nn);
    AA2 = speye(ntot);
    bb  = zeros(ntot,1);

    % Ogni ciclo costruisce le 4 equazioni associate allo stesso x_i.
    for i = 1:nn
        rj  = 4*(i-1) + 1;
        rcu = rj + 1;
        rsc = rj + 2;
        rhe = rj + 3;

        % ---------------- jacket ----------------
        % conduzione centrata + scambio con rame
        AA1(rj,rj)  = 1 + aPj(i) + beta_j_cu(i);
        AA1(rj,rcu) = -beta_j_cu(i);

        if i > 1
            AA1(rj,rj-4) = -aWj(i);
        end
        if i < nn
            AA1(rj,rj+4) = -aEj(i);
        end

        % ---------------- rame ----------------
        % conduzione centrata + scambio con jacket, SC ed He
        AA1(rcu,rcu) = 1 + aPcu(i) + beta_cu_j(i) + ...
                       beta_cu_sc(i) + beta_cu_he(i);
        AA1(rcu,rj)  = -beta_cu_j(i);
        AA1(rcu,rsc) = -beta_cu_sc(i);
        AA1(rcu,rhe) = -beta_cu_he(i);

        if i > 1
            AA1(rcu,rcu-4) = -aWcu(i);
        end
        if i < nn
            AA1(rcu,rcu+4) = -aEcu(i);
        end

        bb(rcu) = dt*qcu(i)/Ccu(i);

        % ---------------- stack SC ----------------
        % conduzione centrata + scambio con rame
        AA1(rsc,rsc) = 1 + aPsc(i) + beta_sc_cu(i);
        AA1(rsc,rcu) = -beta_sc_cu(i);

        if i > 1
            AA1(rsc,rsc-4) = -aWsc(i);
        end
        if i < nn
            AA1(rsc,rsc+4) = -aEsc(i);
        end

        bb(rsc) = dt*qsc(i)/Csc(i);

        % ---------------- elio ----------------
        % accumulo + upwind all'indietro + scambio con il rame
        AA1(rhe,rhe) = 1 + cc(i) + beta_he_cu(i);
        AA1(rhe,rcu) = -beta_he_cu(i);

        if i > 1
            AA1(rhe,rhe-4) = -cc(i);
        end
    end

    % -------------------------------------------------------------
    % Condizione al contorno dell'elio all'ingresso
    % -------------------------------------------------------------
    % T_He(x=0) = T_in. Come nell'esame, la riga viene rimossa sia
    % dalla matrice del passo nuovo sia da quella del passo precedente.
    AA1(ihe(1),:) = 0;
    AA2(ihe(1),:) = 0;
    AA1(ihe(1),ihe(1)) = 1;
    bb(ihe(1)) = cool.T_in;

    % -------------------------------------------------------------
    % Soluzione implicita
    % -------------------------------------------------------------
    Tmpiu1 = AA1\(AA2*Tm + bb);

    if any(~isfinite(Tmpiu1)) || min(Tmpiu1) < 2 || max(Tmpiu1) > 60
        R.ok = false;
        break
    end

    % Diagnostica, analoga all'esame. Non arresta il transitorio perche'
    % qui la durata fisica e' fissata dalla rampa e dalla coda post-rampa.
    err = norm(Tmpiu1-Tm)/max(norm(Tmpiu1-cool.T_in),eps); %#ok<NASGU>

    Tm = Tmpiu1;

    % Estrazione delle quattro variabili dalla griglia interlacciata.
    Tj  = Tm(1:4:ntot);
    Tcu = Tm(2:4:ntot);
    Tsc = Tm(3:4:ntot);
    The = Tm(4:4:ntot);

    Tj_max  = max(Tj_max,Tj);
    Tcu_max = max(Tcu_max,Tcu);
    Tsc_max = max(Tsc_max,Tsc);
    The_max = max(The_max,The);

    % -------------------------------------------------------------
    % Idraulica quasi-stazionaria
    % -------------------------------------------------------------
    gradp = fD.*mdot^2./(2*rhoHe*cool.A_He^2*cool.D_h);
    dp_now = trapz(xx,gradp);
    dp_max = max(dp_max,dp_now);

    if wantHist
        histHe(iter+1) = max(The);
        histSc(iter+1) = max(Tsc);
    end
end

% -----------------------------------------------------------------
% Proprieta' finali per output e grafici
% -----------------------------------------------------------------
[rhoHe,~,cpHe,~,muHe,kHe] = heProps(The,cool.p_out*ones(nn,1)); %#ok<ASGLU>
[fD,Nu] = tubeCorrelation7(mdot,cool.A_He,cool.D_h,cpHe,muHe,kHe);

hHe = Nu.*kHe/cool.D_h;
v = mdot./(rhoHe*cool.A_He);

gradp = fD.*mdot^2./(2*rhoHe*cool.A_He^2*cool.D_h);
dpCum = cumtrapz(xx,gradp);
dpFinal = dpCum(end);
p = cool.p_out + (dpFinal-dpCum);

R.T_j  = Tj;
R.T_c  = Tcu;
R.T_sc = Tsc;
R.T_He = The;

R.Tj_max_t  = Tj_max;
R.Tcu_max_t = Tcu_max;
R.Tc_max_t  = Tsc_max;
R.THe_max_t = The_max;

R.p = p;
R.v = v;
R.rho = rhoHe;
R.Re = mdot*cool.D_h./(cool.A_He*muHe);
R.hHe = hHe;
R.dp_max = max(dp_max,dpFinal);
R.t_transit = trapz(xx,1./max(v,eps));

if wantHist
    nDone = find(~isnan(histHe),1,'last');
    R.t_hist = tt(1:nDone);
    R.THe_hist = histHe(1:nDone);
    R.Tc_hist = histSc(1:nDone);
else
    R.t_hist = [];
    R.THe_hist = [];
    R.Tc_hist = [];
end
end

% ------------------------------------------------------------------------
function [aP,aW,aE] = axialCoeffExam7(kA,C,dt,dx)
% Coefficienti dimensionless della conduzione assiale dopo aver moltiplicato
% l'equazione energetica per dt/C. Agli estremi si usa mezzo volume di
% controllo, cioe' flusso nullo senza eliminare l'equazione dinamica del nodo.
kA = kA(:);
C = C(:);
N = numel(kA);

aW = zeros(N,1);
aE = zeros(N,1);

if N == 1
    aP = zeros(N,1);
    return
end

kf = 2*kA(1:end-1).*kA(2:end)./max(kA(1:end-1)+kA(2:end),eps);

aE(1) = 2*dt*kf(1)/(C(1)*dx^2);
aW(end) = 2*dt*kf(end)/(C(end)*dx^2);

if N > 2
    aW(2:end-1) = dt*kf(1:end-1)./(C(2:end-1)*dx^2);
    aE(2:end-1) = dt*kf(2:end)./(C(2:end-1)*dx^2);
end

aP = aW + aE;
end

% ------------------------------------------------------------------------
function [fD,Nu] = tubeCorrelation7(mdot,A,Dh,cp,mu,k)
% Tubo liscio: laminare 64/Re e Nu=4.36; turbolento Petukhov/Gnielinski;
% interpolazione lineare 2300-4000 per evitare discontinuita' numeriche.
Re=max(mdot*Dh./(A*mu),1);
Pr=max(cp.*mu./max(k,eps),1e-6);
fLam=64./Re;
NuLam=4.36*ones(size(Re));
fTur=(0.79*log(max(Re,3000))-1.64).^(-2);
NuTur=(fTur/8).*(Re-1000).*Pr ./ ...
      max(1+12.7*sqrt(fTur/8).*(Pr.^(2/3)-1),0.1);
NuTur=max(NuTur,4.36);

fD=fLam; Nu=NuLam;
mTur=Re>=4000;
fD(mTur)=fTur(mTur); Nu(mTur)=NuTur(mTur);
mTr=Re>2300 & Re<4000;
w=(Re(mTr)-2300)/(4000-2300);
fD(mTr)=(1-w).*fLam(mTr)+w.*fTur(mTr);
Nu(mTr)=(1-w).*NuLam(mTr)+w.*NuTur(mTr);
end

% ------------------------------------------------------------------------
function [rhoCp,kst] = stackVolProps7(T,cool)
% Proprietà volumetriche equivalenti del tape stack.
% Ag ~= Cu; REBCO+buffer ~= substrato Hastelloy per cp e k, per mancanza di
% correlazioni criogeniche dedicate nel materiale fornito al progetto.
cpCuT=cpCu(T);
kCuT =kCu(T,cool.MR_Cu,cool.RRR);
cpBase=cpHast(T);
kBase =kHast(T);
rhoCp = cool.f_Cu_tape*8960.*cpCuT + ...
        cool.f_Ag_tape*10490.*cpCuT + ...
        cool.f_base_tape*cool.rho_base_tape.*cpBase;
kst = (cool.f_Cu_tape+cool.f_Ag_tape).*kCuT + cool.f_base_tape.*kBase;
end

% ------------------------------------------------------------------------
function [Tmap, X, Y, info] = sectionThermal2D(geo, kmap, src, hHe, T_He, n)
% Campo di temperatura 2D sulla SEZIONE del cavo, ricostruito a valle del
% modello 1D. Serve a vedere come si distribuisce il calore dentro la sezione
% e a controllare visivamente l'ipotesi 1D: se la mappa e' quasi isoterma
% rispetto al salto lungo il canale, comprimere la sezione in nodi
% concentrati era legittimo.
%
% NON e' una seconda soluzione del problema: e' un POST-PROCESSING coerente
% col modello 1D. Prende dal modello la temperatura di bulk dell'elio e il
% coefficiente di scambio nella stazione scelta, prende dallo Step 6 le
% sorgenti locali, e risolve la conduzione stazionaria nella sezione:
%
%     div(k grad T) + q''' = 0     nei metalli e negli stack
%     -k dT/dn = hHe (T - T_He)    sulla parete del canale   (Robin)
%     dT/dn = 0                    sulla superficie esterna  (adiabatica)
%
% Stazionario e non transitorio perche' nella sezione i tempi di diffusione
% sono di millisecondi contro i 250 s della rampa: la sezione insegue.
%
% ANISOTROPIA DELLO STACK, dichiarata: nel modello 1D il calore esce dai
% BORDI dei nastri (le 59 interfacce in vuoto chiudono il percorso
% trasverso). Qui lo stack riceve la conducibilita' IN-PIANO in entrambe le
% direzioni, quindi la mappa mostra lo stack piu' isotermo del vero. Il
% salto che il modello 1D calcola fra stack e former resta il riferimento
% quantitativo; questa mappa serve alla lettura d'insieme.
if nargin < 6 || isempty(n), n = 241; end
x  = linspace(-geo.r_cab, geo.r_cab, n);
h  = x(2) - x(1);
[X, Y] = ndgrid(x, x);
R  = hypot(X, Y);
TH = atan2(Y, X);

%  mappa delle conducibilita' e delle sorgenti volumetriche [W/m^3]
kf  = zeros(n);   qv = zeros(n);
jac = (R <= geo.r_cab)   & (R > geo.r_form_o);
cu  = (R <= geo.r_form_o) & (R > geo.r_chan);
kf(jac) = kmap.k_SS;    kf(cu) = kmap.k_Cu;
qv(cu)  = src.q_cu/kmap.A_Cu;              % accoppiamento -> former
%  GLI SLOT SONO RETTANGOLI RUOTATI, non settori anulari: e' la geometria
%  vera del cavo, la stessa disegnata in 5.7. Ogni slot ha lato radiale
%  slot_d e lato azimutale slot_w, centrato sul centroide R_core e ruotato
%  dell'angolo dello slot. Usare un settore anulare cambierebbe l'area e la
%  forma dei web, quindi anche il campo termico.
slot = false(n);
for k = 0:geo.N_slots-1
    a  = 2*pi*k/geo.N_slots;
    xc = geo.R_core*cos(a);   yc = geo.R_core*sin(a);
    u  =  (X-xc)*cos(a) + (Y-yc)*sin(a);      % coordinata RADIALE locale
    w  = -(X-xc)*sin(a) + (Y-yc)*cos(a);      % coordinata AZIMUTALE locale
    slot = slot | (abs(u) <= geo.slot_d/2 & abs(w) <= geo.slot_w/2);
end
kf(slot) = kmap.k_st;
qv(slot) = src.q_sc/kmap.A_SC;             % magnetizzazione -> stack
kf(R <= geo.r_chan) = 0;                   % canale: non e' un solido

solid = kf > 0;
N   = nnz(solid);
idx = zeros(n);  idx(solid) = 1:N;
[ii, jj] = find(solid);
lc  = sub2ind([n n], ii, jj);
p   = (1:N).';
shift = [1 0; -1 0; 0 1; 0 -1];

%  PRIMO PASSAGGIO: quante facce guardano il canale. Serve perche' un cerchio
%  discretizzato a scalini ha un perimetro piu' lungo di quello vero (il
%  fattore e' 4/pi per una circonferenza), e senza correzione l'area di
%  scambio risulterebbe sovrastimata del ~27%, con un salto di film
%  altrettanto sottostimato. Si riscala la conduttanza di ogni faccia in modo
%  che la somma delle aree riproduca il perimetro bagnato VERO.
isChanFace = false(N,4);
for d = 1:4
    ai = ii + shift(d,1);  aj = jj + shift(d,2);
    ok = ai >= 1 & ai <= n & aj >= 1 & aj <= n;
    la = ones(N,1);  la(ok) = sub2ind([n n], ai(ok), aj(ok));
    isChanFace(:,d) = ok & ~solid(la) & (hypot(X(la), Y(la)) <= geo.r_chan);
end
nRob = nnz(isChanFace);
if nRob == 0
    error('sect2D:noWall','Nessuna faccia sul canale: griglia troppo grossolana.');
end
fcor = (2*pi*geo.r_chan)/(nRob*h);       % correzione di perimetro a scalini

dg  = zeros(N,1);
b   = -qv(lc)*h^2;                       % sorgente: div(k grad T) = -q'''
I = [];  J = [];  V = [];
for d = 1:4
    ai = ii + shift(d,1);  aj = jj + shift(d,2);
    ok = ai >= 1 & ai <= n & aj >= 1 & aj <= n;
    la = ones(N,1);  la(ok) = sub2ind([n n], ai(ok), aj(ok));
    isS = ok & solid(la);
    %  facce solido-solido: conduttanza di faccia = media armonica di k
    %  (per unita' di lunghezza del cavo vale k, perche' area h / distanza h)
    kP  = kf(lc(isS));   kN = kf(la(isS));
    kfz = 2*kP.*kN./(kP + kN);
    dg(isS) = dg(isS) - kfz;
    I = [I; p(isS)];  J = [J; idx(la(isS))];  V = [V; kfz];        %#ok<AGROW>
    %  facce sul CANALE: Robin, film in serie con la mezza cella di solido.
    %  G = 1/(1/(h_He*dx) + 1/(2k))  [W/(m K)], gia' un'AREA: non va
    %  rimoltiplicata per h. Le facce verso l'esterno restano adiabatiche,
    %  quindi non contribuiscono affatto.
    rob = isChanFace(:,d);
    if any(rob)
        G = fcor ./ (1/(hHe*h) + 1./(2*kf(lc(rob))));
        dg(rob) = dg(rob) - G;
        b(rob)  = b(rob)  - G*T_He;
    end
end
I = [I; p];  J = [J; p];  V = [V; dg];
A = sparse(I, J, V, N, N);
Tv = A\b;

Tmap = nan(n);
Tmap(solid) = Tv;
%  VERIFICA INTERNA: in stazionario tutta la potenza immessa esce dal canale,
%  quindi il salto medio parete-fluido deve valere P/(h_He*P_bagnato). Se il
%  conto non torna, l'errore e' nell'assemblaggio, non nella fisica.
dT_ref = (src.q_cu + src.q_sc)/(hHe*2*pi*geo.r_chan);
dT_num = min(Tv) - T_He;
if abs(dT_num - dT_ref) > 0.1*max(dT_ref, eps)
    warning('sect2D:balance', ...
        ['Salto di film numerico %.4f K contro %.4f K analitico: ' ...
         'scarto oltre il 10%%, controllare la griglia.'], dT_num, dT_ref);
end
info = struct('T_min', min(Tv), 'T_max', max(Tv), 'n', n, ...
              'T_stack', mean(Tv(idx(slot & solid))), ...
              'T_wall',  min(Tv), 'dT_film', dT_num, 'dT_film_ref', dT_ref, ...
              'perim_corr', fcor);
end

% ------------------------------------------------------------------------
function [gf, info] = transverseGeomFactor(geo, n, mode)
% Fattore geometrico del percorso trasverso nella sezione del CICC.
%
% COSA CALCOLA. Le correnti di accoppiamento indotte da dB/dt si chiudono
% ATTRAVERSO il former: non seguono una retta, devono aggirare i sei slot
% occupati dagli stack (trasversalmente isolanti, 59 interfacce in vuoto) e
% il canale centrale dell'elio. Si omogeneizza quindi la sezione: si impone
% un gradiente di potenziale uniforme sul bordo esterno del former,
%     div(sigma grad phi) = 0 ,   phi = E*x  su  r = r_former_out
% si ricava sigma_eff = <J_x>/E mediata sull'area, e si definisce
%     geomFac = sigma_Cu / sigma_eff
% cioe' di quanto la sezione reale e' piu' resistiva di un disco di rame
% pieno dello stesso raggio. Il fattore contiene INSIEME la riduzione di
% area e l'allungamento del percorso, che e' esattamente cio' che moltiplica
% rho_Cu in rho_eff (6.5).
%
% COSA NON CALCOLA. Non c'entra il current sharing: quello e' il passaggio
% di corrente di TRASPORTO dal nastro al rame quando T supera T_cs, e in
% questo progetto non avviene mai. Le correnti di accoppiamento esistono
% anche a margine intatto, perche' sono indotte dalla variazione di campo.
%
% VALIDAZIONE INTERNA: la stessa routine girata su un disco di rame pieno
% deve restituire sigma_eff/sigma_Cu = 1. Lo scarto e' in info.ref_err ed e'
% il controllo che la discretizzazione non stia introducendo un bias.
%
% geo: struct con r_form_o, r_chan, r_slot_in, r_slot_out, slot_w, N_slots
% n  : lato della griglia quadrata (default 201; il risultato converge da
%      sopra, con n = 201 -> 401 il valore si muove dello 0.8%)
%
% mode: 'homog' (default) e' l'omogeneizzazione descritta sopra, la sola che
%       entra in rho_eff. 'pair' e' una DIAGNOSTICA: resistenza fra due slot
%       ADIACENTI presi come elettrodi, rapportata alla stessa sezione a rame
%       pieno. Dice quanto e' diretto il percorso fra due strand vicini, ma
%       non e' una resistivita' e non va sostituita a geomFac.
if nargin < 2 || isempty(n), n = 201; end
if nargin < 3 || isempty(mode), mode = 'homog'; end
sig0 = 1e-9;                       % "isolante": stack e canale elio

x  = linspace(-geo.r_form_o, geo.r_form_o, n);
h  = x(2) - x(1);
[X, Y] = ndgrid(x, x);
R  = hypot(X, Y);
TH = atan2(Y, X);

sigRef = sig0*ones(n);  sigRef(R <= geo.r_form_o) = 1;      % disco pieno
sigRea = sigRef;
sigRea(R <= geo.r_chan) = sig0;                            % canale elio
%  Slot = rettangoli ruotati, come la sezione vera (5.7) e come la mappa
%  termica 2D: settori anulari darebbero un'area e dei web diversi.
for k = 0:geo.N_slots-1
    a  = 2*pi*k/geo.N_slots;
    xc = geo.R_core*cos(a);   yc = geo.R_core*sin(a);
    u  =  (X-xc)*cos(a) + (Y-yc)*sin(a);
    w  = -(X-xc)*sin(a) + (Y-yc)*cos(a);
    sigRea(abs(u) <= geo.slot_d/2 & abs(w) <= geo.slot_w/2) = sig0;
end

if strcmpi(mode, 'pair')
    %  Elettrodi: slot 0 e slot 1 (adiacenti), tenuti a +1 e -1. Gli altri
    %  quattro slot restano isolanti nel caso reale. Gli elettrodi devono
    %  essere conduttori, altrimenti la corrente non riesce a uscirne.
    elec = false(n, n, 2);
    for e = 0:1
        a  = 2*pi*e/geo.N_slots;
        xc = geo.R_core*cos(a);   yc = geo.R_core*sin(a);
        u  =  (X-xc)*cos(a) + (Y-yc)*sin(a);
        w  = -(X-xc)*sin(a) + (Y-yc)*cos(a);
        elec(:,:,e+1) = abs(u) <= geo.slot_d/2 & abs(w) <= geo.slot_w/2;
    end
    sA = sigRea;  sA(elec(:,:,1) | elec(:,:,2)) = 1;
    sB = sigRef;  sB(R <= geo.r_chan) = sig0;
    Rreal = pairResistance(sA, elec, R, h, geo.r_form_o);
    Rfull = pairResistance(sB, elec, R, h, geo.r_form_o);
    gf   = Rreal/Rfull;
    info = struct('R_real', Rreal, 'R_full', Rfull, 'n', n, 'mode', 'pair');
    return
end
sRef = homogSigma(sigRef, X, R, h, geo.r_form_o);
sRea = homogSigma(sigRea, X, R, h, geo.r_form_o);
gf   = sRef/sRea;
info = struct('sigma_ref', sRef, 'sigma_real', sRea, ...
              'ref_err', abs(sRef-1), 'n', n, ...
              'area_frac', nnz(sigRea > 0.5)/max(nnz(R <= geo.r_form_o),1));
end

% ------------------------------------------------------------------------
function Rp = pairResistance(sig, elec, R, h, rOut)
% Resistenza fra i due elettrodi di 'elec' (potenziale +1 e -1) attraverso
% la mappa di conducibilita' 'sig'. R = dV/I con dV = 2.
n     = size(sig, 1);
eA    = elec(:,:,1);  eB = elec(:,:,2);
fixed = eA | eB;
phiF  = zeros(n);  phiF(eA) = 1;  phiF(eB) = -1;
unk   = (R <= rOut) & ~fixed;
N     = nnz(unk);
idx   = zeros(n);  idx(unk) = 1:N;
[ii, jj] = find(unk);
lc = sub2ind([n n], ii, jj);
p  = (1:N).';
dg = zeros(N,1);  b = zeros(N,1);
I = [];  J = [];  V = [];
shift = [1 0; -1 0; 0 1; 0 -1];
for d = 1:4
    ai = ii + shift(d,1);  aj = jj + shift(d,2);
    ok = ai >= 1 & ai <= n & aj >= 1 & aj <= n;
    la = sub2ind([n n], ai(ok), aj(ok));
    sf = 2*sig(lc(ok)).*sig(la)./(sig(lc(ok)) + sig(la));
    dg(ok) = dg(ok) - sf;
    isU = unk(la);
    pk  = p(ok);
    I = [I; pk(isU)];  J = [J; idx(la(isU))];  V = [V; sf(isU)];   %#ok<AGROW>
    bi = pk(~isU);
    b(bi) = b(bi) - sf(~isU).*phiF(la(~isU));
end
I = [I; p];  J = [J; p];  V = [V; dg];
A = sparse(I, J, V, N, N);
phi = phiF;
phi(unk) = A\b;

Itot = 0;
[ia, ja] = find(eA);
for q = 1:numel(ia)
    for d = 1:4
        ai = ia(q) + shift(d,1);  aj = ja(q) + shift(d,2);
        if ai < 1 || ai > n || aj < 1 || aj > n, continue; end
        if eA(ai, aj), continue; end
        sf = 2*sig(ia(q),ja(q))*sig(ai,aj)/(sig(ia(q),ja(q)) + sig(ai,aj));
        Itot = Itot + sf*(phi(ia(q),ja(q)) - phi(ai,aj));
    end
end
Rp = 2/max(Itot, eps);
end

% ------------------------------------------------------------------------
function sEff = homogSigma(sig, X, R, h, rOut)
% Risolve div(sigma grad phi) = 0 con phi = X sul bordo e restituisce la
% conducibilita' efficace <J_x>/E (E = 1). Conduttanze di faccia in media
% ARMONICA: e' la forma corretta a volumi finiti quando sigma salta di ordini
% di grandezza fra celle adiacenti.
n      = size(sig, 1);
inside = R <= rOut - h;               % i nodi esterni sono Dirichlet
N      = nnz(inside);
idx    = zeros(n);  idx(inside) = 1:N;
[ii, jj] = find(inside);
lc     = sub2ind([n n], ii, jj);
p      = (1:N).';
dg     = zeros(N,1);
b      = zeros(N,1);
I = [];  J = [];  V = [];
shift = [1 0; -1 0; 0 1; 0 -1];
for d = 1:4
    la  = sub2ind([n n], ii+shift(d,1), jj+shift(d,2));
    sf  = 2*sig(lc).*sig(la)./(sig(lc)+sig(la));
    dg  = dg - sf;
    isI = inside(la);
    I = [I; p(isI)];  J = [J; idx(la(isI))];  V = [V; sf(isI)];   %#ok<AGROW>
    b(~isI) = b(~isI) - sf(~isI).*X(la(~isI));
end
I = [I; p];  J = [J; p];  V = [V; dg];
A = sparse(I, J, V, N, N);
phi = X;
phi(inside) = A\b;

dpdx = zeros(n);
dpdx(2:end-1,:) = (phi(3:end,:) - phi(1:end-2,:))/(2*h);
dpdx(1,:)       = (phi(2,:)   - phi(1,:))/h;
dpdx(end,:)     = (phi(end,:) - phi(end-1,:))/h;
Jx   = -sig.*dpdx;
mask = R <= rOut;
sEff = -mean(Jx(mask));
end

% ------------------------------------------------------------------------
function c = cpCu(T)
% Calore specifico del rame OFHC (UNS C10100/C10200).
% CORRELAZIONE UFFICIALE NIST Cryogenic Material Properties, pagina
% "OFHC Copper" (rev. 02/03/2010), non un fit nostro:
%     log10(cp) = a + b*x + c*x^2 + ... + i*x^8 ,  x = log10(T)
% Validita' dichiarata 4-300 K; errore del fit rispetto ai dati 10% sotto
% 15 K, 5% sopra. Fuori range si clampa agli estremi (dichiarato in 1.x).
co = [-1.91844, -0.15973, 8.61013, -18.996, 21.9661, ...
      -12.7328, 3.54322, -0.3797, 0];
x = log10(min(max(T, 4), 300));
e = zeros(size(x));
for k = 1:numel(co), e = e + co(k)*x.^(k-1); end
c = 10.^e;
end

% ------------------------------------------------------------------------
function k = kCu(T, MR, RRR)
% Conducibilita' termica del rame OFHC a campo NULLO, dalla correlazione
% ufficiale NIST Cryogenic Material Properties, pagina "OFHC Copper"
% (rev. 02/03/2010):
%     log10(k) = (a + c*T^0.5 + e*T + g*T^1.5 + i*T^2)
%                / (1 + b*T^0.5 + d*T + f*T^1.5 + h*T^2)
% con un set di coefficienti per ciascun RRR tabulato (50/100/150/300/500).
% Errore del fit dichiarato da NIST: 1-2%.
%
% MAGNETORESISTENZA: la correlazione e' a campo nullo. In campo la
% resistivita' cresce e, per Wiedemann-Franz, k scende in proporzione:
% k(B) = k_NIST(T,RRR)/MR. MR resta l'unico fattore assunto della catena.
%
if nargin < 2 || isempty(MR),  MR  = 4.0; end
if nargin < 3 || isempty(RRR), RRR = 100; end
RRRtab = [50 100 150 300 500];
C = [ 1.8743 -0.41538 -0.6018   0.13294  0.26426 -0.0219   -0.051276 0.0014871 0.003723
      2.2154 -0.47461 -0.88068  0.13871  0.29505 -0.02043  -0.04831  0.001281  0.003207
      2.3797 -0.4918  -0.98615  0.13942  0.30475 -0.019713 -0.046897 0.0011969 0.0029988
      1.357   0.3981   2.669   -0.1346  -0.6683   0.01342   0.05773  0.0002147 0
      2.8075 -0.54074 -1.2777   0.15362  0.36444 -0.02105  -0.051727 0.0012226 0.0030964];
[~, j] = min(abs(RRRtab - RRR));
a = C(j,1); b = C(j,2); c = C(j,3); d = C(j,4); e = C(j,5);
f = C(j,6); g = C(j,7); h = C(j,8); ii = C(j,9);
Tc = min(max(T, 4), 300);
num = a + c*Tc.^0.5 + e*Tc + g*Tc.^1.5 + ii*Tc.^2;
den = 1 + b*Tc.^0.5 + d*Tc + f*Tc.^1.5 + h*Tc.^2;
k = 10.^(num./den) ./ MR;
end

% ------------------------------------------------------------------------
function k = kHast(T)
% Conducibilita' termica del SUBSTRATO del nastro (Hastelloy C-276, Ni-Mo-Cr).
%
% NIST NON HA l'Hastelloy nel database criogenico (l'indice contiene solo
% acciai, alluminii, Inconel 718, Invar, Ti-6Al-4V, ...). Si usa quindi la
% correlazione ufficiale NIST per l'INCONEL 718 (UNS N07718), stessa forma
% logaritmica, dati 6-275 K, errore del fit 2%.
%
% PERCHE' INCONEL 718 E NON L'ACCIAIO 304, che era il surrogato precedente:
%  - sono entrambe superleghe a base NICHEL, come il C-276, mentre il 304 e'
%    a base ferro;
%  - a 300 K la correlazione NIST dell'Inconel 718 da' 9.79 W/m/K, che
%    coincide col valore di letteratura dell'Hastelloy C-276 (~9.8 W/m/K);
%    il 304 da' 15.3 W/m/K, cioe' il 56% in piu';
%  - a bassa T il trasporto e' elettronico e va come L0*T/rho: la
%    resistivita' a temperatura ambiente dell'Inconel 718 (~1.25 uOhm*m) e'
%    praticamente quella del C-276 (~1.30), mentre quella del 304 e' circa
%    la meta'. Il surrogato giusto e' quindi quello con la resistivita'
%    simile, non quello con la densita' simile.
% RIFERIMENTO PRIMARIO da citare per il C-276 vero: J. Lu, E. S. Choi,
% H. D. Zhou, "Physical properties of Hastelloy C-276 at cryogenic
% temperatures", J. Appl. Phys. 103, 064908 (2008) - misure 2-200 K su
% nastro 4 mm x 0.1 mm, cioe' esattamente il nostro formato.
co = [-8.28921, 39.4470, -83.4353, 98.1690, -67.2088, ...
      26.7082, -5.7205, 0.51115, 0];
x = log10(min(max(T, 4), 300));
e = zeros(size(x));
for kk = 1:numel(co), e = e + co(kk)*x.^(kk-1); end
k = 10.^e;
end

% ------------------------------------------------------------------------
function c = cpHast(T)
% Calore specifico del substrato del nastro, Hastelloy C-276.
%
% FONTE, con i coefficienti PUBBLICATI (non piu' letti dal grafico):
%   J. Lu, E. S. Choi, H. D. Zhou, "Physical properties of Hastelloy C-276
%   at cryogenic temperatures", J. Appl. Phys. 103, 064908 (2008).
% Gli autori riportano  Cp = gamma*T + A*T^3  con
%   gamma = 1.33e-4 J/(g K^2)   ->  0.133   J/(kg K^2)
%   A     = 5.99e-7 J/(g K^4)   ->  5.99e-4 J/(kg K^4)
% ricavati fittando la regione LINEARE del piano Cp/T contro T^2 fra 15 e
% 40 K. Da A discende una temperatura di Debye di 371 K, con densita' di
% 8.89 g/cm^3: quest'ultima coincide con quella che il nostro calcolo della
% densita' dello stack usa per l'Hastelloy (8890 kg/m^3), quindi le due
% fonti sono coerenti fra loro.
%
% UNITA': il fattore 1000 fra J/(g K) e J/(kg K) e' l'unico passaggio di
% conversione, ed e' il punto dove si sbaglia di tre ordini di grandezza.
% Controllo rapido: a 20 K questa formula da' 7.45 J/(kg K), cioe' 7.45
% mJ/(g K), coerente con la figura del Cp del paper.
%
% DUE TRATTI, PERCHE' UNA FORMA SOLA NON COPRE IL RANGE.
%  - Sopra i 10 K vale la forma a due termini di Lu et al., fittata dagli
%    autori fra 15 e 40 K. La legge T^3 di Debye varrebbe sotto Theta/50 =
%    7.4 K, ma la loro simulazione mostra che lo scostamento resta piccolo
%    fino a 40 K.
%  - Sotto i 10 K quella forma NON vale: gli autori stessi misurano un
%    upturn di Cp/T attribuito a transizione tipo spin-glass e a piccoli
%    cluster ferromagnetici, che i due termini non contengono. Li' si usa
%    quindi il DATO MISURATO letto sulla curva del C-276, che e' lineare da
%    1 a 2 mJ/(g K) fra 1 e 10 K, cioe' da 1 a 2 J/(kg K).
%
% RACCORDO. A 10 K il tratto misurato da' 2.00 J/(kg K) e la forma di Lu
% 1.93: c'e' un gradino del 3.5%, che si accetta invece di ritoccare i
% coefficienti pubblicati. Sull'inerzia termica dello stack pesa meno
% dell'incertezza di lettura del grafico.
%
gam   = 0.133;     % [J/kg/K^2] termine elettronico   [Lu et al. 2008]
Acoef = 5.99e-4;   % [J/kg/K^4] termine di Debye      [Lu et al. 2008]
cp1K  = 1.0;       % [J/kg/K] cp misurato a  1 K      [curva C-276]
cp10K = 2.0;       % [J/kg/K] cp misurato a 10 K      [curva C-276]
T_lin = 10;        % [K] estremo superiore del tratto misurato
Tc = max(T, 1);
c  = gam*Tc + Acoef*Tc.^3;
%  tratto misurato: lineare fra 1 e 10 K
lo = Tc <= T_lin;
if any(lo(:))
    c(lo) = cp1K + (cp10K - cp1K)*(Tc(lo) - 1)/(T_lin - 1);
end
%  Sopra i 40 K la forma a due termini esce dalla regione di fit (il Debye
%  satura). Li' si riusa la FORMA del 304 riscalata per raccordarsi con
%  continuita' a 40 K: serve solo alle stampe di controllo, il modello
%  termico non supera mai i 20 K.
hi = Tc > 40;
if any(hi(:))
    sc = (gam*40 + Acoef*40^3)/cpSS(40);
    c(hi) = sc*cpSS(Tc(hi));
end
end


% ------------------------------------------------------------------------
function c = cpSS(T)
% Calore specifico dell'acciaio austenitico 304 (UNS S30400).
% CORRELAZIONE UFFICIALE NIST Cryogenic Material Properties, pagina
% "304 Stainless": log10(cp) = a + b*x + ... + h*x^7, x = log10(T),
% validita' 4-300 K, errore del fit dichiarato 5%.
% massa: e' un SURROGATO dichiarato, non lo stesso materiale.
% bassa temperatura era sbagliato di parecchio: dava 0.24 J/(kg K) a 4.2 K
% contro i 2.07 della correlazione NIST, cioe' un fattore 8.5, e 2.15 contro
% 5.28 a 10 K. L'effetto e' sulla capacita' termica dei solidi, quantificata
% in 7.7 contro quella dell'elio.
co = [22.0061, -127.5528, 303.647, -381.0098, 274.0328, ...
      -112.9212, 24.7593, -2.239153, 0];
x = log10(min(max(T, 4), 300));
e = zeros(size(x));
for k = 1:numel(co), e = e + co(k)*x.^(k-1); end
c = 10.^e;
end

% ------------------------------------------------------------------------
function k = kSS(T)
% Conducibilita' termica dell'acciaio austenitico 304 (UNS S30400).
% CORRELAZIONE UFFICIALE NIST Cryogenic Material Properties, pagina
% "304 Stainless", stessa forma logaritmica del cp: validita' 4-300 K,
% errore del fit dichiarato 2%.
% progressivamente salendo (0.79 contro 0.90 a 10 K, 1.74 contro 2.17 a 20 K).
co = [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, ...
      0.4256, -0.4658, 0.1650, -0.0199];
x = log10(min(max(T, 4), 300));
e = zeros(size(x));
for kk = 1:numel(co), e = e + co(kk)*x.^(kk-1); end
k = 10.^e;
end

% ------------------------------------------------------------------------
function p = toPascal(p)
% Normalizza una pressione a PASCAL deducendo l'unita' dall'ordine di grandezza.
% Serve perche' dataAcquisition memorizza pin_range in bar mentre lo Step 7
% lavora in SI. Le soglie sono larghissime: per un impianto criogenico a elio
% la pressione di esercizio sta fra ~1 e ~50 bar, cioe' 1e5 - 5e6 Pa, quindi
% i tre casi non si sovrappongono e la deduzione non e' ambigua.
%   valore > 1e4   -> gia' Pascal
%   valore > 1e2   -> kPa
%   altrimenti     -> bar   (il caso di dataAcquisition)
% Se un giorno i dati arrivassero in un'unita' diversa, questa funzione va
% aggiornata: e' preferibile un punto unico e dichiarato a una conversione
% sparsa nel codice.
m = max(abs(p(:)));
if m > 1e4
    return;                     % Pa
elseif m > 1e2
    p = p*1e3;                  % kPa -> Pa
else
    p = p*1e5;                  % bar -> Pa
end
end

% ------------------------------------------------------------------------
function g = heTableGrid(which)
% Griglie della tabella dell'elio (serve solo alle stampe di riepilogo).
[Tg, pg] = heTableData();
if strcmpi(which,'T'), g = Tg; else, g = pg; end
end

% ------------------------------------------------------------------------
function h = heEnthalpy(T, p)
% Entalpia specifica dell'elio [J/kg], dalla stessa tabella di heProps.
% Serve solo per il controllo di bilancio entalpico dello Step 7.4.
[Tg, pg, ~, ~, ~, ~, ~, ~, ~, t_h] = heTableData();
T = min(max(T, Tg(1)), Tg(end));
p = min(max(p, pg(1)), pg(end));
h = interp2(pg, Tg, t_h, p, T, 'linear');
end

% ------------------------------------------------------------------------
function [rho, cv, cp, csnd, mu, kf, phi] = heProps(T, p)
% PROPRIETA' DELL'ELIO SUPERCRITICO, interpolate da tabella.
%
% La tabella e' generata da CoolProp (equazione di stato di Ortiz-Vega et al.,
% la stessa di NIST REFPROP per He-4) sulla griglia
%   p = 5, 10, 15, 20 bar   (la banda ammessa dal brief)
%   T = 4.2 ... 50 K        (fitta sotto 8 K, dove il fluido attraversa la
%                            linea pseudo-critica e cp varia di un fattore 5)
% PERCHE' UNA TABELLA E NON UNA CORRELAZIONE: a queste pressioni l'elio e'
% lontanissimo dal gas ideale (rho passa da 150 a 15 kg/m3 fra 4.5 e 20 K a
% 10 bar) e cp ha un massimo pronunciato sulla pseudo-critica. Qualunque forma
% chiusa semplice sbaglierebbe proprio dove sta il punto di progetto.
%
% phi = parametro di Grueneisen = (1/rho)*(dp/di)_rho = (dp/dT)_rho/(rho*cv).
% Fuori tabella i valori sono SATURATI ai bordi, con un avviso una tantum:
% saturare in silenzio nasconderebbe un progetto che e' uscito dal dominio.
persistent Tg pg t_rho t_cv t_cp t_cs t_mu t_kf t_phi warned
if isempty(Tg)
    [Tg, pg, t_rho, t_cv, t_cp, t_cs, t_mu, t_kf, t_phi] = heTableData();
    warned = false;
end
Tq = min(max(T, Tg(1)), Tg(end));
pq = min(max(p, pg(1)), pg(end));
if ~warned && (any(T(:) > Tg(end)) || any(T(:) < Tg(1)) || ...
               any(p(:) > pg(end)) || any(p(:) < pg(1)))
    warning('heProps:outOfTable', ...
        ['Stato dell''elio fuori tabella (T = %.2f ... %.2f K, p = %.2f ... ' ...
         '%.2f bar): valori SATURATI ai bordi. I risultati oltre il bordo non ' ...
         'sono attendibili.'], min(T(:)), max(T(:)), min(p(:))/1e5, max(p(:))/1e5);
    warned = true;
end
rho  = interp2(pg, Tg, t_rho, pq, Tq, 'linear');
if nargout > 1
    cv   = interp2(pg, Tg, t_cv,  pq, Tq, 'linear');
    cp   = interp2(pg, Tg, t_cp,  pq, Tq, 'linear');
    csnd = interp2(pg, Tg, t_cs,  pq, Tq, 'linear');
    mu   = interp2(pg, Tg, t_mu,  pq, Tq, 'linear');
    kf   = interp2(pg, Tg, t_kf,  pq, Tq, 'linear');
    phi  = interp2(pg, Tg, t_phi, pq, Tq, 'linear');
end
end

% ------------------------------------------------------------------------
function [Tg, pg, t_rho, t_cv, t_cp, t_cs, t_mu, t_kf, t_phi, t_h] = heTableData()
% TABELLA DELLE PROPRIETA' DELL'ELIO - dati, nient'altro.
% Generata con CoolProp (EOS di Ortiz-Vega, la stessa di NIST REFPROP per He-4)
% sulla griglia p = [5 6 7 8 9 10 12 14 16 18 20] bar x T = [4.2 ... 50] K.
% GRIGLIA DI PRESSIONE INFITTITA: con sole 4 colonne (5/10/15/20 bar) l'errore
% di interpolazione bilineare a meta' fra due colonne arrivava al 39% su cp e
% al 12% su rho, perche' vicino alla linea pseudo-critica le proprieta' variano
% rapidamente anche con p. Con 11 colonne l'errore massimo scende al 6% su cp
% e all'1.6% su rho (medio: 0.1%). Nel punto di progetto (10 bar, che e' un
% nodo) l'errore era gia' <1%, ma la banda ammessa dal brief e' 5-20 bar e il
% codice deve reggere qualunque p_out l'utente scelga dentro quella banda.
% Unita' SI: rho [kg/m3], cv e cp [J/kg/K], cs [m/s], mu [Pa s], kf [W/m/K],
% h [J/kg], phi [-] (parametro di Grueneisen).
% Le righe sono le temperature (Tg), le colonne le pressioni (pg).
%
% FORMATO: ogni riga della matrice termina con ';' PRIMA del '...'. Senza il
% punto e virgola il '...' e' una semplice continuazione di riga e MATLAB
% costruisce un vettore 1 x NTOT invece della matrice nT x nP: interp2 fallisce
% poi con "Interpolation requires at least two sample points for each grid
% dimension", che non lascia intuire la vera causa.
Tg = [4.2 4.4 4.6 4.8 5 5.2 5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.2 ...
      7.4 7.6 7.8 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 ...
      15 16 17 18 19 20 22 24 26 28 30 35 40 45 50].';
pg = [5 6 7 8 9 10 12 14 16 18 20]*1e5;
t_rho = [ ...
    1.4030e+02 1.4272e+02 1.4491e+02 1.4691e+02 1.4876e+02 1.5048e+02 1.5361e+02 1.5643e+02 1.5899e+02 1.6135e+02 1.6354e+02; ...
    1.3785e+02 1.4051e+02 1.4288e+02 1.4502e+02 1.4699e+02 1.4882e+02 1.5212e+02 1.5506e+02 1.5773e+02 1.6017e+02 1.6243e+02; ...
    1.3513e+02 1.3808e+02 1.4067e+02 1.4299e+02 1.4510e+02 1.4704e+02 1.5053e+02 1.5361e+02 1.5639e+02 1.5892e+02 1.6126e+02; ...
    1.3212e+02 1.3543e+02 1.3828e+02 1.4080e+02 1.4307e+02 1.4514e+02 1.4884e+02 1.5208e+02 1.5497e+02 1.5761e+02 1.6002e+02; ...
    1.2873e+02 1.3250e+02 1.3567e+02 1.3843e+02 1.4089e+02 1.4312e+02 1.4704e+02 1.5046e+02 1.5349e+02 1.5623e+02 1.5873e+02; ...
    1.2490e+02 1.2926e+02 1.3283e+02 1.3587e+02 1.3855e+02 1.4095e+02 1.4514e+02 1.4875e+02 1.5192e+02 1.5478e+02 1.5738e+02; ...
    1.2052e+02 1.2566e+02 1.2971e+02 1.3310e+02 1.3604e+02 1.3864e+02 1.4313e+02 1.4695e+02 1.5028e+02 1.5327e+02 1.5597e+02; ...
    1.1542e+02 1.2162e+02 1.2629e+02 1.3010e+02 1.3334e+02 1.3618e+02 1.4100e+02 1.4505e+02 1.4856e+02 1.5168e+02 1.5450e+02; ...
    1.0940e+02 1.1706e+02 1.2253e+02 1.2684e+02 1.3044e+02 1.3354e+02 1.3875e+02 1.4306e+02 1.4677e+02 1.5003e+02 1.5297e+02; ...
    1.0222e+02 1.1190e+02 1.1836e+02 1.2329e+02 1.2731e+02 1.3073e+02 1.3637e+02 1.4097e+02 1.4489e+02 1.4832e+02 1.5138e+02; ...
    9.3529e+01 1.0605e+02 1.1376e+02 1.1943e+02 1.2395e+02 1.2772e+02 1.3385e+02 1.3878e+02 1.4292e+02 1.4653e+02 1.4973e+02; ...
    8.3198e+01 9.9439e+01 1.0870e+02 1.1524e+02 1.2033e+02 1.2452e+02 1.3120e+02 1.3648e+02 1.4088e+02 1.4467e+02 1.4801e+02; ...
    7.2822e+01 9.2062e+01 1.0318e+02 1.1073e+02 1.1647e+02 1.2111e+02 1.2840e+02 1.3407e+02 1.3874e+02 1.4274e+02 1.4624e+02; ...
    6.4609e+01 8.4125e+01 9.7229e+01 1.0590e+02 1.1235e+02 1.1750e+02 1.2546e+02 1.3155e+02 1.3652e+02 1.4073e+02 1.4441e+02; ...
    5.8443e+01 7.6376e+01 9.0930e+01 1.0079e+02 1.0801e+02 1.1370e+02 1.2238e+02 1.2893e+02 1.3422e+02 1.3866e+02 1.4251e+02; ...
    5.3650e+01 6.9659e+01 8.4503e+01 9.5470e+01 1.0348e+02 1.0972e+02 1.1916e+02 1.2620e+02 1.3183e+02 1.3652e+02 1.4056e+02; ...
    4.9781e+01 6.4147e+01 7.8331e+01 9.0039e+01 9.8806e+01 1.0562e+02 1.1583e+02 1.2337e+02 1.2935e+02 1.3431e+02 1.3855e+02; ...
    4.6567e+01 5.9615e+01 7.2774e+01 8.4657e+01 9.4060e+01 1.0141e+02 1.1239e+02 1.2046e+02 1.2680e+02 1.3203e+02 1.3648e+02; ...
    4.3836e+01 5.5819e+01 6.7962e+01 7.9527e+01 8.9324e+01 9.7160e+01 1.0888e+02 1.1746e+02 1.2418e+02 1.2969e+02 1.3436e+02; ...
    4.1478e+01 5.2579e+01 6.3830e+01 7.4822e+01 8.4705e+01 9.2915e+01 1.0531e+02 1.1440e+02 1.2149e+02 1.2729e+02 1.3218e+02; ...
    3.6748e+01 4.6170e+01 5.5722e+01 6.5189e+01 7.4344e+01 8.2756e+01 9.6394e+01 1.0659e+02 1.1458e+02 1.2109e+02 1.2655e+02; ...
    3.3158e+01 4.1368e+01 4.9717e+01 5.8008e+01 6.6138e+01 7.3969e+01 8.7840e+01 9.8819e+01 1.0755e+02 1.1470e+02 1.2070e+02; ...
    3.0319e+01 3.7604e+01 4.5037e+01 5.2444e+01 5.9727e+01 6.6830e+01 8.0097e+01 9.1381e+01 1.0063e+02 1.0829e+02 1.1478e+02; ...
    2.8005e+01 3.4562e+01 4.1264e+01 4.7970e+01 5.4584e+01 6.1056e+01 7.3411e+01 8.4517e+01 9.4011e+01 1.0203e+02 1.0889e+02; ...
    2.6074e+01 3.2044e+01 3.8148e+01 4.4278e+01 5.0345e+01 5.6297e+01 6.7753e+01 7.8383e+01 8.7846e+01 9.6055e+01 1.0316e+02; ...
    2.4431e+01 2.9920e+01 3.5526e+01 4.1171e+01 4.6778e+01 5.2295e+01 6.2958e+01 7.3008e+01 8.2228e+01 9.0451e+01 9.7698e+01; ...
    2.3012e+01 2.8099e+01 3.3287e+01 3.8518e+01 4.3729e+01 4.8873e+01 5.8850e+01 6.8328e+01 7.7187e+01 8.5281e+01 9.2552e+01; ...
    2.1771e+01 2.6517e+01 3.1349e+01 3.6223e+01 4.1091e+01 4.5909e+01 5.5287e+01 6.4241e+01 7.2701e+01 8.0568e+01 8.7765e+01; ...
    2.0674e+01 2.5128e+01 2.9653e+01 3.4218e+01 3.8783e+01 4.3313e+01 5.2163e+01 6.0646e+01 6.8712e+01 7.6304e+01 8.3354e+01; ...
    1.9695e+01 2.3895e+01 2.8154e+01 3.2448e+01 3.6747e+01 4.1021e+01 4.9399e+01 5.7459e+01 6.5157e+01 7.2457e+01 7.9314e+01; ...
    1.8814e+01 2.2792e+01 2.6818e+01 3.0874e+01 3.4937e+01 3.8982e+01 4.6934e+01 5.4613e+01 6.1970e+01 6.8984e+01 7.5625e+01; ...
    1.8017e+01 2.1797e+01 2.5618e+01 2.9463e+01 3.3316e+01 3.7155e+01 4.4721e+01 5.2053e+01 5.9099e+01 6.5840e+01 7.2259e+01; ...
    1.6627e+01 2.0073e+01 2.3546e+01 2.7035e+01 3.0529e+01 3.4015e+01 4.0911e+01 4.7634e+01 5.4130e+01 6.0377e+01 6.6369e+01; ...
    1.5453e+01 1.8625e+01 2.1814e+01 2.5014e+01 2.8216e+01 3.1411e+01 3.7746e+01 4.3950e+01 4.9974e+01 5.5794e+01 6.1400e+01; ...
    1.4445e+01 1.7388e+01 2.0341e+01 2.3301e+01 2.6260e+01 2.9212e+01 3.5073e+01 4.0831e+01 4.6445e+01 5.1890e+01 5.7156e+01; ...
    1.3568e+01 1.6317e+01 1.9070e+01 2.1826e+01 2.4580e+01 2.7327e+01 3.2783e+01 3.8155e+01 4.3410e+01 4.8524e+01 5.3486e+01; ...
    1.2798e+01 1.5378e+01 1.7960e+01 2.0542e+01 2.3119e+01 2.5690e+01 3.0797e+01 3.5833e+01 4.0771e+01 4.5592e+01 5.0282e+01; ...
    1.2115e+01 1.4548e+01 1.6981e+01 1.9411e+01 2.1836e+01 2.4254e+01 2.9057e+01 3.3797e+01 3.8455e+01 4.3013e+01 4.7458e+01; ...
    1.0955e+01 1.3143e+01 1.5327e+01 1.7506e+01 1.9679e+01 2.1844e+01 2.6144e+01 3.0393e+01 3.4578e+01 3.8687e+01 4.2711e+01; ...
    1.0005e+01 1.1996e+01 1.3981e+01 1.5960e+01 1.7932e+01 1.9896e+01 2.3796e+01 2.7652e+01 3.1455e+01 3.5197e+01 3.8872e+01; ...
    9.2118e+00 1.1040e+01 1.2862e+01 1.4677e+01 1.6485e+01 1.8285e+01 2.1858e+01 2.5392e+01 2.8880e+01 3.2318e+01 3.5701e+01; ...
    8.5383e+00 1.0230e+01 1.1914e+01 1.3592e+01 1.5263e+01 1.6926e+01 2.0227e+01 2.3492e+01 2.6718e+01 2.9899e+01 3.3034e+01; ...
    7.9588e+00 9.5333e+00 1.1101e+01 1.2663e+01 1.4217e+01 1.5764e+01 1.8834e+01 2.1871e+01 2.4872e+01 2.7835e+01 3.0758e+01; ...
    6.8097e+00 8.1545e+00 9.4933e+00 1.0826e+01 1.2152e+01 1.3472e+01 1.6092e+01 1.8685e+01 2.1250e+01 2.3785e+01 2.6289e+01; ...
    5.9545e+00 7.1299e+00 8.2999e+00 9.4644e+00 1.0624e+01 1.1777e+01 1.4068e+01 1.6335e+01 1.8579e+01 2.0800e+01 2.2996e+01; ...
    5.2923e+00 6.3370e+00 7.3770e+00 8.4123e+00 9.4428e+00 1.0469e+01 1.2506e+01 1.4524e+01 1.6522e+01 1.8501e+01 2.0460e+01; ...
    4.7638e+00 5.7045e+00 6.6410e+00 7.5735e+00 8.5019e+00 9.4262e+00 1.1262e+01 1.3082e+01 1.4885e+01 1.6671e+01 1.8441e+01];
t_cv = [ ...
    2.2799e+03 2.2592e+03 2.2395e+03 2.2205e+03 2.2022e+03 2.1844e+03 2.1502e+03 2.1178e+03 2.0871e+03 2.0580e+03 2.0302e+03; ...
    2.3373e+03 2.3167e+03 2.2972e+03 2.2786e+03 2.2606e+03 2.2432e+03 2.2098e+03 2.1782e+03 2.1482e+03 2.1197e+03 2.0925e+03; ...
    2.3921e+03 2.3712e+03 2.3518e+03 2.3334e+03 2.3158e+03 2.2987e+03 2.2660e+03 2.2352e+03 2.2059e+03 2.1780e+03 2.1515e+03; ...
    2.4447e+03 2.4231e+03 2.4036e+03 2.3853e+03 2.3679e+03 2.3511e+03 2.3192e+03 2.2891e+03 2.2605e+03 2.2333e+03 2.2075e+03; ...
    2.4959e+03 2.4727e+03 2.4528e+03 2.4345e+03 2.4173e+03 2.4008e+03 2.3695e+03 2.3401e+03 2.3122e+03 2.2858e+03 2.2606e+03; ...
    2.5464e+03 2.5206e+03 2.4997e+03 2.4812e+03 2.4641e+03 2.4478e+03 2.4172e+03 2.3884e+03 2.3613e+03 2.3355e+03 2.3110e+03; ...
    2.5978e+03 2.5672e+03 2.5447e+03 2.5257e+03 2.5085e+03 2.4924e+03 2.4623e+03 2.4342e+03 2.4078e+03 2.3827e+03 2.3589e+03; ...
    2.6525e+03 2.6133e+03 2.5881e+03 2.5681e+03 2.5506e+03 2.5346e+03 2.5049e+03 2.4775e+03 2.4518e+03 2.4274e+03 2.4043e+03; ...
    2.7142e+03 2.6602e+03 2.6303e+03 2.6087e+03 2.5907e+03 2.5745e+03 2.5453e+03 2.5185e+03 2.4935e+03 2.4698e+03 2.4474e+03; ...
    2.7878e+03 2.7093e+03 2.6721e+03 2.6478e+03 2.6289e+03 2.6124e+03 2.5834e+03 2.5572e+03 2.5329e+03 2.5099e+03 2.4882e+03; ...
    2.8762e+03 2.7623e+03 2.7141e+03 2.6859e+03 2.6654e+03 2.6484e+03 2.6195e+03 2.5938e+03 2.5701e+03 2.5478e+03 2.5268e+03; ...
    2.9676e+03 2.8199e+03 2.7569e+03 2.7232e+03 2.7005e+03 2.6827e+03 2.6535e+03 2.6283e+03 2.6052e+03 2.5836e+03 2.5632e+03; ...
    3.0242e+03 2.8798e+03 2.8009e+03 2.7600e+03 2.7343e+03 2.7153e+03 2.6856e+03 2.6608e+03 2.6383e+03 2.6174e+03 2.5977e+03; ...
    3.0397e+03 2.9346e+03 2.8453e+03 2.7966e+03 2.7672e+03 2.7465e+03 2.7160e+03 2.6914e+03 2.6695e+03 2.6492e+03 2.6302e+03; ...
    3.0406e+03 2.9735e+03 2.8878e+03 2.8326e+03 2.7992e+03 2.7765e+03 2.7447e+03 2.7202e+03 2.6988e+03 2.6792e+03 2.6608e+03; ...
    3.0405e+03 2.9947e+03 2.9250e+03 2.8671e+03 2.8300e+03 2.8052e+03 2.7719e+03 2.7474e+03 2.7264e+03 2.7074e+03 2.6897e+03; ...
    3.0426e+03 3.0057e+03 2.9537e+03 2.8989e+03 2.8595e+03 2.8326e+03 2.7976e+03 2.7729e+03 2.7524e+03 2.7339e+03 2.7168e+03; ...
    3.0470e+03 3.0132e+03 2.9738e+03 2.9266e+03 2.8870e+03 2.8586e+03 2.8219e+03 2.7969e+03 2.7767e+03 2.7588e+03 2.7423e+03; ...
    3.0529e+03 3.0201e+03 2.9880e+03 2.9492e+03 2.9119e+03 2.8829e+03 2.8448e+03 2.8196e+03 2.7996e+03 2.7822e+03 2.7663e+03; ...
    3.0594e+03 3.0274e+03 2.9989e+03 2.9672e+03 2.9337e+03 2.9053e+03 2.8664e+03 2.8408e+03 2.8211e+03 2.8042e+03 2.7888e+03; ...
    3.0753e+03 3.0470e+03 3.0217e+03 2.9988e+03 2.9753e+03 2.9519e+03 2.9142e+03 2.8884e+03 2.8692e+03 2.8533e+03 2.8392e+03; ...
    3.0879e+03 3.0652e+03 3.0425e+03 3.0225e+03 3.0041e+03 2.9860e+03 2.9531e+03 2.9284e+03 2.9100e+03 2.8950e+03 2.8822e+03; ...
    3.0969e+03 3.0798e+03 3.0608e+03 3.0429e+03 3.0267e+03 3.0118e+03 2.9842e+03 2.9616e+03 2.9443e+03 2.9305e+03 2.9187e+03; ...
    3.1031e+03 3.0905e+03 3.0755e+03 3.0602e+03 3.0458e+03 3.0327e+03 3.0090e+03 2.9889e+03 2.9731e+03 2.9604e+03 2.9497e+03; ...
    3.1072e+03 3.0982e+03 3.0868e+03 3.0742e+03 3.0618e+03 3.0502e+03 3.0294e+03 3.0115e+03 2.9971e+03 2.9856e+03 2.9760e+03; ...
    3.1101e+03 3.1036e+03 3.0951e+03 3.0852e+03 3.0749e+03 3.0648e+03 3.0464e+03 3.0305e+03 3.0173e+03 3.0068e+03 2.9983e+03; ...
    3.1121e+03 3.1075e+03 3.1012e+03 3.0936e+03 3.0853e+03 3.0768e+03 3.0607e+03 3.0465e+03 3.0345e+03 3.0249e+03 3.0172e+03; ...
    3.1135e+03 3.1102e+03 3.1056e+03 3.0999e+03 3.0933e+03 3.0864e+03 3.0726e+03 3.0601e+03 3.0493e+03 3.0404e+03 3.0334e+03; ...
    3.1146e+03 3.1122e+03 3.1089e+03 3.1046e+03 3.0996e+03 3.0941e+03 3.0825e+03 3.0716e+03 3.0619e+03 3.0539e+03 3.0474e+03; ...
    3.1154e+03 3.1138e+03 3.1114e+03 3.1082e+03 3.1044e+03 3.1001e+03 3.0907e+03 3.0813e+03 3.0728e+03 3.0655e+03 3.0596e+03; ...
    3.1161e+03 3.1150e+03 3.1133e+03 3.1110e+03 3.1081e+03 3.1048e+03 3.0973e+03 3.0895e+03 3.0821e+03 3.0756e+03 3.0702e+03; ...
    3.1166e+03 3.1159e+03 3.1147e+03 3.1131e+03 3.1111e+03 3.1086e+03 3.1027e+03 3.0963e+03 3.0900e+03 3.0843e+03 3.0795e+03; ...
    3.1175e+03 3.1173e+03 3.1169e+03 3.1162e+03 3.1152e+03 3.1140e+03 3.1107e+03 3.1067e+03 3.1025e+03 3.0984e+03 3.0948e+03; ...
    3.1182e+03 3.1184e+03 3.1184e+03 3.1183e+03 3.1180e+03 3.1175e+03 3.1160e+03 3.1139e+03 3.1114e+03 3.1088e+03 3.1063e+03; ...
    3.1188e+03 3.1192e+03 3.1195e+03 3.1198e+03 3.1199e+03 3.1200e+03 3.1197e+03 3.1189e+03 3.1178e+03 3.1164e+03 3.1150e+03; ...
    3.1193e+03 3.1199e+03 3.1204e+03 3.1209e+03 3.1213e+03 3.1217e+03 3.1222e+03 3.1225e+03 3.1224e+03 3.1220e+03 3.1216e+03; ...
    3.1198e+03 3.1204e+03 3.1211e+03 3.1218e+03 3.1224e+03 3.1230e+03 3.1241e+03 3.1250e+03 3.1257e+03 3.1261e+03 3.1264e+03; ...
    3.1201e+03 3.1209e+03 3.1217e+03 3.1225e+03 3.1233e+03 3.1241e+03 3.1255e+03 3.1269e+03 3.1281e+03 3.1292e+03 3.1301e+03; ...
    3.1207e+03 3.1217e+03 3.1226e+03 3.1236e+03 3.1246e+03 3.1255e+03 3.1275e+03 3.1294e+03 3.1313e+03 3.1331e+03 3.1349e+03; ...
    3.1212e+03 3.1222e+03 3.1233e+03 3.1243e+03 3.1254e+03 3.1265e+03 3.1287e+03 3.1310e+03 3.1332e+03 3.1354e+03 3.1376e+03; ...
    3.1215e+03 3.1226e+03 3.1237e+03 3.1248e+03 3.1260e+03 3.1272e+03 3.1295e+03 3.1319e+03 3.1343e+03 3.1368e+03 3.1392e+03; ...
    3.1217e+03 3.1228e+03 3.1240e+03 3.1252e+03 3.1264e+03 3.1276e+03 3.1300e+03 3.1325e+03 3.1350e+03 3.1376e+03 3.1401e+03; ...
    3.1218e+03 3.1230e+03 3.1242e+03 3.1254e+03 3.1266e+03 3.1278e+03 3.1303e+03 3.1329e+03 3.1354e+03 3.1380e+03 3.1406e+03; ...
    3.1219e+03 3.1231e+03 3.1244e+03 3.1256e+03 3.1268e+03 3.1280e+03 3.1305e+03 3.1331e+03 3.1356e+03 3.1382e+03 3.1408e+03; ...
    3.1219e+03 3.1230e+03 3.1242e+03 3.1254e+03 3.1267e+03 3.1279e+03 3.1303e+03 3.1327e+03 3.1352e+03 3.1377e+03 3.1402e+03; ...
    3.1217e+03 3.1228e+03 3.1240e+03 3.1252e+03 3.1263e+03 3.1275e+03 3.1298e+03 3.1322e+03 3.1346e+03 3.1370e+03 3.1393e+03; ...
    3.1215e+03 3.1226e+03 3.1237e+03 3.1248e+03 3.1259e+03 3.1270e+03 3.1293e+03 3.1315e+03 3.1338e+03 3.1361e+03 3.1384e+03];
t_cp = [ ...
    3.4125e+03 3.2683e+03 3.1520e+03 3.0550e+03 2.9721e+03 2.9001e+03 2.7795e+03 2.6813e+03 2.5986e+03 2.5274e+03 2.4649e+03; ...
    3.6914e+03 3.5114e+03 3.3701e+03 3.2548e+03 3.1581e+03 3.0750e+03 2.9381e+03 2.8285e+03 2.7373e+03 2.6595e+03 2.5918e+03; ...
    4.0105e+03 3.7814e+03 3.6076e+03 3.4692e+03 3.3551e+03 3.2587e+03 3.1025e+03 2.9795e+03 2.8787e+03 2.7935e+03 2.7200e+03; ...
    4.3836e+03 4.0861e+03 3.8691e+03 3.7011e+03 3.5655e+03 3.4527e+03 3.2735e+03 3.1351e+03 3.0231e+03 2.9297e+03 2.8497e+03; ...
    4.8306e+03 4.4355e+03 4.1605e+03 3.9542e+03 3.7917e+03 3.6590e+03 3.4522e+03 3.2958e+03 3.1712e+03 3.0683e+03 2.9811e+03; ...
    5.3808e+03 4.8428e+03 4.4888e+03 4.2330e+03 4.0366e+03 3.8794e+03 3.6397e+03 3.4623e+03 3.3232e+03 3.2097e+03 3.1144e+03; ...
    6.0769e+03 5.3255e+03 4.8631e+03 4.5425e+03 4.3034e+03 4.1161e+03 3.8371e+03 3.6352e+03 3.4796e+03 3.3542e+03 3.2499e+03; ...
    6.9810e+03 5.9054e+03 5.2938e+03 4.8886e+03 4.5958e+03 4.3715e+03 4.0455e+03 3.8153e+03 3.6408e+03 3.5020e+03 3.3878e+03; ...
    8.1828e+03 6.6078e+03 5.7918e+03 5.2774e+03 4.9173e+03 4.6480e+03 4.2661e+03 4.0030e+03 3.8071e+03 3.6534e+03 3.5281e+03; ...
    9.8272e+03 7.4577e+03 6.3670e+03 5.7139e+03 5.2711e+03 4.9476e+03 4.5000e+03 4.1991e+03 3.9790e+03 3.8085e+03 3.6710e+03; ...
    1.2154e+04 8.4789e+03 7.0240e+03 6.2007e+03 5.6593e+03 5.2719e+03 4.7479e+03 4.4040e+03 4.1566e+03 3.9676e+03 3.8166e+03; ...
    1.4724e+04 9.6973e+03 7.7606e+03 6.7357e+03 6.0809e+03 5.6208e+03 5.0102e+03 4.6178e+03 4.3402e+03 4.1306e+03 3.9651e+03; ...
    1.4686e+04 1.1084e+04 8.5679e+03 7.3100e+03 6.5311e+03 5.9919e+03 5.2864e+03 4.8407e+03 4.5297e+03 4.2978e+03 4.1163e+03; ...
    1.2902e+04 1.2230e+04 9.4253e+03 7.9101e+03 7.0007e+03 6.3797e+03 5.5748e+03 5.0719e+03 4.7250e+03 4.4689e+03 4.2703e+03; ...
    1.1402e+04 1.2334e+04 1.0249e+04 8.5180e+03 7.4767e+03 6.7758e+03 5.8722e+03 5.3103e+03 4.9255e+03 4.6438e+03 4.4269e+03; ...
    1.0339e+04 1.1559e+04 1.0810e+04 9.1023e+03 7.9450e+03 7.1693e+03 6.1737e+03 5.5539e+03 5.1304e+03 4.8219e+03 4.5858e+03; ...
    9.5769e+03 1.0655e+04 1.0862e+04 9.5934e+03 8.3882e+03 7.5491e+03 6.4731e+03 5.7996e+03 5.3381e+03 5.0025e+03 4.7467e+03; ...
    9.0109e+03 9.9020e+03 1.0478e+04 9.8810e+03 8.7790e+03 7.9040e+03 6.7632e+03 6.0438e+03 5.5469e+03 5.1846e+03 4.9089e+03; ...
    8.5752e+03 9.3150e+03 9.9376e+03 9.8909e+03 9.0734e+03 8.2203e+03 7.0374e+03 6.2819e+03 5.7542e+03 5.3669e+03 5.0717e+03; ...
    8.2281e+03 8.8568e+03 9.4200e+03 9.6679e+03 9.2233e+03 8.4789e+03 7.2897e+03 6.5094e+03 5.9571e+03 5.5478e+03 5.2342e+03; ...
    7.5991e+03 8.0674e+03 8.4593e+03 8.8065e+03 8.9562e+03 8.7321e+03 7.7909e+03 7.0064e+03 6.4260e+03 5.9804e+03 5.6305e+03; ...
    7.1684e+03 7.5574e+03 7.8581e+03 8.1195e+03 8.3451e+03 8.4419e+03 8.0368e+03 7.3686e+03 6.8083e+03 6.3596e+03 5.9946e+03; ...
    6.8520e+03 7.1897e+03 7.4483e+03 7.6526e+03 7.8359e+03 7.9886e+03 7.9884e+03 7.5692e+03 7.0831e+03 6.6614e+03 6.3048e+03; ...
    6.6099e+03 6.9070e+03 7.1421e+03 7.3195e+03 7.4662e+03 7.5992e+03 7.7537e+03 7.5981e+03 7.2446e+03 6.8771e+03 6.5478e+03; ...
    6.4197e+03 6.6819e+03 6.8990e+03 7.0639e+03 7.1915e+03 7.3013e+03 7.4780e+03 7.4908e+03 7.2946e+03 7.0084e+03 6.7212e+03; ...
    6.2673e+03 6.4993e+03 6.6993e+03 6.8565e+03 6.9758e+03 7.0718e+03 7.2320e+03 7.3168e+03 7.2506e+03 7.0615e+03 6.8295e+03; ...
    6.1430e+03 6.3489e+03 6.5324e+03 6.6824e+03 6.7980e+03 6.8876e+03 7.0293e+03 7.1321e+03 7.1458e+03 7.0477e+03 6.8797e+03; ...
    6.0403e+03 6.2238e+03 6.3915e+03 6.5336e+03 6.6465e+03 6.7338e+03 6.8634e+03 6.9635e+03 7.0153e+03 6.9854e+03 6.8805e+03; ...
    5.9544e+03 6.1186e+03 6.2716e+03 6.4053e+03 6.5150e+03 6.6013e+03 6.7251e+03 6.8178e+03 6.8828e+03 6.8953e+03 6.8434e+03; ...
    5.8816e+03 6.0293e+03 6.1690e+03 6.2939e+03 6.3996e+03 6.4849e+03 6.6068e+03 6.6933e+03 6.7601e+03 6.7946e+03 6.7812e+03; ...
    5.8195e+03 5.9530e+03 6.0805e+03 6.1969e+03 6.2979e+03 6.3817e+03 6.5031e+03 6.5861e+03 6.6506e+03 6.6947e+03 6.7056e+03; ...
    5.7660e+03 5.8871e+03 6.0039e+03 6.1121e+03 6.2080e+03 6.2895e+03 6.4107e+03 6.4923e+03 6.5540e+03 6.6010e+03 6.6258e+03; ...
    5.6787e+03 5.7800e+03 5.8786e+03 5.9719e+03 6.0573e+03 6.1330e+03 6.2519e+03 6.3338e+03 6.3926e+03 6.4386e+03 6.4733e+03; ...
    5.6110e+03 5.6970e+03 5.7813e+03 5.8621e+03 5.9376e+03 6.0064e+03 6.1202e+03 6.2027e+03 6.2617e+03 6.3063e+03 6.3419e+03; ...
    5.5574e+03 5.6314e+03 5.7042e+03 5.7746e+03 5.8414e+03 5.9034e+03 6.0100e+03 6.0916e+03 6.1516e+03 6.1963e+03 6.2314e+03; ...
    5.5139e+03 5.5785e+03 5.6421e+03 5.7039e+03 5.7631e+03 5.8188e+03 5.9175e+03 5.9965e+03 6.0569e+03 6.1024e+03 6.1377e+03; ...
    5.4782e+03 5.5351e+03 5.5911e+03 5.6458e+03 5.6986e+03 5.7487e+03 5.8394e+03 5.9147e+03 5.9745e+03 6.0207e+03 6.0565e+03; ...
    5.4484e+03 5.4990e+03 5.5488e+03 5.5976e+03 5.6448e+03 5.6901e+03 5.7731e+03 5.8442e+03 5.9025e+03 5.9487e+03 5.9852e+03; ...
    5.4018e+03 5.4427e+03 5.4830e+03 5.5225e+03 5.5610e+03 5.5982e+03 5.6680e+03 5.7302e+03 5.7837e+03 5.8285e+03 5.8653e+03; ...
    5.3673e+03 5.4012e+03 5.4345e+03 5.4673e+03 5.4993e+03 5.5304e+03 5.5895e+03 5.6434e+03 5.6915e+03 5.7332e+03 5.7689e+03; ...
    5.3410e+03 5.3695e+03 5.3977e+03 5.4253e+03 5.4524e+03 5.4788e+03 5.5293e+03 5.5762e+03 5.6189e+03 5.6571e+03 5.6906e+03; ...
    5.3203e+03 5.3447e+03 5.3688e+03 5.3925e+03 5.4158e+03 5.4385e+03 5.4822e+03 5.5232e+03 5.5611e+03 5.5956e+03 5.6267e+03; ...
    5.3037e+03 5.3250e+03 5.3459e+03 5.3664e+03 5.3866e+03 5.4064e+03 5.4445e+03 5.4806e+03 5.5143e+03 5.5455e+03 5.5740e+03; ...
    5.2742e+03 5.2898e+03 5.3051e+03 5.3201e+03 5.3349e+03 5.3494e+03 5.3776e+03 5.4046e+03 5.4303e+03 5.4545e+03 5.4772e+03; ...
    5.2551e+03 5.2670e+03 5.2787e+03 5.2902e+03 5.3016e+03 5.3127e+03 5.3345e+03 5.3554e+03 5.3755e+03 5.3947e+03 5.4129e+03; ...
    5.2420e+03 5.2514e+03 5.2606e+03 5.2698e+03 5.2787e+03 5.2876e+03 5.3048e+03 5.3215e+03 5.3376e+03 5.3531e+03 5.3680e+03; ...
    5.2326e+03 5.2402e+03 5.2477e+03 5.2550e+03 5.2623e+03 5.2695e+03 5.2835e+03 5.2972e+03 5.3104e+03 5.3231e+03 5.3354e+03];
t_cs = [ ...
    2.4155e+02 2.5110e+02 2.5979e+02 2.6781e+02 2.7530e+02 2.8237e+02 2.9545e+02 3.0743e+02 3.1854e+02 3.2895e+02 3.3876e+02; ...
    2.3594e+02 2.4619e+02 2.5541e+02 2.6386e+02 2.7171e+02 2.7907e+02 2.9263e+02 3.0499e+02 3.1641e+02 3.2706e+02 3.3708e+02; ...
    2.2978e+02 2.4085e+02 2.5068e+02 2.5960e+02 2.6784e+02 2.7552e+02 2.8960e+02 3.0236e+02 3.1410e+02 3.2501e+02 3.3525e+02; ...
    2.2302e+02 2.3506e+02 2.4558e+02 2.5503e+02 2.6369e+02 2.7173e+02 2.8636e+02 2.9954e+02 3.1161e+02 3.2280e+02 3.3327e+02; ...
    2.1559e+02 2.2878e+02 2.4010e+02 2.5015e+02 2.5928e+02 2.6769e+02 2.8292e+02 2.9654e+02 3.0896e+02 3.2044e+02 3.3114e+02; ...
    2.0742e+02 2.2201e+02 2.3423e+02 2.4494e+02 2.5458e+02 2.6341e+02 2.7927e+02 2.9336e+02 3.0615e+02 3.1792e+02 3.2887e+02; ...
    1.9846e+02 2.1471e+02 2.2798e+02 2.3943e+02 2.4963e+02 2.5890e+02 2.7543e+02 2.9001e+02 3.0318e+02 3.1526e+02 3.2646e+02; ...
    1.8867e+02 2.0692e+02 2.2137e+02 2.3362e+02 2.4442e+02 2.5416e+02 2.7140e+02 2.8650e+02 3.0007e+02 3.1247e+02 3.2393e+02; ...
    1.7810e+02 1.9867e+02 2.1443e+02 2.2755e+02 2.3899e+02 2.4923e+02 2.6721e+02 2.8285e+02 2.9683e+02 3.0955e+02 3.2129e+02; ...
    1.6692e+02 1.9010e+02 2.0725e+02 2.2128e+02 2.3337e+02 2.4413e+02 2.6287e+02 2.7906e+02 2.9346e+02 3.0652e+02 3.1853e+02; ...
    1.5568e+02 1.8139e+02 1.9996e+02 2.1489e+02 2.2764e+02 2.3890e+02 2.5841e+02 2.7516e+02 2.8999e+02 3.0339e+02 3.1568e+02; ...
    1.4645e+02 1.7278e+02 1.9272e+02 2.0851e+02 2.2187e+02 2.3361e+02 2.5386e+02 2.7116e+02 2.8642e+02 3.0017e+02 3.1274e+02; ...
    1.4290e+02 1.6476e+02 1.8569e+02 2.0225e+02 2.1616e+02 2.2834e+02 2.4927e+02 2.6710e+02 2.8278e+02 2.9687e+02 3.0973e+02; ...
    1.4350e+02 1.5848e+02 1.7907e+02 1.9626e+02 2.1062e+02 2.2317e+02 2.4469e+02 2.6301e+02 2.7909e+02 2.9352e+02 3.0665e+02; ...
    1.4527e+02 1.5533e+02 1.7321e+02 1.9063e+02 2.0537e+02 2.1819e+02 2.4019e+02 2.5892e+02 2.7538e+02 2.9013e+02 3.0353e+02; ...
    1.4739e+02 1.5499e+02 1.6873e+02 1.8552e+02 2.0046e+02 2.1348e+02 2.3582e+02 2.5489e+02 2.7167e+02 2.8671e+02 3.0038e+02; ...
    1.4968e+02 1.5604e+02 1.6623e+02 1.8116e+02 1.9597e+02 2.0911e+02 2.3165e+02 2.5096e+02 2.6801e+02 2.8331e+02 2.9722e+02; ...
    1.5210e+02 1.5763e+02 1.6558e+02 1.7789e+02 1.9201e+02 2.0512e+02 2.2774e+02 2.4718e+02 2.6442e+02 2.7994e+02 2.9406e+02; ...
    1.5458e+02 1.5946e+02 1.6613e+02 1.7596e+02 1.8872e+02 2.0155e+02 2.2413e+02 2.4359e+02 2.6094e+02 2.7663e+02 2.9094e+02; ...
    1.5711e+02 1.6141e+02 1.6732e+02 1.7529e+02 1.8629e+02 1.9846e+02 2.2083e+02 2.4024e+02 2.5762e+02 2.7341e+02 2.8786e+02; ...
    1.6352e+02 1.6662e+02 1.7135e+02 1.7706e+02 1.8428e+02 1.9351e+02 2.1407e+02 2.3302e+02 2.5016e+02 2.6593e+02 2.8053e+02; ...
    1.6985e+02 1.7211e+02 1.7596e+02 1.8072e+02 1.8619e+02 1.9279e+02 2.0968e+02 2.2755e+02 2.4415e+02 2.5956e+02 2.7401e+02; ...
    1.7600e+02 1.7770e+02 1.8082e+02 1.8491e+02 1.8958e+02 1.9481e+02 2.0800e+02 2.2387e+02 2.3964e+02 2.5449e+02 2.6854e+02; ...
    1.8191e+02 1.8327e+02 1.8582e+02 1.8935e+02 1.9348e+02 1.9802e+02 2.0869e+02 2.2205e+02 2.3656e+02 2.5072e+02 2.6423e+02; ...
    1.8760e+02 1.8874e+02 1.9088e+02 1.9392e+02 1.9760e+02 2.0168e+02 2.1087e+02 2.2203e+02 2.3489e+02 2.4814e+02 2.6104e+02; ...
    1.9305e+02 1.9407e+02 1.9591e+02 1.9856e+02 2.0184e+02 2.0556e+02 2.1382e+02 2.2340e+02 2.3458e+02 2.4669e+02 2.5887e+02; ...
    1.9829e+02 1.9925e+02 2.0089e+02 2.0322e+02 2.0616e+02 2.0954e+02 2.1714e+02 2.2566e+02 2.3543e+02 2.4631e+02 2.5764e+02; ...
    2.0334e+02 2.0427e+02 2.0577e+02 2.0786e+02 2.1050e+02 2.1360e+02 2.2064e+02 2.2845e+02 2.3715e+02 2.4687e+02 2.5727e+02; ...
    2.0822e+02 2.0914e+02 2.1054e+02 2.1245e+02 2.1485e+02 2.1769e+02 2.2425e+02 2.3152e+02 2.3945e+02 2.4820e+02 2.5769e+02; ...
    2.1293e+02 2.1386e+02 2.1520e+02 2.1697e+02 2.1918e+02 2.2179e+02 2.2791e+02 2.3475e+02 2.4212e+02 2.5011e+02 2.5877e+02; ...
    2.1750e+02 2.1845e+02 2.1974e+02 2.2141e+02 2.2346e+02 2.2589e+02 2.3161e+02 2.3807e+02 2.4501e+02 2.5241e+02 2.6037e+02; ...
    2.2194e+02 2.2291e+02 2.2417e+02 2.2576e+02 2.2769e+02 2.2996e+02 2.3532e+02 2.4145e+02 2.4803e+02 2.5498e+02 2.6238e+02; ...
    2.3046e+02 2.3147e+02 2.3271e+02 2.3420e+02 2.3595e+02 2.3798e+02 2.4275e+02 2.4828e+02 2.5427e+02 2.6058e+02 2.6716e+02; ...
    2.3857e+02 2.3963e+02 2.4086e+02 2.4229e+02 2.4393e+02 2.4579e+02 2.5012e+02 2.5514e+02 2.6065e+02 2.6646e+02 2.7249e+02; ...
    2.4633e+02 2.4743e+02 2.4867e+02 2.5006e+02 2.5163e+02 2.5337e+02 2.5736e+02 2.6198e+02 2.6706e+02 2.7247e+02 2.7808e+02; ...
    2.5378e+02 2.5491e+02 2.5616e+02 2.5754e+02 2.5905e+02 2.6071e+02 2.6446e+02 2.6874e+02 2.7347e+02 2.7852e+02 2.8378e+02; ...
    2.6097e+02 2.6213e+02 2.6338e+02 2.6475e+02 2.6623e+02 2.6782e+02 2.7138e+02 2.7540e+02 2.7982e+02 2.8456e+02 2.8952e+02; ...
    2.6791e+02 2.6910e+02 2.7036e+02 2.7172e+02 2.7317e+02 2.7472e+02 2.7813e+02 2.8194e+02 2.8611e+02 2.9057e+02 2.9526e+02; ...
    2.8118e+02 2.8240e+02 2.8368e+02 2.8502e+02 2.8644e+02 2.8793e+02 2.9113e+02 2.9462e+02 2.9840e+02 3.0243e+02 3.0667e+02; ...
    2.9374e+02 2.9498e+02 2.9627e+02 2.9760e+02 2.9899e+02 3.0044e+02 3.0349e+02 3.0677e+02 3.1027e+02 3.1398e+02 3.1787e+02; ...
    3.0571e+02 3.0696e+02 3.0824e+02 3.0957e+02 3.1094e+02 3.1235e+02 3.1530e+02 3.1843e+02 3.2172e+02 3.2518e+02 3.2880e+02; ...
    3.1716e+02 3.1841e+02 3.1970e+02 3.2101e+02 3.2236e+02 3.2375e+02 3.2661e+02 3.2962e+02 3.3276e+02 3.3603e+02 3.3943e+02; ...
    3.2816e+02 3.2942e+02 3.3070e+02 3.3200e+02 3.3333e+02 3.3469e+02 3.3749e+02 3.4039e+02 3.4341e+02 3.4653e+02 3.4976e+02; ...
    3.5404e+02 3.5528e+02 3.5653e+02 3.5780e+02 3.5909e+02 3.6039e+02 3.6305e+02 3.6577e+02 3.6855e+02 3.7141e+02 3.7433e+02; ...
    3.7802e+02 3.7923e+02 3.8045e+02 3.8168e+02 3.8292e+02 3.8418e+02 3.8671e+02 3.8930e+02 3.9192e+02 3.9458e+02 3.9729e+02; ...
    4.0048e+02 4.0166e+02 4.0285e+02 4.0404e+02 4.0524e+02 4.0644e+02 4.0888e+02 4.1134e+02 4.1383e+02 4.1635e+02 4.1890e+02; ...
    4.2170e+02 4.2284e+02 4.2399e+02 4.2514e+02 4.2630e+02 4.2746e+02 4.2981e+02 4.3217e+02 4.3455e+02 4.3695e+02 4.3937e+02];
t_mu = [ ...
    3.9984e-06 4.1682e-06 4.3324e-06 4.4925e-06 4.6495e-06 4.8043e-06 5.1091e-06 5.4102e-06 5.7097e-06 6.0092e-06 6.3097e-06; ...
    3.8906e-06 4.0612e-06 4.2248e-06 4.3832e-06 4.5379e-06 4.6897e-06 4.9874e-06 5.2800e-06 5.5700e-06 5.8590e-06 6.1483e-06; ...
    3.7821e-06 3.9548e-06 4.1186e-06 4.2762e-06 4.4292e-06 4.5787e-06 4.8703e-06 5.1555e-06 5.4369e-06 5.7165e-06 5.9956e-06; ...
    3.6726e-06 3.8488e-06 4.0139e-06 4.1714e-06 4.3234e-06 4.4711e-06 4.7577e-06 5.0364e-06 5.3103e-06 5.5814e-06 5.8513e-06; ...
    3.5612e-06 3.7429e-06 3.9106e-06 4.0689e-06 4.2205e-06 4.3670e-06 4.6496e-06 4.9227e-06 5.1898e-06 5.4534e-06 5.7150e-06; ...
    3.4473e-06 3.6368e-06 3.8084e-06 3.9684e-06 4.1203e-06 4.2662e-06 4.5457e-06 4.8141e-06 5.0753e-06 5.3321e-06 5.5862e-06; ...
    3.3298e-06 3.5300e-06 3.7070e-06 3.8697e-06 4.0227e-06 4.1686e-06 4.4460e-06 4.7103e-06 4.9665e-06 5.2172e-06 5.4646e-06; ...
    3.2075e-06 3.4221e-06 3.6063e-06 3.7727e-06 3.9275e-06 4.0740e-06 4.3501e-06 4.6113e-06 4.8629e-06 5.1083e-06 5.3496e-06; ...
    3.0794e-06 3.3128e-06 3.5060e-06 3.6773e-06 3.8346e-06 3.9822e-06 4.2580e-06 4.5166e-06 4.7644e-06 5.0050e-06 5.2408e-06; ...
    2.9450e-06 3.2022e-06 3.4063e-06 3.5835e-06 3.7440e-06 3.8933e-06 4.1694e-06 4.4261e-06 4.6707e-06 4.9071e-06 5.1380e-06; ...
    2.8040e-06 3.0911e-06 3.3077e-06 3.4914e-06 3.6557e-06 3.8070e-06 4.0842e-06 4.3397e-06 4.5815e-06 4.8142e-06 5.0406e-06; ...
    2.6605e-06 2.9809e-06 3.2109e-06 3.4016e-06 3.5700e-06 3.7236e-06 4.0024e-06 4.2571e-06 4.4966e-06 4.7260e-06 4.9485e-06; ...
    2.5417e-06 2.8733e-06 3.1173e-06 3.3149e-06 3.4873e-06 3.6433e-06 3.9240e-06 4.1782e-06 4.4158e-06 4.6423e-06 4.8613e-06; ...
    2.4691e-06 2.7730e-06 3.0283e-06 3.2323e-06 3.4082e-06 3.5664e-06 3.8491e-06 4.1030e-06 4.3390e-06 4.5630e-06 4.7787e-06; ...
    2.4296e-06 2.6903e-06 2.9456e-06 3.1548e-06 3.3336e-06 3.4935e-06 3.7777e-06 4.0315e-06 4.2660e-06 4.4878e-06 4.7005e-06; ...
    2.4094e-06 2.6325e-06 2.8720e-06 3.0834e-06 3.2642e-06 3.4252e-06 3.7102e-06 3.9636e-06 4.1969e-06 4.4165e-06 4.6266e-06; ...
    2.4012e-06 2.5966e-06 2.8118e-06 3.0192e-06 3.2006e-06 3.3619e-06 3.6468e-06 3.8995e-06 4.1314e-06 4.3491e-06 4.5567e-06; ...
    2.4006e-06 2.5759e-06 2.7673e-06 2.9638e-06 3.1434e-06 3.3042e-06 3.5878e-06 3.8393e-06 4.0698e-06 4.2855e-06 4.4907e-06; ...
    2.4055e-06 2.5656e-06 2.7374e-06 2.9186e-06 3.0932e-06 3.2523e-06 3.5336e-06 3.7832e-06 4.0118e-06 4.2256e-06 4.4284e-06; ...
    2.4143e-06 2.5624e-06 2.7188e-06 2.8844e-06 3.0506e-06 3.2065e-06 3.4843e-06 3.7313e-06 3.9577e-06 4.1693e-06 4.3699e-06; ...
    2.4478e-06 2.5744e-06 2.7048e-06 2.8397e-06 2.9794e-06 3.1198e-06 3.3831e-06 3.6205e-06 3.8394e-06 4.0446e-06 4.2391e-06; ...
    2.4916e-06 2.6034e-06 2.7174e-06 2.8332e-06 2.9514e-06 3.0723e-06 3.3124e-06 3.5369e-06 3.7457e-06 3.9426e-06 4.1300e-06; ...
    2.5415e-06 2.6422e-06 2.7445e-06 2.8475e-06 2.9512e-06 3.0565e-06 3.2701e-06 3.4785e-06 3.6756e-06 3.8628e-06 4.0419e-06; ...
    2.5952e-06 2.6871e-06 2.7804e-06 2.8740e-06 2.9676e-06 3.0615e-06 3.2519e-06 3.4424e-06 3.6268e-06 3.8035e-06 3.9735e-06; ...
    2.6512e-06 2.7360e-06 2.8220e-06 2.9083e-06 2.9942e-06 3.0799e-06 3.2520e-06 3.4253e-06 3.5965e-06 3.7625e-06 3.9232e-06; ...
    2.7087e-06 2.7877e-06 2.8676e-06 2.9478e-06 3.0276e-06 3.1070e-06 3.2650e-06 3.4235e-06 3.5817e-06 3.7372e-06 3.8888e-06; ...
    2.7672e-06 2.8412e-06 2.9160e-06 2.9910e-06 3.0657e-06 3.1400e-06 3.2869e-06 3.4334e-06 3.5799e-06 3.7251e-06 3.8678e-06; ...
    2.8261e-06 2.8959e-06 2.9663e-06 3.0369e-06 3.1073e-06 3.1772e-06 3.3152e-06 3.4518e-06 3.5882e-06 3.7240e-06 3.8583e-06; ...
    2.8853e-06 2.9514e-06 3.0180e-06 3.0848e-06 3.1514e-06 3.2175e-06 3.3480e-06 3.4766e-06 3.6045e-06 3.7319e-06 3.8584e-06; ...
    2.9445e-06 3.0074e-06 3.0707e-06 3.1341e-06 3.1973e-06 3.2602e-06 3.3842e-06 3.5061e-06 3.6268e-06 3.7468e-06 3.8663e-06; ...
    3.0036e-06 3.0637e-06 3.1240e-06 3.1844e-06 3.2447e-06 3.3046e-06 3.4230e-06 3.5392e-06 3.6538e-06 3.7675e-06 3.8806e-06; ...
    3.0625e-06 3.1201e-06 3.1778e-06 3.2355e-06 3.2931e-06 3.3504e-06 3.4637e-06 3.5749e-06 3.6843e-06 3.7926e-06 3.9002e-06; ...
    3.1795e-06 3.2327e-06 3.2860e-06 3.3391e-06 3.3922e-06 3.4450e-06 3.5495e-06 3.6523e-06 3.7532e-06 3.8526e-06 3.9511e-06; ...
    3.2951e-06 3.3447e-06 3.3942e-06 3.4436e-06 3.4929e-06 3.5420e-06 3.6393e-06 3.7350e-06 3.8291e-06 3.9217e-06 4.0130e-06; ...
    3.4091e-06 3.4556e-06 3.5020e-06 3.5483e-06 3.5944e-06 3.6403e-06 3.7315e-06 3.8213e-06 3.9097e-06 3.9966e-06 4.0823e-06; ...
    3.5214e-06 3.5652e-06 3.6089e-06 3.6525e-06 3.6960e-06 3.7392e-06 3.8251e-06 3.9098e-06 3.9933e-06 4.0755e-06 4.1564e-06; ...
    3.6319e-06 3.6734e-06 3.7148e-06 3.7561e-06 3.7972e-06 3.8382e-06 3.9194e-06 3.9998e-06 4.0790e-06 4.1570e-06 4.2339e-06; ...
    3.7407e-06 3.7802e-06 3.8196e-06 3.8588e-06 3.8979e-06 3.9368e-06 4.0141e-06 4.0905e-06 4.1660e-06 4.2404e-06 4.3138e-06; ...
    3.9533e-06 3.9894e-06 4.0254e-06 4.0612e-06 4.0969e-06 4.1324e-06 4.2030e-06 4.2729e-06 4.3421e-06 4.4104e-06 4.4779e-06; ...
    4.1596e-06 4.1929e-06 4.2261e-06 4.2592e-06 4.2922e-06 4.3250e-06 4.3903e-06 4.4549e-06 4.5189e-06 4.5823e-06 4.6450e-06; ...
    4.3599e-06 4.3910e-06 4.4220e-06 4.4528e-06 4.4835e-06 4.5141e-06 4.5750e-06 4.6353e-06 4.6951e-06 4.7544e-06 4.8130e-06; ...
    4.5548e-06 4.5840e-06 4.6131e-06 4.6420e-06 4.6709e-06 4.6996e-06 4.7568e-06 4.8135e-06 4.8697e-06 4.9255e-06 4.9808e-06; ...
    4.7447e-06 4.7722e-06 4.7997e-06 4.8270e-06 4.8543e-06 4.8815e-06 4.9355e-06 4.9891e-06 5.0424e-06 5.0952e-06 5.1476e-06; ...
    5.1998e-06 5.2242e-06 5.2485e-06 5.2728e-06 5.2969e-06 5.3210e-06 5.3690e-06 5.4166e-06 5.4640e-06 5.5110e-06 5.5577e-06; ...
    5.6310e-06 5.6531e-06 5.6751e-06 5.6971e-06 5.7190e-06 5.7409e-06 5.7844e-06 5.8277e-06 5.8708e-06 5.9135e-06 5.9561e-06; ...
    6.0423e-06 6.0626e-06 6.0830e-06 6.1032e-06 6.1235e-06 6.1436e-06 6.1838e-06 6.2237e-06 6.2635e-06 6.3030e-06 6.3423e-06; ...
    6.4368e-06 6.4558e-06 6.4748e-06 6.4937e-06 6.5126e-06 6.5314e-06 6.5689e-06 6.6062e-06 6.6434e-06 6.6803e-06 6.7171e-06];
t_kf = [ ...
    2.1120e-02 2.1575e-02 2.1999e-02 2.2396e-02 2.2771e-02 2.3129e-02 2.3801e-02 2.4425e-02 2.5012e-02 2.5568e-02 2.6098e-02; ...
    2.1399e-02 2.1899e-02 2.2361e-02 2.2792e-02 2.3198e-02 2.3583e-02 2.4302e-02 2.4969e-02 2.5593e-02 2.6184e-02 2.6746e-02; ...
    2.1605e-02 2.2155e-02 2.2658e-02 2.3126e-02 2.3563e-02 2.3977e-02 2.4747e-02 2.5457e-02 2.6120e-02 2.6746e-02 2.7341e-02; ...
    2.1741e-02 2.2343e-02 2.2891e-02 2.3397e-02 2.3869e-02 2.4313e-02 2.5136e-02 2.5891e-02 2.6594e-02 2.7255e-02 2.7883e-02; ...
    2.1811e-02 2.2467e-02 2.3061e-02 2.3608e-02 2.4116e-02 2.4592e-02 2.5470e-02 2.6271e-02 2.7015e-02 2.7713e-02 2.8373e-02; ...
    2.1826e-02 2.2531e-02 2.3171e-02 2.3759e-02 2.4304e-02 2.4814e-02 2.5750e-02 2.6599e-02 2.7385e-02 2.8120e-02 2.8814e-02; ...
    2.1792e-02 2.2542e-02 2.3225e-02 2.3855e-02 2.4438e-02 2.4982e-02 2.5977e-02 2.6877e-02 2.7705e-02 2.8477e-02 2.9206e-02; ...
    2.1705e-02 2.2502e-02 2.3228e-02 2.3897e-02 2.4518e-02 2.5097e-02 2.6154e-02 2.7104e-02 2.7976e-02 2.8788e-02 2.9550e-02; ...
    2.1524e-02 2.2404e-02 2.3178e-02 2.3888e-02 2.4548e-02 2.5162e-02 2.6281e-02 2.7285e-02 2.8202e-02 2.9052e-02 2.9850e-02; ...
    2.1170e-02 2.2218e-02 2.3067e-02 2.3827e-02 2.4528e-02 2.5179e-02 2.6363e-02 2.7420e-02 2.8383e-02 2.9273e-02 3.0106e-02; ...
    2.0562e-02 2.1906e-02 2.2878e-02 2.3707e-02 2.4457e-02 2.5149e-02 2.6399e-02 2.7512e-02 2.8521e-02 2.9452e-02 3.0320e-02; ...
    1.9688e-02 2.1445e-02 2.2596e-02 2.3521e-02 2.4333e-02 2.5071e-02 2.6394e-02 2.7563e-02 2.8620e-02 2.9591e-02 3.0495e-02; ...
    1.8762e-02 2.0854e-02 2.2220e-02 2.3265e-02 2.4153e-02 2.4946e-02 2.6347e-02 2.7576e-02 2.8681e-02 2.9693e-02 3.0632e-02; ...
    1.8125e-02 2.0203e-02 2.1774e-02 2.2950e-02 2.3923e-02 2.4777e-02 2.6263e-02 2.7553e-02 2.8706e-02 2.9759e-02 3.0735e-02; ...
    1.7770e-02 1.9616e-02 2.1300e-02 2.2596e-02 2.3654e-02 2.4570e-02 2.6144e-02 2.7496e-02 2.8699e-02 2.9794e-02 3.0805e-02; ...
    1.7596e-02 1.9201e-02 2.0850e-02 2.2232e-02 2.3364e-02 2.4338e-02 2.5997e-02 2.7411e-02 2.8663e-02 2.9798e-02 3.0844e-02; ...
    1.7534e-02 1.8953e-02 2.0475e-02 2.1884e-02 2.3071e-02 2.4094e-02 2.5829e-02 2.7302e-02 2.8601e-02 2.9775e-02 3.0856e-02; ...
    1.7542e-02 1.8826e-02 2.0204e-02 2.1576e-02 2.2793e-02 2.3851e-02 2.5648e-02 2.7173e-02 2.8516e-02 2.9729e-02 3.0843e-02; ...
    1.7597e-02 1.8780e-02 2.0033e-02 2.1326e-02 2.2541e-02 2.3619e-02 2.5463e-02 2.7031e-02 2.8414e-02 2.9662e-02 3.0808e-02; ...
    1.7683e-02 1.8786e-02 1.9941e-02 2.1143e-02 2.2324e-02 2.3407e-02 2.5280e-02 2.6881e-02 2.8298e-02 2.9578e-02 3.0753e-02; ...
    1.7979e-02 1.8938e-02 1.9927e-02 2.0939e-02 2.1972e-02 2.2994e-02 2.4867e-02 2.6507e-02 2.7977e-02 2.9316e-02 3.0551e-02; ...
    1.8337e-02 1.9198e-02 2.0081e-02 2.0972e-02 2.1874e-02 2.2784e-02 2.4558e-02 2.6181e-02 2.7660e-02 2.9025e-02 3.0296e-02; ...
    1.8729e-02 1.9513e-02 2.0321e-02 2.1132e-02 2.1943e-02 2.2757e-02 2.4383e-02 2.5942e-02 2.7394e-02 2.8752e-02 3.0031e-02; ...
    1.9138e-02 1.9863e-02 2.0610e-02 2.1362e-02 2.2110e-02 2.2854e-02 2.4339e-02 2.5802e-02 2.7202e-02 2.8529e-02 2.9792e-02; ...
    1.9557e-02 2.0233e-02 2.0931e-02 2.1635e-02 2.2335e-02 2.3029e-02 2.4402e-02 2.5763e-02 2.7093e-02 2.8372e-02 2.9602e-02; ...
    1.9980e-02 2.0615e-02 2.1272e-02 2.1936e-02 2.2598e-02 2.3253e-02 2.4540e-02 2.5811e-02 2.7064e-02 2.8287e-02 2.9474e-02; ...
    2.0405e-02 2.1005e-02 2.1627e-02 2.2257e-02 2.2886e-02 2.3509e-02 2.4730e-02 2.5926e-02 2.7107e-02 2.8270e-02 2.9407e-02; ...
    2.0830e-02 2.1400e-02 2.1991e-02 2.2591e-02 2.3192e-02 2.3788e-02 2.4955e-02 2.6091e-02 2.7209e-02 2.8313e-02 2.9400e-02; ...
    2.1254e-02 2.1798e-02 2.2361e-02 2.2935e-02 2.3511e-02 2.4083e-02 2.5204e-02 2.6292e-02 2.7357e-02 2.8408e-02 2.9446e-02; ...
    2.1676e-02 2.2196e-02 2.2735e-02 2.3285e-02 2.3839e-02 2.4390e-02 2.5471e-02 2.6520e-02 2.7541e-02 2.8545e-02 2.9537e-02; ...
    2.2095e-02 2.2594e-02 2.3111e-02 2.3640e-02 2.4173e-02 2.4705e-02 2.5752e-02 2.6766e-02 2.7751e-02 2.8716e-02 2.9667e-02; ...
    2.2511e-02 2.2991e-02 2.3489e-02 2.3999e-02 2.4513e-02 2.5028e-02 2.6042e-02 2.7027e-02 2.7981e-02 2.8912e-02 2.9828e-02; ...
    2.3333e-02 2.3781e-02 2.4245e-02 2.4721e-02 2.5203e-02 2.5686e-02 2.6644e-02 2.7579e-02 2.8483e-02 2.9361e-02 3.0220e-02; ...
    2.4144e-02 2.4563e-02 2.4999e-02 2.5446e-02 2.5900e-02 2.6357e-02 2.7266e-02 2.8157e-02 2.9022e-02 2.9860e-02 3.0675e-02; ...
    2.4942e-02 2.5337e-02 2.5748e-02 2.6171e-02 2.6601e-02 2.7034e-02 2.7900e-02 2.8753e-02 2.9583e-02 3.0389e-02 3.1171e-02; ...
    2.5728e-02 2.6102e-02 2.6492e-02 2.6893e-02 2.7301e-02 2.7714e-02 2.8542e-02 2.9360e-02 3.0160e-02 3.0938e-02 3.1694e-02; ...
    2.6504e-02 2.6859e-02 2.7230e-02 2.7611e-02 2.8001e-02 2.8395e-02 2.9188e-02 2.9975e-02 3.0748e-02 3.1501e-02 3.2233e-02; ...
    2.7269e-02 2.7607e-02 2.7961e-02 2.8325e-02 2.8698e-02 2.9075e-02 2.9837e-02 3.0596e-02 3.1343e-02 3.2074e-02 3.2785e-02; ...
    2.8770e-02 2.9079e-02 2.9404e-02 2.9739e-02 3.0082e-02 3.0430e-02 3.1137e-02 3.1845e-02 3.2547e-02 3.3237e-02 3.3912e-02; ...
    3.0235e-02 3.0521e-02 3.0821e-02 3.1132e-02 3.1450e-02 3.1774e-02 3.2434e-02 3.3099e-02 3.3761e-02 3.4415e-02 3.5058e-02; ...
    3.1670e-02 3.1936e-02 3.2215e-02 3.2504e-02 3.2802e-02 3.3105e-02 3.3724e-02 3.4351e-02 3.4978e-02 3.5600e-02 3.6213e-02; ...
    3.3076e-02 3.3325e-02 3.3586e-02 3.3857e-02 3.4136e-02 3.4422e-02 3.5005e-02 3.5598e-02 3.6193e-02 3.6786e-02 3.7374e-02; ...
    3.4456e-02 3.4690e-02 3.4936e-02 3.5191e-02 3.5454e-02 3.5724e-02 3.6276e-02 3.6838e-02 3.7406e-02 3.7972e-02 3.8535e-02; ...
    3.7806e-02 3.8010e-02 3.8225e-02 3.8449e-02 3.8680e-02 3.8917e-02 3.9404e-02 3.9904e-02 4.0411e-02 4.0920e-02 4.1430e-02; ...
    4.1032e-02 4.1214e-02 4.1406e-02 4.1606e-02 4.1812e-02 4.2024e-02 4.2461e-02 4.2911e-02 4.3370e-02 4.3833e-02 4.4299e-02; ...
    4.4153e-02 4.4318e-02 4.4492e-02 4.4673e-02 4.4860e-02 4.5053e-02 4.5450e-02 4.5860e-02 4.6279e-02 4.6704e-02 4.7133e-02; ...
    4.7181e-02 4.7334e-02 4.7494e-02 4.7660e-02 4.7832e-02 4.8008e-02 4.8373e-02 4.8751e-02 4.9137e-02 4.9529e-02 4.9926e-02];
t_h = [ ...
    1.3313e+03 1.8048e+03 2.2937e+03 2.7936e+03 3.3013e+03 3.8149e+03 4.8536e+03 5.9016e+03 6.9541e+03 8.0082e+03 9.0619e+03; ...
    2.0411e+03 2.4824e+03 2.9457e+03 3.4244e+03 3.9142e+03 4.4122e+03 5.4253e+03 6.4525e+03 7.4877e+03 8.5269e+03 9.5675e+03; ...
    2.8105e+03 3.2111e+03 3.6431e+03 4.0965e+03 4.5653e+03 5.0454e+03 6.0292e+03 7.0332e+03 8.0493e+03 9.0722e+03 1.0099e+04; ...
    3.6489e+03 3.9972e+03 4.3903e+03 4.8132e+03 5.2571e+03 5.7164e+03 6.6667e+03 7.6446e+03 8.6394e+03 9.6444e+03 1.0656e+04; ...
    4.5688e+03 4.8486e+03 5.1927e+03 5.5783e+03 5.9926e+03 6.4273e+03 7.3391e+03 8.2876e+03 9.2588e+03 1.0244e+04 1.1239e+04; ...
    5.5879e+03 5.7753e+03 6.0570e+03 6.3966e+03 6.7750e+03 7.1809e+03 8.0482e+03 8.9633e+03 9.9081e+03 1.0872e+04 1.1848e+04; ...
    6.7308e+03 6.7907e+03 6.9913e+03 7.2736e+03 7.6086e+03 7.9802e+03 8.7957e+03 9.6730e+03 1.0588e+04 1.1528e+04 1.2485e+04; ...
    8.0325e+03 7.9119e+03 8.0060e+03 8.2161e+03 8.4981e+03 8.8286e+03 9.5837e+03 1.0418e+04 1.1300e+04 1.2214e+04 1.3148e+04; ...
    9.5430e+03 9.1610e+03 9.1133e+03 9.2319e+03 9.4489e+03 9.7302e+03 1.0415e+04 1.1200e+04 1.2045e+04 1.2929e+04 1.3840e+04; ...
    1.1335e+04 1.0565e+04 1.0328e+04 1.0330e+04 1.0467e+04 1.0689e+04 1.1291e+04 1.2020e+04 1.2824e+04 1.3675e+04 1.4560e+04; ...
    1.3520e+04 1.2155e+04 1.1666e+04 1.1521e+04 1.1560e+04 1.1711e+04 1.2216e+04 1.2880e+04 1.3637e+04 1.4453e+04 1.5308e+04; ...
    1.6222e+04 1.3970e+04 1.3143e+04 1.2814e+04 1.2733e+04 1.2800e+04 1.3191e+04 1.3782e+04 1.4487e+04 1.5263e+04 1.6087e+04; ...
    1.9220e+04 1.6046e+04 1.4775e+04 1.4218e+04 1.3994e+04 1.3961e+04 1.4221e+04 1.4728e+04 1.5373e+04 1.6106e+04 1.6895e+04; ...
    2.1982e+04 1.8389e+04 1.6573e+04 1.5739e+04 1.5347e+04 1.5198e+04 1.5307e+04 1.5719e+04 1.6299e+04 1.6982e+04 1.7733e+04; ...
    2.4404e+04 2.0866e+04 1.8543e+04 1.7382e+04 1.6795e+04 1.6513e+04 1.6451e+04 1.6757e+04 1.7264e+04 1.7893e+04 1.8603e+04; ...
    2.6572e+04 2.3263e+04 2.0656e+04 1.9145e+04 1.8337e+04 1.7908e+04 1.7656e+04 1.7843e+04 1.8269e+04 1.8840e+04 1.9504e+04; ...
    2.8560e+04 2.5482e+04 2.2832e+04 2.1017e+04 1.9971e+04 1.9380e+04 1.8921e+04 1.8978e+04 1.9316e+04 1.9822e+04 2.0437e+04; ...
    3.0416e+04 2.7535e+04 2.4971e+04 2.2969e+04 2.1689e+04 2.0926e+04 2.0244e+04 2.0163e+04 2.0405e+04 2.0841e+04 2.1403e+04; ...
    3.2173e+04 2.9454e+04 2.7013e+04 2.4951e+04 2.3476e+04 2.2539e+04 2.1625e+04 2.1396e+04 2.1535e+04 2.1896e+04 2.2401e+04; ...
    3.3852e+04 3.1270e+04 2.8948e+04 2.6910e+04 2.5309e+04 2.4210e+04 2.3058e+04 2.2675e+04 2.2706e+04 2.2988e+04 2.3432e+04; ...
    3.7798e+04 3.5484e+04 3.3398e+04 3.1530e+04 2.9885e+04 2.8539e+04 2.6837e+04 2.6059e+04 2.5805e+04 2.5871e+04 2.6149e+04; ...
    4.1483e+04 3.9383e+04 3.7467e+04 3.5750e+04 3.4210e+04 3.2848e+04 3.0806e+04 2.9659e+04 2.9117e+04 2.8959e+04 2.9057e+04; ...
    4.4985e+04 4.3065e+04 4.1288e+04 3.9686e+04 3.8249e+04 3.6955e+04 3.4823e+04 3.3400e+04 3.2595e+04 3.2218e+04 3.2134e+04; ...
    4.8348e+04 4.6586e+04 4.4932e+04 4.3425e+04 4.2069e+04 4.0848e+04 3.8763e+04 3.7199e+04 3.6182e+04 3.5606e+04 3.5350e+04; ...
    5.1603e+04 4.9982e+04 4.8440e+04 4.7018e+04 4.5731e+04 4.4570e+04 4.2571e+04 4.0976e+04 3.9821e+04 3.9081e+04 3.8670e+04; ...
    5.4774e+04 5.3275e+04 5.1838e+04 5.0497e+04 4.9271e+04 4.8161e+04 4.6247e+04 4.4679e+04 4.3460e+04 4.2601e+04 4.2061e+04; ...
    5.7875e+04 5.6486e+04 5.5145e+04 5.3880e+04 5.2713e+04 5.1649e+04 4.9810e+04 4.8291e+04 4.7061e+04 4.6131e+04 4.5490e+04; ...
    6.0920e+04 5.9628e+04 5.8375e+04 5.7183e+04 5.6073e+04 5.5054e+04 5.3282e+04 5.1814e+04 5.0602e+04 4.9641e+04 4.8932e+04; ...
    6.3918e+04 6.2713e+04 6.1540e+04 6.0417e+04 5.9363e+04 5.8387e+04 5.6678e+04 5.5258e+04 5.4076e+04 5.3112e+04 5.2364e+04; ...
    6.6877e+04 6.5750e+04 6.4649e+04 6.3592e+04 6.2591e+04 6.1658e+04 6.0010e+04 5.8635e+04 5.7487e+04 5.6535e+04 5.5771e+04; ...
    6.9802e+04 6.8745e+04 6.7711e+04 6.6714e+04 6.5765e+04 6.4874e+04 6.3287e+04 6.1954e+04 6.0839e+04 5.9907e+04 5.9143e+04; ...
    7.2698e+04 7.1704e+04 7.0732e+04 6.9791e+04 6.8891e+04 6.8041e+04 6.6515e+04 6.5223e+04 6.4139e+04 6.3230e+04 6.2476e+04; ...
    7.8418e+04 7.7536e+04 7.6670e+04 7.5830e+04 7.5020e+04 7.4250e+04 7.2844e+04 7.1634e+04 7.0610e+04 6.9747e+04 6.9024e+04; ...
    8.4062e+04 8.3272e+04 8.2498e+04 8.1744e+04 8.1016e+04 8.0317e+04 7.9028e+04 7.7900e+04 7.6935e+04 7.6118e+04 7.5430e+04; ...
    8.9645e+04 8.8935e+04 8.8240e+04 8.7561e+04 8.6903e+04 8.6270e+04 8.5092e+04 8.4046e+04 8.3140e+04 8.2367e+04 8.1715e+04; ...
    9.5180e+04 9.4539e+04 9.3912e+04 9.3299e+04 9.2704e+04 9.2130e+04 9.1054e+04 9.0089e+04 8.9243e+04 8.8516e+04 8.7898e+04; ...
    1.0068e+05 1.0010e+05 9.9528e+04 9.8973e+04 9.8434e+04 9.7913e+04 9.6931e+04 9.6043e+04 9.5258e+04 9.4576e+04 9.3995e+04; ...
    1.0614e+05 1.0561e+05 1.0510e+05 1.0459e+05 1.0410e+05 1.0363e+05 1.0274e+05 1.0192e+05 1.0120e+05 1.0056e+05 1.0001e+05; ...
    1.1699e+05 1.1655e+05 1.1613e+05 1.1571e+05 1.1531e+05 1.1491e+05 1.1417e+05 1.1349e+05 1.1288e+05 1.1233e+05 1.1186e+05; ...
    1.2775e+05 1.2739e+05 1.2704e+05 1.2670e+05 1.2636e+05 1.2604e+05 1.2543e+05 1.2486e+05 1.2435e+05 1.2389e+05 1.2349e+05; ...
    1.3846e+05 1.3816e+05 1.3787e+05 1.3759e+05 1.3731e+05 1.3705e+05 1.3654e+05 1.3608e+05 1.3566e+05 1.3528e+05 1.3495e+05; ...
    1.4912e+05 1.4888e+05 1.4864e+05 1.4840e+05 1.4818e+05 1.4796e+05 1.4755e+05 1.4717e+05 1.4683e+05 1.4653e+05 1.4626e+05; ...
    1.5974e+05 1.5954e+05 1.5935e+05 1.5916e+05 1.5898e+05 1.5881e+05 1.5848e+05 1.5818e+05 1.5791e+05 1.5767e+05 1.5746e+05; ...
    1.8618e+05 1.8607e+05 1.8597e+05 1.8587e+05 1.8577e+05 1.8568e+05 1.8552e+05 1.8538e+05 1.8525e+05 1.8515e+05 1.8507e+05; ...
    2.1250e+05 2.1246e+05 2.1242e+05 2.1239e+05 2.1236e+05 2.1233e+05 2.1229e+05 2.1227e+05 2.1226e+05 2.1226e+05 2.1229e+05; ...
    2.3875e+05 2.3876e+05 2.3877e+05 2.3879e+05 2.3881e+05 2.3883e+05 2.3889e+05 2.3895e+05 2.3904e+05 2.3913e+05 2.3923e+05; ...
    2.6493e+05 2.6498e+05 2.6504e+05 2.6510e+05 2.6516e+05 2.6522e+05 2.6536e+05 2.6550e+05 2.6565e+05 2.6582e+05 2.6599e+05];
t_phi = [ ...
    1.4221e+00 1.4324e+00 1.4412e+00 1.4493e+00 1.4570e+00 1.4645e+00 1.4793e+00 1.4943e+00 1.5094e+00 1.5248e+00 1.5405e+00; ...
    1.4091e+00 1.4223e+00 1.4334e+00 1.4432e+00 1.4523e+00 1.4610e+00 1.4776e+00 1.4938e+00 1.5099e+00 1.5259e+00 1.5419e+00; ...
    1.3915e+00 1.4083e+00 1.4220e+00 1.4337e+00 1.4444e+00 1.4543e+00 1.4728e+00 1.4904e+00 1.5075e+00 1.5242e+00 1.5406e+00; ...
    1.3692e+00 1.3905e+00 1.4071e+00 1.4211e+00 1.4335e+00 1.4448e+00 1.4654e+00 1.4845e+00 1.5026e+00 1.5200e+00 1.5369e+00; ...
    1.3417e+00 1.3687e+00 1.3890e+00 1.4056e+00 1.4199e+00 1.4327e+00 1.4556e+00 1.4762e+00 1.4955e+00 1.5137e+00 1.5312e+00; ...
    1.3083e+00 1.3428e+00 1.3676e+00 1.3872e+00 1.4037e+00 1.4183e+00 1.4437e+00 1.4660e+00 1.4864e+00 1.5056e+00 1.5237e+00; ...
    1.2678e+00 1.3124e+00 1.3428e+00 1.3661e+00 1.3852e+00 1.4017e+00 1.4298e+00 1.4540e+00 1.4757e+00 1.4958e+00 1.5146e+00; ...
    1.2190e+00 1.2771e+00 1.3146e+00 1.3422e+00 1.3643e+00 1.3829e+00 1.4141e+00 1.4403e+00 1.4635e+00 1.4846e+00 1.5041e+00; ...
    1.1604e+00 1.2363e+00 1.2826e+00 1.3155e+00 1.3411e+00 1.3622e+00 1.3968e+00 1.4252e+00 1.4499e+00 1.4721e+00 1.4925e+00; ...
    1.0923e+00 1.1897e+00 1.2469e+00 1.2860e+00 1.3156e+00 1.3396e+00 1.3779e+00 1.4087e+00 1.4351e+00 1.4585e+00 1.4798e+00; ...
    1.0185e+00 1.1381e+00 1.2075e+00 1.2538e+00 1.2880e+00 1.3152e+00 1.3577e+00 1.3911e+00 1.4193e+00 1.4440e+00 1.4662e+00; ...
    9.4959e-01 1.0831e+00 1.1650e+00 1.2190e+00 1.2583e+00 1.2890e+00 1.3361e+00 1.3723e+00 1.4024e+00 1.4286e+00 1.4518e+00; ...
    9.0133e-01 1.0282e+00 1.1205e+00 1.1822e+00 1.2268e+00 1.2613e+00 1.3133e+00 1.3526e+00 1.3848e+00 1.4124e+00 1.4367e+00; ...
    8.7264e-01 9.7804e-01 1.0756e+00 1.1442e+00 1.1940e+00 1.2323e+00 1.2894e+00 1.3319e+00 1.3663e+00 1.3955e+00 1.4210e+00; ...
    8.5270e-01 9.3791e-01 1.0324e+00 1.1060e+00 1.1604e+00 1.2024e+00 1.2646e+00 1.3104e+00 1.3471e+00 1.3779e+00 1.4047e+00; ...
    8.3697e-01 9.0854e-01 9.9303e-01 1.0687e+00 1.1268e+00 1.1720e+00 1.2391e+00 1.2883e+00 1.3273e+00 1.3599e+00 1.3880e+00; ...
    8.2399e-01 8.8651e-01 9.5937e-01 1.0332e+00 1.0938e+00 1.1417e+00 1.2132e+00 1.2656e+00 1.3070e+00 1.3413e+00 1.3708e+00; ...
    8.1312e-01 8.6885e-01 9.3205e-01 1.0007e+00 1.0620e+00 1.1119e+00 1.1872e+00 1.2426e+00 1.2863e+00 1.3224e+00 1.3532e+00; ...
    8.0389e-01 8.5407e-01 9.1003e-01 9.7189e-01 1.0319e+00 1.0830e+00 1.1613e+00 1.2194e+00 1.2653e+00 1.3032e+00 1.3354e+00; ...
    7.9597e-01 8.4143e-01 8.9187e-01 9.4721e-01 1.0042e+00 1.0554e+00 1.1360e+00 1.1963e+00 1.2442e+00 1.2837e+00 1.3173e+00; ...
    7.8034e-01 8.1675e-01 8.5723e-01 9.0061e-01 9.4695e-01 9.9394e-01 1.0761e+00 1.1401e+00 1.1918e+00 1.2348e+00 1.2715e+00; ...
    7.6869e-01 7.9891e-01 8.3239e-01 8.6814e-01 9.0588e-01 9.4545e-01 1.0229e+00 1.0881e+00 1.1416e+00 1.1869e+00 1.2259e+00; ...
    7.5959e-01 7.8547e-01 8.1385e-01 8.4412e-01 8.7586e-01 9.0904e-01 9.7776e-01 1.0413e+00 1.0953e+00 1.1415e+00 1.1819e+00; ...
    7.5220e-01 7.7492e-01 7.9954e-01 8.2569e-01 8.5303e-01 8.8144e-01 9.4111e-01 1.0004e+00 1.0535e+00 1.0997e+00 1.1405e+00; ...
    7.4604e-01 7.6635e-01 7.8814e-01 8.1114e-01 8.3512e-01 8.5993e-01 9.1194e-01 9.6553e-01 1.0163e+00 1.0618e+00 1.1023e+00; ...
    7.4080e-01 7.5921e-01 7.7879e-01 7.9936e-01 8.2070e-01 8.4271e-01 8.8861e-01 9.3650e-01 9.8388e-01 1.0279e+00 1.0676e+00; ...
    7.3626e-01 7.5313e-01 7.7095e-01 7.8957e-01 8.0883e-01 8.2862e-01 8.6965e-01 9.1250e-01 9.5595e-01 9.9773e-01 1.0364e+00; ...
    7.3229e-01 7.4787e-01 7.6425e-01 7.8128e-01 7.9885e-01 8.1685e-01 8.5397e-01 8.9257e-01 9.3216e-01 9.7125e-01 1.0083e+00; ...
    7.2876e-01 7.4326e-01 7.5842e-01 7.7414e-01 7.9031e-01 8.0683e-01 8.4076e-01 8.7585e-01 9.1194e-01 9.4819e-01 9.8332e-01; ...
    7.2562e-01 7.3918e-01 7.5330e-01 7.6790e-01 7.8289e-01 7.9818e-01 8.2947e-01 8.6164e-01 8.9468e-01 9.2817e-01 9.6120e-01; ...
    7.2278e-01 7.3553e-01 7.4875e-01 7.6239e-01 7.7637e-01 7.9061e-01 8.1967e-01 8.4940e-01 8.7984e-01 9.1078e-01 9.4168e-01; ...
    7.2021e-01 7.3224e-01 7.4468e-01 7.5748e-01 7.7058e-01 7.8391e-01 8.1106e-01 8.3873e-01 8.6694e-01 8.9563e-01 9.2446e-01; ...
    7.1572e-01 7.2653e-01 7.3766e-01 7.4907e-01 7.6072e-01 7.7254e-01 7.9658e-01 8.2095e-01 8.4561e-01 8.7059e-01 8.9580e-01; ...
    7.1191e-01 7.2174e-01 7.3182e-01 7.4211e-01 7.5259e-01 7.6323e-01 7.8481e-01 8.0663e-01 8.2861e-01 8.5075e-01 8.7305e-01; ...
    7.0865e-01 7.1766e-01 7.2686e-01 7.3624e-01 7.4577e-01 7.5543e-01 7.7500e-01 7.9478e-01 8.1465e-01 8.3457e-01 8.5458e-01; ...
    7.0581e-01 7.1412e-01 7.2260e-01 7.3121e-01 7.3994e-01 7.4878e-01 7.6669e-01 7.8478e-01 8.0292e-01 8.2107e-01 8.3924e-01; ...
    7.0331e-01 7.1103e-01 7.1888e-01 7.2684e-01 7.3490e-01 7.4305e-01 7.5955e-01 7.7620e-01 7.9289e-01 8.0958e-01 8.2625e-01; ...
    7.0110e-01 7.0831e-01 7.1561e-01 7.2301e-01 7.3050e-01 7.3805e-01 7.5333e-01 7.6875e-01 7.8421e-01 7.9966e-01 8.1507e-01; ...
    6.9735e-01 7.0371e-01 7.1013e-01 7.1661e-01 7.2315e-01 7.2974e-01 7.4305e-01 7.5645e-01 7.6991e-01 7.8335e-01 7.9674e-01; ...
    6.9430e-01 6.9997e-01 7.0569e-01 7.1146e-01 7.1726e-01 7.2310e-01 7.3487e-01 7.4671e-01 7.5859e-01 7.7046e-01 7.8230e-01; ...
    6.9176e-01 6.9688e-01 7.0204e-01 7.0722e-01 7.1243e-01 7.1767e-01 7.2820e-01 7.3879e-01 7.4941e-01 7.6003e-01 7.7062e-01; ...
    6.8962e-01 6.9428e-01 6.9896e-01 7.0367e-01 7.0840e-01 7.1314e-01 7.2267e-01 7.3223e-01 7.4182e-01 7.5141e-01 7.6097e-01; ...
    6.8779e-01 6.9206e-01 6.9635e-01 7.0066e-01 7.0497e-01 7.0930e-01 7.1799e-01 7.2671e-01 7.3544e-01 7.4417e-01 7.5287e-01; ...
    6.8419e-01 6.8772e-01 6.9125e-01 6.9479e-01 6.9833e-01 7.0188e-01 7.0898e-01 7.1609e-01 7.2320e-01 7.3030e-01 7.3738e-01; ...
    6.8157e-01 6.8456e-01 6.8755e-01 6.9054e-01 6.9353e-01 6.9652e-01 7.0250e-01 7.0848e-01 7.1445e-01 7.2041e-01 7.2636e-01; ...
    6.7957e-01 6.8216e-01 6.8474e-01 6.8732e-01 6.8989e-01 6.9247e-01 6.9762e-01 7.0276e-01 7.0790e-01 7.1302e-01 7.1813e-01; ...
    6.7801e-01 6.8027e-01 6.8254e-01 6.8480e-01 6.8706e-01 6.8932e-01 6.9382e-01 6.9832e-01 7.0281e-01 7.0729e-01 7.1175e-01];

% --- controllo di forma: intercetta subito un errore di formattazione -----
nT = numel(Tg);  nP = numel(pg);
tabs = {t_rho, t_cv, t_cp, t_cs, t_mu, t_kf, t_phi, t_h};
nms  = {'t_rho','t_cv','t_cp','t_cs','t_mu','t_kf','t_phi','t_h'};
for kk = 1:numel(tabs)
    if ~isequal(size(tabs{kk}), [nT nP])
        error('heTableData:badShape', ...
            ['%s ha dimensione [%s] invece di [%d %d]. Causa tipica: manca il ' ...
             ''';'' a fine riga nel blocco di dati, quindi le righe si fondono ' ...
             'in un unico vettore.'], nms{kk}, num2str(size(tabs{kk})), nT, nP);
    end
end
end