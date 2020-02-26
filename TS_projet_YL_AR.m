%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fichier script crée par Alexandre RIVIERE et Yann LACROIX %%%
%%%                  Traitement du signal                     %%%
%%%               Projet SN DVB-RCS 2019-2020                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3.2.1 Modulation bande base

%Chargement des données
load donnees1.mat; load donnees2.mat;

%Déclaration des constantes
fp1 = 0;        %Porteuse allouée pour l'utilisateur 1 (Hz)
fp2 = 46000;    %Porteuse allouée pour l'utilisateur 2 (Hz)
fe = 128000;    %Fréquence d'échantillonnage (Hz)
Te = 1/fe;      %Intervalle temporel entre les échantillons (s)
T = 0.040;      %Intervalle pour les slots (s)
Ns = 10;        %Nombre d'échantillons pour le signal type NRZ
Ts = Ns*Te;     %Durée allouée pour chaque bit dans notre signal



%Création des messages bruts à partir des données
m1 = kron(2 * bits_utilisateur1 - 1, ones(1,Ns));
m2 = kron(2 * bits_utilisateur2 - 1, ones(1,Ns));


%Axe des temps
t = linspace(0, Ts*Ns, length(m1));

%Tracé des deux signaux bruts
figure;
plot(t, m1, t, m2);
xlabel('Temps (s)');
ylabel('Signal généré amplitude (V)');
title('Figure 1 - Signal NRZ généré pour les deux utilisateurs');
legend('m1', 'm2');


%FFT des deux signaux bruts
M1 = fft(m1);
M2 = fft(m2);

%Densités spectrales associées
DSP_M1 = abs(M1.^2)/length(m1);
DSP_M2 = abs(M2.^2)/length(m2);

%Axe en fréquence
f = linspace(-fe/2,fe/2, length(m1));

%Tracé des deux densités spectrales
figure;
semilogy(f, fftshift(DSP_M1), f, fftshift(DSP_M2));
xlabel('Fréquence (Hz)');
ylabel('Densité spectrale (V.V.s)');
title('Figure 2 - Densités spectrales pour les deux utilisateurs');
legend('m1', 'm2');

%% 3.2.2 Construction du signal MF-TDMA

%Creation des 5 slots pour l'utilisateur 1
signal_slots_1 = zeros(length(m1),5);
%Attribution du slot 2 pour le message de l'utilisateur 1
signal_slots_1(:,2) = m1;
signal_slots_vectorise_1 = signal_slots_1(:);

%Creation des slots pour l'utilisateur 2
signal_slots_2 = zeros(length(m2),5);
%Attribution du slot 5 pour le message de l'utilisateur 2
signal_slots_2(:,5) = m2;
signal_slots_vectorise_2 = signal_slots_2(:);

%Axe des temps
t_slots = linspace(0, 5*T, length(signal_slots_vectorise_1));
%Construction des signaux sur la fréquence de leur porteuse
x1 = signal_slots_vectorise_1; %Pour l'utilisateur 1 la fréquence de porteuse allouée est 0!
x2 = signal_slots_vectorise_2 .* cos(2*pi*fp2.*t_slots');

%Sommation des deux signaux x1 et x2, on obtient le signal x qui va être
%transmis
x = x1 + x2;

%Tracé de ce signal
figure;
subplot(1,2,1);
plot(t_slots, x);
xlabel('Temps (s)');
ylabel('Signal amplitude (V)');
title('Signal MF-TDMA non bruité');

%Paramètres du bruit gaussien
RSB = 75; %Rapport Signal/Bruit en dB (on met 75dB car à 100dB le bruit ne se voyait quasiment pas sur les tracés)
puissance_s = sum(abs(fft(x).^2)/length(x)); %Puissance du signal x
puissance_b = puissance_s*10^(-RSB/10); %Formule établie dans la partie théorique

%Génération du bruit gaussien
BG = sqrt(puissance_b).*randn(length(x),1);

%Ajout du bruit gaussien au signal, on obtient le signal bruité
x_bruite = x + BG;

%Tracé du signal bruité
subplot(1,2,2);
plot(t_slots, x_bruite);
xlabel('Temps (s)');
ylabel('Signal amplitude (V)');
title('Signal MF-TDMA bruité');

%DSP du signal bruité
tdf_xb = fft(x_bruite);
DSP_b = abs(tdf_xb.^2)/length(tdf_xb);

%Axe des fréquences
f = linspace(-fe/2,fe/2, length(DSP_b));
%Tracé de cette DSP
figure;
semilogy(f,fftshift(DSP_b));
xlabel('Fréquence (Hz)');
ylabel('Amplitude de la DSP (V.V.s)');
title('DSP du signal MF-TDMA bruité');

%% 4.1.1 Synthèse du filtre passe-bas

N0 = 301; %Ordre du filtre
fc1 = 20000; % On prend 20kHz comme fréquence de coupure

%Génération de la réponse impulsionnelle
k = -(N0-1)/2 : 1 : (N0-1)/2;
rep_imp_pb = 2*fc1/fe*sinc(2*fc1/fe*k);

%Tracé de la réponse impulsionnelle
figure;
subplot(1,2,1);
plot(k, rep_imp_pb);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Filtre passe-bas - Réponse impulsionnelle');


%Transformée de fourier de la réponse impulsionnelle
TF_pb = fft(rep_imp_pb);
k_f = linspace(-fe/2, fe/2, N0);

%Tracé de la réponse impulsionnelle
subplot(1,2,2);
semilogy(k_f, fftshift(abs(TF_pb)));
xlabel('Fréquence (Hz)');
ylabel('Amplitude de la TF');
title('Filtre passe-bas - TF de la réponse impulsionnelle');


%Tracé de la DSP et de la TF de notre filtre
figure
semilogy(f,fftshift(DSP_b), k_f, fftshift(abs(TF_pb)));
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('Superposition de la DSP du signal reçu et de la TF du filtre passe-bas');
legend('DSP signal reçu', 'TF du filtre passe-bas');

%% 4.1.2 Synthèse du filtre passe-haut

%On commence à générer un passe-bas
N2 = 301; %Ordre du filtre
fc2 = 3*fp2/4; %Fréquence de coupure

%Génération de la réponse impulsionnelle du passe-bas
k2 = -(N2-1)/2 : 1 : (N2-1)/2;
rep_imp_pb2 = 2*fc2/fe*sinc(2*fc2/fe*k2);

%On en déduit la réponse impulsionnelle du passe-haut
rep_imp_ph = - rep_imp_pb2;
rep_imp_ph((N2+1)/2) = rep_imp_ph((N2+1)/2) + 1;

%Transformée de Fourier de notre passe-haut
TF_ph = fft(rep_imp_ph);
k_f2 = linspace(-fe/2, fe/2, N2);

%Tracé de la réponse impulsionnelle de notre passe-haut...
figure
subplot(1,2,1);
plot(k2, rep_imp_ph);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Filtre passe-haut - Réponse impulsionnelle');
%...Et de sa Transformée de Fourier
subplot(1,2,2);
semilogy(k_f2, fftshift(abs(TF_ph)));
xlabel('Fréquence (Hz)');
ylabel('Amplitude de la TF');
title('Filtre passe-haut - TF de la réponse impulsionnelle');

%Tracé de la DSP du filtre passe-haut
figure
semilogy(f,fftshift(DSP_b), k_f2, fftshift(abs(TF_ph)));
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('Superposition de la DSP du signal reçu et de la TF du filtre passe-haut');
legend('DSP signal reçu', 'TF du filtre passe-haut');

%% 4.1.3 Filtrage

signal_bf = conv(x_bruite, rep_imp_pb, 'same'); %Filtrage du signal porté à basse fréquence
signal_hf = conv(x_bruite, rep_imp_ph, 'same'); %Filtrage du signal porté à haute fréquence

%Tracé des deux signaux obtenu après filtrage par respectivement
%un passe-bas et un passe-haut
figure
subplot(1,2,1);
plot(t_slots, signal_bf);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal retrouvé par filtrage passe-bas');

subplot(1,2,2);
plot(t_slots, signal_hf);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal retrouvé par filtrage passe-haut');

%% 4.2 Retour en bande de base

%On multiplie par le cosinus qui a servi à moduler le signal pour retourner
%en bande de base, pour le signal 1 c'ets inutile puisque nous n'avions pas
%modulé le signal (fp1 = 0Hz !)
signal_hf2 = signal_hf.*cos(2*pi*fp2.*t_slots');

%On fait un nouveau filtre passe-bas pour le signal haute fréquence
N3 = 301; %Ordre du filtre
fc3 = 25000; %On prend 25kHz comme fréquence de coupure

%Génération de la réponse impulsionnelle de notre filtre passe-bas
k3 = -(N3-1)/2 : 1 : (N3-1)/2;
rep_imp_pb3 = 2*fc3/fe*sinc(2*fc3/fe*k3);

%On filtre avec un passe bas le signal retrouvé après retour en ban
signal_bf_bande_base = signal_bf;   %Ici pas besoin de filtrer puisqu'on n'a
                                    %pas eu besoin de faire le retour en
                                    %bande de base (pas de modulation)
signal_hf_bande_base = conv(signal_hf2, rep_imp_pb3, 'same');

%Tracé des signaux récupéré après retour en bande de base et filtrage par
%un passe-bas
figure
subplot(1,2,1);
plot(t_slots, signal_bf_bande_base);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal BF Bande de base');

subplot(1,2,2);
plot(t_slots, signal_hf_bande_base);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal HF Bande de base');

%% 4.3 Détection slot utile

L=length(signal_bf_bande_base)/5; %Nombre d'échantillons dans un slot
Mod_signal_bf_bb = abs(signal_bf_bande_base.^2)/length(signal_bf_bande_base);

Energie_max = 0;
slot_ut1 = 1;

%On regarde un à un les slots...
for i = 1:5
    energie = sum(Mod_signal_bf_bb((i-1)*L+1:i*L))/L;
    
    if (energie > Energie_max)
        slot_ut1 = i;
        Energie_max = energie;
    end
end
%...Et on garde celui avec l'énergie maximale (ie le signal du slot utile)
X_bf=signal_bf_bande_base((slot_ut1-1)*L+1:slot_ut1*L);



L=length(signal_hf_bande_base)/5;
Mod_signal_hf_bb = abs(signal_hf_bande_base.^2)/length(signal_hf_bande_base);
Energie_max = 0;
slot_ut2 = 1;

%Même opération ici
for i = 1:5
    energie = sum(Mod_signal_hf_bb((i-1)*L+1:i*L))/L;
    
    if (energie > Energie_max)
        slot_ut2 = i;
        Energie_max = energie;
    end
end
%On récupére le signal du slot utile
X_hf=signal_hf_bande_base((slot_ut2-1)*L+1:slot_ut2*L);

%Axe des temps
t_slots_2 = linspace(0, T, length(X_hf));
%Tracé des deux signaux isolés après détection du slot utile
figure
subplot(1,2,1);
plot(t_slots_2, X_bf);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal slot utile BF Bande de base');

subplot(1,2,2);
plot(t_slots_2, X_hf);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal slot utile HF Bande de base');

%% 4.4 Démodulation bande de base

%Démodulation bande de base pour l'utilisateur 1
SignalFiltre1=conv(X_bf,ones(1,Ns),'same');
SignalEchantillonne1=SignalFiltre1(1 :Ns :end);
BitsRecup1=(sign(SignalEchantillonne1)+1)/2;

%Démodulation bande de base pour l'utilisateur 2
SignalFiltre2=conv(X_hf,ones(1,Ns),'same');
SignalEchantillonne2=SignalFiltre2(1 :Ns :end);
BitsRecup2=(sign(SignalEchantillonne2)+1)/2;

%Conversion binaire->string
texte1 = bin2str(BitsRecup1);
texte2 = bin2str(BitsRecup2);

%Affiche du texte des messages dans la console
disp(texte1);
disp(texte2);