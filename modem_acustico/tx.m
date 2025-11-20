%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg = "Engenharia"; % mensagem (dados)
RB = 100; % taxa de bits (bps)
Fp = 2000; % frequencia da portadora
Fa = 8000; % frequencia ou taxa de amostragem

OCTAVE = 1;

if OCTAVE == 1
  pkg load signal;
end

[ytx,bits] = transmissor(msg, RB, Fp, Fa);

disp(['Tamanho do quadro: ' num2str(length(bits)) ' bits'])
disp(['Tamanho do sinal: ' num2str(length(ytx)) ' amostras'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmissao do sinal usando a placa de som        %
if OCTAVE == 1
    soundsc(ytx,Fa);
else % MATLAB
    p = audioplayer(ytx, Fa);
    play(p);
end


% grafico do sinal transmitido
t = (0:length(ytx) - 1)/Fa;
stem(t,ytx);

% -------------------------------------------------------------------------
% GRAFICO DO ESPECTRO DE FREQUENCIA (PSD)
% -------------------------------------------------------------------------
L = length(ytx);
NFFT = 2^nextpow2(L); % Próxima potência de 2

% Calcula a FFT
Y = fft(ytx, NFFT);

% Vetor de frequência de 0 até Fa/2
f = Fa/NFFT*(0:NFFT/2);

% Calcula o espectro de magnitude
P1 = abs(Y/L);
P1 = P1(1:NFFT/2+1);
P1(2:end-1) = 2*P1(2:end-1); % Para espectro de um lado

% Plota em escala logarítmica (dB)
figure;
plot(f, 20*log10(P1));
title('Espectro de Frequência do Sinal Transmitido');
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
grid on;