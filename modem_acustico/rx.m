clear all; clc;

OCTAVE = 1;

if OCTAVE == 1
pkg load communications;
end

RB = 100; % taxa de bits (bps)
Fp = 2000; 
Fa = 8000;

on = 1; % controle se opera ou nao em tempo real

filename = 'sinal/sinal_rx_teste.wav';

if on == 1
% captura do sinal de audio
NS = 0;             % nivel de sinal
NS_min = 0.1;       % nivel de sensibilidade - ajustar conforme necessario
t_captura = 4;      % tempo de captura do sinal
while NS < NS_min
    y = record(t_captura, Fa);
    if ~isempty(y)
      NS = max(y)
    end
end
disp(['Valor de deteccao do sinal: ' num2str(NS)]);

% salva arquivo com sinal para usar no receptor
audiowrite (filename, y, Fa);

else

[y, fs] = audioread(filename);
plot(y);

endif

% retorna a mensagem em bits e o tamanho dela
[m,l] = receptor_2(y, RB, Fp, Fa);

% conversao de bits para texto
x = bi2de(reshape(m,8,l)')';
msg = char(x);
disp(["Tamanho da mensagem recebida: " num2str(l) " bytes"]);
disp(["Mensagem recebida: " msg]);


