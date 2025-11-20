function [y, msg_encoded] = transmissor(msg, RB, Fp, Fa)

%%%%%%%%%%%%%%%%%%%%%%  Entrada do transmissor   %%%%%%%%%%%%%%%%
%   msg           -- String (texto) a ser transmitida
%   RB            -- Taxa de bits (Bits/s)
%   Fp            -- Frequencia da portadora (Hz)
%   Fa            -- Frequencia de amostragem da placa de som (Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saida:
%   y             -- Sinal modulado pronto para transmissão
%   msg_encoded   -- Vetor de bits final após enquadramento e codificação
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detectar ambiente (Octave ou MATLAB)
% Assume-se Octave se não for MATLAB (para fins de carregamento de pacotes)
IS_OCTAVE = exist('OCTAVE_VERSION', 'builtin') ~= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pacotes de software
if IS_OCTAVE
    pkg load communications;
    pkg load signal; % Para rcosine
end

%% Definições de Enquadramento
SFD = [1 1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 0 0 1 1 1 0 0 0 1 1 0 0 1 0 1 0]; % 32 bits
% PREAMBLE: 10 unsynchronized bits, followed by 15 symbol intervals of 1010...
PRE_UNSYNC = ones(1, 10);
PRE_SYNC_PATTERN = upsample([1 0 1 0 1 0 1 0 1 0 1 0 1 0 1], 2); % 30 bits (15 símbolos x 2)

% Se o preâmbulo é usado apenas para a detecção de energia e sincronização de clock,
% o original pode ser simplificado (original: [ones(1,10) upsample(ones(1,15),2)]).
% Adotando a sequência original, mas renomeando para clareza:
PREAMBLE = [ones(1, 10) upsample(ones(1, 15), 2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Conversão de Texto para Bits e Codificação
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(msg) && ~isempty(msg)
    l = length(msg);                  % Comprimento da mensagem em bytes
    msg_dec = double(msg);            % Converte string para ASCII (decimal)

    % Monta o quadro de dados (Byte de Comprimento + Dados ASCII)
    % A função de2bi/dec2bin é usada para converter para bits (8 por caractere).
    % Assume-se a ordem de bits necessária para o receptor (geralmente LSB ou MSB no último).
    % Usando 'right-msb' (MSB da direita/primeiro) para compatibilidade:
    
    % Byte de Comprimento (l)
    l_bits = de2bi(l, 8, 'right-msb'); % Byte de comprimento como vetor [1x8]

    % Dados da Mensagem (msg_dec)
    msg_data_bits = reshape(de2bi(msg_dec, 8, 'right-msb')', 1, []);

    % Vetor de bits antes da codificação: [Byte_Comprimento | Dados]
    msg_bits_raw = [l_bits msg_data_bits];

    % Codificação de canal (Hamming)
    msg_encoded = hamming_encode(msg_bits_raw);
else
    % Caso de teste ou input não-string (mantém a compatibilidade com o caso MATLAB original)
    if isnumeric(msg)
       msg_encoded = msg;
    else
       error('Input "msg" deve ser uma string de texto.');
    end
end

disp(['Tamanho do payload (bits) antes do enquadramento: ' num2str(length(msg_encoded))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Enquadramento                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PREAMBLE | SFD | PAYLOAD | PREAMBLE]
msg_framed = [PREAMBLE SFD msg_encoded PREAMBLE];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Codificacao Polar (Banda Base)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Codificação polar: 0 -> -1, 1 -> +1
s = [-1, 1];
y_polar = s(msg_framed + 1);

disp(['Tamanho do quadro: ' num2str(length(y_polar)) ' símbolos']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Filtragem para Formatação de Pulso (RC/RRC) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0.5;                             % Fator de decaimento (roll-off)
sps = floor(Fa / RB);                % Amostras por Símbolo

if IS_OCTAVE
    % Ajuste da taxa de bits para garantir um número inteiro de amostras/símbolo
    RB_f = Fa / sps;
    
    % Projeto do filtro Raised-Cosine (Octave: usa rcosine)
    num = rcosine(RB_f, Fa, 'default', r); % Coeficientes do filtro
    
    % Filtragem e Oversampling
    y_filtered = rcosflt(y_polar, RB_f, Fa, 'filter', num)'; % rcosflt já faz upsampling/filtragem
else % MATLAB
    % Projeto do filtro Raised-Cosine (MATLAB: usa rcosdesign)
    % 6 é o atraso em símbolos
    h = rcosdesign(r, 6, sps);
    
    % Upsampling e Filtragem FIR (upfirdn)
    y_filtered = upfirdn(y_polar, h, sps);
end

% Ajuste de comprimento (opcional, dependendo do rcosflt/upfirdn)
% Se necessário, o sinal pode ser normalizado ou truncado aqui.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Modulacao em Banda Passante (BPSK)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (0:length(y_filtered) - 1)/Fa;
y = y_filtered .* cos(2 * pi * Fp * t); % Modulação BPSK

% Retorna o vetor de bits codificado e enquadrado para inspeção
msg_encoded = msg_framed;

end


function coded_bits = hamming_encode(msg_bits)
% hamming_encode - Codifica um vetor de bits usando Hamming (7,4)k = 4;  % Número de bits de dados
n = 7;  % Tamanho do bloco codificado
len_msg = length(msg_bits);
num_blocks = ceil(len_msg / k);

% Padding (preenche com zeros se necessário para completar o último bloco)
padded_bits = [msg_bits zeros(1, num_blocks * k - len_msg)];

coded_bits = [];

for i = 1:num_blocks
    data_block = padded_bits((i-1)*k + 1 : i*k);
    % A função 'encode' no Octave pode retornar um vetor coluna, por isso a transposição (').
    coded_block = encode(data_block, n, k, 'hamming')';
    coded_bits = [coded_bits coded_block];
end

end