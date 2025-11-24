function [x,l] = receptor_2(yrx, TB, Fp, Fa)
%%%%%%%%%%%%%%%%%%%%%%  Entrada do receptor   %%%%%%%%%%%%%%%%
%   yrx -- sinal de audio capturado
%   TB  -- Taxa de bits (Bits/s)
%   Fp  -- Frequencia da portadora (Hz)
%   Fa  -- Frequencia de amostragem (Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saida:
%   x   -- Vetor de bits da mensagem (payload)
%   l   -- Comprimento da mensagem em bytes (excluindo o byte de comprimento)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pacotes de software
pkg load signal;
pkg load communications;

%% Inicialização
l = 0;
x = [];

if isempty(yrx)
    warning('receptor_2:emptyInput', 'Sinal de entrada vazio.');
    return;
end

len_yrx = length(yrx);
t = (0:len_yrx - 1)'/Fa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Sincronizacao com Costas Loop   	    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fle = 64;
h = fir1(fle, 0.001);
mu = 0.03;
theta = zeros(len_yrx, 1);
theta(1) = 0;

zs_buffer = zeros(fle+1, 1);
zc_buffer = zeros(fle+1, 1);
h_flip = flipud(h(:));

for k = 1:len_yrx - 1
  I_k = 2 * yrx(k) * cos(2 * pi * Fp * t(k) + theta(k));
  Q_k = 2 * yrx(k) * sin(2 * pi * Fp * t(k) + theta(k));

  zs_buffer = [zs_buffer(2:end); Q_k];
  zc_buffer = [zc_buffer(2:end); I_k];

  lpfs = h_flip' * zs_buffer;
  lpfc = h_flip' * zc_buffer;

  theta(k+1) = theta(k) - mu * lpfs * lpfc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Demodulacao Coerente             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_demod = yrx .* cos(2 * pi * Fp * t + theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Sincronizacao de Simbolo (Early-Late Gate)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yf = y_demod;
st = round(Fa / TB);
d = 1;
err = 0.01;
si = [];
t_samp = st;
max_idx = len_yrx - st + 1;

while t_samp <= max_idx
    if t_samp > d && (t_samp + d) <= len_yrx
        dif = abs(yf(t_samp - d)) - abs(yf(t_samp + d));

        if dif > err
            t_samp = t_samp - 1;
        elseif dif < -err
            t_samp = t_samp + 1;
        end

        si = [si; t_samp];
    end
    t_samp = t_samp + st;
end

if isempty(si)
    warning('receptor_2:symbolSync', 'Falha na sincronização de símbolo.');
    return;
end

y_simbolo = yf(si);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Decodificacao de Bits (BPSK)               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_demod = (y_simbolo > 0)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Sincronizacao do Quadro (SFD)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SFD = [1 1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 0 0 1 1 1 0 0 0 1 1 0 0 1 0 1 0];
SFD_bipolar = SFD * 2 - 1;
x_bipolar = x_demod * 2 - 1;

xc = xcorr(SFD_bipolar, x_bipolar);
[a_max, b_idx] = max(abs(xc));

if xc(b_idx) < 0
    x_demod = ~x_demod;
end

start_idx_x = (length(x_demod) - length(SFD)) - (b_idx - length(SFD)) + 1;
start_payload_idx = start_idx_x + length(SFD);

if start_payload_idx > length(x_demod)
    warning('receptor_2:frameSync', 'SFD encontrado muito no final.');
    return;
end

bitstream_com_hamming = x_demod(start_payload_idx:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Mensagem SEM correção (removendo ECC)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 4;      % dados
n = 7;      % bloco Hamming

num_blocks_raw = floor(length(bitstream_com_hamming) / n);
raw_payload_bits = [];

for i = 1:num_blocks_raw
    bloco = bitstream_com_hamming((i-1)*n + 1 : i*n);

    % Extrair somente os dados sem corrigir
    % Ordem do Hamming(7,4): [p1 p2 d1 p3 d2 d3 d4]
    d1 = bloco(3);
    d2 = bloco(5);
    d3 = bloco(6);
    d4 = bloco(7);

    raw_payload_bits = [raw_payload_bits d1 d2 d3 d4];
end

disp('--- Mensagem interpretada ANTES da correcao (sem ECC)---');

length_byte_raw = raw_payload_bits(1:8);
l_raw = double(bi2de(length_byte_raw, 'right-msb'));
fprintf('Byte de comprimento sem correcao: %d\n', l_raw);

fprintf('Mensagem sem correção (bits): ');
disp(num2str(raw_payload_bits));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Remoção do PREAMBLE_FINAL                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_bits_data = (l_raw + 1) * 8;

num_blocks_valid = ceil(num_bits_data / k);
num_bits_hamming_valid = num_blocks_valid * n;

num_bits_hamming_valid = min(num_bits_hamming_valid, length(bitstream_com_hamming));
bitstream_valido = bitstream_com_hamming(1:num_bits_hamming_valid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Decodificacao Hamming (7,4)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_decodificado = hamming_decode(bitstream_valido);

disp(['Número de bits após Hamming: ' num2str(length(x_decodificado))]);
disp(['Primeiros 16 bits: ' num2str(x_decodificado(1:min(16,end)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Remoção de Cabeçalho e Payload               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(x_decodificado) < 8
    warning('receptor_2:short', 'Bitstream muito curto.');
    return;
end

length_byte = x_decodificado(1:8);
l = double(bi2de(length_byte, 'right-msb'));

l_max_bytes = floor(length(x_decodificado)/8) - 1;

if l < 0 || l > l_max_bytes
    warning('receptor_2:badLength', 'Comprimento inválido.');
    return;
end

expected_bits = 8 * (l + 1);
expected_bits = min(expected_bits, length(x_decodificado));

if expected_bits < 8
    x = [];
    l = 0;
    return;
end

x = x_decodificado(9:expected_bits);

endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Função de Decodificação Hamming (7,4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decoded_bits = hamming_decode(msg_bits)
    k = 4;
    n = 7;
    len_msg = length(msg_bits);
    fprintf('Comprimento da Mensagem %d\n: ', len_msg)
    num_blocks = floor(len_msg / n);
    decoded_bits = [];
    fprintf('Numero de blocos %d\n: ', num_blocks);

    for i = 1:num_blocks
        block = msg_bits((i-1)*n + 1 : i*n);

        % Extrair os 4 bits de dados do bloco original
        dados_recebidos = block([3 5 6 7]);   % [d1 d2 d3 d4]

        % Decodificar
        data_block = decode(block(:), n, k, 'hamming')';

        % Mostrar bloco recebido
        fprintf('Bloco %d recebido : %s\n', i, num2str(block));

        % Testar se houve correção
        if isequal(dados_recebidos, data_block)
            fprintf('Bloco %d corrigido: não houve correção\n', i);
        else
            fprintf('Bloco %d corrigido: %s\n', i, num2str(data_block));
        end

        decoded_bits = [decoded_bits data_block];
    end
    
end

