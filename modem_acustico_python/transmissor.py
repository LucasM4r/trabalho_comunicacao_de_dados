import numpy as np

# ====================== Funções auxiliares ======================

def bits_from_int(value, width=8):
    """
    Converte inteiro em lista de bits LSB-first
    (equivalente ao de2bi(..., 'right-msb') do MATLAB).
    """
    return [(value >> i) & 1 for i in range(width)]


# ====================== Hamming (7,4) ======================

def hamming74_encode(msg_bits):
    """
    Hamming (7,4) encoder.
    Entrada: vetor 1D de bits (0/1).
    Saída: numpy array de bits codificados.
    Ordem em cada bloco: [p1 p2 d1 p3 d2 d3 d4]
    """
    msg_bits = np.array(msg_bits, dtype=int).flatten()
    k = 4
    n = 7
    len_msg = len(msg_bits)
    num_blocks = int(np.ceil(len_msg / k))

    # Padding com zeros se necessário
    padded = np.concatenate([msg_bits,
                             np.zeros(num_blocks * k - len_msg, dtype=int)])

    coded_blocks = []
    for i in range(num_blocks):
        d1, d2, d3, d4 = padded[i*k:(i+1)*k]
        # Paridades para [p1 p2 d1 p3 d2 d3 d4]
        p1 = (d1 ^ d2 ^ d4) & 1      # posições 1,3,5,7
        p2 = (d1 ^ d3 ^ d4) & 1      # posições 2,3,6,7
        p3 = (d2 ^ d3 ^ d4) & 1      # posições 4,5,6,7
        block = np.array([p1, p2, d1, p3, d2, d3, d4], dtype=int)
        coded_blocks.append(block)

    return np.concatenate(coded_blocks)


# ====================== Filtro RRC (equivalente ao rcosdesign) ======================

def rrcosfilter(span, sps, rolloff):
    """
    Gera um filtro raiz cosseno levantado (RRC) similar ao rcosdesign do MATLAB.
    span: atraso em símbolos (lado inteiro, total = span símbolos)
    sps: amostras por símbolo
    rolloff: fator de roll-off (0 < rolloff <= 1)
    Retorna: (t, h)
    """
    N = span * sps
    t = np.arange(-N/2, N/2 + 1) / float(sps)  # em tempos de símbolo
    h = np.zeros_like(t, dtype=float)

    for i, ti in enumerate(t):
        if abs(ti) < 1e-8:
            h[i] = 1.0 - rolloff + (4 * rolloff / np.pi)
        elif abs(abs(4 * rolloff * ti) - 1.0) < 1e-8:
            # ti = ± Ts / (4*rolloff)
            h[i] = (rolloff / np.sqrt(2)) * (
                ((1 + 2/np.pi) * np.sin(np.pi / (4 * rolloff))) +
                ((1 - 2/np.pi) * np.cos(np.pi / (4 * rolloff)))
            )
        else:
            num = (np.sin(np.pi * ti * (1 - rolloff)) +
                   4 * rolloff * ti * np.cos(np.pi * ti * (1 + rolloff)))
            den = (np.pi * ti * (1 - (4 * rolloff * ti)**2))
            h[i] = num / den

    # Normaliza energia
    h = h / np.sqrt(np.sum(h**2))
    return t, h


# ====================== TRANSMISSOR ======================

def transmissor(msg, RB, Fp, Fa):
    """
    Transmissor BPSK com:
      - Byte de comprimento
      - Hamming(7,4)
      - PREAMBLE + SFD + payload + PREAMBLE
      - Pulso RRC
    msg : string de texto
    RB  : taxa de bits (bits/s)
    Fp  : frequência da portadora (Hz)
    Fa  : frequência de amostragem (Hz)
    Retorna:
      y          : sinal modulado (float)
      msg_framed : bits após enquadramento
    """

    # SFD (32 bits)
    SFD = np.array([
        1,1,0,0,1,1,1,0,0,0,1,1,1,1,0,0,
        0,0,1,1,1,0,0,0,1,1,0,0,1,0,1,0
    ], dtype=int)

    # PREAMBLE original: [ones(1,10) upsample(ones(1,15),2)]
    pre_sync = np.zeros(15*2, dtype=int)
    pre_sync[::2] = 1          # [1 0 1 0 ...]
    PREAMBLE = np.concatenate([np.ones(10, dtype=int), pre_sync])

    # 1. Texto -> bits
    if isinstance(msg, str) and len(msg) > 0:
        l = len(msg)
        msg_dec = [ord(c) for c in msg]

        # Byte de comprimento
        l_bits = np.array(bits_from_int(l, 8), dtype=int)

        # Bits dos dados (8 bits por caractere)
        msg_data_bits = []
        for val in msg_dec:
            msg_data_bits.extend(bits_from_int(val, 8))
        msg_data_bits = np.array(msg_data_bits, dtype=int)

        # [Byte_Comprimento | Dados]
        msg_bits_raw = np.concatenate([l_bits, msg_data_bits])

        # Hamming(7,4)
        msg_encoded = hamming74_encode(msg_bits_raw)
    elif isinstance(msg, (list, np.ndarray)):
        msg_encoded = np.array(msg, dtype=int).flatten()
    else:
        raise ValueError('Input "engenharia" deve ser uma string ou vetor de bits.')

    print(f"Tamanho do payload (bits) antes do enquadramento: {len(msg_encoded)}")

    # 2. Enquadramento
    msg_framed = np.concatenate([PREAMBLE, SFD, msg_encoded, PREAMBLE])

    # 3. Polar: 0 -> -1, 1 -> +1
    y_polar = 2 * msg_framed - 1

    print(f"Tamanho do quadro: {len(y_polar)} símbolos")

    # 4. Pulso RRC
    r = 0.5
    sps = int(np.floor(Fa / RB))  # samples per symbol
    if sps < 1:
        raise ValueError("Fa/RB precisa ser >= 1 para ao menos 1 amostra por símbolo")

    span = 6
    _, h = rrcosfilter(span, sps, r)

    # Upsampling
    up = np.zeros(len(y_polar) * sps)
    up[::sps] = y_polar

    # Filtragem
    y_filtered = np.convolve(up, h, mode='full')

    # 5. Modulação BPSK
    t = np.arange(len(y_filtered)) / float(Fa)
    y = y_filtered * np.cos(2 * np.pi * Fp * t)

    return y, msg_framed
