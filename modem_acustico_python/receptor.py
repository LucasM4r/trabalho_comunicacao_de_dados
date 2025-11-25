import numpy as np
from scipy import signal

# ====================== Funções auxiliares ======================

def int_from_bits(bits):
    """
    Converte lista/vetor de bits LSB-first em inteiro
    (equivalente ao bi2de(..., 'right-msb')).
    """
    val = 0
    for i, b in enumerate(bits):
        val |= (int(b) & 1) << i
    return val


# ====================== Hamming (7,4) ======================

def hamming74_decode(msg_bits, verbose=True):
    """
    Hamming (7,4) decoder (corrige 1 erro por bloco).
    Entrada: bits 0/1, comprimento múltiplo de 7 (resto é descartado).
    Saída: bits de dados decodificados.
    """
    msg_bits = np.array(msg_bits, dtype=int).flatten()
    n = 7
    k = 4
    num_blocks = len(msg_bits) // n
    decoded = []

    for i in range(num_blocks):
        block = msg_bits[i*n:(i+1)*n].copy()
        b1, b2, b3, b4, b5, b6, b7 = block

        # Síndrome
        s1 = (b1 ^ b3 ^ b5 ^ b7) & 1   # paridade 1,3,5,7
        s2 = (b2 ^ b3 ^ b6 ^ b7) & 1   # paridade 2,3,6,7
        s3 = (b4 ^ b5 ^ b6 ^ b7) & 1   # paridade 4,5,6,7

        error_pos = s1 * 1 + s2 * 2 + s3 * 4  # índice 1-based

        if error_pos != 0:
            idx = error_pos - 1  # 0-based
            if idx < n:
                block[idx] ^= 1

        # Dados [d1 d2 d3 d4] = posições [3,5,6,7]
        d1, d2, d3, d4 = block[2], block[4], block[5], block[6]
        decoded.extend([d1, d2, d3, d4])

        if verbose:
            print(f"Bloco {i+1} recebido : {msg_bits[i*n:(i+1)*n]}")
            dados_recebidos = msg_bits[i*n:(i+1)*n][[2,4,5,6]]
            corrected = np.array([d1, d2, d3, d4], dtype=int)
            if np.array_equal(dados_recebidos, corrected):
                print(f"Bloco {i+1} corrigido: não houve correção")
            else:
                print(f"Bloco {i+1} corrigido: {corrected}")

    return np.array(decoded, dtype=int)


# ====================== RECEPTOR ======================

def receptor_2(yrx, TB, Fp, Fa):
    """
    Receptor BPSK correspondente ao transmissor().
    yrx : sinal recebido
    TB  : taxa de bits (bits/s) (igual a RB)
    Fp  : frequência da portadora (Hz)
    Fa  : frequência de amostragem (Hz)
    Retorna:
      x : bits do payload (sem o byte de comprimento)
      l : comprimento da mensagem em bytes
    """
    yrx = np.array(yrx, dtype=float).flatten()
    if yrx.size == 0:
        print("Aviso: sinal de entrada vazio.")
        return np.array([], dtype=int), 0

    len_yrx = len(yrx)
    t = np.arange(len_yrx) / float(Fa)

    # 1. Costas loop (sincronização de fase)
    fle = 64
    h_lp = signal.firwin(fle + 1, 0.001)
    mu = 0.03
    theta = np.zeros(len_yrx)
    zs_buffer = np.zeros(fle + 1)
    zc_buffer = np.zeros(fle + 1)
    h_flip = h_lp[::-1]

    for k in range(len_yrx - 1):
        Ik = 2 * yrx[k] * np.cos(2 * np.pi * Fp * t[k] + theta[k])
        Qk = 2 * yrx[k] * np.sin(2 * np.pi * Fp * t[k] + theta[k])

        zs_buffer = np.concatenate([zs_buffer[1:], [Qk]])
        zc_buffer = np.concatenate([zc_buffer[1:], [Ik]])

        lpfs = np.dot(h_flip, zs_buffer)
        lpfc = np.dot(h_flip, zc_buffer)

        theta[k+1] = theta[k] - mu * lpfs * lpfc

    # 2. Demodulação coerente
    y_demod = yrx * np.cos(2 * np.pi * Fp * t + theta)

    # 3. Sincronização de símbolo (Early-Late)
    yf = y_demod
    st = int(round(Fa / float(TB)))
    d = 1
    err = 0.01
    si = []
    t_samp = st
    max_idx = len_yrx - st

    while t_samp <= max_idx:
        if t_samp > d and (t_samp + d) < len_yrx:
            dif = abs(yf[t_samp - d]) - abs(yf[t_samp + d])

            if dif > err:
                t_samp -= 1
            elif dif < -err:
                t_samp += 1

            si.append(t_samp)
        t_samp += st

    if len(si) == 0:
        print("Aviso: falha na sincronização de símbolo.")
        return np.array([], dtype=int), 0

    y_symbol = yf[si]

    # 4. Decisão BPSK
    x_demod = (y_symbol > 0).astype(int)

    # 5. Sincronização de quadro (busca SFD)
    SFD = np.array([
        1,1,0,0,1,1,1,0,0,0,1,1,1,1,0,0,
        0,0,1,1,1,0,0,0,1,1,0,0,1,0,1,0
    ], dtype=int)
    SFD_bipolar = 2 * SFD - 1
    x_bipolar = 2 * x_demod - 1

    best_corr = -np.inf
    best_idx = 0
    best_sign = 1

    for i in range(0, len(x_demod) - len(SFD) + 1):
        seg = x_bipolar[i:i+len(SFD)]
        corr = np.sum(SFD_bipolar * seg)
        if abs(corr) > best_corr:
            best_corr = abs(corr)
            best_idx = i
            best_sign = 1 if corr >= 0 else -1

    if best_sign < 0:
        x_demod = 1 - x_demod  # inverte todos os bits

    start_payload_idx = best_idx + len(SFD)
    if start_payload_idx >= len(x_demod):
        print("Aviso: SFD encontrado muito no final.")
        return np.array([], dtype=int), 0

    bitstream_com_hamming = x_demod[start_payload_idx:]
    n = 7
    max_blocks = len(bitstream_com_hamming) // n

    if max_blocks < 2:
        print("Aviso: quadro muito curto para extrair comprimento.")
        return np.array([], dtype=int), 0

    # Decodifica os 2 primeiros blocos para obter o byte de comprimento
    length_bits_dec = hamming74_decode(bitstream_com_hamming[:2 * n], verbose=False)
    if len(length_bits_dec) < 8:
        print("Aviso: não foi possível extrair o comprimento.")
        return np.array([], dtype=int), 0

    length_byte = length_bits_dec[:8]
    l = int_from_bits(length_byte)

    expected_bits = 8 * (l + 1)
    blocks_needed = int(np.ceil(expected_bits / 4))

    if blocks_needed > max_blocks:
        print("Aviso: quadro recebido menor que o esperado; truncando dados.")
        blocks_needed = max_blocks
        expected_bits = min(expected_bits, blocks_needed * 4)

    bitstream_util = bitstream_com_hamming[:blocks_needed * n]

    # 6. Payload bruto (sem correção, só debug)
    num_blocks_raw = len(bitstream_util) // n
    raw_payload_bits = []

    for i in range(num_blocks_raw):
        bloco = bitstream_util[i*n:(i+1)*n]
        d1 = bloco[2]
        d2 = bloco[4]
        d3 = bloco[5]
        d4 = bloco[6]
        raw_payload_bits.extend([d1, d2, d3, d4])

    if len(raw_payload_bits) >= 8:
        length_byte_raw = raw_payload_bits[:8]
        l_raw = int_from_bits(length_byte_raw)
        print('--- Mensagem interpretada ANTES da correcao (sem ECC)---')
        print(f'Byte de comprimento sem correcao: {l_raw}')
        print('Mensagem sem correção (bits):',
              ''.join(str(b) for b in raw_payload_bits))

    # 7. Decodificação Hamming(7,4)
    x_decodificado = hamming74_decode(bitstream_util, verbose=True)
    x_decodificado = x_decodificado[:expected_bits]
    print(f'Número de bits após Hamming (úteis): {len(x_decodificado)}')
    print('Primeiros 16 bits:',
          x_decodificado[:min(16, len(x_decodificado))])

    # 8. Remoção do cabeçalho (byte de comprimento) e extração do payload
    if len(x_decodificado) < 8:
        print("Aviso: bitstream muito curto após decodificação.")
        return np.array([], dtype=int), 0

    if expected_bits < 8:
        return np.array([], dtype=int), 0

    # Bits só do payload (sem o byte de comprimento)
    x = x_decodificado[8:expected_bits]
    return x, l
