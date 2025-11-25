from receptor import int_from_bits, receptor_2
from transmissor import transmissor


def bits_para_string(bits, comprimento):
    if len(bits) >= comprimento * 8 and comprimento > 0:
        chars = []
        for i in range(comprimento):
            bits_byte = bits[i*8:(i+1)*8]
            val = int_from_bits(bits_byte)
            chars.append(chr(val))
        return ''.join(chars)
    return ""


def executar_demo():
    msg = "Engenharia"
    RB = 100      # taxa de bits
    Fp = 2000     # portadora
    Fa = 8000     # amostragem

    print("Mensagem transmitida:", msg)
    ytx, _ = transmissor(msg, RB, Fp, Fa)
    m, l = receptor_2(ytx, RB, Fp, Fa)

    msg_rx = bits_para_string(m, l)
    print("Mensagem recebida:", msg_rx)


if __name__ == "__main__":
    executar_demo()
