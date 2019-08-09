NEXUS = '#NEXUS\n' + \
        'begin data;\n' + \
        'dimensions ntax={num_taxa} nchar={num_char};\n' + \
        'format datatype={type} interleave=no gap=-;\n' + \
        'matrix\n\n'


def nexus_fmt(num_seq, seq_len, data_type='protein'):
    return NEXUS.format(num_taxa=num_seq, num_char=seq_len, type=data_type)
