NEXUS = '\n'.join(['#NEXUS', 'begin data;',
                   'dimensions ntax={num_taxa} nchar={num_char};',
                   'format datatype={type} interleave=no gap=-;',
                   'matrix', '', ''])


def nexus_fmt(num_seq, seq_len, data_type='protein'):
    return NEXUS.format(num_taxa=num_seq, num_char=seq_len, type=data_type)
