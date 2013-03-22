#encoding=utf8
import re
patter = re.compile('\d+[MIDNSHP]')
patterPysam = re.compile('[MIDNSHP]\d+')


def parse_CIGAR(cigar):
    '''
    parse sam CIGAR    
    '''
    #  注意这里的cigar是sam格式中的cigar eg: 20S100M2I21S
    aDict = {'M': [], 'I':[], 'D': [], 'N': [], 'S': [], 'H': [], 'P': []}
    for match in patter.findall(cigar):
        aDict[match[-1]].append(int(match[0:-1]))
    return aDict

def sam_read_len(cigar):
    #  注意这里的cigar是sam格式中的cigar
    cigarDict = parse_CIGAR(cigar)
    total = sum(cigarDict['M']) + sum(cigarDict['S']) + sum(cigarDict['I']) 
    return total

def read_simulate_ref(read, cigarstring):
    #  注意这里的cigarstring是pysam的cigarstring eg: S20M100I2S21
    replaceLetters = {'N': 'N', 'D': '-', 'I': '', 'S': '', 'H': '', 'P': ''}
    readSimulateRef = ''
    flag = 0
    for match in patterPysam.findall(cigarstring):
        cigarLater = match[0]
        value = int(match[1:])
        if cigarLater == 'M':
            readSimulateRef = ''.join((readSimulateRef, read[flag:(flag+value)]))
            flag += value
        else:
            readSimulateRef = ''.join((readSimulateRef, replaceLetters[cigarLater]*value))
            if cigarLater  in 'IS':  #  这里要注意 即只有当cigarstring是MSI所对应的值和才是read的长度，所以只有当cigarstring是MSI是才能加flag 
                flag += value
    return readSimulateRef

def pysam_cigarstring_to_cigar(cigarstring):
    cigar = ''
    for match in patterPysam.findall(cigarstring):
        cigar = ''.join((cigar,match[1:], match[0]))
    return cigar
# print read_simulate_ref('Iaaatttssssii', '1I3M1I2M3N4S2I2D')


