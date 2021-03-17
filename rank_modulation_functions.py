import numpy as np
import itertools
import os
import bchlib

############ Reading/Writing to Flash cells ####################
def read_cells(Cs, Qrs):
    '''
    Cs: numpy array of cell charge levels
    Qrs: list of charge sensing levels
    
    returns mp: numpy array of estimated charge levels
    '''
    mp = np.array([1+np.argwhere(Qrs < c).squeeze().max() if c > Qrs[0] else 0 
               for c in Cs], dtype=int)
    return mp

def program_cells(m, Ls, sig):
    '''
    m: list representation of symbols to program
    Ls: list of nominal charge levels for each symbol (ascending order)
    sig: standard deviation of Gaussian noise about nominal charge level
    
    returns Cs: numpy array of cell charge levels
    '''
    n = len(m)
    Del = (Ls[1]-Ls[0])/2 - 1
    write_errs = np.clip(sig*np.random.randn(n), a_min=-Del, a_max=Del)
    Cs = np.array([np.round(Ls[mi] + (Ls[mi]>0)*write_errs[i])
                   for i, mi in enumerate(m)], dtype=int)

    return Cs

def leak_cells(Cs, leak, t):
    '''
    Cs: numpy array of cell charge levels
    leak: scalar value for leak
    t: scalar value for retention interval
    
    returns Cs_corr: numpy array of cell charge levels after leak
    '''
    n = len(Cs)
    leakage = np.random.poisson(leak*t, n)
    Cs_corr = np.clip(Cs-leakage, a_min=0, a_max=None)
    return Cs_corr


################# Programming via rank modulation #########################
def rank_modulate(pi, Ls, sig):
    '''
    pi: numpy array of permutation to encode
    Ls: list of nominal charge levels for each symbol in permutation (ascending order)
    sig: standard deviation of Gaussian noise about nominal charge level
    
    returns Cs: numpy array of cell charge levels
    '''
    n = len(pi)
    Del = (Ls[1]-Ls[0])/2 - 1
    write_errs = np.clip(sig*np.random.randn(n), a_min=-Del, a_max=Del)
    Cs = np.zeros((n,), dtype=int)
    for i in range(1,n):
        Cs[pi[i]] = np.round(Ls[i] + (Ls[i]>0)*write_errs[i])

    return Cs

################### Construction A #####################
def pi2lehmer(pi):
    '''
    pi: numpy array of permutation
    
    returns inv: inversion vector for pi from lehmer code
    '''
    n = len(pi)
    inv = np.copy(pi-1)
    for i in range(n-1):
        for j in range(i+1, n):
            if inv[j] > inv[i]:
                inv[j] -= 1
    inv = inv[:n-1]
    
    return inv

def lehmer2pi(inv):
    '''
    inv: numpy array of inversion vector (from lehmer code)
    
    returns pip: numpy array of permutation from inv
    '''
    n = len(inv)
    pip = np.zeros(n+1, dtype=int)
    pip[:n] = inv
    for i in range(n-1,-1,-1):
        for j in range(i+1, n+1):
            if pip[j] >= pip[i]:
                pip[j] += 1
    pip += 1
    return pip

##################### Construction B #########################
def gray2dec(n): 
    inv = 0 
    while(n): 
        inv = inv ^ n
        n = n >> 1
    return inv

def dec2gray(n):
    return n ^ (n >> 1)

def code2inv(binary_code):
    '''
    binary_code: numpy array of binary codeword
    
    returns inv_vec: inversion vector for given codeword
    '''
    remaining = len(binary_code)
    pos = 0
    i = 2
    inv_vec = []
    while remaining > 0:
        mi = np.floor(np.log2(i)).astype(int)
        inv_vec.append(gray2dec(int(binary_code[pos:pos+mi], 2)))
        pos = pos + mi
        i += 1
        remaining = remaining - mi
    
    return inv_vec

def inv2pi(inv_vec):
    '''
    inv_vec: numpy array of inversion vector
    
    returns pi: permutation form of given inversion vector
    '''
    n = len(inv_vec)
    pi = np.zeros(n+1, dtype=int)
    pos = list(np.arange(len(pi)))
    for i in range(1, n+1):
        ki = pos[inv_vec[n-i]]
        pi[ki] = i
        pos.remove(ki)
    pi[pos] = n+1
    
    return pi

def pi2inv(pi):
    '''
    pi: numpy array of permutation vector
    
    returns inv_vec: inversion vector form of given permutation
    '''
    n = len(pi)-1
    inv = np.ones(n)
    for i in range(n):
        ki = np.argwhere(pi == i+1)[0][0]
        inv[i] = np.sum(pi[:ki]>i+1)
    inv = np.flip(inv)
    y = np.copy(inv)
    for i in range(1, n):
        mi = np.floor(np.log2(i+1))
        if inv[i-1] > 2**mi - 1:
            y[i] = 2**mi - 1
    
    return y

def inv2code(inv_vec):
    '''
    inv_vec: numpy array of inversion vector
    
    returns binstring: 0/1 string representation of the given inversion vector
    '''
    binstring = ''
    for i, val in enumerate(inv_vec):
        mi = int(np.floor(np.log2(i+2)))
        bstr = bin(dec2gray(int(val))).replace("0b", "")
        binstring += (mi-len(bstr))*'0' + bstr
    return binstring

#################### Testing functions ###########################
def test_RBER(n, Ls, sig, Qrs, leak, T):
    '''
    n: integer number of cells
    Ls: nominal charge levels per symbol (ascending order)
    sig: standard deviation of Gaussian noise about nominal charge level
    Qrs: list of charge sensing levels
    leak: scalar charge leakage rate
    T: retention interval to test
    
    returns RBER: raw bit error rate
    '''
    
    m = np.random.choice([0,1,2,3], size=n)
    Cs = program_cells(m, Ls, sig=sig)
    mp = read_cells(leak_cells(Cs, leak, T), Qrs)
    RBER = np.sum(m != mp)/n
    return RBER
    
def test_constr_A(N, bch, Ls, sig, leak, T):
    '''
    n: integer number of cells
    bch: BCH encoder/decoder object
    Ls: nominal charge levels per symbol (ascending order)
    sig: standard deviation of Gaussian noise about nominal charge level
    leak: scalar charge leakage rate
    T: retention interval to test
    
    returns avg_BER: average bit error rate
    '''
    n = bch.n
    h = bch.ecc_bits
    k = n - h
    
    num_chunks = int(np.ceil(N // k))
    errs = 0.0
    for chunk in range(num_chunks):
        m = os.urandom(k//8)
        data = bytearray(m)
        ecc = bch.encode(data)
        codeword = ''.join(format(byte, '08b') for byte in (data+ecc))[:n]
        inv_vec = [i+1 if int(ci)>0 else 0 for i, ci in enumerate(codeword)]
        pi = lehmer2pi(np.flip(inv_vec))
        
        Cs = rank_modulate(pi-1, Ls, sig)
        Csp = leak_cells(Cs, leak, T)
        pi_corrupted = 1 + np.argsort(Csp)
        
        inv_vec_c = np.flip(pi2lehmer(pi_corrupted))
        code_c = ''.join(['1' if xi > np.floor((i+1)/2) else '0' for i, xi in enumerate(inv_vec_c)])
        binstring = code_c[:n]
        data_p = binstring[:-h]
        data_p = bytearray(int(data_p, 2).to_bytes((len(data_p) + 7) // 8, byteorder='big'))
        ecc_p = binstring[-h:] + int(np.ceil(1+h/8)*8 - bch.ecc_bits)*'0'
        ecc_p = bytearray(int(ecc_p, 2).to_bytes((len(ecc_p) + 7) // 8, byteorder='big'))#[:int((h + 7) // 8)]
        bitflips, data_pp, ecc_pp = bch.decode(data_p, ecc_p)
        
        codeword_new = ''.join(format(byte, '08b') for byte in (data_pp+ecc_pp))[:n]
        errs += sum(c1 != c2 for c1, c2 in zip(codeword, codeword_new))
    
    avg_BER = errs/(n*num_chunks)
    return avg_BER

def test_constr_B(N, bch, Ls, sig, leak, T):
    '''
    n: integer number of cells
    bch: BCH encoder/decoder object
    Ls: nominal charge levels per symbol (ascending order)
    sig: standard deviation of Gaussian noise about nominal charge level
    leak: scalar charge leakage rate
    T: retention interval to test
    
    returns avg_BER: average bit error rate
    '''
    n = bch.n
    h = bch.ecc_bits
    k = n - h
    
    num_chunks = int(np.ceil(N // k))
    errs = 0.0
    for chunk in range(num_chunks):
        m = os.urandom(k//8)
        data = bytearray(m)
        ecc = bch.encode(data)
        codeword = ''.join(format(byte, '08b') for byte in (data+ecc))[:n]
        
        inv_vec = code2inv(codeword)
        pi = inv2pi(inv_vec)
        Ls = [40*i for i in range(len(pi))]
        
        Cs = rank_modulate(pi-1, Ls, sig)
        Csp = leak_cells(Cs, leak, T)
        pi_corrupted = 1 + np.argsort(Csp)
        
        inv_vec_c = pi2inv(pi_corrupted)
        code_c = inv2code(inv_vec_c)
        binstring = code_c[:n]
        data_p = binstring[:-h]
        data_p = bytearray(int(data_p, 2).to_bytes((len(data_p) + 7) // 8, byteorder='big'))
        ecc_p = binstring[-h:] + int(np.ceil(1+h/8)*8 - bch.ecc_bits)*'0'
        ecc_p = bytearray(int(ecc_p, 2).to_bytes((len(ecc_p) + 7) // 8, byteorder='big'))#[:int((h + 7) // 8)]
        bitflips, data_pp, ecc_pp = bch.decode(data_p, ecc_p)
        
        codeword_new = ''.join(format(byte, '08b') for byte in (data_pp+ecc_pp))[:n]
        errs += sum(c1 != c2 for c1, c2 in zip(codeword, codeword_new))
    avg_BER = errs/(n*num_chunks)
    return avg_BER