from os.path import dirname
from os.path import realpath, dirname
from ctypes import CDLL, c_double, c_uint, c_ulong, c_bool, Structure, POINTER

libc = CDLL("libc.so.6")
libsimulat = CDLL(dirname(dirname(__file__)) + "/lib/libsimulat.so")

libc.srand(libc.time(None))


def gen_uniform_random(count):
    c_gen_uniform_random = libsimulat.gen_uniform_random
    c_gen_uniform_random.restype = c_double
    for i in range(count):
        yield c_gen_uniform_random()

def gen_rayleigh_random(sigma, count):
    c_gen_rayleigh_random = libsimulat.gen_rayleigh_random
    c_gen_rayleigh_random.argtypes = [c_double]
    c_gen_rayleigh_random.restype = c_double
    for i in range(count):
        yield c_gen_rayleigh_random(sigma)

def gen_standard_normal_random(count):
    c_gen_standard_normal_random = libsimulat.gen_standard_normal_random
    c_gen_standard_normal_random.restype = c_double
    for i in range(count):
        yield c_gen_standard_normal_random()

def gen_normal_random(mean, sigma, count):
    c_gen_normal_random = libsimulat.gen_normal_random
    c_gen_normal_random.argtypes = [c_double, c_double]
    c_gen_normal_random.restype = c_double
    for i in range(count):
        yield c_gen_normal_random(sigma, mean)

def gen_binomial_random(p, count):
    c_gen_binomial_random = libsimulat.gen_binomial_random
    c_gen_binomial_random.argtypes = [c_double]
    for i in range(count):
        yield c_gen_binomial_random(p)

def simulat_BPSK(N, DB_MAX, fading):
    class Point(Structure):
        _fields_ = [('SNR_db', c_double),
                    ('BER', c_double)]
    c_simulat_BPSK = libsimulat.simulat_BPSK
    c_simulat_BPSK.argtypes = [c_ulong, c_uint, c_bool, POINTER(Point)]

    snr_ber = (Point * DB_MAX)()
    c_simulat_BPSK(N, DB_MAX, fading, snr_ber)
    rev = dict()
    for point in snr_ber:
        rev[point.SNR_db] = point.BER
    return rev
