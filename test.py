"""
测试均匀, 高斯, 和瑞利分布
"""

from ctypes import CDLL, c_double, c_uint, c_ulong, c_bool, Structure
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

class TestRandom(object):
    """
    用fit()方法对随机变量取样序列进行拟合
    比较随机变量X的概率密度函数和对数组x进行直方图统计的结果
    http://hyry.dip.jp/tech/book/page/scipy/scipy_stats.html
    """
    def __init__(self, N):
        libc = CDLL("libc.so.6")
        libsimulat = CDLL("./libsimulat.so")

        libc.srand(libc.time(None))

        self.gen_uniform_random = libsimulat.gen_uniform_random
        self.gen_uniform_random.restype = c_double

        self.gen_rayleigh_random = libsimulat.gen_rayleigh_random
        self.gen_rayleigh_random.argtypes = [c_double]
        self.gen_rayleigh_random.restype = c_double

        self.gen_standard_normal_random = libsimulat.gen_standard_normal_random
        self.gen_standard_normal_random.restype = c_double

        self.gen_normal_random = libsimulat.gen_normal_random
        self.gen_normal_random.argtypes = [c_double, c_double]
        self.gen_normal_random.restype = c_double

        self.gen_binomial_random = libsimulat.gen_binomial_random
        self.gen_binomial_random.argtypes = [c_double]

        self.N = N


    #plot混在里面了, 抽空改了
    def test_uniform(self):
        """
        测试 0-1 均匀分布
        """
        randoms = [ self.gen_uniform_random() for i in range(self.N) ]
        print("==========均匀分布==========\n", stats.uniform.fit(randoms), "\n")
        t1 = np.linspace(-1, 2, 201)
        p, t2 = np.histogram(randoms, bins=50, normed=True) #区间和区间概率
        t2 = (t2[:-1] + t2[1:]) / 2
        plt.plot(t2, p, "ro-", t1, stats.uniform(loc=0, scale=1).pdf(t1))

    def test_rayleigh_random(self, sigma):
        """
        测试瑞利分布
        """
        randoms = [ self.gen_rayleigh_random(sigma) for i in range(self.N) ]
        print("==========瑞利分布==========\n", stats.rayleigh.fit(randoms), "\n")
        t1 = np.linspace(0, 10, 101)
        p, t2 = np.histogram(randoms, bins=100, normed=True)
        t2 = (t2[:-1] + t2[1:]) / 2
        plt.plot(t2, p, "ro-", t1, stats.rayleigh(loc=0, scale=sigma).pdf(t1))

    def test_standard_normal_random(self):
        """
        测试标准正态分布
        """
        randoms = [ self.gen_standard_normal_random() for i in range(self.N) ]
        print("==========标准正态分布==========\n", stats.norm.fit(randoms), "\n")
        t1 = np.linspace(-10, 10, 201)
        p, t2 = np.histogram(randoms, bins=100, normed=True)
        t2 = (t2[:-1] + t2[1:]) / 2
        plt.plot(t2, p, "ro-", t1, stats.norm(loc=0, scale=1).pdf(t1))

    def test_normal_random(self, mean, sigma):
        """
        测试正态分布
        """
        randoms = [ self.gen_normal_random(mean, sigma) for i in range(self.N) ]
        print("==========正态分布==========\n", stats.norm.fit(randoms), "\n")
        t1 = np.linspace(-10, 10, 201)
        p, t2 = np.histogram(randoms, bins=100, normed=True)
        t2 = (t2[:-1] + t2[1:]) / 2
        plt.plot(t2, p, "ro-", t1, stats.norm(loc=mean, scale=sigma).pdf(t1))


class TestBPSK(object):
    def __init__(self, N, DB_MAX, fading):
        libsimulat = CDLL("./libsimulat.so")

        self.N = N
        self.DB_MAX = DB_MAX
        self.fading = fading

        self.type_SNR_BER = self.Point * self.DB_MAX #创建结构数组数据类型, 为了引用内部类, 必须使用self.innerclass
        self.snr_ber = self.type_SNR_BER() #创建数组存放返回值

        self.simulat_BPSK = libsimulat.simulat_BPSK
        self.simulat_BPSK.argtypes = [c_ulong, c_uint, c_bool, self.type_SNR_BER]

    def test_simulat_BPSK(self):
        self.simulat_BPSK(self.N, self.DB_MAX, self.fading, self.snr_ber)
        SNR_db, BER = [], []
        for point in self.snr_ber:
            SNR_db.append(point.SNR_db)
            BER.append(point.BER)
        return SNR_db, BER

    class Point(Structure):
        _fields_ = [("SNR_db", c_double),
                    ("BER", c_double)]



if __name__ == '__main__':
    """
    反注释, 然后与wiki上的pdf图像比较, 结果一致
    """
    #t = TestRandom(1000000)
    #t.test_uniform()
    #t.test_rayleigh_random(1)
    #t.test_standard_normal_random()
    #t.test_normal_random(1, 2)

    #t = TestBPSK(10000000, 41, True)
    #SNR_db_fading, BER_fading = t.test_simulat_BPSK()
    #t.fading = False
    #SNR_db_no_fading, BER_no_fading = t.test_simulat_BPSK()

    #plt.title("red for fading and noise, blue for noise only")
    #plt.yscale("log")
    #plt.plot(SNR_db_fading, BER_fading, 'r^-',
    #         SNR_db_no_fading, BER_no_fading, 'bo-')
    #plt.grid(True)
    #plt.show()
