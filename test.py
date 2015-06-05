"""
测试均匀, 高斯, 和瑞利分布
"""

from ctypes import CDLL, c_double
from scipy import stats
import numpy as np
import matplotlib.pyplot as pl

class Test(object):
    """
    用fit()方法对随机变量取样序列进行拟合
    比较随机变量X的概率密度函数和对数组x进行直方图统计的结果
    http://hyry.dip.jp/tech/book/page/scipy/scipy_stats.html
    """
    def __init__(self, N):
        libc = CDLL("libc.so.6")
        librandom = CDLL("./librandom.so")

        libc.srand(libc.time(None))

        self.gen_uniform_random = librandom.gen_uniform_random
        self.gen_uniform_random.restype = c_double

        self.gen_rayleigh_random = librandom.gen_rayleigh_random
        self.gen_rayleigh_random.argtypes = [c_double]
        self.gen_rayleigh_random.restype = c_double

        self.gen_standard_normal_random = librandom.gen_standard_normal_random
        self.gen_standard_normal_random.restype = c_double

        self.gen_normal_random = librandom.gen_normal_random
        self.gen_normal_random.argtypes = [c_double, c_double]
        self.gen_normal_random.restype = c_double

        self.gen_binomial_random = librandom.gen_binomial_random
        self.gen_binomial_random.argtypes = [c_double]

        self.N = N



    def test_uniform(self):
        """
        测试 0-1 均匀分布
        """
        randoms = [ self.gen_uniform_random() for i in range(self.N) ]
        print("==========均匀分布==========\n", stats.uniform.fit(randoms), "\n")
        t1 = np.linspace(-1, 2, 201)
        p, t2 = np.histogram(randoms, bins=50, normed=True) #区间和区间概率
        t2 = (t2[:-1] + t2[1:]) / 2
        pl.plot(t2, p, "ro-", t1, stats.uniform(loc=0, scale=1).pdf(t1))

    def test_rayleigh_random(self, sigma):
        """
        测试瑞利分布
        """
        randoms = [ self.gen_rayleigh_random(sigma) for i in range(self.N) ]
        print("==========瑞利分布==========\n", stats.rayleigh.fit(randoms), "\n")
        t1 = np.linspace(0, 10, 101)
        p, t2 = np.histogram(randoms, bins=100, normed=True)
        t2 = (t2[:-1] + t2[1:]) / 2
        pl.plot(t2, p, "ro-", t1, stats.rayleigh(loc=0, scale=sigma).pdf(t1))

    def test_standard_normal_random(self):
        """
        测试标准正态分布
        """
        randoms = [ self.gen_standard_normal_random() for i in range(self.N) ]
        print("==========标准正态分布==========\n", stats.norm.fit(randoms), "\n")
        t1 = np.linspace(-10, 10, 201)
        p, t2 = np.histogram(randoms, bins=100, normed=True)
        t2 = (t2[:-1] + t2[1:]) / 2
        pl.plot(t2, p, "ro-", t1, stats.norm(loc=0, scale=1).pdf(t1))

    def test_normal_random(self, mean, sigma):
        """
        测试正态分布
        """
        randoms = [ self.gen_normal_random(mean, sigma) for i in range(self.N) ]
        print("==========正态分布==========\n", stats.norm.fit(randoms), "\n")
        t1 = np.linspace(-10, 10, 201)
        p, t2 = np.histogram(randoms, bins=100, normed=True)
        t2 = (t2[:-1] + t2[1:]) / 2
        pl.plot(t2, p, "ro-", t1, stats.norm(loc=mean, scale=sigma).pdf(t1))




if __name__ == '__main__':
    """
    反注释, 然后与wiki上的pdf图像比较, 结果一致
    """
    t = Test(1000000)
    #t.test_uniform()
    #t.test_rayleigh_random(1)
    #t.test_standard_normal_random()
    #t.test_normal_random(1, 2)
    pl.show()
