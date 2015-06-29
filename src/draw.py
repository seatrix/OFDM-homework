"""
绘出误码率/信噪比曲线
"""

import sys
import csv
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from collections import namedtuple

def get_x_y(result_file):
    """
    绘出误码率/信噪比曲线
    """
    with open(file=result_file, mode='r') as data:
        data_csv = csv.reader(data)
        headings = next(data)
        SNR_BER = namedtuple("SNR_BER", headings)
        SNR_db, BER = [], []
        for line in data_csv:
            row = SNR_BER(*line)
            SNR_db.append(float(row.SNR_db))
            BER.append(float(row.BER))
    return SNR_db, BER



if __name__ == '__main__':
    files = sys.argv[1:]
    for file in files:
        SNR_db, BER = get_x_y(file)
        plt.title("")
        plt.yscale("log")
        plt.plot(SNR_db, BER, 'r^-')
    plt.grid(True)
    plt.show()
