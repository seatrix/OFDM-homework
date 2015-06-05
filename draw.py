"""
绘出误码率/信噪比曲线
"""

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
    SNR_db_fading, BER_fading = get_x_y("fading.csv")
    SNR_db_no_fading, BER_no_fading = get_x_y("no_fading.csv")

    plt.title("red for fading and noise, blue for noise only")
    plt.yscale("log")
    plt.plot(SNR_db_fading, BER_fading, 'r^-',
             SNR_db_no_fading, BER_no_fading, 'bo-')
    plt.grid(True)
    plt.show()
