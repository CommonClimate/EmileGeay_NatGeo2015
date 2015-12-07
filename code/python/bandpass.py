'''
This script computes a 2-7 year bandpass filtered coral series and plot it above the original one
Yuxin Zhou
'''
from scipy.signal import butter, lfilter, filtfilt

def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandpass')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

# name = 'CCSM4_historical_r1i1p1_185001-200512_Indo_Pacific_coral.nc'
# coralVar = DJF_mean.averager(name, 240, 300)
# coral = coralVar.getValue()

# fs = 1680
# lowcut = 20
# highcut = 70
# f_order = 3
# y = butter_bandpass_filter(coral, lowcut, highcut, fs, f_order)

# plt.clf()
# f1, = plt.plot(coral, alpha=0.3)
# f2, = plt.plot(y)
# plt.legend([f1,f2],['Original','filtered'])
# plt.title('monthly data, order = %s' % f_order)
# plt.savefig('CCSM4_coral.pdf')
