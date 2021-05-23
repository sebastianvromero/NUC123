import os
import numpy  as np
import pandas as pd

import matplotlib.pyplot as plt

def readData(fname = 'BBN_data.csv'):

    if not os.path.isfile(fname):

        dataNuc = np.loadtxt('nuc123.dat', skiprows = 11 , max_rows = 141)
        dataNHe = np.loadtxt('nuc123.dat', skiprows = 9  , max_rows = 1, dtype = str)

        dataCos = np.loadtxt('nuc123.dat', skiprows = 155, max_rows = 141)
        dataCos = np.delete(dataCos, 0, axis = 1)
        dataCHe = np.loadtxt('nuc123.dat', skiprows = 153, max_rows = 1, dtype = str)
        dataCHe = dataCHe[1 : ]

        data = np.hstack((dataNuc, dataCos))
        head = np.hstack((dataNHe, dataCHe))

        df = pd.DataFrame(data, columns = head)
        df.to_csv(fname, index = False)

    else:

        df = pd.read_csv(fname)
    
    return df

df = readData()

plt.loglog(df['T'], df['N/H'], label = '$n$')
plt.loglog(df['T'], df['P'], label = '$p$')
plt.loglog(df['T'], df['D/H'], label = '$d$')
plt.loglog(df['T'], df['T/H'], label = '$t$')
plt.loglog(df['T'], df['He3/H'], label = '$^3$He')
plt.loglog(df['T'], df['He4'], label = '$^4$He')
plt.loglog(df['T'], df['Li6/H'], label = '$^6$Li')
plt.loglog(df['T'], df['Li7/H'], label = '$^7$Li')
plt.loglog(df['T'], df['Be7/H'], label = '$^7$Be')
plt.loglog(df['T'], df['Li8/H&up'], label = '$^8$Li and above')
plt.xlabel('$t$ (s)')
plt.ylabel('Relative abundance')
plt.grid(which = 'both')
plt.xlim(min(df['T']), max(df['T']))
plt.yticks([1e-27, 1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1e0])
plt.ylim(1e-26, 1e1)
plt.legend()
plt.show()