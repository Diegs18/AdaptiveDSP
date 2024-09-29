
import numpy as np

import matplotlib.pyplot as plt
import scipy.fft as fft
import scipy.signal as sig



def freq_resp(b = [1], a = [1, -.9]):
    
    # Calculate the frequency response
    w, h = sig.freqz(b, a)

    # Plot the magnitude and phase response
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(w/np.pi, 20 * np.log10(abs(h)))
    plt.xlabel('Normalized Frequency (rad/sample)')
    plt.ylabel('Magnitude (dB)')
    plt.title('Magnitude Response')

    plt.subplot(2, 1, 2)
    plt.plot(w, np.angle(h, deg=True))
    plt.xlabel('Normalized Frequency (rad/sample)')
    plt.ylabel('Phase (degrees)')
    plt.title('Phase Response')

    plt.tight_layout()
    plt.show()


def imp_resp(b = [1], a = [1, -.9]):
    
    system = (b, a, 1)
    # Calculate the frequency response
    t, y = sig.dimpulse(system)
    y=np.squeeze(y)
    # Plot the magnitude and phase response
    plt.figure(figsize=(10, 6))
    plt.stem(t, y)
    plt.xlabel('Time (N samples)')
    plt.ylabel('Response')
    plt.title('Impulse Response')

    plt.tight_layout()
    plt.show()


def calc_autocorr(x):
    #x is the signal
    #N is the number of points, optional
    total = 0
    N=len(x)
    rxx = np.zeros(N)
    for m in range(0, N):
        total = 0
        for n in range(0, N-m):
            total += x[n+m]*x[n]
        rxx[m] = total/N

    plt.plot(rxx)
    plt.xlabel('m')
    plt.ylabel('Autocorrelation')
    return rxx



import numpy as np
def LMS(d, x, L, mu, N):
    #d  = desired signal
    #x  = input siganl with noise
    #L  = filter order, num of coefficients - 1
    #mu = step size
    #N  = number of iterations
    # return the filter estimation
    y = np.zeros_like(d)
    e = np.zeros_like(d)
    bhat = np.zeros([1, L+1])
    for n in range(L,N):
        xv = np.flip(x[n-L:n+1])
        y = np.dot(bhat, xv)
        e[n] = d[n]-y
        bhat = bhat + 2*mu*e[n]*xv

    freq_resp(bhat.reshape(L+1,1), [1])
    plt.plot(10*np.log10(e*e+.0000001))
    plt.ylabel('MSE (dB)')
    plt.xlabel('Iteration')
    string = 'Mu = ' + str(mu)
    plt.title(string)
    plt.show()
    return bhat


import numpy as np
def LMS2(d, x, L, mu, N):
    #similar to LMS function but dont plot and return the error instead of filter
    #d  = desired signal
    #x  = input siganl with noise
    #L  = filter order, num of coefficients - 1
    #mu = step size
    #N  = number of iterations
    # return the filter estimation
    y = np.zeros_like(d)
    e = np.zeros_like(d)
    bhat = np.zeros([1, L+1])
    for n in range(L,N):
        xv = np.flip(x[n-L:n+1])
        y = np.dot(bhat, xv)
        e[n] = d[n]-y
        bhat = bhat + 2*mu*e[n]*xv
    return e





