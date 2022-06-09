# longer simulation, fewer points runng on 17th of may

from QTT import *

def mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=10000):

    # input variables

    seed = 19037840

    env = 'z'
    meas_basis = 'x'
    dt = simtime/N
    
    class_instance = QTT(env, meas_basis, theta=theta, seed=seed)

    result = class_instance.freqSimulationResult(S, N, burnin, psi0, simtime, ncpu, delta, epsilon, res)


    path = '../data/multi_eps_data/delta_eps/'
    filename = path + 'theta' + str(theta).replace('.', '') + \
        '_eps' + str(epsilon).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S_' + str(S) + '_N_' + str(N) + '_freqN_' + str(delta.size) + '_burninN_' + str(burnin) 

    np.save(filename, result)

    return 0

def plotDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon):

    dt = simtime/N

    path = '../data/freq/'
    loadname = path + 'theta' + str(theta).replace('.', '') + \
        '_eps' + str(epsilon).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S_' + str(S) + '_N_' + str(N) + '_freqN_' + str(delta.size) + '_burninN_' + str(burnin) 
    
    result = np.load(loadname + '.npy')


    measured_freq = []

    ### OLD METHOD
    """
    for MC in result:

        class_instance = QTT(env='z', meas_basis='x', theta=theta, seed=seed)
        #tmp = class_instance.fft_freq(MC, dt)
        tmp = class_instance.old_freq(MC, simtime)
        measured_freq.append(tmp)

    Omega = np.asarray(measured_freq)

    plt.style.use('ggplot')
    fig = plt.figure()
    
    ax1 = fig.add_subplot(121)
    ax1.plot(delta, Omega-delta)
    ax1.scatter(delta, Omega-delta, s=20.0, color='blue')
    ax1.set_xlabel(r'$\Delta$')
    ax1.set_ylabel(r'$\Omega - \Delta$')

    ax2 = fig.add_subplot(122)
    ax2.plot(delta, Omega)
    ax2.scatter(delta, Omega, s=20.0, color='blue')
    ax2.set_xlabel(r'$\Delta$')
    ax2.set_ylabel(r'$\Omega$')
    
    fig.suptitle(' [Old method] Synchronization plot, with epsilon = ' + str(epsilon))
    fig.set_size_inches(16,9)

    path = '../figure/freq/freq4/'
    figname = path + 'eps' + str(epsilon).replace('.', '') + \
        '_Nfreq' + str(delta.size) + '_oldmethod' +  \
        '_deltavals_' + str(delta[0]).replace('.', '') + '_' + str(delta[-1]).replace('.', '') + '.png'

    fig.savefig(figname, dpi=400)
    """
    
    ### FFT METHOD
    
    for MC in result:

        class_instance = QTT(env='z', meas_basis='x', theta=theta, seed=seed)
        tmp = class_instance.fft_freq(MC, dt)
        measured_freq.append(tmp)
    
    
    
    
    Omega = np.asarray(measured_freq)

    plt.style.use('ggplot')
    fig = plt.figure()
    
    ax1 = fig.add_subplot(121)
    ax1.plot(delta, Omega-delta)
    ax1.scatter(delta, Omega-delta, s=20.0, color='blue')
    ax1.set_xlabel(r'$\Delta$')
    ax1.set_ylabel(r'$\Omega - \Delta$')

    ax2 = fig.add_subplot(122)
    ax2.plot(delta, Omega)
    ax2.scatter(delta, Omega, s=20.0, color='blue')
    ax2.set_xlabel(r'$\Delta$')
    ax2.set_ylabel(r'$\Omega$')
    
    fig.suptitle(' [fft method] Synchronization plot, with epsilon = ' + str(epsilon))
    fig.set_size_inches(16,9)

    path = '../figure/freq/freq4/'
    figname = path + 'eps' + str(epsilon).replace('.', '') + \
        '_Nfreq' + str(delta.size) + '_fftmethod' + \
        '_deltavals_' + str(delta[0]).replace('.', '') + '_' + str(delta[-1]).replace('.', '') + '.png'

    fig.savefig(figname, dpi=400)
    
    




    """
    Omega = np.concatenate([  Omega[:10]  ,  -1.0 * Omega[10:]  ])

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(delta, Omega-delta)
    ax.scatter(delta, Omega-delta, s=20.0, color='blue')
    ax.set_xlabel(r'$\Delta$')
    ax.set_ylabel(r'$\Omega - \Delta$')
    
    fig.suptitle('[TESTING taking negative frequencies] Synchronization plot, with epsilon = ' + str(epsilon))
    fig.set_size_inches(16,9)

    path = '../figure/freq/freq4/'
    figname = path + 'eps' + str(epsilon).replace('.', '') + \
        '_Nfreq' + str(delta.size) + '_alt' + '.png'

    fig.savefig(figname, dpi=400)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(delta, Omega)
    ax.scatter(delta, Omega, s=20.0, color='blue')
    ax.set_xlabel(r'$\Delta$')
    ax.set_ylabel(r'$\Omega$')
    
    fig.suptitle('[TESTING taking negative frequencies] Synchronization plot, with epsilon = ' + str(epsilon))
    fig.set_size_inches(16,9)

    path = '../figure/freq/freq4/'
    figname = path + 'eps' + str(epsilon).replace('.', '') + \
        '_Nfreq' + str(delta.size) + '_alt2' + '.png'

    fig.savefig(figname, dpi=400)

    """

    return 0















# some run ideas. Try running for 10* longer, but only 8 points at eps, 0.5eps, 0.3ep, 0.1eps




""" THIS IS A GOOD RUN """
### test test test
S = 4
N = 100
burnin = 100

epsilon = 0.01
delta = np.array([-0.1, -0.04, -0.02, -0.01, -0.005, -0.002, -0.001])
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus


#mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=100)
###

"""

S = 64
N = 400000
burnin = 100000

epsilon = 0.01
delta = np.array([-0.1, -0.04, -0.02, -0.01, -0.005, -0.002, -0.001])
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus


mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=40000)



S = 64
N = 400000
burnin = 100000

epsilon = 0.02
delta = np.array([-0.1, -0.04, -0.02, -0.01, -0.005, -0.002, -0.001])
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus


mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=40000)



S = 64
N = 400000
burnin = 100000

epsilon = 0.04
delta = np.array([-0.1, -0.04, -0.02, -0.01, -0.005, -0.002, -0.001])
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus


mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=40000)



S = 64
N = 400000
burnin = 100000

epsilon = 0.06
delta = np.array([-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, -0.005, -0.002, -0.001, \
                    0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06])
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus


mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=40000)

"""



"""

S = 80
N = 300000
burnin = 100000

epsilon = 0.005 # times 2 = real eps value
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-0.0025, 0.0025, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=60000)


S = 80
N = 300000
burnin = 100000

epsilon = 0.075
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-0.0025, 0.0025, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=60000)


S = 80
N = 300000
burnin = 100000

epsilon = 0.01
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-0.0025, 0.0025, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=60000)



S = 80
N = 300000
burnin = 100000

epsilon = 0.015
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-0.0025, 0.0025, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=60000)


S = 80
N = 300000
burnin = 100000

epsilon = 0.02
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-0.0025, 0.0025, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=60000)

"""




S = 64
N = 400000
burnin = 130000

epsilon = 0.04
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-1*epsilon, epsilon, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=40000)



S = 64
N = 400000
burnin = 130000

epsilon = 0.06
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-1*epsilon, epsilon, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=40000)


S = 64
N = 400000
burnin = 130000

epsilon = 0.08
print('simulating for epsilon : ', epsilon)
delta = np.linspace(-1*epsilon, epsilon, 24)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

mainDeltaArray(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, res=40000)



### remember to change seed if simulating over same parameters

### end