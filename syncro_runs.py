# Simulating trajectories using QTT on @munin.uio.no
# Use the lwQTT class, with is a light weight class only using numpy

from lwQTT import *

# first need to fix lwQTT so it runs purely on numpy

### ! Fixed system hamiltonian to include 0.5 factor and wrap_freq to return angular frequency

def measured_frequency_result(result, simtime):

    measured_freq = []
    std_freq_result = []

    for MC in result:

        class_instance = QTT(env='z', meas_basis='x')
        tmp = class_instance.wrap_freq(MC, simtime)
        measured_freq.append(tmp)
        tmp = class_instance.std_freq(MC, simtime)
        std_freq_result.append(tmp)
    
    Omega   = np.asarray(measured_freq)
    stdfreq = np.asarray(std_freq_result)

    return Omega, stdfreq


def main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=10000):

    seed = time.time()

    env = 'z'
    meas_basis = 'x'
    dt = simtime/N
    
    ### main simulation

    class_instance = QTT(env, meas_basis, theta=theta, seed=seed)

    result = class_instance.freqSimulationResult(S, N, burnin, psi0, simtime, ncpu, delta, epsilon, res)

    datapath = 'data/'
    filename = path + datapath + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin) 

    np.save(filename, result)

    ### calculating frequencies and standard deviation

    measfreqres, stdfreq = measured_frequency_result(result, simtime)

    freqpath = 'dtune_smalleps/'
    savename = path + freqpath + 'OMEGA_' + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin)

    np.save(savename, measfreqres)

    savenamestd = path + freqpath + 'STD_' + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin)

    np.save(savenamestd, stdfreq)


    return 0


def mainT01(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=10000):

    seed = time.time()

    env = 'z'
    meas_basis = 'x'
    dt = simtime/N

    ### main simulation

    class_instance = QTT(env, meas_basis, theta=theta, temperature=0.1, seed=seed)

    result = class_instance.freqSimulationResult(S, N, burnin, psi0, simtime, ncpu, delta, epsilon, res)

    datapath = 'data/'
    filename = path + datapath + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin) + '_T01'

    np.save(filename, result)

    ### calculating frequencies and standard deviation

    measfreqres, stdfreq = measured_frequency_result(result, simtime)

    freqpath = 'dtune_smalleps_T01/'
    savename = path + freqpath + 'OMEGA_' + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin) + '_T01'

    np.save(savename, measfreqres)

    savenamestd = path + freqpath + 'STD_' + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin) + '_T01'

    np.save(savenamestd, stdfreq)


    return 0


def mainZY(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=10000):

    seed = time.time()

    env = 'z'
    meas_basis = 'y'
    dt = simtime/N

    ### main simulation

    class_instance = QTT(env, meas_basis, theta=theta, temperature=0.5, seed=seed)

    result = class_instance.freqSimulationResult(S, N, burnin, psi0, simtime, ncpu, delta, epsilon, res)

    datapath = 'data/'
    filename = path + datapath + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin) + '_T01'

    np.save(filename, result)

    ### calculating frequencies and standard deviation

    measfreqres, stdfreq = measured_frequency_result(result, simtime)

    freqpath = 'dtune_envZmeasY/'
    savename = path + freqpath + 'OMEGA_' + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin) + '_T01'

    np.save(savename, measfreqres)

    savenamestd = path + freqpath + 'STD_' + 'eps' + str(epsilon).replace('.', '') + \
        '_theta' + str(theta).replace('.', '') + '_dt' + str(dt).replace('.', '') + \
        '_S' + str(S) + '_N' + str(N) + '_freqN' + str(delta.size) + '_burninN' + str(burnin) + '_T01'

    np.save(savenamestd, stdfreq)


    return 0




path = '../../../njord/erlenalo/'

### test test test
"""
S = 4
N = 100
burnin = 100
res = 10

epsilon = 0.01
ndelta = 32
delta = np.linspace(-2.5*epsilon, 2.5*epsilon, ndelta )
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 4
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)

"""

""" # eps 0005-002
epsilon = 0.005

S = 512
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 36
delta = np.linspace(-3.0*epsilon, 3.0*epsilon, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 256
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)



epsilon = 0.01

S = 512
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 36
delta = np.linspace(-3.0*epsilon, 3.0*epsilon, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 256
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)



epsilon = 0.015

S = 512
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 36
delta = np.linspace(-3.0*epsilon, 3.0*epsilon, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 256
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)




epsilon = 0.02

S = 512
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 36
delta = np.linspace(-3.0*epsilon, 3.0*epsilon, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
seed = 1947571
ncpu = 256
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)

"""

# eps 0007-0001
"""
epsilon = 0.001

S = 256
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 48
delta = np.linspace(-0.03, 0.03, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
ncpu = 256
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)


epsilon = 0.0035

S = 256
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 48
delta = np.linspace(-0.03, 0.03, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
ncpu = 256
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)


epsilon = 0.007

S = 256
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 48
delta = np.linspace(-0.03, 0.03, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
ncpu = 256
psi0 = xplus

main(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)
"""


# eps 0007-0001, temperature 0.1
"""
epsilon = 0.001

S = 256
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 48
delta = np.linspace(-0.03, 0.03, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
ncpu = 256
psi0 = xplus

mainT01(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)


epsilon = 0.0035

S = 256
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 48
delta = np.linspace(-0.03, 0.03, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
ncpu = 256
psi0 = xplus

mainT01(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)


epsilon = 0.007

S = 256
N = 400000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 48
delta = np.linspace(-0.03, 0.03, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
ncpu = 256
psi0 = xplus

mainT01(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)
"""


epsilon = 0.003

S = 128
N = 200000
burnin = 150000
res = 10000

print('simulating for epsilon : ', epsilon)
ndelta = 24
delta = np.linspace(-0.015, 0.015, ndelta)
theta = 0.01

dt = 0.01
simtime = N*dt
ncpu = 128
psi0 = xplus


mainZY(S, N, theta, simtime, psi0, ncpu, burnin, delta, epsilon, path, res=res)

# end