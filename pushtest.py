# testing pushing programs to gitclust and syncing with munin.uio.no

from lwQTT import *

# test read write to njord

# idé til en annen gang -> hvis du er usikker på saving loading paths, lag en test slik som denne
# med et enkelt array, skriv inn path fra command line, kjør testen, sett prgm på vent og gi
# resultat fra testen og så prompt continue :)

savearray = np.arange(9)
print('saved array=', savearray)

#path = 'data/' # local path

path = '../../../njord/erlenalo/data/'

savename = path + 'savearray'
np.save(savename, savearray)

loadname = savename + '.npy'
loadarray= np.load(loadname)

print('load array=', loadarray)

print('test success')


