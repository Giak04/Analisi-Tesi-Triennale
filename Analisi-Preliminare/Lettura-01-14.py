import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
#import gammapy.utils.regions
#from gammapy.data import EventList
#D0114cc00_00 = '$Analisi-Gammmapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc00.lv3.fits'
#events = EventList.read(D0114cc00_00)
#events.peek()

#Lettura preliminare del livello cc00

histra = np.array([])
histdec = np.array([])
gammaness = np.array([])


with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc00['RA'])
    histdec = np.append(histdec, d0114_000cc00['DEC'])
    gammaness = np.append(gammaness, d0114_000cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p090_00001507_R_003541_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc00['RA'])
    histdec = np.append(histdec, d0114_090cc00['DEC'])
    gammaness = np.append(gammaness, d0114_090cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p180_00001508_R_003542_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_180cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_180cc00['RA'])
    histdec = np.append(histdec, d0114_180cc00['DEC'])
    gammaness = np.append(gammaness, d0114_180cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p270_00001509_R_003543_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_270cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_270cc00['RA'])
    histdec = np.append(histdec, d0114_270cc00['DEC'])
    gammaness = np.append(gammaness, d0114_270cc00['GAMMANESS'])
    hdul.close()

print(d0114_000cc00.columns)
#print(histra)
#print(histdec)

plt.hist2d(histra, histdec, bins=500)
label = plt.xlabel('Ra (°)')
label = plt.ylabel('Dec (°)')
plt.title('Grafico RA-DEC taglio cc00 14-01')
plt.savefig('./Analisi-Preliminare/Stacked-01-14/Imm_cc00.png')
plt.close()

plt.hist(gammaness, bins=50, color='lightblue', edgecolor='black')
plt.xlabel('Gammaness')
plt.ylabel('Conteggio')
plt.title('Grafico Gammaness taglio cc00 14-01')
plt.semilogy()
plt.savefig('./Analisi-Preliminare/Stacked-01-14/Gammaness_cc00')
plt.close()


#Lettura preliminare del livello cc11

histra = np.array([])
histdec = np.array([])
gammaness = np.array([])


with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc11/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc11['RA'])
    histdec = np.append(histdec, d0114_000cc11['DEC'])
    gammaness = np.append(gammaness, d0114_000cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc11/20260114_MA_Mrk421_W0.50p090_00001507_R_003541_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc11['RA'])
    histdec = np.append(histdec, d0114_090cc11['DEC'])
    gammaness = np.append(gammaness, d0114_090cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc11/20260114_MA_Mrk421_W0.50p180_00001508_R_003542_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_180cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_180cc11['RA'])
    histdec = np.append(histdec, d0114_180cc11['DEC'])
    gammaness = np.append(gammaness, d0114_180cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc11/20260114_MA_Mrk421_W0.50p270_00001509_R_003543_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_270cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_270cc11['RA'])
    histdec = np.append(histdec, d0114_270cc11['DEC'])
    gammaness = np.append(gammaness, d0114_270cc11['GAMMANESS'])
    hdul.close()

print(d0114_000cc11.columns)


plt.hist2d(histra, histdec, bins=500)
label = plt.xlabel('Ra (°)')
label = plt.ylabel('Dec (°)')
plt.title('Grafico RA-DEC taglio cc11 14-01')
plt.savefig('./Analisi-Preliminare/Stacked-01-14/Imm_cc11.png')
plt.close()

plt.hist(gammaness, bins=50, color='lightblue', edgecolor='black')
plt.xlabel('Gammaness')
plt.ylabel('Conteggio')
plt.title('Grafico Gammaness taglio cc11 14-01')
plt.semilogy()
plt.savefig('./Analisi-Preliminare/Stacked-01-14/Gammaness_cc11')
plt.close()


#Lettura preliminare del livello cc17

histra = np.array([])
histdec = np.array([])
gammaness = np.array([])


with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc17/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc17['RA'])
    histdec = np.append(histdec, d0114_000cc17['DEC'])
    gammaness = np.append(gammaness, d0114_000cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc17/20260114_MA_Mrk421_W0.50p090_00001507_R_003541_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc17['RA'])
    histdec = np.append(histdec, d0114_090cc17['DEC'])
    gammaness = np.append(gammaness, d0114_090cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc17/20260114_MA_Mrk421_W0.50p180_00001508_R_003542_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_180cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_180cc17['RA'])
    histdec = np.append(histdec, d0114_180cc17['DEC'])
    gammaness = np.append(gammaness, d0114_180cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc17/20260114_MA_Mrk421_W0.50p270_00001509_R_003543_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_270cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_270cc17['RA'])
    histdec = np.append(histdec, d0114_270cc17['DEC'])
    gammaness = np.append(gammaness, d0114_270cc17['GAMMANESS'])
    hdul.close()

print(d0114_000cc17.columns)


plt.hist2d(histra, histdec, bins=500)
label = plt.xlabel('Ra (°)')
label = plt.ylabel('Dec (°)')
plt.title('Grafico RA-DEC taglio cc17 14-01')
plt.savefig('./Analisi-Preliminare/Stacked-01-14/Imm_cc17.png')
plt.close()

plt.hist(gammaness, bins=50, color='lightblue', edgecolor='black')
plt.xlabel('Gammaness')
plt.ylabel('Conteggio')
plt.semilogy()
plt.title('Grafico Gammaness taglio cc17 14-01')
plt.savefig('./Analisi-Preliminare/Stacked-01-14/Gammaness_cc17')
plt.close()


'''
#usa nome.field(i) o nome['nome colonna'] per le colonne, non serve list

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc00/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    #hdul[1].data #mostra la prima tabella (per ASTRI sono gli eventi), così restituisce una lista di liste
    d0114_00cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    hdul.close()

#guarda le funzioni disponibili per Table, in genere i dati più interessanti sono nomi delle colonne, tipi delle colonne e lunghezza
print(d0114_00cc00.columns)
#print(d0114_00cc00[4][5])
#print(list(zip(*d0114_00cc00))[3])

plt.hist2d(list(zip(*d0114_00cc00))[3], list(zip(*d0114_00cc00))[4], bins=500)
label = plt.xlabel('Ra (°)')
label = plt.ylabel('Dec (°)')
plt.savefig('Imm_cc00.png')

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc10/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc10.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    #hdul[1].data #mostra la prima tabella (per ASTRI sono gli eventi), così restituisce una lista di liste
    d0114_00cc10 = Table(hdul[1].data) #mostra la prima tabella come tabella
    hdul.close()

plt.hist2d(list(zip(*d0114_00cc10))[3], list(zip(*d0114_00cc10))[4], bins=500)
label = plt.xlabel('Ra (°)')
label = plt.ylabel('Dec (°)')
plt.savefig('Imm_cc10.png')

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-14/cc17/20260114_MA_Mrk421_W0.50p000_00001506_R_003540_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    #hdul[1].data #mostra la prima tabella (per ASTRI sono gli eventi), così restituisce una lista di liste
    d0114_00cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    hdul.close()

plt.hist2d(list(zip(*d0114_00cc17))[3], list(zip(*d0114_00cc17))[4], bins=500)
label = plt.xlabel('Ra (°)')
label = plt.ylabel('Dec (°)')
plt.savefig('Imm_cc17.png')
'''

