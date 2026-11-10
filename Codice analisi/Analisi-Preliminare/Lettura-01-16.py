import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np

histra = np.array([])
histdec = np.array([])
gammaness = np.array([])


with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc00/20260116_MA_Mrk421_W0.50p000_00001532_R_003567_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc00['RA'])
    histdec = np.append(histdec, d0114_000cc00['DEC'])
    gammaness = np.append(gammaness, d0114_000cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc00/20260116_MA_Mrk421_W0.50p000_00001536_R_003571_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc00['RA'])
    histdec = np.append(histdec, d0114_000cc00['DEC'])
    gammaness = np.append(gammaness, d0114_000cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc00/20260116_MA_Mrk421_W0.50p090_00001533_R_003568_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc00['RA'])
    histdec = np.append(histdec, d0114_090cc00['DEC'])
    gammaness = np.append(gammaness, d0114_090cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc00/20260116_MA_Mrk421_W0.50p090_00001537_R_003572_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc00['RA'])
    histdec = np.append(histdec, d0114_090cc00['DEC'])
    gammaness = np.append(gammaness, d0114_090cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc00/20260116_MA_Mrk421_W0.50p180_00001534_R_003569_0201_cc00.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_180cc00 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_180cc00['RA'])
    histdec = np.append(histdec, d0114_180cc00['DEC'])
    gammaness = np.append(gammaness, d0114_180cc00['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc00/20260116_MA_Mrk421_W0.50p270_00001535_R_003570_0201_cc00.lv3.fits.gz') as hdul:
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
plt.title('Grafico RA-DEC taglio cc00 16-01')
plt.savefig('./Analisi-Preliminare/Stacked-01-16/Imm_cc00.png')
plt.close()

plt.hist(gammaness, bins=50, color='lightblue', edgecolor='black')
plt.xlabel('Gammaness')
plt.ylabel('Conteggio')
plt.title('Grafico Gammaness taglio cc00 16-01')
plt.semilogy()
plt.savefig('./Analisi-Preliminare/Stacked-01-16/Gammaness_cc00')
plt.close()


#Lettura preliminare del livello cc11

histra = np.array([])
histdec = np.array([])
gammaness = np.array([])


with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc11/20260116_MA_Mrk421_W0.50p000_00001532_R_003567_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc11['RA'])
    histdec = np.append(histdec, d0114_000cc11['DEC'])
    gammaness = np.append(gammaness, d0114_000cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc11/20260116_MA_Mrk421_W0.50p000_00001536_R_003571_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc11['RA'])
    histdec = np.append(histdec, d0114_000cc11['DEC'])
    gammaness = np.append(gammaness, d0114_000cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc11/20260116_MA_Mrk421_W0.50p090_00001533_R_003568_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc11['RA'])
    histdec = np.append(histdec, d0114_090cc11['DEC'])
    gammaness = np.append(gammaness, d0114_090cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc11/20260116_MA_Mrk421_W0.50p090_00001537_R_003572_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc11['RA'])
    histdec = np.append(histdec, d0114_090cc11['DEC'])
    gammaness = np.append(gammaness, d0114_090cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc11/20260116_MA_Mrk421_W0.50p180_00001534_R_003569_0201_cc11.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_180cc11 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_180cc11['RA'])
    histdec = np.append(histdec, d0114_180cc11['DEC'])
    gammaness = np.append(gammaness, d0114_180cc11['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc11/20260116_MA_Mrk421_W0.50p270_00001535_R_003570_0201_cc11.lv3.fits.gz') as hdul:
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
plt.title('Grafico RA-DEC taglio cc11 16-01')
plt.savefig('./Analisi-Preliminare/Stacked-01-16/Imm_cc11.png')
plt.close()

plt.hist(gammaness, bins=50, color='lightblue', edgecolor='black')
plt.xlabel('Gammaness')
plt.ylabel('Conteggio')
plt.title('Grafico Gammaness taglio cc11 16-01')
plt.semilogy()
plt.savefig('./Analisi-Preliminare/Stacked-01-16/Gammaness_cc11')
plt.close()


#Lettura preliminare del livello cc17

histra = np.array([])
histdec = np.array([])
gammaness = np.array([])


with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc17/20260116_MA_Mrk421_W0.50p000_00001532_R_003567_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc17['RA'])
    histdec = np.append(histdec, d0114_000cc17['DEC'])
    gammaness = np.append(gammaness, d0114_000cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc17/20260116_MA_Mrk421_W0.50p000_00001536_R_003571_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_000cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_000cc17['RA'])
    histdec = np.append(histdec, d0114_000cc17['DEC'])
    gammaness = np.append(gammaness, d0114_000cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc17/20260116_MA_Mrk421_W0.50p090_00001533_R_003568_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc17['RA'])
    histdec = np.append(histdec, d0114_090cc17['DEC'])
    gammaness = np.append(gammaness, d0114_090cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc17/20260116_MA_Mrk421_W0.50p090_00001537_R_003572_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_090cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_090cc17['RA'])
    histdec = np.append(histdec, d0114_090cc17['DEC'])
    gammaness = np.append(gammaness, d0114_090cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc17/20260116_MA_Mrk421_W0.50p180_00001534_R_003569_0201_cc17.lv3.fits.gz') as hdul:
    hdul[0].header
    hdul.info() #informazioni sul contenuto
    d0114_180cc17 = Table(hdul[1].data) #mostra la prima tabella come tabella
    histra = np.append(histra, d0114_180cc17['RA'])
    histdec = np.append(histdec, d0114_180cc17['DEC'])
    gammaness = np.append(gammaness, d0114_180cc17['GAMMANESS'])
    hdul.close()

with fits.open('./Analisi-Gammapy/Mrk421/data/2026-01-16/cc17/20260116_MA_Mrk421_W0.50p270_00001535_R_003570_0201_cc17.lv3.fits.gz') as hdul:
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
plt.title('Grafico RA-DEC taglio cc17 16-01')
plt.savefig('./Analisi-Preliminare/Stacked-01-16/Imm_cc17.png')
plt.close()

plt.hist(gammaness, bins=50, color='lightblue', edgecolor='black')
plt.xlabel('Gammaness')
plt.ylabel('Conteggio')
plt.title('Grafico Gammaness taglio cc17 16-01')
plt.semilogy()
plt.savefig('./Analisi-Preliminare/Stacked-01-16/Gammaness_cc17')
plt.close()
