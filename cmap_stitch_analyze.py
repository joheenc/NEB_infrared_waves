from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import detrend
from scipy.fft import fft
from astropy.timeseries import LombScargle

def cmap_stitch_analyze(infile, outfile="out.png", nbins=12, cutoff=0.15, gradient=6):
    observations = []
    fitsfiles = [x.strip() for x in open(infile, 'r').readlines()]
    for file in fitsfiles:
        obs = np.flip(fits.open(file)[0].data, 0)
        observations.append(obs)

    nebStart = 135
    nebEnd = 156
    binstarts = [int(x) for x in np.linspace(0, 720, nbins+1)]
    keepbins = []
    lmap = np.full((360, 720), 0.)
    cmaps = []

    #choose the maps with the highest average brightness over the bin
    for bin in range(len(binstarts)-1):
        bestObs = None
        maxBrightness = 0
        for obs in observations:
            avgBrightness = np.mean(obs[:, binstarts[bin]:binstarts[bin+1]])
            if avgBrightness > maxBrightness:
                bestObs = obs
                maxBrightness = avgBrightness
        if maxBrightness < cutoff:
            lmap[:, binstarts[bin]:binstarts[bin+1]] = np.full((len(lmap), binstarts[bin+1]-binstarts[bin]), -1)
        else:
            lmap[:, binstarts[bin]:binstarts[bin+1]] = bestObs[:, binstarts[bin]:binstarts[bin+1]]
            keepbins.append(bin)
        cmaps.append(bestObs)

    #stitch together at the bin edges
    for idx in range(len(keepbins)-1):
        bin = keepbins[idx]
        nextBin = keepbins[idx+1]
        if bin + 1 == nextBin:
            stitchStart = binstarts[nextBin]-gradient//2
            stitchEnd = binstarts[nextBin]+gradient//2
            lmap[:, stitchStart:stitchEnd] = (cmaps[bin][:, stitchStart:stitchEnd] + cmaps[nextBin][:, stitchStart:stitchEnd]) / 2.0
    lmap = lmap[:, np.argwhere(np.mean(lmap, 0)!=-1).flatten()]
    neb = lmap[nebStart:nebEnd, :]
    nebAvg = np.mean(neb, 0)
    #detrend with median filter, bins=separate cmap regions
    avg_detrended = detrend(nebAvg, type='constant', bp=[int(i) for i in np.linspace(0, len(lmap), len(keepbins))])
    fig, axes = plt.subplots(2,1, figsize=(20,20))
    axes[0].imshow(lmap, cmap='gist_heat')
    axes[1].imshow(neb, cmap='gist_heat')
    plt.savefig(outfile)

    #fourier analysis
    fourier = fft(nebAvg)
    xf = np.arange(len(nebAvg)//2)
    fftRaw = (720 / len(nebAvg) * np.flip(np.argsort(fourier[0:len(nebAvg)//2]))[:10])
    fourier_detrended = fft(avg_detrended)
    fftDetrended = (720 / len(nebAvg) * np.flip(np.argsort(fourier_detrended[0:len(nebAvg)//2]))[:10])

    #Lomb-Scargle periodogram
    lon = np.arange(len(nebAvg))
    lon_normed = (lon-np.min(lon))/len(lon)
    frequency, power = LombScargle(lon_normed,nebAvg).autopower()
    lsRaw = (720 / len(lon_normed) * frequency[np.flip(np.argsort(power[0:len(frequency)//20]))][:10])
    frequency, power_detrended = LombScargle(lon_normed,avg_detrended).autopower()
    lsDetrended = (720 / len(lon_normed) * frequency[np.flip(np.argsort(power_detrended[0:len(frequency)//20]))][:10])

    return fftRaw, fftDetrended, lsRaw, lsDetrended
