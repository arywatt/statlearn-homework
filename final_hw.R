# *********************************************************************** #
#                         ** Statistical Learning **                      #
#                 ** Final Homework / Some hints with R **                #
# *********************************************************************** #

#' ---
#' title:  "Statistical Learning"
#' author: "Pierpaolo Brutti"
#' date:   "Final Homework / Some hints with R"
#' ---

# Packages ----------------------------------------------------------------

suppressMessages(require(signal, quietly = T))
# library(help = signal)
suppressMessages(require(audio, quietly = T)) 
# library(help = audio)

suppressMessages(require(wrassp,  quietly = T))
# library(help = wrassp)

# library(warbleR)
# library(help = warbleR)
# library(tuneR)
# library(help = tuneR)
# library(help = audiolyzR)

# Read --------------------------------------------------------------------

library(help = wrassp)
?read.AsspDataObj

# My path to the data
auPath <- "data_example"
# List the .au files
auFiles <- list.files(auPath, pattern=glob2rx('*.au'), full.names=TRUE)

# Load an audio file, e.g. the first one in the list above
x <- read.AsspDataObj(auFiles)
x
# Take a look
names(x)
dim(x$audio)
head(x$audio)   # actual sound wave

attributes(x)
attributes(x)$sampleRate

# Number of records
numRecs.AsspDataObj(x)

# and we can of course also plot these samples 
# (only plot every 10th element to accelerate plotting)
plot(seq(0,numRecs.AsspDataObj(x) - 1, 10) / rate.AsspDataObj(x), 
     x$audio[c(TRUE, rep(FALSE,9))], 
     type='l', 
     xlab='time (s)', 
     ylab='Audio samples')

# Play --------------------------------------------------------------------

suppressMessages(require(tuneR, quietly = T))
?Wave
?WaveMC
?getWavPlayer
?play

# Transform
xwv <- Wave( as.numeric(x$audio), samp.rate = rate.AsspDataObj(x), bit = 16)
setWavPlayer("afplay") # needed for MAC OS not for Windoze
# play(xwv) # do not run in Markdown
plot(xwv)

xw.dwn = downsample(xwv, samp.rate = 11025)
xwv
xw.dwn@samp.rate

xw.dwn

?attributes

str(xwv)
  
# Powerspectrum -----------------------------------------------------------

?powspec
out = powspec( as.numeric(x$audio), sr = rate.AsspDataObj(x))
# The output is a matrix, where each column represents a power spectrum 
# for a given time frame and each row represents a frequency.
dim(out)
# 
image(out)

# Spetral info ------------------------------------------------------------

# calculate the fundamental frequency contour
f0vals = ksvF0("data_example/f1.au", toFile=F)
# plot the fundamental frequency contour
plot(seq(0,numRecs.AsspDataObj(f0vals) - 1) / rate.AsspDataObj(f0vals) +
       attr(f0vals, 'startTime'),
     f0vals$F0, 
     type='l', 
     xlab='time (s)', 
     ylab='F0 frequency (Hz)')


# STFT --------------------------------------------------------------------

suppressMessages(require(signal, quietly = T))
library(help = signal)

# Short Time Fourier Transform (default values)
# fs       <- rate.AsspDataObj(x)   # sampling rate
fs       <- xw.dwn@samp.rate
winsize  <- 2048
nfft     <- 2048
hopsize  <- 512
noverlap <- winsize - hopsize
# sp  <- specgram(x = x$audio, n = nfft, Fs = fs, window = winsize, overlap = noverlap)
sp  <- specgram(x = xw.dwn@left, n = nfft, Fs = fs, window = winsize, overlap = noverlap)
sp
names(sp)
class(sp)

# Setup -------------------------------------------------------------------

# STFT
winsize  <- 2048
nfft     <- 2048
hopsize  <- 512
noverlap <- winsize - hopsize

# Frequency bands selection
nb   <- 2^3
lowB <- 100
eps  <- .Machine$double.eps

# Number of seconds of the analyzed window
corrtime     <- 15

# Analysis ----------------------------------------------------------------

# Sampling rate
# fs <- rate.AsspDataObj(x)
# Short-time fourier transform
# sp       <- specgram(x = x$audio, n = nfft, Fs = fs, window = winsize, overlap = noverlap)
ntm      <- ncol(sp$S)  # number of (overlapping) time segments

# Energy of bands
fco    <- round( c(0, lowB*(fs/2/lowB)^((0:(nb-1))/(nb-1)))/fs*nfft )
energy <- matrix(0, nb, ntm)
for (tm in 1:ntm){
  for (i in 1:nb){
    lower_bound <- 1 + fco[i]
    upper_bound <- min( c( 1 + fco[i + 1], nrow(sp$S) ) )
    energy[i, tm] <- sum( abs(sp$S[ lower_bound:upper_bound, tm ])^2 )
  }
}
energy[energy < eps] <- eps
energy = 10*log10(energy)

dim(energy)

# Init pdf-plot
# pdf(file = paste(gsub(".wav","",nomefl),"-plot.pdf", sep = ""), height = 11, width = 8)
par(mfrow = c(2,1))

# Take a look
matplot(sp$t, t(energy), type = "l", lty = 1, 
        main = "Energy Band Envelopes",
        xlab = "Time (secs)", ylab = "Energy", xaxt = "n",
        col = viridis::viridis(nb+1, .3), lwd = 1.5)
axis(1, seq(min(sp$t), max(sp$t), length.out = 15),
     round(seq(0,30,length.out = 15)))
legend( "bottom", paste("Band-",1:nb,sep=""), lwd = 8, col =  viridis::viridis(nb+1, .3),
        horiz = TRUE, bty = "n", cex = .7 )

# 6sec zoom
matplot(sp$t, t(energy), type = "l", lty = 1, 
        main = "Energy Band Envelopes / 6s-Zoom",
        xlab = "Time (secs)", ylab = "Energy", xaxt = "n", xlim = c(0, 6.5),
        col = viridis::viridis(nb+1, .3), lwd = 1.5)
axis(1, seq(min(sp$t), max(sp$t), length.out = 15), 
     round(seq(0,30,length.out = 15)))
legend( "bottom", paste("Band-",1:nb,sep=""), lwd = 8, col =  viridis::viridis(nb+1, .3),
        horiz = TRUE, bty = "n", cex = .7 )
