# Acoustic RADAR in Julia
Just small experiments in Julia

## Description
Small experiments in Julia with audio signals and RADAR. Implemented the chirp generation and saving to WAV. Implemented the reading from WAV and range compression, that is, the correlation between the received signal and the transmitted signal. Implemented the cross correlation using FFT. Did some final automatic analysis and plots.

## How to use it
First install Julia and then install the dependencies. Then...

``` shell
# To generate the chirp WAV file, you can use the following:

$ julia gen_chirp_wav.jl

# You can use a bash shell with Sox, the audio Suisse army nice,
# on Linux to play the Chirp WAV file:
$ play chirp_signal.wav

# You can use a bash shell on Linux to record the echos into a WAV
# file:
$ sox -d echos_3.wav 

# To analyze the echos WAV file, you can use the following:
$ julia audio_radar.jl
```

## License
MIT Open Source License

## Have fun
Best regards, <br>
João Carvalho
