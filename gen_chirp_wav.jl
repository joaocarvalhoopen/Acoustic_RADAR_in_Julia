using WAV
using Plots

# pyplot()  # Set PyPlot as the backend

function chirp_wav_write(
      filename:: String,
      chirp_duration_sec:: Float64,
      silence_duration_sec:: Float64,
      num_repeat:: Int )

    # Parameters
    fs       = 48000  # Sample rate in Hz
    duration = chirp_duration_sec  # Duration of the chirp in seconds
    f_start  = 1000.0  # Start frequency in Hz
    f_end    = 15000.0  # End frequency in Hz

    # Generate the time vector for the chirp signal
    t = 0:1/fs:duration-1/fs

    # Generate the chirp signal
    k = (f_end - f_start) / duration
    chirp_signal = sin.(2π * (f_start .* t .+ 0.5 .* k .* t .^ 2))

    final_signal = zeros( 0 )

    for i = 1:num_repeat
        append!( final_signal, chirp_signal )
        append!( final_signal, zeros( Int( silence_duration_sec * fs ) ) )
    end


    # Generate the time vector for the chirp signal
    t_final = 1:1/fs:length( final_signal ) /fs


    # Plot the chirp signal
    p = plot(t_final[1:500], final_signal[1:500], xlabel="Time (seconds)", ylabel="Amplitude", title="Chirp Signal", legend=false)
    display( p )


    # Ask user to press a key to continue
    println("Press Enter to continue...")
    readline()



    # # Prepare stereo data: chirp on the left channel and zeros on the right channel
    stereo_data = hcat( final_signal, zeros( length( final_signal ) ) )

    # # Save to WAV file
    wavwrite(stereo_data, filename, Fs=fs)

end

chirp_duration_sec  = 0.3 # seconds
silent_duration_sec = 0.7 # seconds
total_duration_sec  = 30 # seconds
chirp_wav_write( "chirp_signal.wav",
                 chirp_duration_sec,
                 silent_duration_sec,
                 total_duration_sec )