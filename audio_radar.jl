# Name        : Acoustic RADAR in Julia
#
# Description : Small experiments in Julia with audio signals and RADAR.
#               Implemented the chirp generation and saving to WAV.
#               Implemented the reading from WAV and range compression, that
#               is, the correlation between the received signal and the
#               transmitted signal. Implemented the cross correlation using
#               FFT. Did some final automatic analysis and plots.
#
# Author      : João Carvalho 
# Data        : 2023.10.28
# License     : MIT Open Source License
#
# Notes       : To generate the chirp WAV file, you can use the following:
#               $ julia gen_chirp_wav.jl
#
#               You can use a bash shell with Sox, the audio Suisse army nice,
#               on Linux to play the Chirp WAV file:
#               $ play chirp_signal.wav
#
#               You can use a bash shell on Linux to record the echos into a WAV
#               file:
#               $ sox -d echos_3.wav 
#
#               To analyze the echos WAV file, you can use the following:
#               $ julia audio_radar.jl


using WAV
using LinearAlgebra
using Plots
using FFTW

using DSP  # Digital signal processing library


# using Statistics

using StatsBase


# load_filename = "./echos.wav"
load_filename = "./echos_2.wav"


VELOCITY_OF_SOUND = 343.0  # m/s

function chirp_signal( chirp_duration_sec :: Float64 ) :: Vector{ Float64 }

    # Parameters
    fs       = 48000    # Sample rate in Hz
    duration = chirp_duration_sec  # Duration of the chirp in seconds
    f_start  = 1000.0   # Start frequency in Hz
    f_end    = 15000.0  # End frequency in Hz

    # Generate the time vector for the chirp signal
    t = 0:1/fs:duration-1/fs

    # Generate the chirp signal
    k = (f_end - f_start) / duration
    chirp_signal = sin.(2π * (f_start .* t .+ 0.5 .* k .* t .^ 2))
    chirp_signal
end

function echos_read_wav_into_mono( filename :: String ) :: Tuple{ Vector{ Float64 }, Int64 }

    # Read the Wave File
    data, fs, nbits, opt = wavread( filename, format="double" )

    # data = vec(data)  # Convert multi-dimensional array to a vector

    # Just process the 1st channel
    echos = data[ : , 1 ] 

    println("echos size: ", size( echos ) )


    println( typeof( echos ) )

    return ( echos, fs )
end

function cross_correlation_fft( x :: Vector{ Float64 },
                                y :: Vector{ Float64 } ) :: Vector{ Float64 }
    # Zero-padding
    n = length(x) + length(y) - 1
    x_padded = vcat(x, zeros(n - length(x)))
    y_padded = vcat(y, zeros(n - length(y)))

    # FFT of both signals
    X = fft(x_padded)
    Y = fft(y_padded)

    # Multiply FFT of x with complex conjugate of FFT of y
    result = X .* conj.(Y)

    # Inverse FFT to get cross-correlation
    cross_corr = real(ifft(result))

    return cross_corr
end

function radar_range_compression( received_signal :: Vector{ Float64 },
                                  transmitted_chirp :: Vector{ Float64 } )  :: Vector{ Float64 }
    # # Time reverse the transmitted chirp for matched filtering
    # filter_signal = reverse( transmitted_chirp )
    
    # # Convolve the received signal with the matched filter
    # compressed_signal = DSP.conv( received_signal, filter_signal ) # , mode="valid" )
    

    println( "received_signal: ", size( received_signal ) )
    println( "transmitted_chirp: ", size( transmitted_chirp ) )

    println( "received_signal length: ",  length( received_signal ) )
    println( "transmitted_chirp length: ",  length( transmitted_chirp ) )


    # transmitted_chirp = copy(transmitted_chirp)

    if length( received_signal ) > length( transmitted_chirp )
        pad_length = length( received_signal ) - length( transmitted_chirp )
        transmitted_chirp = vcat( transmitted_chirp, zeros( pad_length ) )
    else
        pad_length = length( transmitted_chirp ) - length( received_signal )
        received_signal = vcat( received_signal, zeros( pad_length ) )
    end

    compressed_signal = crosscor( received_signal, transmitted_chirp )

    compressed_signal
end


function sec_2_samples( sec :: Float64, fs :: Int64 = 48_000 ) :: Int64
    Int64( floor( sec * fs ) )
end

function samples_2_sec( sample_number :: Int64, fs :: Int64 = 48_000 ) :: Float64
    Float64( sample_number * (1.0 / fs ) ) 
end

function samples_2_distance( sample_number :: Int64, fs :: Int64 = 48_000 ) :: Float64
    Float64( sample_number * (1.0 / fs ) * VELOCITY_OF_SOUND ) 
end


# echos, fs = echos_read_wav_into_mono( "./echos.wav" )
# echos, fs = echos_read_wav_into_mono( "./echos_2.wav" )
echos, fs = echos_read_wav_into_mono( load_filename )


# Clean the first 20_000 samples, from the starting value.
# echos[ 1 : 20_000 ] .= 0.0

echos[ 1 : 40_000 ] .= 0.0


chirp_duration_sec  = 0.3  # seconds
transmitted_chirp = chirp_signal( chirp_duration_sec )
received_signal = echos
# compressed_range = radar_range_compression( received_signal, transmitted_chirp )

compressed_range = cross_correlation_fft( received_signal, transmitted_chirp )


# Generate the time vector for the chirp signal
t = 0 : 1 / fs : ( length( echos ) - 1 ) / fs


# Plot the chirp signal
p1 = plot( t[ 1 : 200_000 ], echos[ 1 : 200_000 ], xlabel="Time (seconds)", ylabel="Amplitude", title="Chirp Signal", legend=false)
display( p1 )


readline()


println( "compressed_range: ", size( compressed_range ) )

println( "compressed_range length: ",  length( compressed_range ) )





# Generate the time vector for the chirp signal
t = 0 : 1 / fs : ( length( compressed_range ) - 1 ) / fs


# Plot the chirp signal
p2 = plot( t, compressed_range, xlabel="Time (seconds)", ylabel="Amplitude", title="Compressed Range Signal", legend=false)
display( p2 )


readline()






# Plot the chirp signal
# p = plot( t[ 1 : 200_000 ], compressed_range[ 1 : 200_000 ], xlabel="Time (seconds)", ylabel="Amplitude", title="Chirp Signal", legend=false)



# Generate the time vector for the chirp signal
t = 0 : 1 / fs : ( length( compressed_range ) - 1 ) / fs



# p3 = plot( t[ sample_f0 : sample_f1 ],
#          compressed_range[ sample_f0 : sample_f1 ],
#          xlabel="Time (seconds)", ylabel="Amplitude", title="Chirp Signal", legend=false )

p3 = plot( t,
         compressed_range,
         xlabel="Time (seconds)", ylabel="Amplitude", title="Compressed Range Signal", legend=false )



# Find indices of elements greater than 6.3
# indices = findall( x -> x > 6.5, compressed_range )
max_val = maximum( compressed_range )
# indices = findall( x -> x > max_val - 0.3, compressed_range )

# indices = findall( x -> x > ( max_val - 2.3 ), compressed_range )

indices = findall( x -> x > ( max_val - 5.3 ), compressed_range )

for _ in 1:length( indices)
    print( indices )
end


# indices = indices .- ( sample_f0 )

# Get the actual values using the indices
values_greater_than_6_5 = compressed_range[ indices ] 


x_points = ( indices .- 1 ) ./ fs
y_points = values_greater_than_6_5
# y_points = zeros( length( x_points ) )
# display( p3 )
scatter!( x_points, y_points, label="Marked Points", color=:red, marker=:circle )

display( p3 )


readline()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test
res = sec_2_samples( 1.0 )
println( "sec_2_samples( 1 ) res: ", res )

res = samples_2_sec( 48_000 )
println( "samples_2_sec( 48_000 ) res: ", res )



# Align to first echo.
first_echo_index = indices[ 1 ]

distance_start_window_to_first_echo_seconds = 0.0015
distance_start_window_to_first_echo_samples = sec_2_samples( distance_start_window_to_first_echo_seconds ) 
before_first_echo_instant = samples_2_sec( first_echo_index ) - distance_start_window_to_first_echo_seconds  # 0.025
before_first_echo_index   = sec_2_samples( before_first_echo_instant )  

delta_first_echo_seconds = 0.003  # 0.008
after_first_echo_instant = before_first_echo_instant + distance_start_window_to_first_echo_seconds + delta_first_echo_seconds
after_first_echo_index   = sec_2_samples( after_first_echo_instant )  

# compressed_range[ before_first_echo_index : after_first_echo_index ]




# sample_f0 = sec_2_samples( 1.130 ) # 1.425 )
# # sample_f1 = sec_2_samples( 1.2 ) # 1.5 )
# sample_f1 = sec_2_samples( 1.135 )

# max_val_2 = maximum( compressed_range[ sample_f0 + 120 : sample_f0 + 150 ] )
# indices_2 = findall( x -> x > max_val_2 - 0.01, compressed_range[ sample_f0 + 120 : sample_f0 + 200 ] )

# println( "indices maximum cabeça (1.130 + 120 sample,  1.130 + 150 samples) : ", indices_2 )

compressed_range_window = compressed_range[ before_first_echo_index : after_first_echo_index ]

distance_start_window_to_first_echo = distance_start_window_to_first_echo_samples 

println( "distance_start_window_to_fist_echo : ", distance_start_window_to_first_echo ) 

println( "length compressed_range_window : ", length( compressed_range_window ) )

compressed_range_window_sub_view = compressed_range_window[ distance_start_window_to_first_echo_samples + 10 : distance_start_window_to_first_echo_samples + 65 ] # 45 ] 


println( "length compressed_range_window_sub_view : ", length( compressed_range_window_sub_view ) )

max_val_2 = maximum( compressed_range_window_sub_view )
# indices_2 = findall( x -> x > max_val_2 - 0.1, compressed_range_window_sub_view )
indices_2 = findall( x -> x >= max_val_2 , compressed_range_window_sub_view )

println( "indices maximum in range sub window : ", indices )

# Generate the time vector for the chirp signal
t = 0 : ( 1 / fs ) * VELOCITY_OF_SOUND : (( length( compressed_range_window ) - 1 ) / fs) * VELOCITY_OF_SOUND



println("")

indexs_final = [ indices[ 1 ] - before_first_echo_index + 1  ]

push!( indexs_final, indices_2[1] + distance_start_window_to_first_echo + 10 )        
println( "indexs final maximum in range window distance xpto cm : ", indexs_final )

peak_values = [ compressed_range_window[ indexs_final[ 1 ] ], 
                compressed_range_window[ indexs_final[ 2 ] ] ]

distance_between_points = round( samples_2_distance( indexs_final[ 2 ] - indexs_final[ 1 ] ) * 100, digits=2 )
println( "Peak values ", distance_between_points, " cm : ", peak_values )



range_cm_samples_28 = 0

p4 = plot( t,
        # t[ before_first_echo_index : after_first_echo_index ],
        compressed_range_window,
        xlabel="meters", ylabel="Amplitude", title= "Compressed Range Microphone Testa " * string( distance_between_points) * " cm", legend=false )


display( p4 )


readline()
        

# indexs_final = [ indices[ 1 ] - before_first_echo_index + 1  ]

# push!( indexs_final, indices_2[1] + distance_start_window_to_first_echo + 10 )        
# println( "indexs final maximum in range window distance ", range_cm_samples_28, " cm : ", indexs_final )

# peak_values = [ compressed_range_window[ indexs_final[ 1 ] ], 
#                 compressed_range_window[ indexs_final[ 2 ] ] ]

# distance_between_points = round( samples_2_distance( indexs_final[ 2 ] - indexs_final[ 1 ] ) * 100, digits=2 )
# println( "Peak values ", distance_between_points, " cm : ", peak_values )


# x_points = ( ( indices .+ before_first_echo_index .- 1 ) ./ fs ) .* VELOCITY_OF_SOUND
x_points = ( ( 1 / fs )  .* ( indexs_final .- 1 )  ) .* VELOCITY_OF_SOUND
y_points = peak_values
# y_points = zeros( length( x_points ) )

scatter!( x_points, y_points, label="Marked Points", color=:red, marker=:circle )
        
display( p4 )


readline()


exit()




















####################################################################################
####################################################################################
####################################################################################



# VELOCITY_OF_SOUND = 343.0  # m/s

sample_f0 = sec_2_samples( 1.130 ) # 1.425 )
# sample_f1 = sec_2_samples( 1.2 ) # 1.5 )
sample_f1 = sec_2_samples( 1.135 )

max_val_2 = maximum( compressed_range[ sample_f0 + 120 : sample_f0 + 150 ] )
indices_2 = findall( x -> x > max_val_2 - 0.01, compressed_range[ sample_f0 + 120 : sample_f0 + 200 ] )

println( "indices maximum cabeça (1.130 + 120 sample,  1.130 + 150 samples) : ", indices_2 )


max_val = maximum( compressed_range[ sample_f0 : sample_f1 ] )
indices = findall( x -> x > max_val - 0.1, compressed_range[ sample_f0 : sample_f1 ] )

println( "indices maximum in range 1.1 sec, 1.2 sec : ", indices )

# Generate the time vector for the chirp signal
t = 0 : ( 1 / fs ) * velocity_of_sound : (( length( compressed_range ) - 1 ) / fs) * velocity_of_sound



println("")

distance_30cm = 0.3 # 30cm

num_samples_range_30cm = Int( floor( distance_30cm * ( fs / velocity_of_sound ) ) )

delta_samples_mic_testa = 28

range_cm_samples_28 = ( delta_samples_mic_testa * ( 1 / fs ) ) * velocity_of_sound * 100 




p4 = plot( t[ sample_f0 : sample_f1 ],
        compressed_range[ sample_f0 : sample_f1 ],
        xlabel="meters", ylabel="Amplitude", title= "Compressed Range Microphone Testa " * string( round( range_cm_samples_28, digits=3 ) ) , legend=false )



# println("")

# distance_30cm = 0.3 # 30cm

# num_samples_range_30cm = Int( floor( distance_30cm * ( fs / velocity_of_sound ) ) )

# delta_samples_mic_testa = 28

# range_cm_samples_28 = ( delta_samples_mic_testa * ( 1 / fs ) ) * velocity_of_sound * 100 



push!( indices, indices[1] + delta_samples_mic_testa ) # ) num_samples_range_30cm )        
println( "indices maximum in range 1.1 sec, 1.2 sec distance ", range_cm_samples_28, " cm : ", indices )

values_greater_than_6_5 = values_greater_than_6_5[ 1 : 1 ]

push!( values_greater_than_6_5, compressed_range[ sample_f0 + indices[2] ]  ) 
println( "values_distance 30 cm : ", values_greater_than_6_5 )



x_points = ( ( indices .+ sample_f0 .- 1 ) ./ fs ) .* velocity_of_sound
y_points = values_greater_than_6_5
# y_points = zeros( length( x_points ) )
# display( p3 )

# scatter!( x_points, y_points, label="Marked Points", color=:red, marker=:circle )
        

display( p4 )


readline()


exit()




