# -*- coding: utf-8 -*-
"""
Created on Sun May 20 18:13:46 2018

@author: ravi
"""

import parselmouth
import numpy as np
import pylab

def __get_formants__(audio_file):
    my_sound = parselmouth.Sound(audio_file)
    my_formants = my_sound.to_formant_burg(time_step=0.001)
    num_frames = my_formants.get_number_of_frames()
    formants_array = np.zeros((num_frames,5))
    start_time = my_formants.start_time
    time_step = my_formants.get_time_step()
    
    for i in range(num_frames):
        formants_array[i,0] = my_formants.get_value_at_time(1, start_time + i*time_step)
        formants_array[i,1] = my_formants.get_value_at_time(2, start_time + i*time_step)
        formants_array[i,2] = my_formants.get_value_at_time(3, start_time + i*time_step)
        formants_array[i,3] = my_formants.get_value_at_time(4, start_time + i*time_step)
        formants_array[i,4] = my_formants.get_value_at_time(5, start_time + i*time_step)
    
    return formants_array
    
if __name__ == '__main__':
    file_loc_angry = '/home/ravi/Desktop/Spectrogram_Registration/angry.wav'
    file_loc_neutral = '/home/ravi/Desktop/Spectrogram_Registration/neutral.wav'
    
    formants_angry = __get_formants__(file_loc_angry)
    formants_neutral = __get_formants__(file_loc_neutral)
    
    pylab.figure()
    pylab.subplot(211), pylab.plot(formants_angry[:,0], label='Angry'), pylab.legend(), pylab.title('F1')
    pylab.subplot(212), pylab.plot(formants_neutral[:,0], label='Neutral'), pylab.legend()
    
    pylab.figure()    
    pylab.subplot(211), pylab.plot(formants_angry[:,1], label='Angry'), pylab.legend(), pylab.title('F2')
    pylab.subplot(212), pylab.plot(formants_neutral[:,1], label='Neutral'), pylab.legend()
    
    pylab.figure()    
    pylab.subplot(211), pylab.plot(formants_angry[:,2], label='Angry'), pylab.legend(), pylab.title('F3')
    pylab.subplot(212), pylab.plot(formants_neutral[:,2], label='Neutral'), pylab.legend()

    pylab.figure()    
    pylab.subplot(211), pylab.plot(formants_angry[:,3], label='Angry'), pylab.legend(), pylab.title('F4')
    pylab.subplot(212), pylab.plot(formants_neutral[:,3], label='Neutral'), pylab.legend()    

    pylab.figure()    
    pylab.subplot(211), pylab.plot(formants_angry[:,4], label='Angry'), pylab.legend(), pylab.title('F5')
    pylab.subplot(212), pylab.plot(formants_neutral[:,4], label='Neutral'), pylab.legend()
