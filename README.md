# VSR
Velocity Selective Recording using a TDNN.

Velocity Selective Recording (VSR) has been proposed as a technique of extracting more information from implantable neural interfaces. Re- cent research has shown that a Time Delay Neural Network (TDNN) has the potential to outperform traditional VSR techniques, but requires a time consuming optimisation process. This project explores a series of analyses performed in order to understand the operation of this model in VSR, in order to eliminate the need for optimisation and potentially improve its performance for low Signal-to-Noise Ratio (SNR) inputs.

I focused on time response analysis of the TDNN, which revealed that it converges to a highly sophisticated version of a traditional method, called delay-and-sum. The key to the performance of the model lies in the timings of the neural signals from multiple channels, creating constructive and destructive interference patterns. It was also shown that the reason that the model fails with real data extracted from a rat, is that these recordings suffer from DC distortion and neural pulse interference, which severely affect the performance of this specific model.

Finally, the theories formed during the analysis of the results were validated by designing FIR filters to simulate the neural signal processing that the TDNN accomplishes. The designs did not reach the level of performance that the TDNN achieved, but important steps were made towards the final goal and future work is proposed.
