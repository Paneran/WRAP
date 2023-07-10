## WRAP
Wireless RF Analog Project (WRAP) is a wireless communications project taught at the IEEE at UCLA branch. Students learn to wirelessly send a string of text over 100m of air, impmlementing ideas from all across EE, including circuit design, controls, communication theory, and signal processing. All ideas are intertwined and can be divided into two main components: The circuit design (hardware) and digital signal processing (software). 

# System Design and This Repository
A microcontroller encodes a bit string as a 1MHz signal and outputs it to the transmitter hardware, which upconverts the signal from 1-27MHz signal. The 27MHz transmits across the channel to the receiver board, which downconverts the signal back to 1MHz and another microcontroller interprets the signal and recovers the transmitted bits. This repository contains the R&D for the microcontroller component of this project. The work here was then translated into lectures and assignments found in this [Google Drive](https://drive.google.com/drive/folders/1hFtTVYYO03o8BIRFq2AxVKLdrDG3jyCn?usp=sharing).

First students simulated the link in MATLAB (MATLAB folder). The two main files of interest are "link_sim.m" and "receiver_scope_data.m". link_sim is a simulation from transmitter to receiver, including non-idealities like noise and attenuation. Receiver_scope_data is the receiver processing data collected from an oscilloscope from a transmitter and receiver board. This is importantn for tuning gain values for Costas Loop and Timing Recovery. 

After students were able to transmit bits in simulation, they moved to programming it in C (C folder). C_testing contains programs that are the MATLAB code but in C. The final C files used the built in STM32 libraries for convolutions, rather than the ones I wrote in tools.c to reduce overhead because tools.c is too slow. The wrap_receiver_2/ folder contains all the necessary files for the STM32 CubeIDE to run it on the MCU. 

# Demo
https://github.com/Paneran/Wireless_RF_Analog_Project/assets/32348542/9b494989-9275-4c86-b272-cef50629c1d1

