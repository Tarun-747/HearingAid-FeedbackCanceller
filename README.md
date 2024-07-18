# Delayless Multiband-Structured Subband Adaptive Feedback Canceller for Digital Hearing Aids

## Introduction
This project focuses on improving acoustic feedback cancellation in digital hearing aids using a proportionate delayless multiband-structured subband adaptive feedback canceller. Acoustic feedback is a common issue in behind-the-ear (BTE) digital hearing aids, which can significantly degrade the performance of the device. This project implements an adaptive filter to address this issue effectively.

## Key Features
- **Adaptive Feedback Cancellation**: Utilizes an adaptive filter to model and cancel the acoustic feedback path, enhancing the performance of the hearing aid.
- **Prediction Error Method (PEM)**: Applies PEM to reduce bias effects caused by the finite correlation between the microphone input signal and the loudspeaker input signal.
- **Subband Implementation**: Implements the adaptive filter in the subband domain to improve convergence rates and reduce computational load.
- **Delayless Multiband-Structured Subband**: Incorporates a delayless subband approach to minimize aliasing, band-edge effects, and signal path delay.

## System Overview
The system is designed to continuously estimate and cancel the acoustic feedback signal in a digital hearing aid. It comprises several components, including:
- **Microphone**: Captures the sound to be amplified.
- **Amplifier**: Amplifies the sound captured by the microphone.
- **Loudspeaker**: Outputs the amplified sound.
- **Adaptive Filter**: Models and cancels the feedback signal using the PEM-based approach.

## Implementation Details
### PEM-Based Feedback Cancellation
- The PEM-based feedback cancellation scheme whitens the input and error signals using a prediction error filter, reducing the correlation and improving the feedback canceller's performance.

### Delayless Subband Implementation
- The adaptive filter operates in the subband domain, splitting the input signal into multiple subbands using an analysis filter bank. The subband signals are then processed and combined using a synthesis filter bank.

### Improved Pro
