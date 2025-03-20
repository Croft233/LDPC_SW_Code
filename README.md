# Slepian-Wolf with LDPC Codes

This repository presents an enhanced study where Low-Density Parity-Check (LDPC) codes are integrated into the Slepian-Wolf coding framework. This approach leverages LDPC's powerful error-correcting capabilities to further improve compression efficiency and approach the theoretical Slepian-Wolf bound.

## Overview

- **Integration of LDPC Codes:**  
  LDPC codes, known for their sparse parity-check matrices and efficient decoding (via the sum-product algorithm), are incorporated into the Slepian-Wolf framework to boost performance.

- **Objectives:**  
  - To implement and evaluate an LDPC-based Slepian-Wolf coding scheme.
  - To analyze how the integration of LDPC codes helps in approaching the theoretical Slepian-Wolf bound with minimal error margins.

- **Key Components:**  
  - **LDPC Code Theory:** Explanation of approximate lower triangular (ALT) encoding and sum-product decoding.
  - **Combined Scheme:** Details on how the parity-check matrix is split and used to encode correlated sources separately and decode jointly.
  - **Simulation Results:** MATLAB-based simulations demonstrating the performance improvements and closeness to the Slepian-Wolf bound.

## Repository Contents

- **Theory & Implementation:**  
  Documentation on LDPC encoding/decoding and its integration with Slepian-Wolf coding.
- **Simulation Scripts:**  
  MATLAB scripts that simulate the LDPC-enhanced Slepian-Wolf coding process, including BER performance analysis.

## How to Run the Simulations

**Requirements:**  
   MATLAB or GNU Octave.


## References

1. Slepian, D., & Wolf, J. K. (1973). Noiseless coding of correlated information sources.
2. Liveris, A. D., Xiong, Z., & Georghiades, C. N. (2002). Compression of binary sources with side information at the decoder using LDPC codes.
3. Additional references as cited in the research paper.
