# Phase change models in Basilisk solver #

## Sharp interface models
This phase change model includes:  
- mass flux calcualtion  
    - [reference](https://github.com/GabrielGLK/thesis-cases/blob/master/reference-papers/phase-change/mass-transfer-model/sharp-interface-model/simpified-sharp-model/Sun%20%E7%AD%89%E3%80%82%20-%202012%20-%20Development%20of%20a%20vapor%E2%80%93liquid%20phase%20change%20model%20f.pdf)
    - [code](https://github.com/GabrielGLK/phase-change/blob/master/phase-change/phase_change_models/sharp-interface-model/film-boiling/phase_change_code/phase-change.h#L143)
    - Feature:  
The method to calcualte temperature gradient is different, we use similar method by [Quentin](http://basilisk.fr/sandbox/qmagdelaine/phase_change/elementary_body.h).